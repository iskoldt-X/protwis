"""Import Schrodinger Engine 1 (small-molecule) interactions into GPCRdb.

The schrodinger_interaction pipeline writes one YAML per ligand instance:

    {data_dir}/{PDB}/{HET}_{chain}_{resnum}{icode}/{HET}_{chain}_{resnum}{icode}.yaml

with the rows under ``doc["result"]["interactions"]``. This module turns those
rows into ResidueFragmentInteraction (RFI) rows hung on the curated
StructureLigandInteraction (SLI) anchors.

Two layers:

* a pure layer (no database access) that selects the instance(s) for an
  anchor and plans the rows to write, with an exact account of every row
  that is not written and why;
* a database layer that replaces the RFI rows of one structure inside one
  transaction.

Instance selection and the receptor chain come from two maps built offline by
build_schrodinger_chain_map (interaction/schrodinger_chain_map.py): product
chain names are mmCIF author names, GPCRdb's are its own stored PDB-format
names, and the two differ for a few dozen structures.

Replacement semantics (ADR-089, ADR-091): every in-scope anchor loses all its
existing rows. Anchors with a product instance get the Schrodinger rows;
anchors without one (the map says no_product) are left empty and reported
with the map's reason. Anchors outside Engine 1 scope (peptide, protein and
placeholder ligands) are never touched. Fragments left unreferenced are
deleted, and so is their PdbData text when nothing else references it.
"""

import collections
import glob
import os
import re

import csv

import yaml
from django.db import transaction

from interaction import schrodinger_chain_map as chain_map

from interaction.models import (
    ResidueFragmentInteraction,
    ResidueFragmentInteractionType,
    StructureLigandInteraction,
)
from residue.models import Residue
from structure.models import Fragment, PdbData, Rotamer


# ---------------------------------------------------------------------------
# Scope
# ---------------------------------------------------------------------------

# Ligand types whose anchors Engine 1 serves. Peptide and protein ligands
# belong to Engine 2, which is frozen; their rows are left as they are.
IN_SCOPE_LIGAND_TYPES = frozenset({"small-molecule", "lipid"})

# pdb_reference values that name no chemical component (lower case in the
# 2026-09 database, upper case in older dumps).
PLACEHOLDER_REFERENCES = frozenset({"PEP", "APO"})

# Families never written to the database. Wat-HBond: the bridging water is not
# recorded in the product, so no row of this family can be checked downstream.
EXCLUDED_FAMILIES = frozenset({"Wat-HBond"})

STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

# Structure types this import serves (team convention: experiment only).
STRUCTURE_ORIGIN = "experiment"


# ---------------------------------------------------------------------------
# Instance discovery and selection (pure)
# ---------------------------------------------------------------------------

INSTANCE_DIR_RE = re.compile(
    r"^(?P<het>[^_]+)_(?P<chain>[^_]+)_(?P<resnum>-?\d+)(?P<icode>[A-Za-z]?)$"
)

def instance_yaml_paths(data_dir, pdb_code):
    """Map instance name -> YAML path for one PDB (flat product layout).

    Only directories whose name parses as an instance and that hold a YAML of
    the same name are returned. A missing PDB directory returns {}.
    """
    top = os.path.join(data_dir, pdb_code.upper())
    found = {}
    for inst_dir in sorted(glob.glob(os.path.join(top, "*"))):
        name = os.path.basename(inst_dir)
        if not os.path.isdir(inst_dir) or not INSTANCE_DIR_RE.match(name):
            continue
        path = os.path.join(inst_dir, name + ".yaml")
        if os.path.isfile(path):
            found[name] = path
    return found


class MapMismatch(ValueError):
    """The chain maps do not describe the anchors in the database."""


class UnresolvedAnchor(ValueError):
    """The anchor map could not decide which product instance an anchor is."""


def _read_map(path):
    header, body = {}, []
    with open(path, newline="") as fh:
        for line in fh:
            # Header lines only precede the column row. A body field may hold a
            # quoted newline whose continuation begins with "# "; without this
            # guard such a line would be read as a header key, and the header
            # is where the schema and the receptor row live.
            if not body and line.startswith("# "):
                key, _, value = line[2:].rstrip("\n").partition("\t")
                header[key] = value
            else:
                body.append(line)
    return header, list(csv.DictReader(body, delimiter="\t"))


def load_anchor_map(path):
    """(header, {(pdb, HET, token): row}) from anchor_instance_map.tsv."""
    header, rows = _read_map(path)
    table = {}
    for r in rows:
        key = (r["pdb"].upper(), r["het"].upper(), r["token"])
        if key in table:
            raise MapMismatch("duplicate anchor map key {}".format(key))
        table[key] = r
    return header, table


def load_receptor_map(path):
    """(header, {pdb: row}) from receptor_chain_map.tsv."""
    header, rows = _read_map(path)
    table = {}
    for r in rows:
        key = r["pdb"].upper()
        if key in table:
            raise MapMismatch("duplicate receptor map key {}".format(key))
        table[key] = r
    return header, table


# ---------------------------------------------------------------------------
# Per-PDB chain maps
# ---------------------------------------------------------------------------

CHAINMAP_NAME = "chainmap.tsv"

# The only chainmap layout this importer reads. A file naming another one is
# refused, never guessed at: the header carries the receptor row, so a layout
# change can move a field without changing any column name.
CHAINMAP_SCHEMA = "engine1-chainmap/1"

# The single receptor row is per structure, so it rides in the header under
# this prefix instead of being a one-row table. Both the builder and the
# reader take it from here; two literals would let a rename pass unnoticed.
CHAINMAP_RECEPTOR_PREFIX = "receptor."

# The producer writes one of these per structure it processed, so its presence
# is the evidence that a run happened. Without it an empty directory cannot be
# told from one the builder created for a structure nobody ran.
PRODUCT_SUMMARY_NAME = "summary.yaml"

# Header keys that say where a chainmap came from. They are expected to differ
# across a tree merged from several production runs; the importer reports the
# distinct values rather than insisting on one.
PROVENANCE_KEYS = ("annotation_commit", "ligands_sha256", "structures_sha256",
                   "builder_sha256")


def load_chainmap(path):
    """(pdb, {(pdb, HET, token): row}, receptor_row, provenance) from one file."""
    header, rows = _read_map(path)
    schema = header.get("schema")
    if schema != CHAINMAP_SCHEMA:
        raise MapMismatch("{}: chainmap schema {!r}, this importer reads {!r}".format(
            path, schema, CHAINMAP_SCHEMA))
    pdb = (header.get("pdb") or "").strip().upper()
    if not pdb:
        raise MapMismatch("{}: the header names no pdb".format(path))
    receptor = {"pdb": pdb}
    for column in chain_map.RECEPTOR_COLUMNS:
        if column == "pdb":
            continue
        key = CHAINMAP_RECEPTOR_PREFIX + column
        if key not in header:
            raise MapMismatch("{}: the header has no {}".format(path, key))
        receptor[column] = header[key]
    anchors = {}
    for r in rows:
        if (r.get("pdb") or "").strip().upper() != pdb:
            raise MapMismatch("{}: a row for {!r} in the chainmap of {}".format(
                path, r.get("pdb"), pdb))
        key = (pdb, (r.get("het") or "").upper(), r.get("token") or "")
        if key in anchors:
            raise MapMismatch("{}: duplicate anchor key {}".format(path, key))
        anchors[key] = r
    provenance = {k: header.get(k, "") for k in PROVENANCE_KEYS}
    return pdb, anchors, receptor, provenance


def load_chainmap_dir(data_dir, pdb_codes):
    """Per-PDB chain maps from {data_dir}/{PDB}/chainmap.tsv.

    Returns (anchor_table, receptor_table, missing, provenance). The two tables
    are shaped exactly like load_anchor_map / load_receptor_map, so
    import_structure does not know which format it was given. `missing` lists
    the PDB codes with no chainmap.tsv: they are reported and failed one at a
    time rather than aborting the run, so one absent file names itself instead
    of hiding every other. `provenance` counts the distinct values of each
    PROVENANCE_KEY, which a tree merged from several runs will show as more
    than one.
    """
    anchors, receptors, missing = {}, {}, []
    provenance = dict((k, {}) for k in PROVENANCE_KEYS)
    for pdb in pdb_codes:
        pdb = pdb.upper()
        path = os.path.join(data_dir, pdb, CHAINMAP_NAME)
        if not os.path.isfile(path):
            missing.append(pdb)
            continue
        named, rows, receptor, prov = load_chainmap(path)
        if named != pdb:
            raise MapMismatch("{}: names pdb {} but sits in the directory of {}".format(
                path, named, pdb))
        anchors.update(rows)
        receptors[pdb] = receptor
        for key, value in prov.items():
            provenance[key][value] = provenance[key].get(value, 0) + 1
    return anchors, receptors, missing, provenance


def product_pdb_codes(data_dir):
    """Every structure directory of a delivered product tree, sorted.

    The tree holds one directory per structure that was run, whether or not it
    yielded a ligand, so this is the corpus of a delivery.
    """
    return sorted(name.upper() for name in os.listdir(data_dir)
                  if os.path.isdir(os.path.join(data_dir, name)))


# Map statuses that name a product instance to import.
IMPORT_STATUSES = frozenset({"ok", "errata"})


def anchor_instances(pdb, het, chain_res, anchor_map, instance_names):
    """Product instances for one anchor, from the anchor map.

    Returns (names, mode, notes). mode is ``mapped`` (every named copy has an
    instance), ``mapped_partial`` (some copies have none), ``all_copies``
    (chain_res names no residue; every copy of the HET) or ``no_product``.
    Raises MapMismatch when a copy is missing from the map, or when the map
    says no_product although ``instance_names`` (the product tree being
    imported) holds a copy of the HET -- a map built against another tree
    would otherwise clear anchors that do have a product. Raises
    UnresolvedAnchor when the map could not decide.
    """
    pdb, het = pdb.upper(), het.upper()
    has_copy = any(n.split("_", 1)[0].upper() == het for n in instance_names)
    tokens = chain_map.split_tokens(chain_res) or [""]
    names, notes, missing = [], [], 0
    for tok in tokens:
        row = anchor_map.get((pdb, het, tok))
        if row is None:
            raise MapMismatch("{} {} {!r} is not in the anchor map".format(pdb, het, tok))
        status = row["status"]
        if status == "unresolved":
            raise UnresolvedAnchor("{} {} {!r}: {}".format(pdb, het, tok, row["note"]))
        if status in IMPORT_STATUSES:
            names.append(row["instance"])
            if status == "errata":
                notes.append("errata: " + row["note"])
        elif status == "all_copies":
            names.extend(n for n in row["instance"].split(";") if n)
        elif status == "no_product":
            if has_copy:
                raise MapMismatch("{} {}: map says no_product but the product tree has "
                                  "a copy; the map was built against another tree".format(pdb, het))
            missing += 1
            notes.append(row["note"])
        else:
            raise MapMismatch("{} {} {!r}: unknown map status {!r}".format(pdb, het, tok, status))
    names = sorted(set(names))
    if tokens == [""]:
        return names, ("all_copies" if names else "no_product"), notes
    if not names:
        return [], "no_product", notes
    return names, ("mapped" if not missing else "mapped_partial"), notes


def instance_chains(pdb, het, chain_res, anchor_map, names):
    """GPCRdb chain for each selected instance of one anchor.

    Named copies take the chain of their chain_res token (GPCRdb naming; the
    token pattern allows one character only). Copies selected as all_copies
    keep their product chain, which must then be a single character; otherwise
    MapMismatch.
    """
    pdb, het = pdb.upper(), het.upper()
    out = {}
    for tok in chain_map.split_tokens(chain_res) or [""]:
        row = anchor_map.get((pdb, het, tok))
        if row is not None and row["status"] in IMPORT_STATUSES and tok:
            out[row["instance"]] = tok.split(":", 1)[0]
    for name in names:
        if name not in out:
            product_chain = name.split("_")[1]
            if len(product_chain) != 1:
                raise MapMismatch("{} {}: no GPCRdb chain for {} (multi-character product "
                                  "chain and no chain_res token)".format(pdb, het, name))
            out[name] = product_chain
    return out


def standard_ligand_block(block, instance, gpcrdb_chain):
    """Rewrite every atom line of a producer ligand block.

    Returns (text, capped_lines): the rewritten block and the output lines
    whose B-factor was capped.
    """
    m = INSTANCE_DIR_RE.match(instance)
    if not m:
        raise MalformedProduct("not an instance name: {!r}".format(instance))
    out, capped = [], []
    for line in (block or "").splitlines():
        if not line.strip():
            continue
        new, was_capped = standard_ligand_line(line, m.group("het"), m.group("chain"),
                                               m.group("resnum"), m.group("icode"), gpcrdb_chain)
        out.append(new)
        if was_capped:
            capped.append(new)
    return "\n".join(out), capped


def receptor_chain(pdb, receptor_map):
    """The product chain that is GPCRdb's preferred chain, from the receptor map."""
    row = receptor_map.get(pdb.upper())
    if row is None:
        raise MapMismatch("{} is not in the receptor map".format(pdb))
    if row["status"] != "ok":
        raise UnresolvedAnchor("{} receptor chain {}: {}".format(pdb, row["status"], row["note"]))
    if not row["auth_chain"]:
        raise MapMismatch("{}: receptor map row is ok but names no chain".format(pdb))
    return row["auth_chain"]


def check_fingerprints(pdb, receptor_map, gpcrdb_text, instance_names):
    """Refuse a structure whose stored text or product instances changed since the build."""
    row = receptor_map.get(pdb.upper())
    if row is None:
        raise MapMismatch("{} is not in the receptor map".format(pdb))
    if row.get("gpcrdb_text_sha256") != chain_map.text_sha256(gpcrdb_text):
        raise MapMismatch("{}: GPCRdb structure text differs from the one the map was built "
                          "from (new dump?); rebuild the maps".format(pdb))
    if row.get("product_instances_sha256") != chain_map.instances_sha256(instance_names):
        raise MapMismatch("{}: product instances differ from the tree the map was built "
                          "from; check --data-dir or rebuild the maps".format(pdb))


# ---------------------------------------------------------------------------
# Row routing (pure)
# ---------------------------------------------------------------------------

_TYPE_MAP = None


def load_type_map():
    """(feature_family, direction) -> slug, from interaction_type_map.yaml."""
    global _TYPE_MAP
    if _TYPE_MAP is None:
        path = os.path.join(os.path.dirname(__file__), "interaction_type_map.yaml")
        with open(path) as fh:
            doc = yaml.safe_load(fh)
        _TYPE_MAP = {
            (rule["family"], rule["direction"]): rule["slug"]
            for rule in doc.get("rules") or []
        }
    return _TYPE_MAP


def required_slugs():
    """Slugs this import can write (map targets minus excluded families)."""
    return frozenset(
        slug for (family, _), slug in load_type_map().items()
        if family not in EXCLUDED_FAMILIES
    ) | {"polar_backbone"}


class UnroutableRow(ValueError):
    """A product row whose (family, direction) has no rule in the map."""


def resolve_slug(feature_family, direction):
    try:
        return load_type_map()[(feature_family, direction or "")]
    except KeyError:
        raise UnroutableRow(
            "no interaction_type_map rule for family={!r} direction={!r}".format(
                feature_family, direction))


_BACKBONE_PROMOTABLE = frozenset({"polar_donor_protein", "polar_acceptor_protein"})
_BACKBONE_ATOMS = frozenset({"N", "O"})


def apply_backbone_override(slug, receptor_atom_name):
    """H-bond to a main-chain N or O becomes ``polar_backbone``."""
    if slug in _BACKBONE_PROMOTABLE and (receptor_atom_name or "").strip() in _BACKBONE_ATOMS:
        return "polar_backbone"
    return slug


def plan_rows(interactions, receptor_chain_name):
    """Plan the RFI rows for the product rows of one anchor.

    ``receptor_chain_name`` is the product (mmCIF author) chain that GPCRdb
    stores as the receptor; '' keeps every chain.

    Returns (records, counts, other_chain_by_chain). ``records`` are dicts
    with ``sequence_number``, ``amino_acid``, ``slug`` and ``ligand_lines``
    (the ligand atom lines of every product row collapsed into the record, in
    first-seen order), one per distinct (sequence_number, slug).
    ``counts`` accounts for every input row:

        rows_in == excluded_family + nonstandard_residue + other_chain
                   + duplicate + len(records)

    ``other_chain`` rows are those on a chain other than the receptor chain:
    GPCRdb residues carry no chain, so they cannot be stored.
    Raises UnroutableRow for a (family, direction) the map does not know.
    """
    counts = collections.Counter()
    counts["rows_in"] = len(interactions)
    other_chain_by_chain = collections.Counter()
    records = []
    by_key = {}
    for row in interactions:
        family = row["feature_family"]
        if family in EXCLUDED_FAMILIES:
            counts["excluded_family"] += 1
            continue
        res = row["receptor_residue"]
        amino_acid = str(res["name_1_letter"]).upper()
        if amino_acid not in STANDARD_AMINO_ACIDS:
            counts["nonstandard_residue"] += 1
            continue
        chain = str(res["chain_id"])
        if receptor_chain_name and chain != receptor_chain_name:
            counts["other_chain"] += 1
            other_chain_by_chain[chain] += 1
            continue
        slug = apply_backbone_override(
            resolve_slug(family, row.get("direction")), row.get("receptor_atom_name"))
        seq = int(res["pdb_residue_number"])
        ligand_lines = [line for line in (row.get("ligand_pdb_block") or "").splitlines()
                        if line.strip()]
        record = by_key.get((seq, slug))
        if record is not None:
            counts["duplicate"] += 1
            for line in ligand_lines:
                if line not in record["ligand_lines"]:
                    record["ligand_lines"].append(line)
            continue
        record = {
            "sequence_number": seq,
            "amino_acid": amino_acid,
            "slug": slug,
            "ligand_lines": [],
        }
        for line in ligand_lines:
            if line not in record["ligand_lines"]:
                record["ligand_lines"].append(line)
        by_key[(seq, slug)] = record
        records.append(record)
    counts["planned"] = len(records)
    return records, counts, dict(other_chain_by_chain)


class MalformedProduct(ValueError):
    """A product YAML that cannot be read or lacks result.interactions."""


def read_instance_rows(path):
    """Return the interaction rows of one instance YAML."""
    try:
        with open(path) as fh:
            doc = yaml.safe_load(fh)
    except (OSError, yaml.YAMLError) as exc:
        raise MalformedProduct("{}: {}".format(path, exc))
    result = (doc or {}).get("result") if isinstance(doc, dict) else None
    if not isinstance(result, dict) or not isinstance(result.get("interactions"), list):
        raise MalformedProduct("{}: no result.interactions list".format(path))
    return result["interactions"]


# ---------------------------------------------------------------------------
# Anchors (database)
# ---------------------------------------------------------------------------

def is_in_scope(sli):
    """True iff this SLI anchor is served by Engine 1."""
    reference = (sli.pdb_reference or "").strip().upper()
    if not reference or reference in PLACEHOLDER_REFERENCES:
        return False
    return sli.ligand.ligand_type.slug in IN_SCOPE_LIGAND_TYPES


# ---------------------------------------------------------------------------
# Import of one structure (database)
# ---------------------------------------------------------------------------

class AnchorOutcome(object):
    """What happened to one in-scope anchor."""

    __slots__ = ("sli_id", "het", "mode", "instances", "notes", "deleted", "written",
                 "counts", "other_chain_by_chain", "dropped", "fragments_created")

    def __init__(self, sli_id, het):
        self.sli_id = sli_id
        self.het = het
        self.mode = ""
        self.instances = []
        self.notes = []
        self.deleted = 0
        self.written = 0
        self.counts = collections.Counter()
        self.other_chain_by_chain = {}
        self.dropped = collections.Counter()
        self.fragments_created = 0


class UnexpectedCascade(RuntimeError):
    """A delete removed objects of a model it was not meant to touch."""


def _only_deleted(deleted_by_model, allowed):
    extra = {k: v for k, v in deleted_by_model.items() if v and k not in allowed}
    if extra:
        raise UnexpectedCascade("delete also removed {}".format(extra))


# ---------------------------------------------------------------------------
# Ligand atom lines -> standard PDB columns (pure)
# ---------------------------------------------------------------------------

# The producer writes ligand atom lines one column short of the PDB standard
# (no altloc column) and widens them for five-character CCD codes and
# multi-character chains. The numeric fields always carry fixed decimals
# (coordinates 3, occupancy and B-factor 2), so they can be read from the
# tail of the line even where widened fields run into each other.
_LIGAND_TAIL_RE = re.compile(
    r"^\s*(?P<x>-?\d+\.\d{3})\s*(?P<y>-?\d+\.\d{3})\s*(?P<z>-?\d+\.\d{3})"
    r"\s*(?P<occ>-?\d+\.\d{2})\s*(?P<b>-?\d+\.\d{2})\s+(?P<element>[A-Za-z]{1,2})\s*$")

# Largest B-factor the standard 6-column field can hold.
_MAX_PDB_B = 999.99
_PDB_ATOM_LINE_WIDTH = 78


class MalformedLigandLine(MalformedProduct):
    """A ligand atom line that does not have the producer's layout."""


def _pdb_atom_name(name, element):
    """Columns 13-16: one-letter elements with short names start in column 14."""
    if len(name) < 4 and len(element) == 1:
        return " " + name.ljust(3)
    return name.ljust(4)


def standard_ligand_line(line, het, product_chain, resnum, icode, gpcrdb_chain):
    """Rewrite one producer ligand atom line in standard PDB v3.3 columns.

    The producer's residue name, chain and residue number are checked against
    the instance the line came from; any disagreement raises
    MalformedLigandLine instead of guessing. In the output the residue name is
    cut to three characters and the chain is GPCRdb's, as in GPCRdb's own
    stored structure text. Returns (standard_line, b_factor_was_capped).
    """
    record = line[:6].strip()
    if record not in ("ATOM", "HETATM"):
        raise MalformedLigandLine("not an atom record: {!r}".format(line[:30]))
    serial = line[6:11].strip()
    name = line[12:16].strip()
    width = max(3, len(het))
    pos = 16
    resname = line[pos:pos + width].strip()
    pos += width + 1
    chain = line[pos:pos + len(product_chain)]
    pos += len(product_chain)
    num_width = max(4, len(str(resnum)))
    number = line[pos:pos + num_width].strip()
    pos += num_width
    ins = line[pos:pos + 1].strip()
    pos += 1
    tail = _LIGAND_TAIL_RE.match(line[pos:])
    if (not serial.isdigit() or not name or resname.upper() != het.upper()
            or chain != product_chain or number != str(resnum) or ins != (icode or "")
            or tail is None):
        raise MalformedLigandLine("{} {}_{}_{}{}: cannot read {!r}".format(
            record, het, product_chain, resnum, icode, line))
    b = float(tail.group("b"))
    capped = b > _MAX_PDB_B
    element = tail.group("element").upper()
    out = "{:<6}{:>5} {} {:>3} {:1}{:>4}{:1}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2}".format(
        record, int(serial) % 100000, _pdb_atom_name(name, element), het.upper()[:3],
        gpcrdb_chain, resnum, icode or "",
        float(tail.group("x")), float(tail.group("y")), float(tail.group("z")),
        float(tail.group("occ")), min(b, _MAX_PDB_B), element)
    # A value too wide for its field (residue number beyond 4 digits, coordinate
    # beyond 8.3f, negative B below -99.99) would shift every later column.
    if len(out) != _PDB_ATOM_LINE_WIDTH or len(gpcrdb_chain) != 1:
        raise MalformedLigandLine("{} {}_{}_{}{}: a field does not fit standard PDB columns: {!r}".format(
            record, het, product_chain, resnum, icode, line))
    return out, capped


def fragment_text(ligand_lines):
    """Stored fragment text: the ligand atoms that make the contact."""
    return "\n".join(ligand_lines) + "\n" if ligand_lines else ""


def delete_orphan_fragments(structure):
    """Delete this structure's fragments that no interaction row references.

    A fragment's PdbData text is deleted too, but only when no other row of
    any model references it: PdbData is shared by rotamers, structures,
    anchors and others, all with on_delete=CASCADE. The referencing models
    are read from Django's model metadata, not from a hand-written list.
    Returns a Counter with fragments_deleted, pdbdata_deleted and
    pdbdata_kept_referenced.
    """
    out = collections.Counter()
    referenced = set(ResidueFragmentInteraction.objects
                     .filter(fragment__structure=structure)
                     .values_list("fragment_id", flat=True))
    orphans = [(f_id, pd_id) for f_id, pd_id in
               Fragment.objects.filter(structure=structure).values_list("id", "pdbdata_id")
               if f_id not in referenced]
    if not orphans:
        return out
    _, deleted = Fragment.objects.filter(id__in=[f for f, _ in orphans]).delete()
    _only_deleted(deleted, {"structure.Fragment"})
    out["fragments_deleted"] = deleted.get("structure.Fragment", 0)

    candidates = {pd for _, pd in orphans}
    still_used = set()
    for rel in PdbData._meta.related_objects:
        still_used.update(rel.related_model._base_manager
                          .filter(**{rel.field.name + "__in": candidates})
                          .values_list(rel.field.name, flat=True))
    free = candidates - still_used
    out["pdbdata_kept_referenced"] = len(candidates & still_used)
    if free:
        _, deleted = PdbData.objects.filter(id__in=free).delete()
        _only_deleted(deleted, {"structure.PdbData"})
        out["pdbdata_deleted"] = deleted.get("structure.PdbData", 0)
    return out


def check_map_covers(pdb, slis, anchor_map):
    """Every copy the database has must be listed; extra copies are allowed.

    Returns the copies the map lists that this database cannot use, sorted, so
    the caller can report them. They are never looked up -- anchor_instances
    only asks for the tokens of the chain_res the database stores -- but they
    are the visible cost of the subset rule, and silence about them would hide
    a ligand copy that was computed and then not imported.

    A map built by build_schrodinger_chain_map lists exactly the database's
    copies by construction: it walks the same StructureLigandInteraction rows
    through the same split_tokens() call. (It also filters experiment-origin
    structures, which the import command does separately.) The subset rule is
    for a map built from the annotation instead (ligands.tsv), where extras are
    normal: the annotation splits a ligand into physical copies, while
    StructureLigandInteraction is keyed on (structure, ligand, ligand_role,
    annotated) and has no copy dimension, so build_structures keeps one row per
    (structure, ligand, role).

    Measured on dump 20260917_phase2 against ligands.tsv at gpcrdb_data
    9fe1875, over the 1,672 in-scope structures: the file side is never short
    of a copy the database has (0 missing) and lists 373 extra copies across
    239 structures. Classifying each extra by whether its chain is one the
    database already uses for the same HET: 328 on another chain, 44 further
    copies on a chain the database does use (35 of them the calcium ion), and
    one whole HET -- 7IPG A1CS8, whose SMILES normalises to the same
    stereochemistry-stripped InChIKey as A1CQL, and Ligand.clean_inchikey is
    unique, so the two annotation rows end up on one ligand and one SLI.

    The cost of allowing extras: set-equality was the only per-structure
    witness that the SLI rows still match the ones the map was built from.
    check_fingerprints does not cover them -- it compares the stored structure
    text and the product instance names, not chain_res, pdb_reference or the
    ligand type. What still fails loud is any (het, token) the database holds
    that the map does not list. What now passes silently is drift that moves a
    copy onto a pair the map already lists: a chain rename onto a listed
    second-chain copy, or two HETs collapsing onto one ligand. Both then import
    what the database still has an anchor for -- in the collapsed case that
    means one real chemical entity is never imported at all -- and both show up
    in the returned unused copies rather than as an exception.
    """
    pdb = pdb.upper()
    wanted = set()
    for sli in slis:
        het = sli.pdb_reference.upper()
        for tok in chain_map.split_tokens(sli.chain_res) or [""]:
            wanted.add((het, tok))
    listed = {(het, tok) for (p, het, tok) in anchor_map if p == pdb}
    missing = wanted - listed
    if missing:
        raise MapMismatch("{}: the database has {} copies, the anchor map lists {} of them "
                          "(plus {} it cannot use); first absent: {}".format(
                              pdb, len(wanted), len(wanted) - len(missing),
                              len(listed - wanted), sorted(missing)[:5]))
    return sorted(listed - wanted)


def import_structure(structure, data_dir, anchor_map, receptor_map):
    """Replace the Engine 1 RFI rows of one structure.

    Runs in one transaction: any exception leaves the structure exactly as it
    was. Returns (outcomes, out_of_scope_count, cleanup_counter, unused_copies),
    where unused_copies are the map's (HET, token) pairs this database has no
    anchor for. Rows planned but not written because the database has no
    matching residue or rotamer are counted in ``outcome.dropped``; for every
    anchor

        counts["planned"] == written + sum(dropped.values())
    """
    pdb_code = structure.pdb_code.index.upper()
    instances = instance_yaml_paths(data_dir, pdb_code)
    types = {t.slug: t for t in ResidueFragmentInteractionType.objects.all()}

    outcomes = []
    out_of_scope = 0
    unused = []
    with transaction.atomic():
        slis = list(StructureLigandInteraction.objects
                    .filter(structure=structure)
                    .select_related("ligand__ligand_type")
                    .order_by("id"))
        in_scope = [sli for sli in slis if is_in_scope(sli)]
        out_of_scope = len(slis) - len(in_scope)
        if in_scope:
            unused = check_map_covers(pdb_code, in_scope, anchor_map)
            check_fingerprints(pdb_code, receptor_map,
                               structure.pdb_data.pdb if structure.pdb_data_id else "", instances)
            chain = receptor_chain(pdb_code, receptor_map)
        for sli in in_scope:
            outcome = AnchorOutcome(sli.id, sli.pdb_reference.upper())
            names, outcome.mode, outcome.notes = anchor_instances(
                pdb_code, sli.pdb_reference, sli.chain_res, anchor_map, instances)
            outcome.instances = names

            rows = []
            capped = set()
            gchains = instance_chains(pdb_code, sli.pdb_reference, sli.chain_res, anchor_map, names)
            for name in names:
                if name not in instances:
                    raise MalformedProduct("{}: the anchor map names {} but the product tree "
                                           "under {} has no such instance".format(pdb_code, name, data_dir))
                for row in read_instance_rows(instances[name]):
                    # Ligand atoms are stored in standard PDB columns (ADR-095).
                    row["ligand_pdb_block"], n = standard_ligand_block(
                        row.get("ligand_pdb_block"), name, gchains[name])
                    # The same block repeats on every interaction row of an
                    # instance; count each capped atom once.
                    capped.update((name, line) for line in n)
                    rows.append(row)
            records, outcome.counts, outcome.other_chain_by_chain = plan_rows(rows, chain)
            outcome.counts["ligand_lines_bfactor_capped"] = len(capped)

            _, deleted_by_model = ResidueFragmentInteraction.objects.filter(
                structure_ligand_pair=sli).delete()
            _only_deleted(deleted_by_model, {"interaction.ResidueFragmentInteraction"})
            outcome.deleted = deleted_by_model.get("interaction.ResidueFragmentInteraction", 0)

            for rec in records:
                residues = list(Residue.objects.filter(
                    protein_conformation=structure.protein_conformation,
                    sequence_number=rec["sequence_number"])[:2])
                if not residues:
                    outcome.dropped["residue_not_found"] += 1
                    continue
                if len(residues) > 1:
                    outcome.dropped["residue_ambiguous"] += 1
                    continue
                residue = residues[0]
                if residue.amino_acid != rec["amino_acid"]:
                    outcome.dropped["amino_acid_mismatch"] += 1
                    continue
                rotamers = list(Rotamer.objects.filter(structure=structure, residue=residue)[:2])
                if len(rotamers) != 1:
                    outcome.dropped["rotamer_not_found" if not rotamers else "rotamer_ambiguous"] += 1
                    continue
                # The fragment holds the ligand atoms of this contact. Reuse only
                # a fragment with exactly this text (one this import wrote
                # earlier); legacy fragments hold other content.
                text = fragment_text(rec["ligand_lines"])
                fragment = (Fragment.objects
                            .filter(ligand=sli.ligand, structure=structure, residue=residue,
                                    pdbdata__pdb=text)
                            .order_by("id").first())
                if fragment is None:
                    fragment = Fragment.objects.create(
                        ligand=sli.ligand, structure=structure, residue=residue,
                        pdbdata=PdbData.objects.create(pdb=text))
                    outcome.fragments_created += 1
                ResidueFragmentInteraction.objects.create(
                    structure_ligand_pair=sli,
                    rotamer=rotamers[0],
                    fragment=fragment,
                    interaction_type=types[rec["slug"]],
                )
                outcome.written += 1
            outcomes.append(outcome)
        cleanup = delete_orphan_fragments(structure) if in_scope else collections.Counter()
    return outcomes, out_of_scope, cleanup, unused
