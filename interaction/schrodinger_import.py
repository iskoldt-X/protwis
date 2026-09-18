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

Replacement semantics: for every in-scope anchor the existing RFI rows are
deleted and the Schrodinger rows written in their place, even when the
Schrodinger result is empty. Anchors outside Engine 1 scope (peptide, protein
and placeholder ligands) are never touched.
"""

import collections
import glob
import os
import re

import yaml
from django.db import transaction

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

CHAIN_RES_RE = re.compile(r"^(?P<chain>[A-Za-z0-9]+):(?P<resnum>-?\d+)(?P<icode>[A-Za-z]?)$")


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


def parse_chain_res(chain_res):
    """Parse an SLI ``chain_res`` of the form ``CHAIN:RESNUM[ICODE]``.

    Returns (chain, resnum_str, icode) or None when the value is empty or
    names a chain only.
    """
    m = CHAIN_RES_RE.match((chain_res or "").strip())
    if not m:
        return None
    return m.group("chain"), m.group("resnum"), m.group("icode")


def select_instances(het_code, chain_res, instance_names):
    """Choose the product instance(s) that belong to one SLI anchor.

    When the anchor names an exact residue (``chain_res`` = ``A:408``), only
    the instance with that chain and residue number is used: the curated
    anchor points at one copy, and other copies of the same HET may be
    crystal contacts or a second site with its own anchor.

    When ``chain_res`` is empty or a bare chain, every instance of the HET is
    used and the preferred-chain filter decides later which rows are kept.

    Returns (names, mode) with mode one of ``exact``, ``exact_missing``,
    ``all_copies``.
    """
    het = het_code.upper()
    copies = sorted(n for n in instance_names if n.split("_", 1)[0].upper() == het)
    parsed = parse_chain_res(chain_res)
    if parsed is None:
        return copies, "all_copies"
    chain, resnum, icode = parsed
    wanted = "{}_{}_{}{}".format(het, chain, resnum, icode)
    if wanted in copies:
        return [wanted], "exact"
    return [], "exact_missing"


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


def preferred_chain_of(structure_preferred_chain):
    """First chain of ``Structure.preferred_chain`` ('' keeps every chain)."""
    return (structure_preferred_chain or "").split(",")[0].strip()


def plan_rows(interactions, preferred_chain):
    """Plan the RFI rows for the product rows of one anchor.

    Returns (records, counts, other_chain_by_chain). ``records`` are dicts
    with ``sequence_number``, ``amino_acid``, ``slug`` and
    ``receptor_pdb_block``, one per distinct (sequence_number, slug).
    ``counts`` accounts for every input row:

        rows_in == excluded_family + nonstandard_residue + other_chain
                   + duplicate + len(records)

    ``other_chain`` rows are those on a chain other than the preferred one:
    GPCRdb residues carry no chain, so they cannot be stored.
    Raises UnroutableRow for a (family, direction) the map does not know.
    """
    counts = collections.Counter()
    counts["rows_in"] = len(interactions)
    other_chain_by_chain = collections.Counter()
    records = []
    seen = set()
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
        if preferred_chain and chain != preferred_chain:
            counts["other_chain"] += 1
            other_chain_by_chain[chain] += 1
            continue
        slug = apply_backbone_override(
            resolve_slug(family, row.get("direction")), row.get("receptor_atom_name"))
        seq = int(res["pdb_residue_number"])
        if (seq, slug) in seen:
            counts["duplicate"] += 1
            continue
        seen.add((seq, slug))
        records.append({
            "sequence_number": seq,
            "amino_acid": amino_acid,
            "slug": slug,
            "receptor_pdb_block": row.get("receptor_pdb_block") or "",
        })
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

    __slots__ = ("sli_id", "het", "mode", "instances", "deleted", "written", "counts",
                 "other_chain_by_chain", "dropped")

    def __init__(self, sli_id, het):
        self.sli_id = sli_id
        self.het = het
        self.mode = ""
        self.instances = []
        self.deleted = 0
        self.written = 0
        self.counts = collections.Counter()
        self.other_chain_by_chain = {}
        self.dropped = collections.Counter()


def import_structure(structure, data_dir):
    """Replace the Engine 1 RFI rows of one structure.

    Runs in one transaction: any exception leaves the structure exactly as it
    was. Returns (outcomes, out_of_scope_count). Rows planned but not written
    because the database has no matching residue or rotamer are counted in
    ``outcome.dropped``; for every anchor

        counts["planned"] == written + sum(dropped.values())
    """
    pdb_code = structure.pdb_code.index.upper()
    preferred = preferred_chain_of(structure.preferred_chain)
    instances = instance_yaml_paths(data_dir, pdb_code)
    types = {t.slug: t for t in ResidueFragmentInteractionType.objects.all()}

    outcomes = []
    out_of_scope = 0
    with transaction.atomic():
        slis = (StructureLigandInteraction.objects
                .filter(structure=structure)
                .select_related("ligand__ligand_type")
                .order_by("id"))
        for sli in slis:
            if not is_in_scope(sli):
                out_of_scope += 1
                continue
            outcome = AnchorOutcome(sli.id, sli.pdb_reference.upper())
            names, outcome.mode = select_instances(sli.pdb_reference, sli.chain_res, instances)
            outcome.instances = names

            rows = []
            for name in names:
                rows.extend(read_instance_rows(instances[name]))
            records, outcome.counts, outcome.other_chain_by_chain = plan_rows(rows, preferred)

            _, deleted_by_model = ResidueFragmentInteraction.objects.filter(
                structure_ligand_pair=sli).delete()
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
                fragment = (Fragment.objects
                            .filter(ligand=sli.ligand, structure=structure, residue=residue)
                            .order_by("id").first())
                if fragment is None:
                    fragment = Fragment.objects.create(
                        ligand=sli.ligand, structure=structure, residue=residue,
                        pdbdata=PdbData.objects.create(pdb=rec["receptor_pdb_block"]))
                ResidueFragmentInteraction.objects.create(
                    structure_ligand_pair=sli,
                    rotamer=rotamers[0],
                    fragment=fragment,
                    interaction_type=types[rec["slug"]],
                )
                outcome.written += 1
            outcomes.append(outcome)
    return outcomes, out_of_scope
