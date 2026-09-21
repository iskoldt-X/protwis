"""Build one chainmap.tsv per PDB, from files only.

Same decisions as build_schrodinger_chain_map: the per-anchor and per-receptor
resolution is interaction.schrodinger_chain_map, imported by both. The branch
for an anchor whose chain_res names no residue is written out in both commands
and has to be kept in step by hand. Every input here comes from files instead of
the database, and the output is one file per structure, written next to that
structure's products so the two can never drift apart::

    python manage.py build_schrodinger_chainmap_files \\
        --gpcrdb-data /app/data/gpcrdb_data --cif-dir /app/data/cif \\
        --data-dir /app/data/schrodinger --annotation-commit 9fe1875

Where each input comes from:

    database                              file
    ------------------------------------  ---------------------------------------
    StructureLigandInteraction            structure_data/annotation/ligands.tsv
      pdb_reference / chain_res             Name / Residue_seq_id
      ligand.ligand_type.slug               Type
    Structure.preferred_chain             structure_data/annotation/structures.tsv
                                            ChainID (column 6)
    Structure.pdb_data.pdb                structure_data/pdbs/<PDB>.pdb
    structure_type.origin == experiment   structures.tsv holds experimental only

The command issues no database query and imports no model of its own; the
importer module it borrows constants and instance discovery from pulls in
Django models at import time, but nothing here touches the ORM.

Writes <out-dir>/<PDB>/chainmap.tsv, one per PDB, each carrying its own
provenance header. No field enumerates what else was in the run, so trees built
at different times can be merged by copying directories; the headers then differ
between files, which is how a merged tree is recognised.

An empty body means the annotation lists no ligand Engine 1 serves for that
structure. It does not mean the engine looked and found none: the body is built
from ligands.tsv, never from the product tree.
"""

import csv
import hashlib
import inspect
import io
import os

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError

from interaction import schrodinger_chain_map as cm
from interaction import schrodinger_import as si

# The format constants live in schrodinger_import, the module that reads these
# files back; one definition, so a bump cannot land on only one side.
SCHEMA = si.CHAINMAP_SCHEMA
RECEPTOR_PREFIX = si.CHAINMAP_RECEPTOR_PREFIX


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _sha256_parts(parts):
    """sha256 over a list, length-prefixed so the parts cannot be re-cut.

    One of the parts is a function body and carries newlines, so joining on a
    separator would not be injective.
    """
    h = hashlib.sha256()
    for part in parts:
        raw = part.encode("utf-8")
        h.update(str(len(raw)).encode("ascii") + b":" + raw)
    return h.hexdigest()


# The header keys that are the same for every structure of a run, plus the two
# that are not, in the order they are written.
HEADER_KEYS = ("annotation_commit", "ligands_sha256", "structures_sha256",
               "gpcrdb_pdb_sha256", "cif_sha256", "builder_sha256")


def has_product_summary(data_dir, pdb):
    """Did the producer leave its per-structure summary in this directory?

    The only evidence in the delivery that Engine 1 ever looked at a structure.
    Read from the product tree, never from the out dir: it is a fact about the
    producer, so building into a separate directory must not change the answer.
    """
    path = os.path.join(data_dir, pdb, si.PRODUCT_SUMMARY_NAME)
    try:
        # An empty file is a truncated copy, not a run. The whole judgement
        # rests on this one answer, so it does not accept a file that says
        # nothing.
        return os.path.isfile(path) and os.path.getsize(path) > 0
    except OSError:
        return False


def chainmap_header(pdb, values, receptor, has_summary):
    """The '# key<TAB>value' lines of one chainmap, in order.

    Every field is a statement about this one structure or about an input that
    was the same for every structure, so merging two runs by copying
    directories leaves every file still telling the truth about itself.
    """
    header = [("schema", SCHEMA), ("pdb", pdb)]
    header += [(key, values[key]) for key in HEADER_KEYS]
    header.append((si.PRODUCT_SUMMARY_KEY, "yes" if has_summary else "no"))
    header += [(RECEPTOR_PREFIX + c, receptor.get(c, "")) for c in cm.RECEPTOR_COLUMNS
               if c != "pdb"]
    return header


def _some(names, limit=20):
    """A list for a log line, honest about what it left out."""
    head = ", ".join(names[:limit])
    return head if len(names) <= limit else "{} (+{} more)".format(head, len(names) - limit)


def _flat(value):
    """A header value may hold free text (note); keep it on one line, one field.

    Only the characters that would break the format are replaced, so the value
    still reads as what it was.
    """
    text = str(value if value is not None else "")
    for bad in ("\t", "\r", "\n"):
        text = text.replace(bad, " ")
    return text


def read_tsv(path):
    """Rows of a gpcrdb_data annotation TSV, keys and values stripped.

    A row wider than the header arrives under the key None as a list: refuse
    the file rather than guess which field went where. A row that is too short
    is padded with empty strings, which drops the anchor from the map and makes
    the importer refuse that structure -- loud, and the safe direction.
    """
    rows = []
    with open(path, newline="") as fh:
        for n, r in enumerate(csv.DictReader(fh, delimiter="\t"), start=2):
            if None in r:
                raise CommandError("{} line {}: more fields than the header has columns".format(
                    path, n))
            rows.append({(k or "").strip(): (v or "").strip() for k, v in r.items()})
    return rows


def preferred_chains(rows):
    """PDB -> preferred chain, as structure.functions.ParseStructureCSV stores it.

    That parser keeps only the first character of a chain id containing a dot,
    which is how a numeric chain id has been seen to arrive; resolve_receptor
    then takes the part before the first comma. Neither branch fires on the
    annotation as it stands -- they are here because the database side applies
    them, and the two must agree whatever the file holds.
    """
    out = {}
    for r in rows:
        pdb = (r.get("PDB") or "").upper()
        chain = r.get("ChainID") or ""
        if "." in chain:
            chain = chain[0]
        if pdb:
            out[pdb] = chain
    return out


def annotation_anchors(rows):
    """PDB -> [(HET, token, chain_res)] in file order, for the anchors Engine 1 serves.

    Mirrors schrodinger_import.is_in_scope: a real chemical component whose
    annotated type is one Engine 1 serves. The (HET, token) pairs are the same
    set check_map_covers compares against the database, so the ligand_role axis
    -- which the annotation does not have -- never enters.
    """
    out = {}
    seen = set()
    for r in rows:
        pdb = (r.get("PDB") or "").upper()
        het = (r.get("Name") or "").strip().upper()
        if not pdb or not het or het in si.PLACEHOLDER_REFERENCES:
            continue
        # The exact annotation spelling, as is_in_scope reads the slug. A new
        # upstream spelling (for instance "small molecule" with a space) drops
        # the anchor here while the database keeps it, which the importer then
        # refuses loudly -- the safe direction.
        if (r.get("Type") or "") not in si.IN_SCOPE_LIGAND_TYPES:
            continue
        chain_res = r.get("Residue_seq_id") or ""
        for tok in cm.split_tokens(chain_res) or [""]:
            if (pdb, het, tok) in seen:
                continue
            seen.add((pdb, het, tok))
            out.setdefault(pdb, []).append((het, tok, chain_res))
    return out


def write_chainmap(path, header, anchor_rows):
    buf = io.StringIO()
    for key, value in header:
        buf.write("# {}\t{}\n".format(key, _flat(value)))
    w = csv.DictWriter(buf, fieldnames=cm.ANCHOR_COLUMNS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(anchor_rows)
    # The out dir is the product tree that gets delivered; never leave a half
    # written map behind if the run dies on the next structure.
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        fh.write(buf.getvalue())
    os.replace(tmp, path)


class Command(BaseCommand):
    help = "Build one per-PDB chainmap.tsv for the Schrodinger importer, from files only."

    def add_arguments(self, parser):
        parser.add_argument("--gpcrdb-data", default=None,
                            help="gpcrdb_data checkout (uses structure_data/annotation and "
                                 "structure_data/pdbs). Defaults to DATA_DIR, so the maps are "
                                 "built from the same annotation the build will run on.")
        parser.add_argument("--cif-dir", required=True, help="Directory of the input mmCIFs, <PDB>.cif.")
        parser.add_argument("--data-dir", required=True, help="Product tree {data_dir}/{PDB}/{instance}/.")
        parser.add_argument("--annotation-commit", required=True,
                            help="gpcrdb_data commit the annotation was read at.")
        parser.add_argument("--out-dir", default=None,
                            help="Where to write {PDB}/chainmap.tsv; defaults to --data-dir, "
                                 "so the map ships with the products.")
        parser.add_argument("--allow-stray", action="store_true",
                            help="Do not refuse product directories that are absent from "
                                 "structures.tsv. They still get no chainmap, and the "
                                 "importer fails on any of them the database still knows, "
                                 "so delete them before delivering the tree.")
        parser.add_argument("--pdb", action="append", default=[],
                            help="Restrict to these PDB codes; default is every structure in "
                                 "structures.tsv.")

    def handle(self, *args, **opt):
        gdata = opt["gpcrdb_data"] or settings.DATA_DIR
        cif_dir, data_dir = opt["cif_dir"], opt["data_dir"]
        out_dir = opt["out_dir"] or data_dir
        for label, path in (("--gpcrdb-data", gdata), ("--cif-dir", cif_dir), ("--data-dir", data_dir)):
            if not os.path.isdir(path):
                raise CommandError("{} {!r} is not a directory".format(label, path))
        # The out dir itself is created; its parent is not, so a typo there
        # would otherwise build a whole tree somewhere nobody looks.
        out_parent = os.path.dirname(os.path.abspath(out_dir))
        if not os.path.isdir(out_parent):
            raise CommandError("--out-dir {!r}: {} does not exist".format(out_dir, out_parent))
        ann = os.path.join(gdata, "structure_data", "annotation")
        ligands_tsv = os.path.join(ann, "ligands.tsv")
        structures_tsv = os.path.join(ann, "structures.tsv")
        pdb_dir = os.path.join(gdata, "structure_data", "pdbs")
        for path in (ligands_tsv, structures_tsv):
            if not os.path.isfile(path):
                raise CommandError("{} is missing under --gpcrdb-data".format(path))
        if not os.path.isdir(pdb_dir):
            raise CommandError("{} is missing under --gpcrdb-data".format(pdb_dir))

        ligand_rows = read_tsv(ligands_tsv)
        labels = cm.annotation_labels(ligand_rows)
        anchors = annotation_anchors(ligand_rows)
        chains = preferred_chains(read_tsv(structures_tsv))

        # The corpus is structures.tsv, the same list the producer runs, not
        # the directories that happen to sit under --data-dir. Driving it off
        # the tree would skip a structure whose run yielded nothing -- exactly
        # the structures whose anchors have to be cleared -- and would write a
        # chainmap into any stray directory, .git included.
        corpus = sorted(chains)
        if not corpus:
            raise CommandError("{} lists no structures".format(structures_tsv))
        if opt["pdb"]:
            pdbs = sorted({p.upper() for p in opt["pdb"]})
            unknown = [p for p in pdbs if p not in chains]
            if unknown:
                raise CommandError("not in {}: {}".format(structures_tsv, ", ".join(unknown)))
        else:
            pdbs = corpus
            # A directory holding products but absent from the corpus would
            # ship with no chainmap.tsv at all. Housekeeping directories (.git,
            # logs) hold no instance and are only worth a quiet note -- keeping
            # the two apart is the point, so that the day a real structure
            # lands in the list it is not read as noise.
            dropped, housekeeping = [], []
            for name in sorted(os.listdir(data_dir)):
                if name.upper() in chains or not os.path.isdir(os.path.join(data_dir, name)):
                    continue
                (dropped if si.instance_yaml_paths(data_dir, name) else housekeeping).append(name)
            if housekeeping:
                self.stdout.write("directories that are not structures, ignored: {}".format(
                    _some(housekeeping)))
            if dropped and not opt["allow_stray"]:
                raise CommandError(
                    "these directories hold Engine 1 products but are not in {}, so they would "
                    "ship without a chainmap: {}. Pass --allow-stray if that is intended, for "
                    "instance after a structure was retired from the annotation.".format(
                        structures_tsv, _some(dropped)))
            if dropped:
                self.stdout.write("product directories with no chainmap (--allow-stray): "
                                  "{}".format(_some(dropped)))

        ligands_sha = _sha256_file(ligands_tsv)
        structures_sha = _sha256_file(structures_tsv)
        # What produced these rows: the algorithm module, this command, and the
        # four things it borrows from the importer module. Hashing the whole of
        # that module instead would move this stamp on every unrelated importer
        # edit, and a merged tree would then report a difference that is not
        # one.
        builder_sha = _sha256_parts([
            _sha256_file(cm.__file__), _sha256_file(os.path.abspath(__file__)),
            repr(sorted(si.IN_SCOPE_LIGAND_TYPES)), repr(sorted(si.PLACEHOLDER_REFERENCES)),
            si.INSTANCE_DIR_RE.pattern, repr(si.INSTANCE_DIR_RE.flags),
            inspect.getsource(si.instance_yaml_paths),
            si.CHAINMAP_SCHEMA, si.CHAINMAP_RECEPTOR_PREFIX, si.PRODUCT_SUMMARY_NAME,
            si.PRODUCT_SUMMARY_KEY])

        counts, rstatus, written, no_products = {}, {}, 0, 0
        for pdb in pdbs:
            has_summary = has_product_summary(data_dir, pdb)
            rows, receptor, note = self.build_one(
                pdb, anchors.get(pdb, []), chains.get(pdb), labels,
                os.path.join(cif_dir, pdb + ".cif"), os.path.join(pdb_dir, pdb + ".pdb"),
                si.instance_yaml_paths(data_dir, pdb), has_summary)
            no_products += receptor["product_instances_sha256"] == cm.instances_sha256([])
            header = chainmap_header(pdb, dict(
                annotation_commit=opt["annotation_commit"], ligands_sha256=ligands_sha,
                structures_sha256=structures_sha, builder_sha256=builder_sha,
                # The bytes on disk. receptor.gpcrdb_text_sha256 is what the
                # importer compares against the database, hashed as decoded
                # text; the two differ only on a CRLF file.
                gpcrdb_pdb_sha256=note["gpcrdb_pdb_sha256"],
                cif_sha256=note["cif_sha256"]), receptor, has_summary)
            target = os.path.join(out_dir, pdb)
            os.makedirs(target, exist_ok=True)
            write_chainmap(os.path.join(target, "chainmap.tsv"), header, rows)
            written += 1
            for r in rows:
                counts[(r["status"], r["source"])] = counts.get((r["status"], r["source"]), 0) + 1
            key = (receptor["status"], receptor["method"])
            rstatus[key] = rstatus.get(key, 0) + 1

        self.stdout.write("out-dir {}".format(out_dir))
        self.stdout.write("annotation_commit {} ligands_sha256 {} structures_sha256 {} "
                          "builder_sha256 {}".format(opt["annotation_commit"], ligands_sha,
                                                     structures_sha, builder_sha))
        self.stdout.write("chainmap.tsv written: {} ({} with no product instance)".format(
            written, no_products))
        self.stdout.write("anchor rows {}: {}".format(sum(counts.values()), sorted(counts.items())))
        self.stdout.write("receptor rows {}: {}".format(sum(rstatus.values()), sorted(rstatus.items())))

    def build_one(self, pdb, anchor_keys, preferred_chain, labels, cif_path, gpcrdb_pdb_path,
                  instances, has_summary):
        """(anchor_rows, receptor_row, provenance) for one structure.

        An unreadable input is not a reason to skip: every anchor is written as
        unresolved with the reason, and the importer refuses the structure. A
        silently missing structure would instead look like one with no anchors.
        """
        note = {"cif_sha256": "", "gpcrdb_pdb_sha256": "", "product_summary": has_summary}
        if preferred_chain is None:
            preferred_chain = ""
        def unresolved(exc):
            reason = "input unreadable: {}: {}".format(type(exc).__name__, exc)[:200]
            rows = [dict({c: "" for c in cm.ANCHOR_COLUMNS}, pdb=pdb, het=het, token=tok,
                         status="unresolved", note=reason)
                    for het, tok, _ in anchor_keys]
            receptor = dict({c: "" for c in cm.RECEPTOR_COLUMNS}, pdb=pdb,
                            preferred_chain=preferred_chain, status="unresolved", note=reason,
                            gpcrdb_text_sha256="",
                            product_instances_sha256=cm.instances_sha256(instances))
            return rows, receptor, note

        try:
            note["cif_sha256"] = _sha256_file(cif_path)
            with open(cif_path) as fh:
                cif_text = fh.read()
            note["gpcrdb_pdb_sha256"] = _sha256_file(gpcrdb_pdb_path)
            with open(gpcrdb_pdb_path) as fh:
                gtext = fh.read()
        except (OSError, UnicodeDecodeError) as exc:
            return unresolved(exc)
        try:
            # Only the parsers get the wide clause: parse_mmcif_atoms raises a
            # bare KeyError for an absent _atom_site column and a ValueError
            # for an unparsable coordinate, and one malformed input must cost
            # one structure rather than the whole run. Anything raised outside
            # these two calls is a bug and is left to surface as one.
            cif_atoms = cm.parse_mmcif_atoms(cif_text)
            gatoms = cm.parse_gpcrdb_pdb(gtext)
        except (cm.ParseError, KeyError, ValueError) as exc:
            return unresolved(exc)

        receptor = cm.resolve_receptor(pdb, preferred_chain, cif_atoms, gatoms)
        receptor["gpcrdb_text_sha256"] = cm.text_sha256(gtext)
        receptor["product_instances_sha256"] = cm.instances_sha256(instances)

        rows = []
        for het, tok, chain_res in anchor_keys:
            if tok:
                rows.append(cm.resolve_anchor(pdb, het, tok, cif_atoms, gatoms, instances,
                                              labels.get((pdb, het, tok))))
                continue
            copies = sorted(n for n in instances if n.split("_", 1)[0].upper() == het)
            rows.append(dict(
                {c: "" for c in cm.ANCHOR_COLUMNS}, pdb=pdb, het=het, token="",
                instance=";".join(copies),
                status="all_copies" if copies else "no_product",
                note=("chain_res {!r} names no residue; every copy used".format(chain_res)
                      if copies else
                      "the product has no instance of {} (chain_res {!r})".format(het, chain_res))))
        return rows, receptor, note
