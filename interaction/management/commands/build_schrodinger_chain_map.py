"""Build the chain-name maps the Schrodinger importer reads.

Usage::

    python manage.py build_schrodinger_chain_map \\
        --cif-dir /app/data/cif --data-dir /app/data/schrodinger \\
        --annotation /runs/ligands.tsv --annotation-commit b2af5d6 \\
        --dump-id 20260917_phase2 --out-dir /runs/chainmap [--pdb 2RH1 ...]

Writes anchor_instance_map.tsv and receptor_chain_map.tsv. Both start with
'#'-prefixed provenance lines (dump id, annotation commit and sha256, input
and product manifests, builder sha256). Read-only with respect to the
database. Rebuild after every new GPCRdb dump or product run, and diff the
result against the previous maps before using it.
"""

import csv
import hashlib
import os

from django.core.management.base import BaseCommand, CommandError

from interaction import schrodinger_chain_map as cm
from interaction import schrodinger_import as si
from interaction.models import StructureLigandInteraction
from structure.models import Structure


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _sha256_lines(lines):
    h = hashlib.sha256()
    for line in lines:
        h.update(line.encode("utf-8") + b"\n")
    return h.hexdigest()


class Command(BaseCommand):
    help = "Build the anchor and receptor chain-name maps for the Schrodinger importer."

    def add_arguments(self, parser):
        parser.add_argument("--cif-dir", required=True, help="Directory of the input mmCIFs, <PDB>.cif.")
        parser.add_argument("--data-dir", required=True, help="Product tree {data_dir}/{PDB}/{instance}/.")
        parser.add_argument("--annotation", required=True, help="Upstream ligands.tsv at a pinned commit.")
        parser.add_argument("--annotation-commit", required=True, help="gpcrdb_data commit of --annotation.")
        parser.add_argument("--dump-id", required=True, help="Name of the GPCRdb dump the database holds.")
        parser.add_argument("--out-dir", required=True)
        parser.add_argument("--pdb", action="append", default=[], help="Restrict to these PDB codes.")

    def handle(self, *args, **opt):
        for key in ("cif_dir", "data_dir"):
            if not os.path.isdir(opt[key]):
                raise CommandError("--{} {!r} is not a directory".format(key.replace("_", "-"), opt[key]))
        os.makedirs(opt["out_dir"], exist_ok=True)
        with open(opt["annotation"], newline="") as fh:
            rows = [{(k or "").strip(): (v or "").strip() for k, v in r.items()}
                    for r in csv.DictReader(fh, delimiter="\t")]
        labels = cm.annotation_labels(rows)

        wanted = {p.upper() for p in opt["pdb"]}
        anchors = {}
        for sli in (StructureLigandInteraction.objects
                    .select_related("structure__pdb_code", "structure__structure_type",
                                    "ligand__ligand_type")
                    .order_by("id")):
            if sli.structure is None or not si.is_in_scope(sli):
                continue
            if sli.structure.structure_type.origin != si.STRUCTURE_ORIGIN:
                continue
            pdb = sli.structure.pdb_code.index.upper()
            if wanted and pdb not in wanted:
                continue
            anchors.setdefault(pdb, []).append(sli)
        if wanted - set(anchors):
            raise CommandError("no in-scope anchors for: {}".format(", ".join(sorted(wanted - set(anchors)))))

        anchor_rows, receptor_rows, cif_manifest, product_manifest = [], [], [], []
        for pdb in sorted(anchors):
            structure = Structure.objects.select_related("pdb_data").get(pk=anchors[pdb][0].structure_id)
            instances = si.instance_yaml_paths(opt["data_dir"], pdb)
            product_manifest.extend("{}\t{}".format(pdb, name) for name in sorted(instances))
            cif_path = os.path.join(opt["cif_dir"], pdb + ".cif")
            gtext = structure.pdb_data.pdb if structure.pdb_data else ""
            try:
                cif_manifest.append("{}\t{}".format(pdb, _sha256_file(cif_path)))
                with open(cif_path) as fh:
                    cif_atoms = cm.parse_mmcif_atoms(fh.read())
                gatoms = cm.parse_gpcrdb_pdb(gtext)
            except (OSError, cm.ParseError) as exc:
                note = "input unreadable: {}".format(exc)[:200]
                for sli in anchors[pdb]:
                    for tok in cm.split_tokens(sli.chain_res) or [""]:
                        anchor_rows.append(dict({c: "" for c in cm.ANCHOR_COLUMNS}, pdb=pdb,
                                                het=sli.pdb_reference.upper(), token=tok,
                                                status="unresolved", note=note))
                receptor_rows.append(dict({c: "" for c in cm.RECEPTOR_COLUMNS}, pdb=pdb,
                                          preferred_chain=structure.preferred_chain or "",
                                          status="unresolved", note=note,
                                          gpcrdb_text_sha256=cm.text_sha256(gtext),
                                          product_instances_sha256=cm.instances_sha256(instances)))
                continue

            receptor_row = cm.resolve_receptor(pdb, structure.preferred_chain, cif_atoms, gatoms)
            # Per-structure fingerprints: the importer refuses a structure whose
            # stored text or product instance list no longer matches the build.
            receptor_row["gpcrdb_text_sha256"] = cm.text_sha256(gtext)
            receptor_row["product_instances_sha256"] = cm.instances_sha256(instances)
            receptor_rows.append(receptor_row)
            seen = set()
            for sli in anchors[pdb]:
                het = sli.pdb_reference.upper()
                tokens = cm.split_tokens(sli.chain_res)
                if not tokens:
                    copies = sorted(n for n in instances if n.split("_", 1)[0].upper() == het)
                    key = (het, "")
                    if key not in seen:
                        seen.add(key)
                        anchor_rows.append(dict({c: "" for c in cm.ANCHOR_COLUMNS}, pdb=pdb, het=het,
                                                token="", instance=";".join(copies),
                                                status="all_copies" if copies else "no_product",
                                                note=("chain_res {!r} names no residue; every copy used".format(
                                                    sli.chain_res or "") if copies else
                                                    "the product has no instance of {} (chain_res {!r})".format(
                                                        het, sli.chain_res or ""))))
                    continue
                for tok in tokens:
                    if (het, tok) in seen:
                        continue
                    seen.add((het, tok))
                    anchor_rows.append(cm.resolve_anchor(
                        pdb, het, tok, cif_atoms, gatoms, instances, labels.get((pdb, het, tok))))

        header = [
            "dump_id\t" + opt["dump_id"],
            "annotation_commit\t" + opt["annotation_commit"],
            "annotation_sha256\t" + _sha256_file(opt["annotation"]),
            "cif_manifest_sha256\t" + _sha256_lines(cif_manifest),
            "product_manifest_sha256\t" + _sha256_lines(product_manifest),
            "builder_sha256\t" + _sha256_lines([_sha256_file(cm.__file__), _sha256_file(__file__)]),
            "structures\t" + str(len(anchors)),
        ]
        for name, columns, out in (("anchor_instance_map.tsv", cm.ANCHOR_COLUMNS, anchor_rows),
                                   ("receptor_chain_map.tsv", cm.RECEPTOR_COLUMNS, receptor_rows)):
            with open(os.path.join(opt["out_dir"], name), "w", newline="") as fh:
                for line in header:
                    fh.write("# " + line + "\n")
                w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
                w.writeheader()
                w.writerows(out)
        counts = {}
        for r in anchor_rows:
            counts[(r["status"], r["source"])] = counts.get((r["status"], r["source"]), 0) + 1
        rcounts = {}
        for r in receptor_rows:
            rcounts[(r["status"], r["method"])] = rcounts.get((r["status"], r["method"]), 0) + 1
        self.stdout.write("anchor rows {}: {}".format(len(anchor_rows), sorted(counts.items())))
        self.stdout.write("receptor rows {}: {}".format(len(receptor_rows), sorted(rcounts.items())))
