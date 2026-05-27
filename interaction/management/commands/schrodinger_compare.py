"""Characterization command for the Schrödinger small-molecule pipeline (Plan A, Step A10).

Runs the *new* Schrödinger interaction processor against real GPCRdb data for one
PDB+ligand, captures the resulting ``ResidueFragmentInteraction`` rows, and diffs
them against a golden baseline (the *old* RDKit pipeline's rows, exported in A0).

By default it is a **dry run**: the processor writes inside a transaction that is
rolled back, so the production database is left untouched — only the diff is
reported. Pass ``--commit`` to persist (e.g. for a real build).

    docker exec gpcrdb-app python manage.py schrodinger_compare
    docker exec gpcrdb-app python manage.py schrodinger_compare --pdb 6LN2 --het 97Y

The diff is the project's first real scientific output: it shows where Schrödinger
(Engine 1) and the legacy RDKit pipeline agree and differ. Differences are expected
and meaningful, not bugs (see Plan A Inbox-3): ``acc`` (accessibility) and
``Van der Waals`` are structural — Engine 1 does not emit them.
"""

import json
import os

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from interaction.models import ResidueFragmentInteraction, StructureLigandInteraction
from interaction.schrodinger_processor import process_schrodinger_sm_interactions

# Slugs Engine 1 cannot emit by construction (Plan A Inbox-3); excluded from the
# "comparable chemistry" overlap metric but still shown in the diff.
STRUCTURAL_ONLY_SLUGS = {"acc", "Van der Waals"}

DEFAULT_YAML_DIR = os.path.join(
    settings.BASE_DIR, "interaction", "tests", "schrodinger_data"
)
DEFAULT_GOLDEN = os.path.join(
    settings.BASE_DIR, "interaction", "tests", "golden_6LN2.json"
)


class Command(BaseCommand):
    help = "Diff the Schrödinger SM pipeline output against the RDKit golden baseline (A10)."

    def add_arguments(self, parser):
        parser.add_argument("--pdb", default="6LN2")
        parser.add_argument("--het", default="97Y")
        parser.add_argument(
            "--yaml-dir",
            default=DEFAULT_YAML_DIR,
            help="SCHRODINGER_INTERACTIONS_DIR override (the Engine 1 nested layout root).",
        )
        parser.add_argument(
            "--golden",
            default=DEFAULT_GOLDEN,
            help="golden_<PDB>.json (old RDKit baseline) to diff against.",
        )
        parser.add_argument(
            "--commit",
            action="store_true",
            default=False,
            help="Persist the new rows (default: dry-run that rolls back).",
        )

    # -- helpers ---------------------------------------------------------------

    def _load_golden_set(self, path):
        with open(path) as fh:
            data = json.load(fh)
        # golden_set is a list of [seq, slug]; normalise to a set of tuples.
        return {(int(seq), slug) for seq, slug in data["golden_set"]}

    def _new_set_for_pair(self, sli):
        qs = ResidueFragmentInteraction.objects.filter(
            structure_ligand_pair=sli
        ).values_list("rotamer__residue__sequence_number", "interaction_type__slug")
        return {(int(seq), slug) for seq, slug in qs}

    def _print_block(self, title, pairs):
        self.stdout.write(self.style.MIGRATE_HEADING(title))
        if not pairs:
            self.stdout.write("    (none)")
            return
        for seq, slug in sorted(pairs):
            self.stdout.write(f"    {seq:>5}  {slug}")

    # -- main ------------------------------------------------------------------

    def handle(self, *args, **opts):
        pdb = opts["pdb"].upper()
        het = opts["het"].upper()

        try:
            sli = StructureLigandInteraction.objects.get(
                structure__pdb_code__index=pdb, pdb_reference=het
            )
        except StructureLigandInteraction.DoesNotExist:
            raise CommandError(
                f"No StructureLigandInteraction for {pdb} / {het}. "
                "The structure+ligand must already exist in the DB."
            )
        except StructureLigandInteraction.MultipleObjectsReturned:
            raise CommandError(f"Multiple SLIs for {pdb} / {het}; ambiguous.")

        golden = self._load_golden_set(opts["golden"])

        # Run the new pipeline. Dry-run by default: capture inside the transaction,
        # then roll back so production data is untouched.
        new_set = set()
        ok = False
        with transaction.atomic():
            ok = process_schrodinger_sm_interactions(
                current_structure_obj=sli.structure,
                current_ligand_db_obj=sli.ligand,
                ligand_pdb_het_code=het,
                pdb_code_str=pdb,
                schrodinger_interactions_dir_override=opts["yaml_dir"],
            )
            new_set = self._new_set_for_pair(sli)
            if not opts["commit"]:
                transaction.set_rollback(True)

        # --- report -----------------------------------------------------------
        both = golden & new_set
        only_old = golden - new_set
        only_new = new_set - golden

        # Comparable-chemistry view: drop slugs Engine 1 cannot emit from the old side.
        golden_comparable = {(s, g) for (s, g) in golden if g not in STRUCTURAL_ONLY_SLUGS}
        comp_both = golden_comparable & new_set
        comp_only_old = golden_comparable - new_set

        mode = "COMMITTED" if opts["commit"] else "DRY-RUN (rolled back)"
        self.stdout.write("")
        self.stdout.write(self.style.MIGRATE_HEADING(
            f"=== Schrödinger vs RDKit diff for {pdb} / {het}  [{mode}] ==="
        ))
        self.stdout.write(
            f"processor returned: {ok}   ligand: {sli.ligand.name}   "
            f"preferred_chain: {sli.structure.preferred_chain}"
        )
        self.stdout.write(
            f"golden (old RDKit): {len(golden)} rows   "
            f"new (Schrödinger): {len(new_set)} rows"
        )
        self.stdout.write("")
        self._print_block("BOTH (identical seq+slug):", both)
        self._print_block("ONLY OLD (RDKit; missing from Schrödinger):", only_old)
        self._print_block("ONLY NEW (Schrödinger; not in RDKit):", only_new)
        self.stdout.write("")
        self.stdout.write(self.style.MIGRATE_HEADING(
            "Comparable-chemistry overlap (golden minus acc/Van der Waals):"
        ))
        denom = len(golden_comparable) or 1
        self.stdout.write(
            f"    {len(comp_both)}/{len(golden_comparable)} "
            f"({100*len(comp_both)//denom}%) of comparable old rows reproduced; "
            f"comparable-only-old: {sorted(comp_only_old)}"
        )

        # Emit machine-readable JSON for downstream report authoring.
        self.stdout.write("")
        self.stdout.write("<<<DIFF_JSON>>>")
        self.stdout.write(json.dumps({
            "pdb": pdb, "het": het, "mode": mode, "processor_ok": ok,
            "golden_count": len(golden), "new_count": len(new_set),
            "both": sorted(both), "only_old": sorted(only_old),
            "only_new": sorted(only_new),
            "comparable_both": sorted(comp_both),
            "comparable_only_old": sorted(comp_only_old),
        }))
        self.stdout.write("<<<END_DIFF_JSON>>>")
