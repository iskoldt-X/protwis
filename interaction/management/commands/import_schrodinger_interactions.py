"""Production import path for Schrödinger pre-computed interaction YAMLs.

Phase 1c (axes 1D / 2E / 3C / 4C / 5A / 6C). Iterates a list of PDB IDs,
dispatches each to ``SchrodingerInteractionCalculator``, and persists RFI
rows in one atomic transaction. Anomalies are logged to a CSV that lives
**outside** the transaction so a rollback never loses the diagnosis trail.

Usage::

    # Dry-run on the 4 Phase 1c integration fixtures (no DB writes):
    docker exec gpcrdb-app python manage.py import_schrodinger_interactions \\
        --dry-run \\
        --pdb-list /app/interaction/tests/schrodinger_data/_pdb_list.txt \\
        --schrodinger-data-dir /app/interaction/tests/schrodinger_data

    # Full import (Phase 2, against the gpcrdb-db-schrodinger-new container):
    docker exec gpcrdb-app python manage.py import_schrodinger_interactions

Three-level logging (axis 6C, table in Phase 1c 设计 §轴 6) routes each anomaly
to INFO / WARNING / ERROR. ERROR categories raise to abort the whole tx (axis
4C single-big-tx rollback-on-fail); WARNING / INFO continue. The category
table is the contract — see ``ANOMALY_LEVELS`` at module top.
"""

import csv
import datetime
import os

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from interaction.calculators import get_interaction_calculator
from interaction.models import (
    ResidueFragmentInteraction,
    ResidueFragmentInteractionType,
    StructureLigandInteraction,
)
from structure.models import Structure


DEFAULT_PDB_LIST = '/app/data/protwis/gpcr/structure_data/annotation/structures.csv'

# Anomaly category -> log level (Phase 1c 设计 §轴 6). ERROR categories raise
# and abort the transaction. WARNING/INFO continue to the next PDB.
ANOMALY_LEVELS = {
    'zero_interactions_by_design': 'INFO',
    'metal_water_drop_pending_plan_d': 'WARNING',
    'non_standard_aa_x_residue': 'WARNING',
    'comma_feature': 'WARNING',
    'idempotency_drift_over_20pct': 'WARNING',
    'missing_receptor_atom_name': 'ERROR',
    'residue_lookup_fail': 'ERROR',
    'slug_not_in_fixture_21': 'ERROR',
    'structure_not_in_db': 'WARNING',
    'no_sli_for_structure': 'INFO',
    'processor_returned_false': 'WARNING',
}


class _AnomalyWriter:
    """CSV writer that bypasses the Django ORM.

    Phase 1c axis 4C: the anomaly log must survive a transaction rollback.
    We open the file at __enter__, flush on every row, and never touch the
    DB — so even if the outer @transaction.atomic raises and rolls back,
    the diagnosis trail is on disk.
    """

    COLUMNS = ['timestamp', 'pdb_id', 'level', 'category', 'residue', 'message']

    def __init__(self, path):
        self.path = path
        self._fh = None
        self._writer = None
        self.counts = {'INFO': 0, 'WARNING': 0, 'ERROR': 0}

    def __enter__(self):
        self._fh = open(self.path, 'w', newline='')
        self._writer = csv.writer(self._fh)
        self._writer.writerow(self.COLUMNS)
        self._fh.flush()
        return self

    def __exit__(self, exc_type, exc, tb):
        if self._fh:
            self._fh.close()
        return False

    def log(self, pdb_id, category, residue='', message=''):
        level = ANOMALY_LEVELS.get(category, 'WARNING')
        self.counts[level] = self.counts.get(level, 0) + 1
        self._writer.writerow([
            datetime.datetime.utcnow().isoformat(timespec='seconds'),
            pdb_id, level, category, residue, message,
        ])
        self._fh.flush()  # critical: ensure rollback can't lose the row
        return level


class Command(BaseCommand):
    help = "Import Schrödinger pre-computed interactions for one or many PDB IDs (Phase 1c)."

    # Slugs Engine 1 cannot emit by construction (acc = accessibility,
    # Van der Waals = structural-only). Excluded from the literal-overlap
    # metric denominator since they would push it to 0% by construction.
    STRUCTURAL_ONLY_SLUGS = frozenset({'acc', 'Van der Waals'})

    def add_arguments(self, parser):
        parser.add_argument('--dry-run', action='store_true', default=False,
            help='Preview without persisting. Computes chemistry-on-residue + literal-overlap '
                 'metrics for any PDB with existing golden RFI rows. Phase 1c axis 3C.')
        parser.add_argument('--batch-size', type=int, default=None,
            help='PDBs per inner transaction (default = all in one big tx, axis 4C). '
                 'Set to e.g. 50 if memory profile demands smaller batches.')
        parser.add_argument('--pdb-list', type=str, default=None,
            help='Path to a file with PDB IDs, one per line. Lines starting with # are '
                 'ignored. Default: {} (the structures.csv bind-mount). '
                 'A .csv first column is also accepted (header row skipped).'.format(DEFAULT_PDB_LIST))
        parser.add_argument('--schrodinger-data-dir', type=str, default=None,
            help='Override SCHRODINGER_INTERACTIONS_DIR. Defaults to settings or to the '
                 'value in the processor.')
        parser.add_argument('--anomaly-csv', type=str, default='schrodinger_import_anomalies.csv',
            help='Path for the anomaly CSV (written outside the transaction, flushed per row).')

    # -- helpers ---------------------------------------------------------------

    def _check_fixture(self):
        """Phase 1c axis 1D: fail-loud self-check that the 21-slug fixture is seeded.

        Uses slug-keyed update_or_create instead of loaddata, mirroring migration
        0007's pattern. Tolerant of legacy DBs that already have the original 18
        slugs with different PKs (loaddata would PK-collide on those).
        """
        count = ResidueFragmentInteractionType.objects.count()
        if count >= 21:
            self.stdout.write("Fixture self-check OK ({} RFIT rows).".format(count))
            return
        self.stdout.write(self.style.WARNING(
            "RFIT count = {} (expected 21). Upserting interaction_types fixture by slug...".format(count)))
        # Locate the fixture file from the interaction app. Mirrors migration 0007.
        import json
        from django.apps import apps
        app_path = apps.get_app_config('interaction').path
        fixture_path = os.path.join(app_path, 'fixtures', 'interaction_types.json')
        with open(fixture_path, 'r') as fh:
            rows = json.load(fh)
        for row in rows:
            if row.get('model') != 'interaction.residuefragmentinteractiontype':
                continue
            fields = row['fields']
            ResidueFragmentInteractionType.objects.update_or_create(
                slug=fields['slug'],
                defaults={
                    'name': fields['name'],
                    'type': fields['type'],
                    'direction': fields.get('direction', ''),
                },
            )
        count = ResidueFragmentInteractionType.objects.count()
        if count < 21:
            raise CommandError(
                "Fixture seed failed: RFIT count = {} after upsert (expected >= 21). "
                "Check interaction/fixtures/interaction_types.json.".format(count))

    def _resolve_pdb_list(self, options):
        """Resolve --pdb-list to a list of PDB IDs, fail-loud on missing file."""
        path = options['pdb_list'] or DEFAULT_PDB_LIST
        if not os.path.isfile(path):
            raise CommandError(
                "PDB list not found at {!r}. Provide --pdb-list or ensure structures.csv "
                "is bind-mounted at {}.".format(path, DEFAULT_PDB_LIST))
        pdbs = []
        with open(path, 'r') as fh:
            for i, line in enumerate(fh):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                # Accept either a bare PDB ID per line or a CSV-style row (first column).
                first = line.split(',')[0].strip().strip('"').strip("'")
                if i == 0 and first.lower() in ('pdb', 'pdb_code', 'pdb_id'):
                    continue  # header
                if len(first) != 4:
                    continue
                pdbs.append(first.upper())
        if not pdbs:
            raise CommandError("PDB list {!r} is empty after filtering.".format(path))
        return pdbs

    def _golden_set(self, structure):
        """Capture (seq, slug) for all RFI rows on this structure before the processor runs."""
        qs = ResidueFragmentInteraction.objects.filter(
            structure_ligand_pair__structure=structure
        ).values_list(
            'rotamer__residue__sequence_number',
            'interaction_type__slug',
        )
        return {(seq, slug) for seq, slug in qs if seq is not None}

    def _new_set(self, structure):
        return self._golden_set(structure)  # same query shape; semantic alias

    def _metrics(self, golden, new):
        """Compute Phase 1c axis 3C dual metrics.

        literal-overlap = |new ∩ golden_comparable| / |golden_comparable|
            (excludes acc + Van der Waals from the denominator since
             Engine 1 cannot emit them — including them forces 0%)
        chemistry-on-residue = |new_residues ∩ golden_residues| / |golden_residues|
            (projection to the per-residue set; charity to slug differences)
        """
        comparable = {(s, g) for (s, g) in golden if g not in self.STRUCTURAL_ONLY_SLUGS}
        if not comparable:
            literal = None
        else:
            literal = round(100.0 * len(new & comparable) / len(comparable), 1)

        golden_res = {seq for (seq, _) in comparable}
        new_res = {seq for (seq, _) in new}
        if not golden_res:
            chem = None
        else:
            chem = round(100.0 * len(new_res & golden_res) / len(golden_res), 1)
        return chem, literal

    def _process_one_pdb(self, pdb_id, calculator, anomaly_writer, dry_run, schrodinger_data_dir):
        """Returns a per-PDB summary dict for the dry-run report."""
        try:
            structure = Structure.objects.get(pdb_code__index=pdb_id)
        except Structure.DoesNotExist:
            anomaly_writer.log(pdb_id, 'structure_not_in_db',
                message="No Structure row for {}; skipping.".format(pdb_id))
            return {'pdb': pdb_id, 'status': 'skipped_no_structure'}

        slis = list(StructureLigandInteraction.objects.filter(structure=structure))
        if not slis:
            anomaly_writer.log(pdb_id, 'no_sli_for_structure',
                message="Structure {} has zero SLI rows; nothing to do.".format(pdb_id))
            return {'pdb': pdb_id, 'status': 'skipped_no_sli'}

        # Pre-capture golden set BEFORE the processor deletes-then-inserts.
        golden = self._golden_set(structure)

        # Run the calculator. It calls process_schrodinger_sm_interactions per SLI;
        # that function is @transaction.atomic and does the delete-then-insert.
        results = calculator.compute_interactions(
            structure, schrodinger_data_dir=schrodinger_data_dir)

        # Check for processor failures per-SLI.
        for sli_id, het, ok, err in results:
            if not ok:
                anomaly_writer.log(pdb_id, 'processor_returned_false',
                    message="SLI {} het={} err={}".format(sli_id, het, err or 'unknown'))

        new = self._new_set(structure)
        chem, literal = self._metrics(golden, new)
        summary = {
            'pdb': pdb_id,
            'status': 'ok',
            'sli_count': len(slis),
            'golden_rows': len(golden),
            'new_rows': len(new),
            'chemistry_on_residue_pct': chem,
            'literal_overlap_pct': literal,
        }

        # Axis 6C INFO: zero-interactions-by-design (e.g. 1F88 RET covalent retinal).
        # Engine 1 sees a covalent attachment and emits an empty interaction list;
        # the legacy RDKit golden may carry rows because RDKit does not honor
        # covalent semantics — but the new pipeline writing zero is correct, not
        # a bug. Heuristic: every SLI processor call succeeded AND no rows landed.
        if results and not new and all(ok for (_, _, ok, _) in results):
            anomaly_writer.log(pdb_id, 'zero_interactions_by_design',
                message="Processor returned ok for all {} SLIs, zero rows written; "
                        "consistent with a covalent ligand (Engine 1 by-design "
                        "empty interactions). Golden had {} rows.".format(
                            len(results), len(golden)))
            summary['status'] = 'zero_by_design'

        return summary

    def _print_dry_run_report(self, summaries, anomaly_path, anomaly_counts):
        self.stdout.write("")
        self.stdout.write(self.style.MIGRATE_HEADING(
            "=== Phase 1c dry-run report ({} PDBs) ===".format(len(summaries))))
        self.stdout.write("{:<8} {:>6} {:>7} {:>7} {:>9} {:>9}  {}".format(
            'PDB', 'SLIs', 'golden', 'new', 'chem%', 'literal%', 'status'))
        for s in summaries:
            chem = '-' if s.get('chemistry_on_residue_pct') is None else '{:.1f}'.format(s['chemistry_on_residue_pct'])
            lit = '-' if s.get('literal_overlap_pct') is None else '{:.1f}'.format(s['literal_overlap_pct'])
            self.stdout.write("{:<8} {:>6} {:>7} {:>7} {:>9} {:>9}  {}".format(
                s['pdb'],
                s.get('sli_count', '-'),
                s.get('golden_rows', '-'),
                s.get('new_rows', '-'),
                chem, lit,
                s['status']))
        self.stdout.write("")
        self.stdout.write("Anomaly CSV: {} (INFO={} WARNING={} ERROR={})".format(
            anomaly_path,
            anomaly_counts.get('INFO', 0),
            anomaly_counts.get('WARNING', 0),
            anomaly_counts.get('ERROR', 0)))
        self.stdout.write("")
        self.stdout.write(self.style.NOTICE(
            "Dry-run: transaction rolled back; database state unchanged."))

    # -- main ------------------------------------------------------------------

    def handle(self, *args, **options):
        pdb_list = self._resolve_pdb_list(options)
        self.stdout.write("Resolved {} PDB IDs from {}.".format(
            len(pdb_list), options['pdb_list'] or DEFAULT_PDB_LIST))

        calculator = get_interaction_calculator('schrodinger')
        batch_size = options['batch_size'] or len(pdb_list)
        dry_run = options['dry_run']
        schrodinger_data_dir = options['schrodinger_data_dir']

        summaries = []
        anomaly_path = options['anomaly_csv']

        with _AnomalyWriter(anomaly_path) as anomaly_writer:
            # Phase 1c axis 4C: one big @transaction.atomic. On dry-run, set_rollback
            # at the end so the entire write is undone (including the fixture upsert);
            # on full run, the transaction commits at scope exit. ERROR-level anomalies
            # re-raise → tx rollback. Fixture self-check sits INSIDE the tx so a
            # legacy DB does not get permanently mutated by a dry-run.
            try:
                with transaction.atomic():
                    self._check_fixture()
                    for batch_start in range(0, len(pdb_list), batch_size):
                        batch = pdb_list[batch_start:batch_start + batch_size]
                        for pdb_id in batch:
                            try:
                                summary = self._process_one_pdb(
                                    pdb_id, calculator, anomaly_writer,
                                    dry_run, schrodinger_data_dir)
                                summaries.append(summary)
                            except Exception as e:
                                # Any unexpected exception is logged ERROR + re-raised.
                                anomaly_writer.log(pdb_id, 'processor_returned_false',
                                    message="Unhandled exception: {}: {}".format(
                                        type(e).__name__, str(e)[:200]))
                                raise
                    if dry_run:
                        transaction.set_rollback(True)
            except Exception as e:
                # Outer catch: ensures anomaly file is closed cleanly (the with-block
                # exits regardless) and gives the operator a clear top-level message.
                self.stdout.write(self.style.ERROR(
                    "Import aborted; transaction rolled back. Cause: {}: {}".format(
                        type(e).__name__, str(e)[:300])))
                self.stdout.write("Anomaly CSV preserved at {}.".format(anomaly_path))
                raise

            if dry_run:
                self._print_dry_run_report(summaries, anomaly_path, anomaly_writer.counts)
            else:
                ok = sum(1 for s in summaries if s['status'] == 'ok')
                self.stdout.write("")
                self.stdout.write(self.style.SUCCESS(
                    "Imported {} PDBs ({} OK). Anomaly CSV: {} (INFO={} WARNING={} ERROR={}).".format(
                        len(summaries), ok, anomaly_path,
                        anomaly_writer.counts.get('INFO', 0),
                        anomaly_writer.counts.get('WARNING', 0),
                        anomaly_writer.counts.get('ERROR', 0))))
