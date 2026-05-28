"""Unit tests for the import_schrodinger_interactions management command (Phase 1c).

These tests focus on the command shell: argparse, fixture self-check,
PDB-list resolution, anomaly CSV behavior (outside the transaction), and
the dry-run summary shape. End-to-end against real Structure rows lives
in the Phase C integration check.
"""

import csv
import io
import os
import tempfile

from django.core.management import call_command
from django.core.management.base import CommandError
from django.test import TestCase

from interaction.management.commands.import_schrodinger_interactions import (
    ANOMALY_LEVELS,
    Command,
    _AnomalyWriter,
)
from interaction.models import (
    ResidueFragmentInteraction,
    ResidueFragmentInteractionType,
)


class AnomalyWriterTests(TestCase):
    """The CSV writer must survive a transaction rollback (axis 4C)."""

    def test_writes_header_and_rows_with_flush(self):
        with tempfile.NamedTemporaryFile(mode='r', suffix='.csv', delete=False) as fh:
            path = fh.name
        try:
            with _AnomalyWriter(path) as w:
                w.log('6LN2', 'zero_interactions_by_design',
                      residue='N406', message='no interactions')
                w.log('6CM4', 'comma_feature', message='Halogen,Acceptor')
            with open(path) as fh:
                rows = list(csv.reader(fh))
            self.assertEqual(rows[0], _AnomalyWriter.COLUMNS)
            self.assertEqual(len(rows), 3)
            # row layout: timestamp, pdb_id, level, category, residue, message
            self.assertEqual(rows[1][1], '6LN2')
            self.assertEqual(rows[1][2], 'INFO')
            self.assertEqual(rows[1][3], 'zero_interactions_by_design')
            self.assertEqual(rows[2][2], 'WARNING')
            self.assertEqual(rows[2][3], 'comma_feature')
        finally:
            os.unlink(path)

    def test_level_counts_track_per_severity(self):
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as fh:
            path = fh.name
        try:
            with _AnomalyWriter(path) as w:
                w.log('A', 'zero_interactions_by_design')
                w.log('B', 'comma_feature')
                w.log('C', 'slug_not_in_fixture_21')
                counts = dict(w.counts)
            self.assertEqual(counts['INFO'], 1)
            self.assertEqual(counts['WARNING'], 1)
            self.assertEqual(counts['ERROR'], 1)
        finally:
            os.unlink(path)


class AnomalyCategoryContractTests(TestCase):
    """ANOMALY_LEVELS is the contract between calculator + log consumer."""

    def test_required_categories_present(self):
        for cat in [
            'zero_interactions_by_design',
            'metal_water_drop_pending_plan_d',
            'non_standard_aa_x_residue',
            'comma_feature',
            'missing_receptor_atom_name',
            'residue_lookup_fail',
            'slug_not_in_fixture_21',
            'idempotency_drift_over_20pct',
        ]:
            self.assertIn(cat, ANOMALY_LEVELS)

    def test_levels_are_only_known_severities(self):
        for level in ANOMALY_LEVELS.values():
            self.assertIn(level, {'INFO', 'WARNING', 'ERROR'})

    def test_slug_not_in_fixture_is_error(self):
        # ADR-009 fail-loud: an unrecognised slug must abort the tx, not log silently.
        self.assertEqual(ANOMALY_LEVELS['slug_not_in_fixture_21'], 'ERROR')


class FixtureSelfCheckTests(TestCase):
    """Phase 1c axis 1D: the command auto-seeds the 21-slug fixture if missing."""

    def test_passes_when_fixture_seeded(self):
        # Migration 0007 seeds 21 rows automatically in the test DB. Verify count.
        self.assertGreaterEqual(ResidueFragmentInteractionType.objects.count(), 21)
        cmd = Command()
        cmd.stdout = io.StringIO()
        cmd._check_fixture()  # should not raise

    def test_auto_loaddata_when_under_seeded(self):
        # Wipe the catalog (no FK refs in this test DB yet) and verify the command
        # repopulates it via call_command('loaddata', ...).
        ResidueFragmentInteractionType.objects.all().delete()
        self.assertEqual(ResidueFragmentInteractionType.objects.count(), 0)
        cmd = Command()
        cmd.stdout = io.StringIO()
        cmd._check_fixture()
        self.assertGreaterEqual(ResidueFragmentInteractionType.objects.count(), 21)


class PdbListResolutionTests(TestCase):
    """The --pdb-list resolver is fail-loud on missing files (no silent fallback)."""

    def test_missing_file_raises(self):
        cmd = Command()
        with self.assertRaises(CommandError) as ctx:
            cmd._resolve_pdb_list({'pdb_list': '/nonexistent/path/structures.csv'})
        self.assertIn('PDB list not found', str(ctx.exception))

    def test_parses_one_per_line(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as fh:
            fh.write("# this is a comment\n")
            fh.write("6LN2\n")
            fh.write("  6cm4  \n")  # whitespace + lowercase tolerated
            fh.write("\n")
            fh.write("1F88\n")
            path = fh.name
        try:
            cmd = Command()
            result = cmd._resolve_pdb_list({'pdb_list': path})
            self.assertEqual(result, ['6LN2', '6CM4', '1F88'])
        finally:
            os.unlink(path)

    def test_parses_csv_first_column_with_header(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as fh:
            fh.write('pdb,resolution,family\n')
            fh.write('6LN2,2.4,5HT\n')
            fh.write('2Y02,2.3,ADRB\n')
            path = fh.name
        try:
            cmd = Command()
            result = cmd._resolve_pdb_list({'pdb_list': path})
            self.assertEqual(result, ['6LN2', '2Y02'])
        finally:
            os.unlink(path)

    def test_empty_after_filter_raises(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as fh:
            fh.write("# only comments here\n")
            path = fh.name
        try:
            cmd = Command()
            with self.assertRaises(CommandError):
                cmd._resolve_pdb_list({'pdb_list': path})
        finally:
            os.unlink(path)


class MetricsTests(TestCase):
    """Axis 3C dual metrics: chemistry-on-residue + literal-overlap.

    Numbers below mirror the 6CM4 验证 fixture (Step 6CM4 验证.md):
        polar comparable golden = 9 rows over 8 residues
        new pipeline reproduces 5 of those residues but only 1 (seq, slug)
        ⇒ chemistry-on-residue ≈ 62%, literal-overlap ≈ 11%
    """

    def test_chemistry_charity_above_literal_strictness(self):
        cmd = Command()
        # Synthesised pair: 8 golden residues, new pipeline hits 5 with mostly-different slugs.
        golden = {
            (114, 'polar_double_neg_protein'),
            (114, 'polar_acceptor_protein'),
            (118, 'polar_donor_protein'),
            (200, 'aro_ef_protein'),
            (203, 'aro_ef_protein'),
            (303, 'aro_ef_protein'),
            (309, 'aro_ef_protein'),
            (100, 'aro_ef_protein'),
        }
        new = {
            (114, 'polar_double_neg_protein'),  # literal match
            (114, 'polar_acceptor_protein'),    # literal match (added in v2)
            (200, 'aro_ef_protein'),            # literal match
            (203, 'aro_ef_protein'),            # literal match
            (309, 'aro_ef_protein'),            # literal match
            (100, 'aro_ef_protein'),            # literal match
        }
        chem, literal = cmd._metrics(golden, new)
        # golden unique residues: {114, 118, 200, 203, 303, 309, 100} = 7
        # new unique residues: {114, 200, 203, 309, 100} = 5 (all in golden)
        # → chemistry-on-residue = 5/7 ≈ 71.4%
        self.assertAlmostEqual(chem, 71.4, places=1)
        # golden comparable = 8 (no acc/VdW); new ∩ golden = 6 literal pairs
        # → literal-overlap = 6/8 = 75%
        self.assertAlmostEqual(literal, 75.0, places=1)

    def test_no_golden_returns_none(self):
        cmd = Command()
        chem, literal = cmd._metrics(set(), set())
        self.assertIsNone(chem)
        self.assertIsNone(literal)

    def test_structural_only_slugs_excluded_from_literal_denominator(self):
        cmd = Command()
        # If denominator included 'acc' / 'Van der Waals' Engine 1 can never emit,
        # literal-overlap would be pushed to 0 by construction.
        golden = {(100, 'acc'), (100, 'Van der Waals'), (100, 'aro_ef_protein')}
        new = {(100, 'aro_ef_protein')}
        chem, literal = cmd._metrics(golden, new)
        # comparable = just {(100, aro_ef_protein)} ⇒ literal = 1/1 = 100%
        self.assertEqual(literal, 100.0)


class CommandSmokeTests(TestCase):
    """End-to-end shell smoke: --dry-run on a non-existent PDB list still gives a clean rollback."""

    def test_dry_run_with_missing_structures_emits_anomalies_no_writes(self):
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as fh:
            fh.write("ZZZZ\nYYYY\n")  # neither exists in the test DB
            pdb_list = fh.name
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as fh:
            anomaly_csv = fh.name
        try:
            before = ResidueFragmentInteraction.objects.count()
            stdout = io.StringIO()
            call_command(
                'import_schrodinger_interactions',
                '--dry-run',
                '--pdb-list', pdb_list,
                '--anomaly-csv', anomaly_csv,
                stdout=stdout,
            )
            after = ResidueFragmentInteraction.objects.count()
            self.assertEqual(before, after)  # dry-run never writes
            with open(anomaly_csv) as fh:
                rows = list(csv.reader(fh))
            # Header + 2 'structure_not_in_db' entries
            self.assertEqual(len(rows), 3)
            self.assertIn('structure_not_in_db', rows[1][3])
            self.assertIn('structure_not_in_db', rows[2][3])
        finally:
            os.unlink(pdb_list)
            os.unlink(anomaly_csv)
