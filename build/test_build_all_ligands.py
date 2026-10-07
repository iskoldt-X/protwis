"""Unit tests for the ligand import steps of build_all and build_all_interactions
(no database: every import command is replaced by a stub).

    python -c "import django; django.setup(); import unittest; \\
        unittest.main(module='build.test_build_all_ligands', argv=['x'])"
"""

import os
import tempfile
import unittest
from unittest import mock

from django.core.management.base import CommandError

from build import ligand_imports
from build.management.commands import build_all


def options(**kw):
    out = {"skip_ligand_import": False, "engine1_data_dir": "/e1", "engine2_data_dir": "/e2",
           "engine1_report_dir": "/r1", "engine2_report_dir": "/r2", "phase": None}
    out.update(kw)
    return out


def no_import(name, **kw):
    raise AssertionError("a real import command was reached: %s" % name)


class LigandImportStepsTests(unittest.TestCase):
    def test_maps_then_both_dry_runs_then_both_imports(self):
        steps = build_all.Command().ligand_import_steps(options())
        self.assertEqual([(c, o.get("dry_run", False), o.get("data_dir")) for c, o in steps], [
            ("build_schrodinger_chainmap_files", False, "/e1"),
            ("build_schrodinger_peptide_maps", False, "/e2"),
            ("import_schrodinger_interactions", True, "/e1"),
            ("import_schrodinger_peptides", True, "/e2"),
            ("import_schrodinger_interactions", False, "/e1"),
            ("import_schrodinger_peptides", False, "/e2"),
            ("remove_schrodinger_maps", False, None),
        ])
        self.assertEqual([o["report_json"] for c, o in steps if c in ligand_imports.COMMANDS],
                         ["/r1/report.dryrun.json", "/r2/report.dryrun.json",
                          "/r1/report.json", "/r2/report.json"])

    def test_the_peptide_maps_read_the_indexes_from_the_engine1_tree(self):
        steps = dict((c, o) for c, o in build_all.Command().ligand_import_steps(options())
                     if c in ligand_imports.MAP_COMMANDS)
        self.assertEqual(steps["build_schrodinger_peptide_maps"],
                         {"data_dir": "/e2", "index_dir": "/e1"})
        # The annotation decides which structures get a map; a product directory
        # it no longer lists is not a reason to stop the build.
        self.assertEqual(steps["build_schrodinger_chainmap_files"],
                         {"data_dir": "/e1", "allow_stray": True})

    def test_split_runs_the_maps_with_the_dry_runs(self):
        before, imports = ligand_imports.split(ligand_imports.steps(options()))
        self.assertEqual([(c, o.get("dry_run", False)) for c, o in before], [
            ("build_schrodinger_chainmap_files", False),
            ("build_schrodinger_peptide_maps", False),
            ("import_schrodinger_interactions", True),
            ("import_schrodinger_peptides", True)])
        self.assertEqual([(c, o.get("dry_run", False)) for c, o in imports], [
            ("import_schrodinger_interactions", False),
            ("import_schrodinger_peptides", False),
            ("remove_schrodinger_maps", False)])

    def test_the_clean_up_names_both_trees(self):
        steps = build_all.Command().ligand_import_steps(options())
        self.assertEqual(steps[-1], ["remove_schrodinger_maps",
                                     {"engine1_dir": "/e1", "engine2_dir": "/e2"}])

    def test_the_peptide_maps_alone_need_the_engine1_delivery(self):
        with tempfile.TemporaryDirectory() as d2:
            with self.assertRaisesRegex(CommandError, "Engine 1 products are missing"):
                ligand_imports.check_deliveries(
                    options(engine1_data_dir="/nonexistent-e1", engine2_data_dir=d2),
                    ["build_schrodinger_peptide_maps"])
            ligand_imports.check_deliveries(
                options(engine1_data_dir="/nonexistent-e1", engine2_data_dir=d2),
                ["import_schrodinger_peptides"])

    def test_skipping_imports_nothing(self):
        self.assertEqual(build_all.Command().ligand_import_steps(options(skip_ligand_import=True)), [])


class BuildAllInteractionsTests(unittest.TestCase):
    """build_all_interactions keeps its job: contact network and the same imports,
    dry-run before the contacts and imported after them."""

    def test_it_takes_the_import_options_and_dry_runs_before_the_contacts(self):
        from tools.management.commands import build_all_interactions as bai
        cmd = bai.Command()
        parser = cmd.create_parser("manage.py", "build_all_interactions")
        opts = vars(parser.parse_args([]))
        for key in ("engine1_data_dir", "engine2_data_dir", "skip_ligand_import"):
            self.assertIn(key, opts)
        calls = []
        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            opts.update(options(engine1_data_dir=d1, engine2_data_dir=d2))
            with mock.patch.object(bai.Command, "prepare_input",
                                   lambda self, proc, pdbs: calls.append("contacts")), \
                    mock.patch.object(ligand_imports, "call_command",
                                      lambda name, **kw: calls.append((name, kw.get("dry_run", False)))):
                cmd.handle(**opts)
        self.assertEqual(calls, [("build_schrodinger_chainmap_files", False),
                                 ("build_schrodinger_peptide_maps", False),
                                 ("import_schrodinger_interactions", True),
                                 ("import_schrodinger_peptides", True),
                                 "contacts",
                                 ("import_schrodinger_interactions", False),
                                 ("import_schrodinger_peptides", False),
                                 ("remove_schrodinger_maps", False)])

    def test_a_failing_dry_run_stops_it_before_the_contacts_and_imports_nothing(self):
        from tools.management.commands import build_all_interactions as bai
        calls = []

        def dry_fails(name, **kw):
            calls.append((name, kw.get("dry_run", False)))
            if kw.get("dry_run") and name == "import_schrodinger_peptides":
                raise CommandError("1 structure(s) failed")

        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            with mock.patch.object(bai.Command, "prepare_input",
                                   lambda self, proc, pdbs: calls.append("contacts")), \
                    mock.patch.object(ligand_imports, "call_command", dry_fails):
                with self.assertRaisesRegex(CommandError, "failed"):
                    bai.Command().handle(**options(engine1_data_dir=d1, engine2_data_dir=d2,
                                                   proc=1))
        self.assertEqual(calls, [("build_schrodinger_chainmap_files", False),
                                 ("build_schrodinger_peptide_maps", False),
                                 ("import_schrodinger_interactions", True),
                                 ("import_schrodinger_peptides", True)])

    def test_a_failing_import_after_the_contacts_makes_the_command_fail(self):
        from tools.management.commands import build_all_interactions as bai
        calls = []

        def real_fails(name, **kw):
            calls.append((name, kw.get("dry_run", False)))
            if name in ligand_imports.COMMANDS and not kw.get("dry_run"):
                raise CommandError("1 structure(s) failed")

        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            with mock.patch.object(bai.Command, "prepare_input",
                                   lambda self, proc, pdbs: calls.append("contacts")), \
                    mock.patch.object(ligand_imports, "call_command", real_fails):
                with self.assertRaisesRegex(CommandError, "failed"):
                    bai.Command().handle(**options(engine1_data_dir=d1, engine2_data_dir=d2,
                                                   proc=1))
        self.assertEqual(calls, [("build_schrodinger_chainmap_files", False),
                                 ("build_schrodinger_peptide_maps", False),
                                 ("import_schrodinger_interactions", True),
                                 ("import_schrodinger_peptides", True),
                                 "contacts",
                                 ("import_schrodinger_interactions", False)])

    def test_a_contact_network_that_raises_does_not_stop_the_imports(self):
        # The contact network is the legacy calculation; its failures are
        # printed and logged as before, and the ligand imports still run.
        from tools.management.commands import build_all_interactions as bai
        calls = []

        def contacts_fail(self, proc, pdbs):
            calls.append("contacts")
            raise RuntimeError("contact network failed")

        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            with mock.patch.object(bai.Command, "prepare_input", contacts_fail), \
                    mock.patch.object(ligand_imports, "call_command",
                                      lambda name, **kw: calls.append((name, kw.get("dry_run", False)))):
                bai.Command().handle(**options(engine1_data_dir=d1, engine2_data_dir=d2, proc=1))
        self.assertEqual(calls, [("build_schrodinger_chainmap_files", False),
                                 ("build_schrodinger_peptide_maps", False),
                                 ("import_schrodinger_interactions", True),
                                 ("import_schrodinger_peptides", True),
                                 "contacts",
                                 ("import_schrodinger_interactions", False),
                                 ("import_schrodinger_peptides", False),
                                 ("remove_schrodinger_maps", False)])

    def test_each_lane_keeps_its_own_accounting_directory_from_one_plan(self):
        from tools.management.commands import build_all_interactions as bai
        paths = []
        with tempfile.TemporaryDirectory() as d1, tempfile.TemporaryDirectory() as d2:
            with mock.patch.object(bai.Command, "prepare_input", lambda self, proc, pdbs: None), \
                    mock.patch.object(ligand_imports, "steps", wraps=ligand_imports.steps) as plan, \
                    mock.patch.object(ligand_imports, "call_command",
                                      lambda name, **kw: paths.append(
                                          (name, kw["report_json"], kw["anomaly_csv"]))
                                      if name in ligand_imports.COMMANDS else None):
                bai.Command().handle(**options(engine1_data_dir=d1, engine2_data_dir=d2,
                                               engine1_report_dir=None, engine2_report_dir=None,
                                               proc=1))
        self.assertEqual(plan.call_count, 1)
        self.assertEqual(len(paths), 4)
        dirs = {}
        for name, report, anomalies in paths:
            # A run's report and anomaly list sit side by side.
            self.assertEqual(os.path.dirname(report), os.path.dirname(anomalies))
            dirs.setdefault(name, set()).add(os.path.dirname(report))
        self.assertEqual(sorted(dirs), ["import_schrodinger_interactions",
                                        "import_schrodinger_peptides"])
        for name, found in dirs.items():
            self.assertEqual(len(found), 1, name)
        # The lanes do not share a directory, and no run overwrites another's file.
        self.assertNotEqual(dirs["import_schrodinger_interactions"],
                            dirs["import_schrodinger_peptides"])
        files = [p for _n, report, anomalies in paths for p in (report, anomalies)]
        self.assertEqual(len(set(files)), 8)

    def test_skip_is_off_by_default_and_skipping_needs_no_delivery(self):
        from tools.management.commands import build_all_interactions as bai
        cmd = bai.Command()
        opts = vars(cmd.create_parser("manage.py", "build_all_interactions").parse_args([]))
        self.assertIs(opts["skip_ligand_import"], False)
        calls = []
        opts.update(options(skip_ligand_import=True, engine1_data_dir="/nonexistent-e1",
                            engine2_data_dir="/nonexistent-e2"))
        with mock.patch.object(bai.Command, "prepare_input",
                               lambda self, proc, pdbs: calls.append("contacts")), \
                mock.patch.object(ligand_imports, "call_command",
                                  lambda name, **kw: calls.append(name)):
            cmd.handle(**opts)
        self.assertEqual(calls, ["contacts"])

    def test_a_missing_engine2_delivery_alone_stops_it_before_the_contacts(self):
        from tools.management.commands import build_all_interactions as bai
        calls = []
        with tempfile.TemporaryDirectory() as d1:
            with mock.patch.object(bai.Command, "prepare_input",
                                   lambda self, proc, pdbs: calls.append("contacts")), \
                    mock.patch.object(ligand_imports, "call_command", no_import):
                with self.assertRaisesRegex(CommandError, "Engine 2 products are missing"):
                    bai.Command().handle(**options(engine1_data_dir=d1,
                                                   engine2_data_dir="/nonexistent-e2", proc=1))
        self.assertEqual(calls, [])

    def test_a_missing_delivery_stops_it_before_the_contacts(self):
        from tools.management.commands import build_all_interactions as bai
        calls = []
        with mock.patch.object(bai.Command, "prepare_input",
                               lambda self, proc, pdbs: calls.append("contacts")), \
                mock.patch.object(ligand_imports, "call_command", no_import):
            with self.assertRaisesRegex(CommandError, "Engine 1 products are missing"):
                bai.Command().handle(**options(engine1_data_dir="/nonexistent-e1",
                                               engine2_data_dir="/nonexistent-e2", proc=1))
        self.assertEqual(calls, [])


if __name__ == "__main__":
    unittest.main()
