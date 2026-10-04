"""Unit tests for the ligand import steps of build_all (no database).

    python -c "import django; django.setup(); import unittest; \\
        unittest.main(module='build.test_build_all_ligands', argv=['x'])"
"""

import unittest

import tempfile
from unittest import mock

from django.core.management.base import CommandError

from build import ligand_imports
from build.management.commands import build_all


def options(**kw):
    out = {"skip_ligand_import": False, "engine1_data_dir": "/e1", "engine2_data_dir": "/e2",
           "engine1_report_dir": "/r1", "engine2_report_dir": "/r2", "phase": None}
    out.update(kw)
    return out


class LigandImportStepsTests(unittest.TestCase):
    def test_both_dry_runs_come_before_either_import(self):
        steps = build_all.Command().ligand_import_steps(options())
        self.assertEqual([(c, o.get("dry_run", False), o["data_dir"]) for c, o in steps], [
            ("import_schrodinger_interactions", True, "/e1"),
            ("import_schrodinger_peptides", True, "/e2"),
            ("import_schrodinger_interactions", False, "/e1"),
            ("import_schrodinger_peptides", False, "/e2"),
        ])
        self.assertEqual([o["report_json"] for _c, o in steps],
                         ["/r1/report.dryrun.json", "/r2/report.dryrun.json",
                          "/r1/report.json", "/r2/report.json"])

    def test_skipping_imports_nothing(self):
        self.assertEqual(build_all.Command().ligand_import_steps(options(skip_ligand_import=True)), [])


class BuildAllInteractionsTests(unittest.TestCase):
    """build_all_interactions keeps its job: contact network, then the same imports."""

    def test_it_takes_the_import_options_and_runs_the_imports_after_the_contacts(self):
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
        self.assertEqual(calls, ["contacts",
                                 ("import_schrodinger_interactions", True),
                                 ("import_schrodinger_peptides", True),
                                 ("import_schrodinger_interactions", False),
                                 ("import_schrodinger_peptides", False)])

    def test_a_missing_delivery_stops_it_before_the_contacts(self):
        from tools.management.commands import build_all_interactions as bai
        calls = []
        with mock.patch.object(bai.Command, "prepare_input",
                               lambda self, proc, pdbs: calls.append("contacts")):
            with self.assertRaisesRegex(CommandError, "Engine 1 products are missing"):
                bai.Command().handle(**options(engine1_data_dir="/nonexistent-e1",
                                               engine2_data_dir="/nonexistent-e2", proc=1))
        self.assertEqual(calls, [])


if __name__ == "__main__":
    unittest.main()
