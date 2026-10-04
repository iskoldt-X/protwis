"""Unit tests for the ligand import steps of build_all (no database).

    python -c "import django; django.setup(); import unittest; \\
        unittest.main(module='build.test_build_all_ligands', argv=['x'])"
"""

import unittest

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


if __name__ == "__main__":
    unittest.main()
