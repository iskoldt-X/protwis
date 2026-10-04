"""Unit tests for interaction.stored_interactions (no database).

    python -c "import django; django.setup(); import unittest; \\
        unittest.main(module='interaction.test_stored_interactions', argv=['x'])"
"""

import unittest

from interaction import stored_interactions as si
from interaction.views import regexaa


class BuildResultsTests(unittest.TestCase):
    ROWS = [
        ("ZMA", "D", 113, "polar_donor_protein", "polar (hydrogen bond)", "polar", "protein"),
        ("ZMA", "F", 290, "hyd", "hydrophobic", "hydrophobic", ""),
        ("CLR", "W", 286, "hyd", "hydrophobic", "hydrophobic", None),
        ("ZMA", "X", 400, "hyd", "hydrophobic", "hydrophobic", ""),
        ("ZMA", "N", 253, "acc", "accessible", "hidden", ""),
    ]

    def test_shape_order_and_residue_names(self):
        out = si.build_results(self.ROWS, "A")
        self.assertEqual(list(out), ["ZMA", "CLR"])
        self.assertEqual(out["ZMA"]["score"], 3)
        self.assertEqual(out["ZMA"]["interactions"][0],
                         ["ASP113A", "", "polar_donor_protein", "polar (hydrogen bond)", "polar", "protein"])
        self.assertEqual(out["CLR"]["interactions"], [["TRP286A", "", "hyd", "hydrophobic", "hydrophobic", ""]])
        # The page splits the residue with regexaa.
        self.assertEqual(regexaa(out["ZMA"]["interactions"][0][0]), ("D", "113", "A"))

    def test_ties_keep_key_order_and_empty_is_empty(self):
        out = si.build_results([("B", "A", 1, "hyd", "h", "hydrophobic", ""),
                                ("A", "A", 2, "hyd", "h", "hydrophobic", "")], "R")
        self.assertEqual(list(out), ["A", "B"])
        self.assertEqual(si.build_results([], "R"), {})

    def test_ligand_key(self):
        self.assertEqual(si.ligand_key("zma", "ZM241385"), "ZMA")
        self.assertEqual(si.ligand_key(" pep ", "DAMGO"), "DAMGO")


if __name__ == "__main__":
    unittest.main()
