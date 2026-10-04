"""Unit tests for interaction.schrodinger_complex (no database).

    python -c "import django; django.setup(); import unittest; \\
        unittest.main(module='interaction.test_schrodinger_complex', argv=['x'])"
"""

import unittest

from interaction import schrodinger_complex as cx


def line(record, serial, name, resname, chain, resnum, x, y, z, element, icode=" "):
    return "{:<6}{:>5} {:<4} {:>3} {}{:>4}{}   {:8.3f}{:8.3f}{:8.3f}  1.00 20.00          {:>2}".format(
        record, serial, name, resname, chain, resnum, icode, x, y, z, element)


TEXT = "\n".join([
    line("ATOM", 1, "N", "ASP", "R", 147, 1.0, 0.0, 0.0, "N"),
    line("ATOM", 2, "OD1", "ASP", "R", 147, 2.0, 0.0, 0.0, "O"),
    line("ATOM", 3, "N", "TRP", "R", 293, 9.0, 0.0, 0.0, "N"),
    line("ATOM", 4, "N", "SER", "R", 300, 20.0, 0.0, 0.0, "N"),
    line("ATOM", 5, "N", "GLY", "R", 147, 30.0, 0.0, 0.0, "N", icode="A"),
    line("HETATM", 6, "C1", "ZMA", "R", 401, 5.000, 5.000, 5.000, "C"),
    line("HETATM", 7, "N2", "ZMA", "R", 401, 6.000, 5.000, 5.000, "N"),
    line("HETATM", 8, "H1", "ZMA", "R", 401, 6.500, 5.000, 5.000, "H"),
    line("HETATM", 9, "C1", "CLR", "R", 402, 40.0, 0.0, 0.0, "C"),
    line("HETATM", 10, "O", "HOH", "R", 501, 5.010, 5.000, 5.000, "O"),
    line("ATOM", 11, "N", "TYR", "P", 1, 3.0, 3.0, 3.0, "N"),
    line("ATOM", 12, "CA", "TYR", "P", 1, 3.5, 3.0, 3.0, "C"),
    line("ATOM", 13, "N", "GLY", "P", 2, 4.0, 3.0, 3.0, "N"),
    "ENDMDL",
    line("HETATM", 14, "C1", "ZMA", "R", 401, 5.000, 5.000, 5.000, "C"),
])


def rows(text):
    return text.splitlines()[:-1]


class ComplexTextTests(unittest.TestCase):
    def test_the_ligand_by_coordinates_and_the_written_residues(self):
        # A product coordinate 0.007 A off, as measured on the corpus.
        out = cx.complex_text(TEXT, "R", {147, 293}, ligand_xyz=[(5.007, 5.0, 5.0)])
        names = [(l[17:20], l[22:26].strip(), l[26]) for l in rows(out)]
        self.assertEqual(names, [("ASP", "147", " "), ("ASP", "147", " "), ("TRP", "293", " "),
                                 ("ZMA", "401", " "), ("ZMA", "401", " "), ("ZMA", "401", " ")])
        self.assertEqual(out.splitlines()[-1], "END")
        # Lines are copied as they are, hydrogen included; water never is.
        self.assertIn(TEXT.splitlines()[7], out.splitlines())
        self.assertNotIn("HOH", out)

    def test_no_ligand_no_file(self):
        self.assertEqual(cx.complex_text(TEXT, "R", {147}, ligand_xyz=[(5.06, 5.0, 5.0)]), "")
        self.assertEqual(cx.complex_text(TEXT, "R", {147}, ligand_xyz=[]), "")
        self.assertEqual(cx.complex_text("", "R", {147}, ligand_xyz=[(5.0, 5.0, 5.0)]), "")

    def test_a_chain_anchor_takes_the_whole_chain(self):
        out = cx.complex_text(TEXT, "R", {300}, ligand_chain="P")
        self.assertEqual([(l[17:20], l[21], l[22:26].strip()) for l in rows(out)],
                         [("SER", "R", "300"), ("TYR", "P", "1"), ("TYR", "P", "1"), ("GLY", "P", "2")])

    def test_only_the_receptor_chain_and_no_insertion_code(self):
        out = cx.complex_text(TEXT, "Q", {147}, ligand_xyz=[(5.0, 5.0, 5.0)])
        self.assertEqual({l[17:20] for l in rows(out)}, {"ZMA"})
        out = cx.complex_text(TEXT, "R", {147}, ligand_xyz=[(5.0, 5.0, 5.0)])
        self.assertNotIn("GLY", out)

    def test_only_the_first_model(self):
        out = cx.complex_text(TEXT, "R", set(), ligand_xyz=[(5.0, 5.0, 5.0)])
        self.assertEqual(len(rows(out)), 3)

    def test_ligand_line_xyz_skips_hydrogens(self):
        self.assertEqual(cx.ligand_line_xyz(TEXT.splitlines()[5:8]),
                         [(5.0, 5.0, 5.0), (6.0, 5.0, 5.0)])


if __name__ == "__main__":
    unittest.main()
