"""Unit tests for the database-free layer of interaction.schrodinger_import.

They need Django configured (the module imports models) but never touch the
database, so they run with plain unittest:

    python -c "import django; django.setup(); import unittest; \
        unittest.main(module='interaction.test_schrodinger_import', argv=['x'])"
"""

import os
import shutil
import tempfile
import types
import unittest

from interaction import schrodinger_import as si


def row(family, direction="", seq=100, aa="F", chain="A", atom="CB", block="ATOM"):
    return {
        "feature_family": family,
        "direction": direction,
        "receptor_atom_name": atom,
        "receptor_pdb_block": block,
        "receptor_residue": {
            "name_1_letter": aa,
            "pdb_residue_number": seq,
            "chain_id": chain,
            "insertion_code": "",
        },
    }


class ParseChainResTests(unittest.TestCase):

    def test_chain_and_number(self):
        self.assertEqual(si.parse_chain_res("A:408"), ("A", "408", ""))

    def test_negative_number_and_insertion_code(self):
        self.assertEqual(si.parse_chain_res("R:-5"), ("R", "-5", ""))
        self.assertEqual(si.parse_chain_res("A:12B"), ("A", "12", "B"))

    def test_chain_only_or_empty_is_none(self):
        for value in ("L", "", None, "A:", ":12"):
            self.assertIsNone(si.parse_chain_res(value), value)


class SelectInstancesTests(unittest.TestCase):

    NAMES = ["CAU_A_408", "CAU_B_408", "CA_A_501", "CLR_A_1203", "CLR_A_1204"]

    def test_exact_copy_only(self):
        self.assertEqual(si.select_instances("CAU", "A:408", self.NAMES),
                         (["CAU_A_408"], "exact"))

    def test_exact_missing_returns_nothing(self):
        self.assertEqual(si.select_instances("CAU", "A:999", self.NAMES), ([], "exact_missing"))

    def test_all_copies_when_chain_res_empty(self):
        self.assertEqual(si.select_instances("CLR", None, self.NAMES),
                         (["CLR_A_1203", "CLR_A_1204"], "all_copies"))

    def test_het_prefix_does_not_leak(self):
        # "CA" must not pick up "CAU_*" instances.
        self.assertEqual(si.select_instances("CA", None, self.NAMES), (["CA_A_501"], "all_copies"))

    def test_het_case_insensitive(self):
        self.assertEqual(si.select_instances("cau", "B:408", self.NAMES),
                         (["CAU_B_408"], "exact"))


class RoutingTests(unittest.TestCase):

    # Every (family, direction) pair in the 2026-09-17 production tree
    # (12,718 instance YAMLs, 124,276 rows).
    PRODUCTION_PAIRS = {
        ("HPhob", ""): "hyd",
        ("Acceptor", "ligand-acceptor"): "polar_donor_protein",
        ("Donor", "ligand-donor"): "polar_acceptor_protein",
        ("Aromatic", "edge-to-face"): "aro_ef_protein",
        ("Aromatic", "face-to-face"): "aro_ff",
        ("NegCharge", "neg-pos"): "polar_double_pos_protein",
        ("PosCharge", "pos-neg"): "polar_double_neg_protein",
        ("Metal", ""): "metal_coordination_protein",
        ("PiCat", "ligand-cation"): "aro_ion_protein",
        ("PiCat", "receptor-cation"): "aro_ion_protein",
        ("Wat-HBond", ""): "water_bridge_protein",
        ("XBond", ""): "halogen_protein",
    }

    def test_every_production_pair_routes(self):
        for (family, direction), slug in self.PRODUCTION_PAIRS.items():
            self.assertEqual(si.resolve_slug(family, direction), slug, (family, direction))

    def test_none_direction_is_empty(self):
        self.assertEqual(si.resolve_slug("HPhob", None), "hyd")

    def test_unknown_pair_raises(self):
        with self.assertRaises(si.UnroutableRow):
            si.resolve_slug("Acceptor", "ligand-donor")

    def test_backbone_override(self):
        self.assertEqual(si.apply_backbone_override("polar_donor_protein", "N"), "polar_backbone")
        self.assertEqual(si.apply_backbone_override("polar_acceptor_protein", " O "),
                         "polar_backbone")
        self.assertEqual(si.apply_backbone_override("polar_donor_protein", "OG"),
                         "polar_donor_protein")
        self.assertEqual(si.apply_backbone_override("hyd", "N"), "hyd")
        self.assertEqual(si.apply_backbone_override("polar_donor_protein", None),
                         "polar_donor_protein")

    def test_only_main_chain_n_and_o_promote(self):
        for atom in ("CA", "C", "CB", "OG", "OG1", "ND2", "NE2", "OH", "H", ""):
            for slug in ("polar_donor_protein", "polar_acceptor_protein"):
                self.assertEqual(si.apply_backbone_override(slug, atom), slug, (slug, atom))

    def test_required_slugs_exclude_water_bridge(self):
        self.assertEqual(si.required_slugs(), frozenset({
            "hyd", "polar_donor_protein", "polar_acceptor_protein", "aro_ef_protein",
            "aro_ff", "polar_double_pos_protein", "polar_double_neg_protein",
            "metal_coordination_protein", "aro_ion_protein", "halogen_protein",
            "polar_backbone",
        }))


class PlanRowsTests(unittest.TestCase):

    def test_accounting_identity_and_each_bucket(self):
        rows = [
            row("HPhob", seq=100),
            row("HPhob", seq=100),                      # duplicate (100, hyd)
            row("Acceptor", "ligand-acceptor", seq=100, aa="F", atom="N"),  # backbone
            row("Wat-HBond", seq=101),                  # excluded family
            row("HPhob", seq=102, aa="X"),              # non-standard residue
            row("HPhob", seq=103, chain="B"),           # other chain
            row("Metal", seq=104, aa="H"),
        ]
        records, counts, by_chain = si.plan_rows(rows, "A")
        self.assertEqual(counts["rows_in"], 7)
        self.assertEqual(counts["duplicate"], 1)
        self.assertEqual(counts["excluded_family"], 1)
        self.assertEqual(counts["nonstandard_residue"], 1)
        self.assertEqual(counts["other_chain"], 1)
        self.assertEqual(by_chain, {"B": 1})
        self.assertEqual(counts["planned"], len(records))
        self.assertEqual(
            counts["rows_in"],
            counts["excluded_family"] + counts["nonstandard_residue"]
            + counts["other_chain"] + counts["duplicate"] + counts["planned"])
        self.assertEqual([(r["sequence_number"], r["slug"]) for r in records],
                         [(100, "hyd"), (100, "polar_backbone"),
                          (104, "metal_coordination_protein")])

    def test_distinct_slugs_on_one_residue_both_survive(self):
        rows = [row("PosCharge", "pos-neg", seq=113, aa="D"),
                row("Donor", "ligand-donor", seq=113, aa="D", atom="OD1")]
        records, counts, _ = si.plan_rows(rows, "A")
        self.assertEqual({r["slug"] for r in records},
                         {"polar_double_neg_protein", "polar_acceptor_protein"})
        self.assertEqual(counts["duplicate"], 0)

    def test_empty_preferred_chain_keeps_every_chain(self):
        records, counts, _ = si.plan_rows([row("HPhob", chain="B")], "")
        self.assertEqual(len(records), 1)
        self.assertEqual(counts["other_chain"], 0)

    def test_unroutable_row_raises(self):
        with self.assertRaises(si.UnroutableRow):
            si.plan_rows([row("Bogus", "x")], "A")

    def test_preferred_chain_takes_first(self):
        self.assertEqual(si.preferred_chain_of("A,B"), "A")
        self.assertEqual(si.preferred_chain_of(None), "")


class ProductFilesTests(unittest.TestCase):

    def setUp(self):
        self.root = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.root)

    def _write(self, rel, text):
        path = os.path.join(self.root, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def test_instance_discovery(self):
        good = self._write("2RH1/CAU_A_408/CAU_A_408.yaml", "result: {interactions: []}\n")
        self._write("2RH1/summary.yaml", "x: 1\n")
        self._write("2RH1/not_an_instance/not_an_instance.yaml", "x: 1\n")
        os.makedirs(os.path.join(self.root, "2RH1", "CLR_A_1"))  # no YAML
        self.assertEqual(si.instance_yaml_paths(self.root, "2rh1"), {"CAU_A_408": good})
        self.assertEqual(si.instance_yaml_paths(self.root, "9ZZZ"), {})

    def test_read_rows(self):
        ok = self._write("a.yaml", "result:\n  interactions: []\n")
        self.assertEqual(si.read_instance_rows(ok), [])
        for text in ("result: {}\n", "- 1\n", "", "result: {interactions: 3}\n", "a: [\n"):
            bad = self._write("b.yaml", text)
            with self.assertRaises(si.MalformedProduct, msg=text):
                si.read_instance_rows(bad)


class ScopeTests(unittest.TestCase):

    @staticmethod
    def sli(reference, ligand_type):
        return types.SimpleNamespace(
            pdb_reference=reference,
            ligand=types.SimpleNamespace(ligand_type=types.SimpleNamespace(slug=ligand_type)))

    def test_scope(self):
        self.assertTrue(si.is_in_scope(self.sli("CAU", "small-molecule")))
        self.assertTrue(si.is_in_scope(self.sli("CLR", "lipid")))
        self.assertFalse(si.is_in_scope(self.sli("pep", "small-molecule")))
        self.assertFalse(si.is_in_scope(self.sli("APO", "none")))
        self.assertFalse(si.is_in_scope(self.sli("XYZ", "peptide")))
        self.assertFalse(si.is_in_scope(self.sli("XYZ", "protein")))
        self.assertFalse(si.is_in_scope(self.sli("", "small-molecule")))
        self.assertFalse(si.is_in_scope(self.sli(None, "small-molecule")))


if __name__ == "__main__":
    unittest.main()
