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


def row(family, direction="", seq=100, aa="F", chain="A", atom="CB", block="ATOM", lig=""):
    return {
        "feature_family": family,
        "direction": direction,
        "receptor_atom_name": atom,
        "receptor_pdb_block": block,
        "ligand_pdb_block": lig,
        "receptor_residue": {
            "name_1_letter": aa,
            "pdb_residue_number": seq,
            "chain_id": chain,
            "insertion_code": "",
        },
    }


ANCHOR_MAP = {
    ("6ZIN", "Q6Q", "A:1000"): {"status": "ok", "instance": "Q6Q_AAA_1000", "note": ""},
    ("6N51", "QUS", "A:903"): {"status": "errata", "instance": "QUS_A_903", "note": "label says B"},
    ("8E0G", "A1A7R", "A:54"): {"status": "no_product", "instance": "", "note": "no instance"},
    ("7V68", "2CU", "R:502"): {"status": "unresolved", "instance": "", "note": "drift"},
    ("9X9X", "U7D", "R:601"): {"status": "ok", "instance": "U7D_R_601", "note": ""},
    ("9X9X", "U7D", "R:602"): {"status": "no_product", "instance": "", "note": "gone"},
    ("7E2X", "CLR", ""): {"status": "all_copies", "instance": "CLR_A_1;CLR_A_2", "note": ""},
    ("6CMO", "RET", ""): {"status": "no_product", "instance": "", "note": "not modelled"},
}


class AnchorInstancesTests(unittest.TestCase):

    def test_renamed_chain(self):
        self.assertEqual(si.anchor_instances("6zin", "q6q", "A:1000", ANCHOR_MAP, []),
                         (["Q6Q_AAA_1000"], "mapped", []))

    def test_errata_still_imports_and_is_reported(self):
        names, mode, notes = si.anchor_instances("6N51", "QUS", "A:903", ANCHOR_MAP, [])
        self.assertEqual((names, mode), (["QUS_A_903"], "mapped"))
        self.assertTrue(notes and notes[0].startswith("errata: "))

    def test_no_product_and_partial(self):
        self.assertEqual(si.anchor_instances("8E0G", "A1A7R", "A:54", ANCHOR_MAP, [])[:2], ([], "no_product"))
        self.assertEqual(si.anchor_instances("9X9X", "U7D", "R:601, R:602", ANCHOR_MAP, [])[:2],
                         (["U7D_R_601"], "mapped_partial"))

    def test_chain_res_without_residue_uses_all_copies_row(self):
        self.assertEqual(si.anchor_instances("7E2X", "CLR", None, ANCHOR_MAP, [])[:2],
                         (["CLR_A_1", "CLR_A_2"], "all_copies"))
        self.assertEqual(si.anchor_instances("6CMO", "RET", "", ANCHOR_MAP, [])[:2], ([], "no_product"))

    def test_unresolved_and_missing_are_loud(self):
        with self.assertRaises(si.UnresolvedAnchor):
            si.anchor_instances("7V68", "2CU", "R:502", ANCHOR_MAP, [])
        with self.assertRaises(si.MapMismatch):
            si.anchor_instances("2RH1", "CAU", "A:408", ANCHOR_MAP, [])

    def test_receptor_chain(self):
        rmap = {"6ZIN": {"status": "ok", "auth_chain": "AAA", "note": ""},
                "7V68": {"status": "unresolved", "auth_chain": "", "note": "x"}}
        self.assertEqual(si.receptor_chain("6zin", rmap), "AAA")
        with self.assertRaises(si.UnresolvedAnchor):
            si.receptor_chain("7V68", rmap)
        with self.assertRaises(si.MapMismatch):
            si.receptor_chain("2RH1", rmap)


class MapGuardTests(unittest.TestCase):

    def test_no_product_with_a_copy_in_the_tree_is_a_map_mismatch(self):
        with self.assertRaises(si.MapMismatch):
            si.anchor_instances("8E0G", "A1A7R", "A:54", ANCHOR_MAP, ["A1A7R_A_54"])
        with self.assertRaises(si.MapMismatch):
            si.anchor_instances("6CMO", "RET", "", ANCHOR_MAP, ["RET_A_1"])
        # a copy of another HET does not count
        self.assertEqual(si.anchor_instances("8E0G", "A1A7R", "A:54", ANCHOR_MAP, ["CLR_A_403"])[:2],
                         ([], "no_product"))

    @staticmethod
    def sli(ref, chain_res):
        return types.SimpleNamespace(pdb_reference=ref, chain_res=chain_res)

    def test_map_must_cover_every_database_copy(self):
        amap = {("9X9X", "U7D", "R:601"): {}, ("9X9X", "U7D", "R:602"): {}, ("9X9X", "CLR", ""): {}}
        si.check_map_covers("9X9X", [self.sli("U7D", "R:601, R:602"), self.sli("CLR", None)], amap)
        with self.assertRaises(si.MapMismatch):   # database has a copy the map does not list
            si.check_map_covers("9X9X", [self.sli("U7D", "R:601, R:602, R:603"), self.sli("CLR", "")], amap)

    def test_map_may_list_copies_the_database_cannot_hold(self):
        """The annotation splits copies by chain; SLI has no copy dimension."""
        amap = {("9X9X", "U7D", "R:601"): {}, ("9X9X", "U7D", "S:601"): {}, ("9X9X", "CLR", ""): {}}
        self.assertEqual(
            si.check_map_covers("9X9X", [self.sli("U7D", "R:601"), self.sli("CLR", "")], amap),
            [("U7D", "S:601")])
        # a whole HET the database does not have is also extra, not an error,
        # but it is still handed back to be reported
        self.assertEqual(si.check_map_covers("9X9X", [self.sli("U7D", "R:601")], amap),
                         [("CLR", ""), ("U7D", "S:601")])

    def test_map_covers_is_case_insensitive_in_the_pdb_code(self):
        amap = {("9X9X", "U7D", "R:601"): {}}
        self.assertEqual(si.check_map_covers("9x9x", [self.sli("U7D", "R:601")], amap), [])

    def test_map_rows_of_another_structure_do_not_count_as_coverage(self):
        amap = {("OTHR", "U7D", "R:601"): {}}
        with self.assertRaises(si.MapMismatch):
            si.check_map_covers("9X9X", [self.sli("U7D", "R:601")], amap)

    def _write(self, text):
        fh = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
        fh.write(text)
        fh.close()
        self.addCleanup(os.unlink, fh.name)
        return fh.name

    def test_map_files_header_and_duplicates(self):
        header = "# dump_id\t20260917_phase2\n# structures\t1\n"
        path = self._write(header + "pdb\thet\ttoken\tinstance\tstatus\n6zin\tq6q\tA:1000\tQ6Q_AAA_1000\tok\n")
        head, table = si.load_anchor_map(path)
        self.assertEqual(head, {"dump_id": "20260917_phase2", "structures": "1"})
        self.assertEqual(table[("6ZIN", "Q6Q", "A:1000")]["instance"], "Q6Q_AAA_1000")
        dup = self._write(header + "pdb\thet\ttoken\n6ZIN\tQ6Q\tA:1000\n6ZIN\tQ6Q\tA:1000\n")
        with self.assertRaises(si.MapMismatch):
            si.load_anchor_map(dup)
        rdup = self._write(header + "pdb\tauth_chain\n6ZIN\tAAA\n6zin\tAAA\n")
        with self.assertRaises(si.MapMismatch):
            si.load_receptor_map(rdup)

    def test_a_quoted_newline_in_a_row_cannot_forge_a_header(self):
        """A note may hold a newline; its continuation must stay in the body."""
        note = 'unreadable:\n# schema\tengine1-chainmap/999'
        path = self._write("# schema\tengine1-chainmap/1\n"
                           "pdb\thet\ttoken\tnote\n"
                           '6ZIN\tQ6Q\tA:1\t"%s"\n' % note)
        head, table = si._read_map(path)
        self.assertEqual(head, {"schema": "engine1-chainmap/1"})
        self.assertEqual(len(table), 1)
        self.assertEqual(table[0]["note"], note)

    def test_fingerprints(self):
        from interaction import schrodinger_chain_map as cm
        rmap = {"6ZIN": {"gpcrdb_text_sha256": cm.text_sha256("ATOM 1\n"),
                         "product_instances_sha256": cm.instances_sha256(["Q6Q_AAA_1000", "CLR_AAA_1"])}}
        si.check_fingerprints("6zin", rmap, "ATOM 1\n", ["CLR_AAA_1", "Q6Q_AAA_1000"])
        with self.assertRaises(si.MapMismatch):
            si.check_fingerprints("6ZIN", rmap, "ATOM 2\n", ["CLR_AAA_1", "Q6Q_AAA_1000"])
        with self.assertRaises(si.MapMismatch):
            si.check_fingerprints("6ZIN", rmap, "ATOM 1\n", ["Q6Q_AAA_1000"])
        with self.assertRaises(si.MapMismatch):
            si.check_fingerprints("2RH1", rmap, "", [])

    def test_ok_receptor_row_without_chain_is_refused(self):
        with self.assertRaises(si.MapMismatch):
            si.receptor_chain("X", {"X": {"status": "ok", "auth_chain": "", "note": ""}})



CHAINMAP_HEADER = "".join(
    "# {}\t{}\n".format(k, v) for k, v in [
        ("schema", "engine1-chainmap/1"), ("pdb", "6ZIN"),
        ("annotation_commit", "9fe1875"), ("ligands_sha256", "aa"),
        ("structures_sha256", "bb"), ("gpcrdb_pdb_sha256", "cc"),
        ("cif_sha256", "dd"), ("builder_sha256", "ee"),
        ("receptor.preferred_chain", "A"), ("receptor.auth_chain", "A"),
        ("receptor.status", "ok"), ("receptor.method", "exact"),
        ("receptor.n_ca_gpcrdb", "10"), ("receptor.n_ca_matched", "10"),
        ("receptor.renumbered", "0"), ("receptor.note", ""),
        ("receptor.gpcrdb_text_sha256", "ff"), ("receptor.product_instances_sha256", "gg"),
    ])
CHAINMAP_BODY = ("pdb\thet\ttoken\tinstance\tstatus\n"
                 "6ZIN\tq6q\tA:1000\tQ6Q_A_1000\tok\n")


class ChainmapDirTests(unittest.TestCase):
    """The per-PDB chain maps the build consumes without being told anything."""

    def setUp(self):
        self.root = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.root)

    def write(self, pdb, text):
        d = os.path.join(self.root, pdb)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, si.CHAINMAP_NAME), "w") as fh:
            fh.write(text)
        return os.path.join(d, si.CHAINMAP_NAME)

    def test_one_file_round_trips(self):
        path = self.write("6ZIN", CHAINMAP_HEADER + CHAINMAP_BODY)
        pdb, anchors, receptor, prov = si.load_chainmap(path)
        self.assertEqual(pdb, "6ZIN")
        self.assertEqual(anchors[("6ZIN", "Q6Q", "A:1000")]["instance"], "Q6Q_A_1000")
        self.assertEqual(receptor["auth_chain"], "A")
        self.assertEqual(receptor["product_instances_sha256"], "gg")
        self.assertEqual(prov["annotation_commit"], "9fe1875")

    def test_an_unknown_schema_is_refused(self):
        for schema in ("engine1-chainmap/2", "", "engine1-chainmap"):
            path = self.write("6ZIN", CHAINMAP_HEADER.replace(
                "engine1-chainmap/1", schema) + CHAINMAP_BODY)
            with self.assertRaises(si.MapMismatch):
                si.load_chainmap(path)
        # a file with no header at all
        path = self.write("6ZIN", CHAINMAP_BODY)
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap(path)

    def test_a_missing_receptor_field_is_refused(self):
        path = self.write("6ZIN", CHAINMAP_HEADER.replace(
            "# receptor.auth_chain\tA\n", "") + CHAINMAP_BODY)
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap(path)

    def test_a_row_of_another_structure_is_refused(self):
        path = self.write("6ZIN", CHAINMAP_HEADER + CHAINMAP_BODY +
                          "2RH1\tCAU\tA:408\tCAU_A_408\tok\n")
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap(path)

    def test_a_duplicate_anchor_key_is_refused(self):
        path = self.write("6ZIN", CHAINMAP_HEADER + CHAINMAP_BODY +
                          "6ZIN\tQ6Q\tA:1000\tQ6Q_A_1000\tok\n")
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap(path)

    def test_a_directory_without_a_map_is_named_not_skipped(self):
        self.write("6ZIN", CHAINMAP_HEADER + CHAINMAP_BODY)
        os.makedirs(os.path.join(self.root, "2RH1"))
        anchors, receptors, missing, prov = si.load_chainmap_dir(self.root, ["6ZIN", "2RH1"])
        self.assertEqual(missing, ["2RH1"])
        self.assertEqual(sorted(receptors), ["6ZIN"])
        self.assertEqual(len(anchors), 1)

    def test_provenance_counts_a_merged_tree(self):
        self.write("6ZIN", CHAINMAP_HEADER + CHAINMAP_BODY)
        self.write("2RH1", CHAINMAP_HEADER.replace("# pdb\t6ZIN", "# pdb\t2RH1")
                   .replace("9fe1875", "deadbee") + "pdb\thet\ttoken\n")
        _, _, _, prov = si.load_chainmap_dir(self.root, ["6ZIN", "2RH1"])
        self.assertEqual(prov["annotation_commit"], {"9fe1875": 1, "deadbee": 1})
        self.assertEqual(prov["ligands_sha256"], {"aa": 2})

    def test_a_file_in_the_wrong_directory_is_refused(self):
        d = os.path.join(self.root, "2RH1")
        os.makedirs(d)
        with open(os.path.join(d, si.CHAINMAP_NAME), "w") as fh:
            fh.write(CHAINMAP_HEADER + CHAINMAP_BODY)      # header says 6ZIN
        with self.assertRaises(si.MapMismatch):
            si.load_chainmap_dir(self.root, ["2RH1"])

    def test_product_pdb_codes_is_directories_only(self):
        self.write("6ZIN", CHAINMAP_HEADER + CHAINMAP_BODY)
        os.makedirs(os.path.join(self.root, "2rh1"))
        with open(os.path.join(self.root, "_import_anomalies.csv"), "w") as fh:
            fh.write("x\n")
        self.assertEqual(si.product_pdb_codes(self.root), ["2RH1", "6ZIN"])

class StandardLigandLineTests(unittest.TestCase):

    # Lines as the producer writes them (no altloc column; widened for long names).
    CAU = "HETATM    1  O17CAU A 408     -33.477  10.957   8.170  1.00 50.96           O"
    FIVE = "HETATM 4817  C9 A1C5S A1202     -20.354 -13.384  28.940  1.00 22.97           C"
    AAA = "HETATM   10  N4 T0B AAA 601      -2.398  28.100  15.000  1.00 30.00           N"
    ABUT = "HETATM    7  C1 LIG A   1    -100.123-200.456-300.789  1.0033084.00           C"
    # A four-character atom name runs straight into the residue name, as the producer writes it.
    HNAME = "HETATM   40 HN12CAU A 408     -30.000  11.000   9.000  1.00 50.96           H"

    def check_standard(self, out, resname, chain, resnum, xyz):
        self.assertEqual(len(out), 78)
        self.assertEqual(out[16], " ")                 # altloc column present and blank
        self.assertEqual(out[17:20].strip(), resname)
        self.assertEqual(out[21], chain)
        self.assertEqual(out[22:26].strip(), resnum)
        self.assertEqual((float(out[30:38]), float(out[38:46]), float(out[46:54])), xyz)

    def test_ordinary_line(self):
        out, capped = si.standard_ligand_line(self.CAU, "CAU", "A", "408", "", "A")
        self.check_standard(out, "CAU", "A", "408", (-33.477, 10.957, 8.17))
        self.assertEqual(out[12:16], " O17")
        self.assertFalse(capped)

    def test_five_char_code_is_cut_to_three(self):
        out, _ = si.standard_ligand_line(self.FIVE, "A1C5S", "A", "1202", "", "A")
        self.check_standard(out, "A1C", "A", "1202", (-20.354, -13.384, 28.94))

    def test_multichar_chain_becomes_gpcrdb_chain(self):
        out, _ = si.standard_ligand_line(self.AAA, "T0B", "AAA", "601", "", "A")
        self.check_standard(out, "T0B", "A", "601", (-2.398, 28.1, 15.0))

    def test_abutting_fields_and_capped_b_factor(self):
        out, capped = si.standard_ligand_line(self.ABUT, "LIG", "A", "1", "", "A")
        self.check_standard(out, "LIG", "A", "1", (-100.123, -200.456, -300.789))
        self.assertTrue(capped)
        self.assertEqual(out[60:66], "999.99")

    def test_four_char_hydrogen_name(self):
        out, _ = si.standard_ligand_line(self.HNAME, "CAU", "A", "408", "", "A")
        self.assertEqual(out[12:16], "HN12")
        self.assertEqual(out[76:78], " H")

    def test_mismatches_raise(self):
        for het, chain, resnum in (("CAZ", "A", "408"), ("CAU", "B", "408"), ("CAU", "A", "409")):
            with self.assertRaises(si.MalformedLigandLine):
                si.standard_ligand_line(self.CAU, het, chain, resnum, "", "A")
        with self.assertRaises(si.MalformedLigandLine):
            si.standard_ligand_line("REMARK nothing here", "CAU", "A", "408", "", "A")
        with self.assertRaises(si.MalformedLigandLine):
            si.standard_ligand_line(self.CAU[:60], "CAU", "A", "408", "", "A")

    def test_block_and_instance_chains(self):
        text, capped = si.standard_ligand_block(self.CAU + "\n\n" + self.HNAME + "\n", "CAU_A_408", "A")
        self.assertEqual(len(text.splitlines()), 2)
        amap = {("6ZIN", "Q6Q", "A:1000"): {"status": "ok", "instance": "Q6Q_AAA_1000", "note": ""},
                ("7E2X", "CLR", ""): {"status": "all_copies", "instance": "CLR_R_602;CLR_R_603", "note": ""},
                ("9X9X", "LIG", ""): {"status": "all_copies", "instance": "LIG_AB_1", "note": ""}}
        self.assertEqual(si.instance_chains("6ZIN", "Q6Q", "A:1000", amap, ["Q6Q_AAA_1000"]), {"Q6Q_AAA_1000": "A"})
        self.assertEqual(si.instance_chains("7E2X", "CLR", "", amap, ["CLR_R_602", "CLR_R_603"]),
                         {"CLR_R_602": "R", "CLR_R_603": "R"})
        with self.assertRaises(si.MapMismatch):
            si.instance_chains("9X9X", "LIG", "", amap, ["LIG_AB_1"])

    def test_capped_lines_are_returned(self):
        _, capped = si.standard_ligand_block(self.CAU, "CAU_A_408", "A")
        self.assertEqual(capped, [])
        _, capped = si.standard_ligand_block(self.ABUT, "LIG_A_1", "A")
        self.assertEqual(len(capped), 1)

    def test_fields_too_wide_raise(self):
        wide = (
            ("HETATM    1  C1 LIG A10000     -1.000   2.000   3.000  1.00 20.00           C", "10000", ""),
            ("HETATM    1  C1 LIG A   1   -1000.500   2.000   3.000  1.00 20.00           C", "1", ""),
            ("HETATM    1  C1 LIG A   1    10000.000   2.000   3.000  1.00 20.00           C", "1", ""),
            ("HETATM    1  C1 LIG A   1      1.000   2.000   3.000  1.00-100.00           C", "1", ""),
            ("HETATM    1  C1 LIG A   1      1.000   2.000   3.0001000.00 20.00           C", "1", ""),
        )
        for line, resnum, icode in wide:
            with self.assertRaises(si.MalformedLigandLine, msg=line):
                si.standard_ligand_line(line, "LIG", "A", resnum, icode, "A")
        with self.assertRaises(si.MalformedLigandLine):
            si.standard_ligand_line(self.CAU, "CAU", "A", "408", "", "")


class CascadeGuardTests(unittest.TestCase):

    def test_only_expected_models_may_be_deleted(self):
        si._only_deleted({"structure.Fragment": 3}, {"structure.Fragment"})
        si._only_deleted({"structure.Fragment": 3, "structure.Rotamer": 0}, {"structure.Fragment"})
        with self.assertRaises(si.UnexpectedCascade):
            si._only_deleted({"structure.PdbData": 1, "structure.Rotamer": 2}, {"structure.PdbData"})

    def test_fragment_text(self):
        self.assertEqual(si.fragment_text(["HETATM 1", "HETATM 2"]), "HETATM 1\nHETATM 2\n")
        self.assertEqual(si.fragment_text([]), "")


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

    def test_receptor_chain_is_the_product_chain(self):
        # The receptor chain is now the product (author) chain from the map.
        records, counts, _ = si.plan_rows([row("HPhob", chain="AAA")], "AAA")
        self.assertEqual((len(records), counts["other_chain"]), (1, 0))


class LigandLinesTests(unittest.TestCase):

    def test_collapsed_rows_merge_ligand_atoms_in_first_seen_order(self):
        rows = [row("HPhob", seq=100, lig="HETATM C1\nHETATM C2\n"),
                row("HPhob", seq=100, lig="HETATM C2\nHETATM C3\n"),
                row("Acceptor", "ligand-acceptor", seq=100, atom="OG", lig="HETATM O1\n")]
        records, counts, _ = si.plan_rows(rows, "A")
        self.assertEqual(counts["duplicate"], 1)
        self.assertEqual([r["ligand_lines"] for r in records],
                         [["HETATM C1", "HETATM C2", "HETATM C3"], ["HETATM O1"]])


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
