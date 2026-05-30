"""Tests for the Engine 2 (peptide / protein-protein) interaction consumer (D3).

Run with::

    docker exec gpcrdb-app python manage.py test interaction

This is the THIRD Strategy implementation. Engine 2 consumes the ADR-014
``engine2/1.0`` schema (receptor-residue ↔ peptide-residue interface) and writes
to the contactnetwork ``InteractingPeptideResiduePair`` / ``InteractionPeptide``
models — the residue-residue interface tables protwis already ships and the
peptide API already reads.

Pure-logic tests (mapping, side resolution, locator) use ``SimpleTestCase``
(no DB; raises if a test touches the DB), faithful to the rule that
parse/map logic must be testable without a database. The single DB-integration
test uses ``TestCase`` and builds a minimal Structure / Residue /
LigandPeptideStructure fixture in-memory.

Fixture data: ``interaction/tests/engine2_data/`` holds the D4 pioneer output
for 7F4F (MC1R / Afamelanotide) and 4GRV (NTSR1 / Neurotensin 8-13), in the
real adapter layout ``{PDB}/{PDB}_{recv}_{lig}/{PDB}_{recv}_{lig}.yaml``.
"""

import os

import yaml
from django.test import SimpleTestCase, TestCase

from interaction.schrodinger_processor import (
    _engine2_atoms_by_side,
    _engine2_pair_partners,
    _engine2_partner_side,
    _three_to_one,
    build_residue_pair_metrics,
    locate_engine2_yamls,
    map_engine2_interaction,
    process_schrodinger_peptide_interactions,
    process_schrodinger_protein_interactions,
)

HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE2_DATA_DIR = os.path.join(HERE, "engine2_data")
YAML_7F4F = os.path.join(ENGINE2_DATA_DIR, "7F4F", "7F4F_R_F", "7F4F_R_F.yaml")
YAML_4GRV = os.path.join(ENGINE2_DATA_DIR, "4GRV", "4GRV_A_B", "4GRV_A_B.yaml")


def _load(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def _interactions(path):
    return _load(path)["interface_interactions"]


class Engine2YamlLocationTests(SimpleTestCase):
    """The locator finds the adapter YAML (engine2/1.0), not the worker sidecar."""

    def test_locate_7f4f(self):
        found = locate_engine2_yamls(ENGINE2_DATA_DIR, "7F4F")
        self.assertEqual(len(found), 1, found)
        self.assertTrue(found[0].endswith("7F4F_R_F/7F4F_R_F.yaml"))

    def test_locate_lowercase_normalised(self):
        found = locate_engine2_yamls(ENGINE2_DATA_DIR, "7f4f")
        self.assertEqual(len(found), 1)

    def test_excludes_worker_sidecar(self):
        # Only the adapter YAML (basename == dir name) is returned; a *_worker.yaml
        # sidecar in the same dir must NOT be picked up.
        found = locate_engine2_yamls(ENGINE2_DATA_DIR, "7F4F")
        self.assertFalse(any("_worker.yaml" in p for p in found))

    def test_missing_pdb_returns_empty(self):
        self.assertEqual(locate_engine2_yamls(ENGINE2_DATA_DIR, "ZZZZ"), [])

    def test_real_yaml_is_engine2_schema(self):
        doc = _load(YAML_7F4F)
        self.assertEqual(doc["schema_version"], "engine2/1.0")
        self.assertEqual(doc["metadata"]["ligand_type"], "peptide")
        self.assertEqual(doc["metadata"]["receptor_chain"], "R")
        self.assertEqual(doc["metadata"]["ligand_chain"], "F")


class PartnerSideTests(SimpleTestCase):
    """selection1 = receptor, selection2 = peptide; every partner is tagged."""

    def test_side_mapping(self):
        self.assertEqual(_engine2_partner_side({"side": "selection1"}), "receptor")
        self.assertEqual(_engine2_partner_side({"side": "selection2"}), "peptide")
        self.assertIsNone(_engine2_partner_side({}))

    def test_pair_partners_hbond_donor_receptor(self):
        # 7F4F entry 0: donor HIE260 (R, selection1=receptor) → acceptor TRP9 (F).
        entry = _interactions(YAML_7F4F)[0]
        recv, pep = _engine2_pair_partners(entry)
        self.assertEqual(recv["chain_id"], "R")
        self.assertEqual(recv["resid"], 260)
        self.assertEqual(pep["chain_id"], "F")
        self.assertEqual(pep["resid"], 9)

    def test_pair_partners_hbond_donor_peptide(self):
        # 7F4F entry 1: donor DPN7 (F, selection2=peptide) → acceptor GLU94 (R).
        entry = _interactions(YAML_7F4F)[1]
        recv, pep = _engine2_pair_partners(entry)
        self.assertEqual(recv["chain_id"], "R")  # receptor always resolved to R
        self.assertEqual(recv["resid"], 94)
        self.assertEqual(pep["chain_id"], "F")
        self.assertEqual(pep["resid"], 7)

    def test_every_interface_entry_resolves_to_opposite_sides(self):
        # The interface worker emits only receptor↔peptide pairs; all should
        # resolve (no None) across both real PDBs.
        for path in (YAML_7F4F, YAML_4GRV):
            for entry in _interactions(path):
                recv, pep = _engine2_pair_partners(entry)
                self.assertIsNotNone(recv, entry)
                self.assertIsNotNone(pep, entry)

    def test_untagged_same_side_returns_none(self):
        # A hand-built intra-chain pair (both selection1) must NOT be written.
        entry = {
            "type": "hydrogen_bond",
            "donor": {"side": "selection1", "resid": 1},
            "acceptor": {"side": "selection1", "resid": 2},
        }
        recv, pep = _engine2_pair_partners(entry)
        self.assertIsNone(recv)
        self.assertIsNone(pep)


class AtomsBySideTests(SimpleTestCase):
    """receptor/peptide atom names are assigned by side, order-independent."""

    def test_receptor_first(self):
        a = {"side": "selection1", "atom_name": "NE2"}
        b = {"side": "selection2", "atom_name": "O"}
        recv_atom, pep_atom = _engine2_atoms_by_side(a, b)
        self.assertEqual(recv_atom, "NE2")
        self.assertEqual(pep_atom, "O")

    def test_peptide_first(self):
        a = {"side": "selection2", "atom_name": "N"}
        b = {"side": "selection1", "atom_name": "OE1"}
        recv_atom, pep_atom = _engine2_atoms_by_side(a, b)
        self.assertEqual(recv_atom, "OE1")
        self.assertEqual(pep_atom, "N")

    def test_missing_atom_name_empty(self):
        a = {"side": "selection1"}
        b = {"side": "selection2"}
        recv_atom, pep_atom = _engine2_atoms_by_side(a, b)
        self.assertEqual(recv_atom, "")
        self.assertEqual(pep_atom, "")


class InteractionMappingTests(SimpleTestCase):
    """Engine 2 type → legacy (interaction_type, specific_type, level)."""

    def test_hydrogen_bond_carries_hbond_class(self):
        entry = _interactions(YAML_7F4F)[0]  # hb_sb
        m = map_engine2_interaction(entry)
        self.assertEqual(m["interaction_type"], "polar")
        self.assertIn("h-bond", m["specific_type"])
        self.assertIn("hb_sb", m["specific_type"])  # 4-class preserved (ADR-014)
        self.assertEqual(m["interaction_level"], 0)

    def test_salt_bridge_is_ionic(self):
        salt = next(e for e in _interactions(YAML_7F4F) if e["type"] == "salt_bridge")
        m = map_engine2_interaction(salt)
        self.assertEqual(m["interaction_type"], "ionic")
        self.assertIn("salt-bridge", m["specific_type"])

    def test_pi_stacking_carries_subtype(self):
        pi = next(e for e in _interactions(YAML_7F4F) if e["type"] == "pi_pi_stacking")
        m = map_engine2_interaction(pi)
        self.assertEqual(m["interaction_type"], "aromatic")
        self.assertIn("Face-to-Face", m["specific_type"])

    def test_pi_cation_is_aromatic(self):
        pc = next(e for e in _interactions(YAML_7F4F) if e["type"] == "pi_cation")
        m = map_engine2_interaction(pc)
        self.assertEqual(m["interaction_type"], "aromatic")
        self.assertIn("pi-cation", m["specific_type"])

    def test_hydrophobic_is_hydrophobic(self):
        h = next(e for e in _interactions(YAML_7F4F) if e["type"] == "hydrophobic_contact")
        m = map_engine2_interaction(h)
        self.assertEqual(m["interaction_type"], "hydrophobic")

    def test_steric_clash_is_mapped(self):
        sc = next(e for e in _interactions(YAML_7F4F) if e["type"] == "steric_clash")
        m = map_engine2_interaction(sc)
        self.assertEqual(m["interaction_type"], "steric-clash")

    def test_unknown_type_fails_loud(self):
        # ADR-009: a new Engine 2 interaction class must surface, not drop.
        with self.assertRaises(ValueError):
            map_engine2_interaction({"type": "halogen_bond"})

    def test_every_real_entry_maps_without_error(self):
        for path in (YAML_7F4F, YAML_4GRV):
            for entry in _interactions(path):
                m = map_engine2_interaction(entry)
                self.assertIn(
                    m["interaction_type"],
                    {"polar", "ionic", "aromatic", "hydrophobic", "steric-clash"},
                )

    def test_hbond_atom_names_resolved(self):
        # 7F4F entry 0: receptor HIE260 NE2 (donor) ↔ peptide TRP9 O (acceptor).
        m = map_engine2_interaction(_interactions(YAML_7F4F)[0])
        self.assertEqual(m["receptor_atom"], "NE2")
        self.assertEqual(m["peptide_atom"], "O")


class ThreeToOneTests(SimpleTestCase):
    """Peptide residues include unnatural / protonation-state aliases."""

    def test_standard(self):
        self.assertEqual(_three_to_one("ARG"), "R")
        self.assertEqual(_three_to_one("trp"), "W")

    def test_protonation_aliases(self):
        # HIE/HID/HIP are His protonation states from PrepWizard.
        self.assertEqual(_three_to_one("HIE"), "H")
        self.assertEqual(_three_to_one("ASH"), "D")

    def test_unknown_falls_back_to_x(self):
        # DPN / NLE are real unnatural residues in 7F4F/4GRV; mapping the full
        # unnatural_amino_acids.yaml is out of D3 scope (flagged in Inbox).
        self.assertEqual(_three_to_one("ZZZ"), "X")

    def test_none_safe(self):
        self.assertEqual(_three_to_one(None), "X")


class SummaryConsistencyTests(SimpleTestCase):
    """The per-PDB summary counts match the instance YAML (cross-check fixture)."""

    def test_7f4f_counts_match_summary(self):
        interactions = _interactions(YAML_7F4F)
        type_counts = {}
        for e in interactions:
            type_counts[e["type"]] = type_counts.get(e["type"], 0) + 1
        summary = _load(os.path.join(ENGINE2_DATA_DIR, "7F4F", "summary.yaml"))
        expected = summary["interface_summary"][0]["counts"]["interaction_type_counts"]
        for itype, n in expected.items():
            self.assertEqual(type_counts.get(itype, 0), n, itype)


class ResiduePairMetricsTests(SimpleTestCase):
    """ADR-015: per-residue-pair buried-SASA / Sc indexed from
    residue_pair_summaries, resolved by side (set_1=receptor, set_2=peptide)."""

    def test_metrics_indexed_by_recv_pep_seq(self):
        summaries = _load(YAML_7F4F)["residue_pair_summaries"]
        metrics = build_residue_pair_metrics(summaries)
        # Pair (recv ASP42, pep HIE6): residue1=selection1=receptor.
        self.assertIn((42, 6), metrics)
        m = metrics[(42, 6)]
        self.assertAlmostEqual(m["buried_sasa_receptor"], 0.23914552299555836)
        self.assertAlmostEqual(m["buried_sasa_peptide"], 0.9303769604153024)
        self.assertAlmostEqual(m["surface_complementarity"], 0.6882730654191149)

    def test_every_summary_pair_has_three_metrics(self):
        for path in (YAML_7F4F, YAML_4GRV):
            summaries = _load(path)["residue_pair_summaries"]
            metrics = build_residue_pair_metrics(summaries)
            self.assertEqual(len(metrics), len(summaries))
            for m in metrics.values():
                self.assertIn("buried_sasa_receptor", m)
                self.assertIn("buried_sasa_peptide", m)
                self.assertIn("surface_complementarity", m)

    def test_side_flip_preserves_receptor_peptide_assignment(self):
        # Hand-built entry where residue1 is the PEPTIDE (selection2): set_1
        # must then be read as the peptide-side SASA, not the receptor's.
        summaries = [{
            "residue1": {"side": "selection2", "resid": 9},  # peptide
            "residue2": {"side": "selection1", "resid": 260},  # receptor
            "properties": {
                "set_1_buried_sasa": 0.11,  # peptide (residue1 = selection2)
                "set_2_buried_sasa": 0.99,  # receptor
                "surface_complementarity": 0.5,
            },
        }]
        metrics = build_residue_pair_metrics(summaries)
        self.assertIn((260, 9), metrics)  # keyed (recv, pep)
        m = metrics[(260, 9)]
        self.assertAlmostEqual(m["buried_sasa_receptor"], 0.99)
        self.assertAlmostEqual(m["buried_sasa_peptide"], 0.11)

    def test_missing_summaries_returns_empty(self):
        self.assertEqual(build_residue_pair_metrics(None), {})
        self.assertEqual(build_residue_pair_metrics([]), {})

    def test_intra_chain_pair_skipped(self):
        # both selection1 — not an interface pair, must be dropped.
        summaries = [{
            "residue1": {"side": "selection1", "resid": 1},
            "residue2": {"side": "selection1", "resid": 2},
            "properties": {"set_1_buried_sasa": 0.5, "set_2_buried_sasa": 0.5,
                           "surface_complementarity": 0.5},
        }]
        self.assertEqual(build_residue_pair_metrics(summaries), {})


class Engine2DbIntegrationTests(TestCase):
    """End-to-end: the processor writes the contactnetwork peptide models.

    Builds a minimal Structure / Residue / LigandPeptideStructure graph for
    7F4F (receptor chain R / peptide chain F), runs the consumer against the
    frozen engine2_data fixture, and checks the persisted
    InteractingPeptideResiduePair / InteractionPeptide rows + idempotency.
    """

    @classmethod
    def setUpTestData(cls):
        from common.models import WebLink, WebResource
        from ligand.models import Ligand, LigandPeptideStructure
        from protein.models import (
            Protein,
            ProteinConformation,
            ProteinFamily,
            ProteinSequenceType,
            ProteinSource,
            ProteinState,
            Species,
        )
        from residue.models import Residue
        from structure.models import Structure, StructureType

        fam = ProteinFamily.objects.create(slug="000_001_001_001", name="MC1R fam")
        species = Species.objects.create(latin_name="Homo sapiens", common_name="Human")
        source = ProteinSource.objects.create(name="SWISSPROT")
        seqtype = ProteinSequenceType.objects.create(slug="wt", name="Wild-type")
        state = ProteinState.objects.create(slug="active", name="Active")

        protein = Protein.objects.create(
            family=fam, species=species, source=source, sequence_type=seqtype,
            entry_name="mc1r_human_test", name="MC1R", sequence="M",
        )
        cls.pconf = ProteinConformation.objects.create(protein=protein, state=state)

        stype = StructureType.objects.create(slug="x-ray", name="X-ray")
        wr = WebResource.objects.create(slug="pdb", name="PDB", url="https://pdb/$index")
        weblink = WebLink.objects.create(index="7F4F", web_resource=wr)
        cls.structure = Structure.objects.create(
            protein_conformation=cls.pconf, structure_type=stype, state=state,
            pdb_code=weblink, preferred_chain="R", publication_date="2021-01-01",
        )

        # Receptor residues that the 7F4F fixture references on chain R.
        # (subset is enough — entries on absent residues are skipped with a warning)
        cls.recv_seqs = {
            42: "D", 45: "F", 94: "E", 117: "D", 118: "N", 121: "D",
            183: "Y", 260: "H",
        }
        for seq, aa in cls.recv_seqs.items():
            Residue.objects.create(
                protein_conformation=cls.pconf, sequence_number=seq, amino_acid=aa,
            )

        lig = Ligand.objects.create(name="Afamelanotide")
        cls.lps = LigandPeptideStructure.objects.create(
            structure=cls.structure, ligand=lig, chain="F",
        )

    def test_writes_peptide_pairs_and_interactions(self):
        from contactnetwork.models import (
            InteractingPeptideResiduePair,
            InteractionPeptide,
        )

        ok = process_schrodinger_peptide_interactions(
            current_structure_obj=self.structure,
            ligand_chain="F",
            pdb_code_str="7F4F",
            schrodinger_interactions_dir_override=ENGINE2_DATA_DIR,
        )
        self.assertTrue(ok)

        pairs = InteractingPeptideResiduePair.objects.filter(peptide=self.lps)
        self.assertGreater(pairs.count(), 0)
        # Every pair's receptor residue is on the receptor conformation.
        for p in pairs:
            self.assertEqual(p.receptor_residue.protein_conformation_id, self.pconf.id)

        rows = InteractionPeptide.objects.filter(
            interacting_peptide_pair__peptide=self.lps
        )
        self.assertGreater(rows.count(), 0)
        # The legacy vocabulary the API understands is produced.
        types = set(rows.values_list("interaction_type", flat=True))
        self.assertTrue(types <= {"polar", "ionic", "aromatic", "hydrophobic", "steric-clash"})
        # 4-class hbond signal preserved in specific_type (ADR-014).
        polar_specifics = set(
            rows.filter(interaction_type="polar").values_list("specific_type", flat=True)
        )
        self.assertTrue(any("hb_" in s for s in polar_specifics), polar_specifics)

    def test_interface_metrics_persisted_on_pair(self):
        # ADR-015: buried-SASA / Sc written onto InteractingPeptideResiduePair.
        from contactnetwork.models import InteractingPeptideResiduePair

        process_schrodinger_peptide_interactions(
            current_structure_obj=self.structure,
            ligand_chain="F",
            pdb_code_str="7F4F",
            schrodinger_interactions_dir_override=ENGINE2_DATA_DIR,
        )
        # Pair receptor ASP42 ↔ peptide HIE6 is in both interface_interactions
        # and residue_pair_summaries, and residue 42 is in the fixture.
        pair = InteractingPeptideResiduePair.objects.get(
            peptide=self.lps,
            receptor_residue__sequence_number=42,
            peptide_sequence_number=6,
        )
        self.assertAlmostEqual(pair.buried_sasa_receptor, 0.23914552299555836, places=6)
        self.assertAlmostEqual(pair.buried_sasa_peptide, 0.9303769604153024, places=6)
        self.assertAlmostEqual(pair.surface_complementarity, 0.6882730654191149, places=6)

        # Every written pair that the summaries listed has non-NULL metrics
        # (all 7F4F interface pairs appear in residue_pair_summaries).
        written = InteractingPeptideResiduePair.objects.filter(peptide=self.lps)
        self.assertGreater(written.count(), 0)
        for p in written:
            self.assertIsNotNone(p.surface_complementarity, p)
            self.assertIsNotNone(p.buried_sasa_receptor, p)
            self.assertIsNotNone(p.buried_sasa_peptide, p)

    def test_idempotent_rerun(self):
        from contactnetwork.models import InteractionPeptide

        kwargs = dict(
            current_structure_obj=self.structure, ligand_chain="F",
            pdb_code_str="7F4F", schrodinger_interactions_dir_override=ENGINE2_DATA_DIR,
        )
        process_schrodinger_peptide_interactions(**kwargs)
        first = InteractionPeptide.objects.filter(
            interacting_peptide_pair__peptide=self.lps).count()
        process_schrodinger_peptide_interactions(**kwargs)
        second = InteractionPeptide.objects.filter(
            interacting_peptide_pair__peptide=self.lps).count()
        self.assertEqual(first, second)  # delete-then-insert, no accumulation

    def test_missing_lps_returns_false(self):
        ok = process_schrodinger_peptide_interactions(
            current_structure_obj=self.structure,
            ligand_chain="Z",  # no LigandPeptideStructure for chain Z
            pdb_code_str="7F4F",
            schrodinger_interactions_dir_override=ENGINE2_DATA_DIR,
        )
        self.assertFalse(ok)

    def test_missing_yaml_returns_false(self):
        ok = process_schrodinger_peptide_interactions(
            current_structure_obj=self.structure,
            ligand_chain="F",
            pdb_code_str="9ZZZ",  # no YAML for this PDB
            schrodinger_interactions_dir_override=ENGINE2_DATA_DIR,
        )
        self.assertFalse(ok)

    def test_protein_path_aliases_peptide_path(self):
        # process_schrodinger_protein_interactions delegates to the peptide path.
        from contactnetwork.models import InteractionPeptide

        ok = process_schrodinger_protein_interactions(
            current_structure_obj=self.structure,
            ligand_chain="F",
            pdb_code_str="7F4F",
            schrodinger_interactions_dir_override=ENGINE2_DATA_DIR,
        )
        self.assertTrue(ok)
        self.assertGreater(
            InteractionPeptide.objects.filter(
                interacting_peptide_pair__peptide=self.lps).count(),
            0,
        )
