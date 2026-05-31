"""Tests for the Engine 2 receptor / G-alpha interface consumer.

Run with::

    docker exec gpcrdb-app python manage.py test interaction

Unlike the peptide path (InteractingPeptideResiduePair, partner = bare seqnum +
three-letter name), the G-protein interface resolves BOTH residues to real
protwis Residue rows and writes the contactnetwork InteractingResiduePair /
Interaction models -- the same target the legacy do_complexes path (cube.py)
writes and that signprot/views.py already consumes.

The bridge under test: the G-alpha residue lives in the structure's "{pdb}_a"
ProteinConformation keyed by PDB author sequence number (both sides store author
numbering -- the WT-replace path was disabled in 2019), reached via
SignprotComplex.alpha. An AA-guard on both sides catches any register offset
(surface, do not silently corrupt).

Fixture data: interaction/tests/engine2_data/ holds the real adapter output for
6CMO (rho / Gi, alpha='A', author approx WT) and 7RYC (OXTR / Gq, alpha='D',
author 1005-1246 high-offset construct -- exercises chain selection + offset).
7RYC also carries a legacy peptide instance (chain L, oxytocin) so the
partner_category routing guard is exercised. The integration fixtures build the
receptor + "_a" conformations IN-MEMORY from the same YAML, so the AA-guard
passes by construction and pair/interaction counts are deterministic.
"""

import os

import yaml
from django.test import SimpleTestCase, TestCase

from interaction.schrodinger_processor import (
    _engine2_pair_partners,
    process_schrodinger_gprotein_interactions,
)

HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE2_DATA_DIR = os.path.join(HERE, "engine2_data")
YAML_6CMO = os.path.join(ENGINE2_DATA_DIR, "6CMO", "6CMO_R_A", "6CMO_R_A.yaml")
YAML_7RYC_GP = os.path.join(ENGINE2_DATA_DIR, "7RYC", "7RYC_O_D", "7RYC_O_D.yaml")
YAML_7RYC_PEP = os.path.join(ENGINE2_DATA_DIR, "7RYC", "7RYC_O_L", "7RYC_O_L.yaml")

# Three-letter to one-letter incl PrepWizard protonation aliases, used only to
# build the in-memory fixture residues so the AA-guard passes by construction.
_T3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "HIE": "H", "HID": "H", "HIP": "H", "ASH": "D",
    "GLH": "E", "LYN": "K", "CYX": "C", "ARN": "R",
}


def _load(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def _pairs_and_residues(doc):
    """From interface_interactions return (recv {seq: resname}, gp {seq:
    resname}, pairs [(recv_seq, gp_seq), ...]) using the same side resolver the
    processor uses."""
    recv, gp, pairs = {}, {}, []
    for entry in doc.get("interface_interactions", []):
        r, p = _engine2_pair_partners(entry)
        if r is None or p is None:
            continue
        recv[int(r["resid"])] = (r.get("resname") or "").upper()
        gp[int(p["resid"])] = (p.get("resname") or "").upper()
        pairs.append((int(r["resid"]), int(p["resid"])))
    return recv, gp, pairs


# ---------------------------------------------------------------------------
# Pure-fixture sanity (no DB)
# ---------------------------------------------------------------------------
class GproteinFixtureShapeTests(SimpleTestCase):
    def test_6cmo_is_galpha_instance(self):
        meta = _load(YAML_6CMO)["metadata"]
        self.assertEqual(meta["partner_category"], "g_protein_alpha")
        self.assertEqual(meta["ligand_chain"], "A")
        self.assertEqual(meta["partner_uniprot"], "gnai1_human")

    def test_7ryc_has_both_galpha_and_legacy_peptide(self):
        # The routing hinge: same PDB, two instances; route on partner_category.
        gp = _load(YAML_7RYC_GP)["metadata"]
        pep = _load(YAML_7RYC_PEP)["metadata"]
        self.assertEqual(gp["partner_category"], "g_protein_alpha")
        self.assertEqual(gp["ligand_chain"], "D")  # not 'A'
        self.assertEqual(pep.get("partner_category") or "", "")
        self.assertEqual(pep["ligand_type"], "peptide")
        self.assertEqual(pep["ligand_chain"], "L")

    def test_6cmo_pair_and_interaction_counts(self):
        recv, gp, pairs = _pairs_and_residues(_load(YAML_6CMO))
        self.assertEqual(len(recv), 17)
        self.assertEqual(len(gp), 14)
        self.assertEqual(len(pairs), 49)          # one Interaction row per entry
        self.assertEqual(len(set(pairs)), 25)     # distinct InteractingResiduePair


# ---------------------------------------------------------------------------
# DB integration -- the in-memory fixture builder
# ---------------------------------------------------------------------------
def _build_complex_fixture(cls, pdb, alpha_chain, recv_entry, alpha_entry, yaml_path):
    """Build a minimal receptor + "_a" graph for one G-protein complex from a
    real fixture YAML, plus its SignprotComplex. Attaches to cls:
    structure, receptor_pc, alpha_pc, and the (recv, gp, pairs) decomposition."""
    from common.models import WebLink, WebResource
    from protein.models import (
        Protein, ProteinConformation, ProteinFamily, ProteinSequenceType,
        ProteinSource, ProteinState, Species,
    )
    from residue.models import Residue, ResidueGenericNumber, ResidueNumberingScheme
    from signprot.models import SignprotComplex
    from structure.models import Structure, StructureType

    fam = ProteinFamily.objects.create(slug=f"000_001_001_{pdb}", name=f"{pdb} fam")
    species = Species.objects.create(latin_name=f"Sp {pdb}", common_name=pdb)
    source = ProteinSource.objects.create(name=f"SRC_{pdb}")
    seqtype = ProteinSequenceType.objects.create(slug=f"wt_{pdb}", name="WT")
    state = ProteinState.objects.create(slug=f"active_{pdb}", name="Active")

    receptor = Protein.objects.create(
        family=fam, species=species, source=source, sequence_type=seqtype,
        entry_name=recv_entry, name=f"{pdb} receptor", sequence="M",
    )
    cls.receptor_pc = ProteinConformation.objects.create(protein=receptor, state=state)

    # The G-alpha WT protein + its "{pdb}_a" structure-bound conformation.
    galpha = Protein.objects.create(
        family=fam, species=species, source=source, sequence_type=seqtype,
        entry_name=alpha_entry, name=f"{pdb} Galpha", sequence="M",
    )
    cls.alpha_pc = ProteinConformation.objects.create(protein=galpha, state=state)

    stype = StructureType.objects.create(slug=f"x-ray_{pdb}", name="X-ray")
    wr = WebResource.objects.create(slug=f"pdb_{pdb}", name="PDB", url="https://pdb/$index")
    weblink = WebLink.objects.create(index=pdb, web_resource=wr)
    cls.structure = Structure.objects.create(
        protein_conformation=cls.receptor_pc, structure_type=stype, state=state,
        pdb_code=weblink, preferred_chain="R", publication_date="2021-01-01",
    )
    SignprotComplex.objects.create(
        protein=galpha, structure=cls.structure, alpha=alpha_chain,
    )

    recv, gp, pairs = _pairs_and_residues(_load(yaml_path))
    cls.recv, cls.gp, cls.pairs = recv, gp, pairs

    for seq, resname in recv.items():
        Residue.objects.create(
            protein_conformation=cls.receptor_pc, sequence_number=seq,
            amino_acid=_T3[resname],
        )

    # G-alpha residues on "_a", each with a generic number (display_generic_number)
    # to prove the consumer linkage comes free (the real "_a" build is fully
    # populated).
    cgn_scheme = ResidueNumberingScheme.objects.create(
        slug=f"cgn_{pdb}", short_name="CGN", name="Common Gprot Numbering",
    )
    for seq, resname in gp.items():
        rgn = ResidueGenericNumber.objects.create(
            scheme=cgn_scheme, label=f"G.X.{seq}",
        )
        Residue.objects.create(
            protein_conformation=cls.alpha_pc, sequence_number=seq,
            amino_acid=_T3[resname], display_generic_number=rgn,
        )


class Engine2Gprotein6CMOTests(TestCase):
    """6CMO (rho / Gi, alpha='A'): the clean, author-approx-WT base case."""

    @classmethod
    def setUpTestData(cls):
        _build_complex_fixture(cls, "6CMO", "A", "rho_6cmo_recv", "6cmo_a", YAML_6CMO)

    def _run(self):
        return process_schrodinger_gprotein_interactions(
            current_structure_obj=self.structure, ligand_chain="A",
            pdb_code_str="6CMO", schrodinger_interactions_dir_override=ENGINE2_DATA_DIR,
        )

    def test_writes_interacting_residue_pairs_not_peptide_model(self):
        from contactnetwork.models import (
            Interaction, InteractingPeptideResiduePair, InteractingResiduePair,
        )

        self.assertTrue(self._run())

        pairs = InteractingResiduePair.objects.filter(referenced_structure=self.structure)
        self.assertEqual(pairs.count(), 25)                       # distinct pairs
        self.assertEqual(
            Interaction.objects.filter(interacting_pair__in=pairs).count(), 49
        )
        # Target model is the residue-residue model -- NOT the peptide model.
        self.assertEqual(InteractingPeptideResiduePair.objects.count(), 0)

    def test_res1_receptor_res2_on_alpha_conformation_with_cgn(self):
        from contactnetwork.models import InteractingResiduePair

        self._run()
        for pair in InteractingResiduePair.objects.filter(referenced_structure=self.structure):
            # res1 = receptor side, res2 = G-alpha side (the bridge target).
            self.assertEqual(pair.res1.protein_conformation_id, self.receptor_pc.id)
            self.assertEqual(pair.res2.protein_conformation_id, self.alpha_pc.id)
            # The generic number comes free off the "_a" residue (consumer reads it).
            self.assertIsNotNone(pair.res2.display_generic_number)

    def test_legacy_vocabulary_and_hbond_class_preserved(self):
        from contactnetwork.models import Interaction, InteractingResiduePair

        self._run()
        pairs = InteractingResiduePair.objects.filter(referenced_structure=self.structure)
        rows = Interaction.objects.filter(interacting_pair__in=pairs)
        types = set(rows.values_list("interaction_type", flat=True))
        self.assertTrue(types <= {"polar", "ionic", "aromatic", "hydrophobic", "steric-clash"})
        polar = set(
            rows.filter(interaction_type="polar").values_list("specific_type", flat=True)
        )
        self.assertTrue(any("hb_" in s for s in polar), polar)  # 4-class preserved

    def test_atom_names_on_correct_side(self):
        # res1=receptor -> atomname_residue1; res2=G-alpha -> atomname_residue2.
        from contactnetwork.models import Interaction, InteractingResiduePair

        self._run()
        # 6CMO entry 0: receptor ARG147 NH1 (donor) / G-alpha ALA31 O (acceptor).
        pair = InteractingResiduePair.objects.get(
            referenced_structure=self.structure,
            res1__sequence_number=147, res2__sequence_number=31,
        )
        row = Interaction.objects.filter(
            interacting_pair=pair, interaction_type="polar"
        ).first()
        self.assertEqual(row.atomname_residue1, "NH1")
        self.assertEqual(row.atomname_residue2, "O")

    def test_idempotent_rerun(self):
        from contactnetwork.models import Interaction, InteractingResiduePair

        self._run()
        pairs1 = InteractingResiduePair.objects.filter(referenced_structure=self.structure).count()
        inter1 = Interaction.objects.filter(
            interacting_pair__referenced_structure=self.structure
        ).count()
        self._run()
        pairs2 = InteractingResiduePair.objects.filter(referenced_structure=self.structure).count()
        inter2 = Interaction.objects.filter(
            interacting_pair__referenced_structure=self.structure
        ).count()
        self.assertEqual((pairs1, inter1), (25, 49))
        self.assertEqual((pairs2, inter2), (25, 49))  # no accumulation

    def test_wrong_chain_returns_false(self):
        # ligand_chain 'L' has no g_protein_alpha YAML for 6CMO -> False.
        self.assertFalse(
            process_schrodinger_gprotein_interactions(
                current_structure_obj=self.structure, ligand_chain="L",
                pdb_code_str="6CMO", schrodinger_interactions_dir_override=ENGINE2_DATA_DIR,
            )
        )


class Engine2Gprotein7RYCTests(TestCase):
    """7RYC (OXTR / Gq, alpha='D'): high-offset construct (author 1005-1246) +
    co-resident legacy peptide instance (chain L). Exercises (a) chain selection
    off SignprotComplex.alpha='D' (not 'A'), (b) the offset join, (c) the
    partner_category routing guard rejecting the peptide instance."""

    @classmethod
    def setUpTestData(cls):
        _build_complex_fixture(cls, "7RYC", "D", "oxtr_7ryc_recv", "7ryc_a", YAML_7RYC_GP)

    def test_high_offset_galpha_joins_on_author_seqnum(self):
        from contactnetwork.models import Interaction, InteractingResiduePair

        ok = process_schrodinger_gprotein_interactions(
            current_structure_obj=self.structure, ligand_chain="D",
            pdb_code_str="7RYC", schrodinger_interactions_dir_override=ENGINE2_DATA_DIR,
        )
        self.assertTrue(ok)
        pairs = InteractingResiduePair.objects.filter(referenced_structure=self.structure)
        self.assertEqual(pairs.count(), 20)
        self.assertEqual(
            Interaction.objects.filter(interacting_pair__in=pairs).count(), 43
        )
        # res2 seqnums are in the high-offset construct range (1005-1246).
        gp_seqs = set(pairs.values_list("res2__sequence_number", flat=True))
        self.assertTrue(all(1000 < s < 1300 for s in gp_seqs), gp_seqs)
        self.assertTrue(gp_seqs <= set(self.gp.keys()))

    def test_legacy_peptide_instance_not_routed_to_gprotein(self):
        # The chain-L oxytocin instance (partner_category='') must be rejected by
        # the G-protein path (it belongs to the peptide path).
        self.assertFalse(
            process_schrodinger_gprotein_interactions(
                current_structure_obj=self.structure, ligand_chain="L",
                pdb_code_str="7RYC", schrodinger_interactions_dir_override=ENGINE2_DATA_DIR,
            )
        )


# ---------------------------------------------------------------------------
# AA-guard + bridge-failure edges (synthetic YAML)
# ---------------------------------------------------------------------------
class Engine2GproteinAaGuardTests(TestCase):
    """Both sides of the interface key on author seqnum, so a register offset
    must be caught. The DB residue's amino acid is cross-checked against the YAML
    resname; on disagreement (or missing residue, or an insertion code that
    cannot be resolved) the entry is SKIPPED with a warning, never written as a
    mis-registered row."""

    PDB = "9GPA"  # synthetic

    @classmethod
    def setUpTestData(cls):
        from common.models import WebLink, WebResource
        from protein.models import (
            Protein, ProteinConformation, ProteinFamily, ProteinSequenceType,
            ProteinSource, ProteinState, Species,
        )
        from residue.models import Residue
        from signprot.models import SignprotComplex
        from structure.models import Structure, StructureType

        fam = ProteinFamily.objects.create(slug="000_001_001_GPA", name="Guard fam")
        species = Species.objects.create(latin_name="Sp guard", common_name="Guard")
        source = ProteinSource.objects.create(name="SRC_GPA")
        seqtype = ProteinSequenceType.objects.create(slug="wt_gpa", name="WT")
        state = ProteinState.objects.create(slug="active_gpa", name="Active")

        receptor = Protein.objects.create(
            family=fam, species=species, source=source, sequence_type=seqtype,
            entry_name="guard_recv", name="GUARD", sequence="M",
        )
        cls.receptor_pc = ProteinConformation.objects.create(protein=receptor, state=state)
        galpha = Protein.objects.create(
            family=fam, species=species, source=source, sequence_type=seqtype,
            entry_name="9gpa_a", name="GUARD Galpha", sequence="M",
        )
        cls.alpha_pc = ProteinConformation.objects.create(protein=galpha, state=state)

        stype = StructureType.objects.create(slug="x-ray_gpa", name="X-ray")
        wr = WebResource.objects.create(slug="pdb_gpa", name="PDB", url="https://pdb/$index")
        weblink = WebLink.objects.create(index=cls.PDB, web_resource=wr)
        cls.structure = Structure.objects.create(
            protein_conformation=cls.receptor_pc, structure_type=stype, state=state,
            pdb_code=weblink, preferred_chain="R", publication_date="2021-01-01",
        )
        SignprotComplex.objects.create(protein=galpha, structure=cls.structure, alpha="A")

        # Receptor residue 100 -> D. G-alpha residues: 200 -> E (DB), 300 -> H (DB).
        Residue.objects.create(protein_conformation=cls.receptor_pc, sequence_number=100, amino_acid="D")
        Residue.objects.create(protein_conformation=cls.alpha_pc, sequence_number=200, amino_acid="E")
        Residue.objects.create(protein_conformation=cls.alpha_pc, sequence_number=300, amino_acid="H")
        # NB: no G-alpha residue 999 (not-found case).

    def _hbond(self, recv_seq, recv_rn, gp_seq, gp_rn, gp_inscode=""):
        return {
            "type": "hydrogen_bond", "hbond_class": "hb_ss",
            "donor": {"chain_id": "R", "resname": recv_rn, "resid": recv_seq,
                      "inscode": "", "atom_name": "N", "side": "selection1"},
            "acceptor": {"chain_id": "A", "resname": gp_rn, "resid": gp_seq,
                         "inscode": gp_inscode, "atom_name": "O", "side": "selection2"},
            "distance_h_acceptor": 1.9, "angle_donor_h_acceptor": 160.0,
        }

    def _write_yaml(self, tmpdir, entries):
        doc = {
            "schema_version": "engine2/1.0",
            "metadata": {
                "ligand_type": "protein", "receptor_chain": "R", "ligand_chain": "A",
                "partner_category": "g_protein_alpha", "partner_uniprot": "test_galpha",
            },
            "interface_interactions": entries,
            "residue_pair_summaries": [],
        }
        inst = os.path.join(tmpdir, self.PDB, f"{self.PDB}_R_A")
        os.makedirs(inst)
        path = os.path.join(inst, f"{self.PDB}_R_A.yaml")
        with open(path, "w") as fh:
            yaml.safe_dump(doc, fh)
        return tmpdir

    def _run(self, base):
        return process_schrodinger_gprotein_interactions(
            current_structure_obj=self.structure, ligand_chain="A",
            pdb_code_str=self.PDB, schrodinger_interactions_dir_override=base,
        )

    def test_galpha_aa_mismatch_skipped(self):
        import tempfile
        from contactnetwork.models import InteractingResiduePair

        with tempfile.TemporaryDirectory() as tmp:
            base = self._write_yaml(tmp, [
                self._hbond(100, "ASP", 200, "GLU"),   # match  (DB E / GLU)
                self._hbond(100, "ASP", 200, "LYS"),   # mismatch on G-alpha (DB E / LYS)
            ])
            with self.assertLogs("interaction.schrodinger_processor", level="WARNING") as cm:
                self.assertTrue(self._run(base))

        pairs = InteractingResiduePair.objects.filter(referenced_structure=self.structure)
        gp_seqs = set(pairs.values_list("res2__sequence_number", flat=True))
        self.assertEqual(gp_seqs, {200})            # only the matching pair written
        # The mismatch row is dropped (1 distinct pair, 1 interaction).
        self.assertEqual(pairs.count(), 1)
        self.assertTrue(any("G-alpha AA mismatch" in m and "seq 200" in m for m in cm.output), cm.output)

    def test_galpha_not_found_skipped(self):
        import tempfile
        from contactnetwork.models import InteractingResiduePair

        with tempfile.TemporaryDirectory() as tmp:
            base = self._write_yaml(tmp, [self._hbond(100, "ASP", 999, "GLY")])
            with self.assertLogs("interaction.schrodinger_processor", level="WARNING") as cm:
                self.assertTrue(self._run(base))
        self.assertEqual(
            InteractingResiduePair.objects.filter(referenced_structure=self.structure).count(), 0
        )
        self.assertTrue(any("G-alpha Residue not found" in m and "seq 999" in m for m in cm.output), cm.output)

    def test_receptor_aa_mismatch_skipped(self):
        import tempfile
        from contactnetwork.models import InteractingResiduePair

        with tempfile.TemporaryDirectory() as tmp:
            # receptor 100 DB=D, YAML says LYS -> mismatch on the receptor side.
            base = self._write_yaml(tmp, [self._hbond(100, "LYS", 200, "GLU")])
            with self.assertLogs("interaction.schrodinger_processor", level="WARNING") as cm:
                self.assertTrue(self._run(base))
        self.assertEqual(
            InteractingResiduePair.objects.filter(referenced_structure=self.structure).count(), 0
        )
        self.assertTrue(any("Receptor AA mismatch" in m and "seq 100" in m for m in cm.output), cm.output)

    def test_protonation_alias_not_a_false_mismatch(self):
        # G-alpha residue 300 DB=H; YAML resname HIE (a His protonation state)
        # must NOT be flagged -- _three_to_one normalises it to H.
        import tempfile
        from contactnetwork.models import InteractingResiduePair

        with tempfile.TemporaryDirectory() as tmp:
            base = self._write_yaml(tmp, [self._hbond(100, "ASP", 300, "HIE")])
            self.assertTrue(self._run(base))
        self.assertTrue(
            InteractingResiduePair.objects.filter(
                referenced_structure=self.structure, res2__sequence_number=300
            ).exists()
        )

    def test_nonempty_insertion_code_skipped(self):
        # The Residue model has no insertion-code field, so a G-alpha partner
        # carrying inscode='A' cannot be resolved unambiguously -- it must skip
        # rather than silently collide with residue 200.
        import tempfile
        from contactnetwork.models import InteractingResiduePair

        with tempfile.TemporaryDirectory() as tmp:
            base = self._write_yaml(tmp, [self._hbond(100, "ASP", 200, "GLU", gp_inscode="A")])
            with self.assertLogs("interaction.schrodinger_processor", level="WARNING") as cm:
                self.assertTrue(self._run(base))
        self.assertEqual(
            InteractingResiduePair.objects.filter(referenced_structure=self.structure).count(), 0
        )
        self.assertTrue(any("insertion code" in m and "seq 200" in m for m in cm.output), cm.output)

    def test_partial_update_replaces_pair_children(self):
        # Per-pair idempotency on a CHANGED re-run: pair (100,200) has 2
        # interactions on run 1, 1 on run 2 -- the stale child must be dropped,
        # not accumulated.
        import tempfile
        from contactnetwork.models import Interaction, InteractingResiduePair

        with tempfile.TemporaryDirectory() as tmp:
            base = self._write_yaml(tmp, [
                self._hbond(100, "ASP", 200, "GLU"),
                self._hbond(100, "ASP", 200, "GLU"),
            ])
            self.assertTrue(self._run(base))
        pair = InteractingResiduePair.objects.get(
            referenced_structure=self.structure, res2__sequence_number=200)
        self.assertEqual(Interaction.objects.filter(interacting_pair=pair).count(), 2)

        with tempfile.TemporaryDirectory() as tmp2:
            base2 = self._write_yaml(tmp2, [self._hbond(100, "ASP", 200, "GLU")])
            self.assertTrue(self._run(base2))
        pair = InteractingResiduePair.objects.get(
            referenced_structure=self.structure, res2__sequence_number=200)
        # Children replaced, not accumulated: 1 pair, 1 interaction.
        self.assertEqual(
            InteractingResiduePair.objects.filter(referenced_structure=self.structure).count(), 1)
        self.assertEqual(Interaction.objects.filter(interacting_pair=pair).count(), 1)


class Engine2DispatchRoutingTests(TestCase):
    """The calculator routes Engine 2 chains by partner_category. 7RYC carries
    BOTH a G-alpha interface (chain D) and a legacy oxytocin peptide (chain L):
    the G-alpha chain -> InteractingResiduePair, the peptide chain -> the peptide
    model, and neither is double-written to the other model."""

    @classmethod
    def setUpTestData(cls):
        from common.models import WebLink, WebResource
        from ligand.models import Ligand, LigandPeptideStructure
        from protein.models import (
            Protein, ProteinConformation, ProteinFamily, ProteinSequenceType,
            ProteinSource, ProteinState, Species,
        )
        from residue.models import Residue, ResidueGenericNumber, ResidueNumberingScheme
        from signprot.models import SignprotComplex
        from structure.models import Structure, StructureType

        fam = ProteinFamily.objects.create(slug="000_001_001_DSP", name="Disp fam")
        species = Species.objects.create(latin_name="Sp disp", common_name="Disp")
        source = ProteinSource.objects.create(name="SRC_DSP")
        seqtype = ProteinSequenceType.objects.create(slug="wt_dsp", name="WT")
        state = ProteinState.objects.create(slug="active_dsp", name="Active")

        receptor = Protein.objects.create(
            family=fam, species=species, source=source, sequence_type=seqtype,
            entry_name="oxtr_disp", name="OXTR", sequence="M",
        )
        cls.receptor_pc = ProteinConformation.objects.create(protein=receptor, state=state)
        galpha = Protein.objects.create(
            family=fam, species=species, source=source, sequence_type=seqtype,
            entry_name="7ryc_a", name="Gq", sequence="M",
        )
        cls.alpha_pc = ProteinConformation.objects.create(protein=galpha, state=state)

        stype = StructureType.objects.create(slug="x-ray_dsp", name="X-ray")
        wr = WebResource.objects.create(slug="pdb_dsp", name="PDB", url="https://pdb/$index")
        weblink = WebLink.objects.create(index="7RYC", web_resource=wr)
        cls.structure = Structure.objects.create(
            protein_conformation=cls.receptor_pc, structure_type=stype, state=state,
            pdb_code=weblink, preferred_chain="O", publication_date="2021-01-01",
        )
        SignprotComplex.objects.create(protein=galpha, structure=cls.structure, alpha="D")

        # Receptor residues = union of chain-O residues across BOTH interfaces.
        recv_gp, gp, _ = _pairs_and_residues(_load(YAML_7RYC_GP))
        recv_pep, _, _ = _pairs_and_residues(_load(YAML_7RYC_PEP))
        recv = {**recv_pep, **recv_gp}
        for seq, resname in recv.items():
            Residue.objects.create(
                protein_conformation=cls.receptor_pc, sequence_number=seq,
                amino_acid=_T3[resname],
            )
        scheme = ResidueNumberingScheme.objects.create(slug="cgn_dsp", short_name="CGN", name="CGN")
        for seq, resname in gp.items():
            rgn = ResidueGenericNumber.objects.create(scheme=scheme, label=f"G.X.{seq}")
            Residue.objects.create(
                protein_conformation=cls.alpha_pc, sequence_number=seq,
                amino_acid=_T3[resname], display_generic_number=rgn,
            )

        # The legacy peptide (chain L) is LPS-keyed; the G-alpha chain (D) is NOT.
        lig = Ligand.objects.create(name="oxytocin")
        cls.lps = LigandPeptideStructure.objects.create(
            structure=cls.structure, ligand=lig, chain="L",
        )

    def test_routes_galpha_and_peptide_to_distinct_models(self):
        from contactnetwork.models import (
            InteractingPeptideResiduePair, InteractingResiduePair,
        )
        from interaction.calculators.schrodinger import SchrodingerInteractionCalculator

        results = SchrodingerInteractionCalculator().compute_interactions(
            self.structure, schrodinger_data_dir=ENGINE2_DATA_DIR,
        )
        labels = {label: ok for (_id, label, ok, _err) in results}
        self.assertEqual(labels.get("g_protein_alpha:D"), True)
        self.assertEqual(labels.get("peptide:L"), True)

        # G-alpha interface -> residue-residue model (res2 on "_a"), 20 pairs.
        gp_pairs = InteractingResiduePair.objects.filter(referenced_structure=self.structure)
        self.assertEqual(gp_pairs.count(), 20)
        self.assertTrue(all(
            p.res2.protein_conformation_id == self.alpha_pc.id for p in gp_pairs))

        # Legacy peptide -> peptide model, keyed off the chain-L LPS.
        pep_pairs = InteractingPeptideResiduePair.objects.filter(peptide=self.lps)
        self.assertGreater(pep_pairs.count(), 0)
        # No G-alpha-chain row leaked into the peptide model.
        self.assertEqual(
            InteractingPeptideResiduePair.objects.exclude(peptide=self.lps).count(), 0)


class Engine2GproteinNoComplexTests(TestCase):
    """A structure with a g_protein_alpha YAML but NO SignprotComplex (or no
    "_a" conformation) must return False gracefully, not raise."""

    PDB = "9NOC"

    @classmethod
    def setUpTestData(cls):
        from common.models import WebLink, WebResource
        from protein.models import (
            Protein, ProteinConformation, ProteinFamily, ProteinSequenceType,
            ProteinSource, ProteinState, Species,
        )
        from structure.models import Structure, StructureType

        fam = ProteinFamily.objects.create(slug="000_001_001_NOC", name="NoC fam")
        species = Species.objects.create(latin_name="Sp noc", common_name="NoC")
        source = ProteinSource.objects.create(name="SRC_NOC")
        seqtype = ProteinSequenceType.objects.create(slug="wt_noc", name="WT")
        state = ProteinState.objects.create(slug="active_noc", name="Active")
        receptor = Protein.objects.create(
            family=fam, species=species, source=source, sequence_type=seqtype,
            entry_name="noc_recv", name="NOC", sequence="M",
        )
        pc = ProteinConformation.objects.create(protein=receptor, state=state)
        stype = StructureType.objects.create(slug="x-ray_noc", name="X-ray")
        wr = WebResource.objects.create(slug="pdb_noc", name="PDB", url="https://pdb/$index")
        weblink = WebLink.objects.create(index=cls.PDB, web_resource=wr)
        cls.structure = Structure.objects.create(
            protein_conformation=pc, structure_type=stype, state=state,
            pdb_code=weblink, preferred_chain="R", publication_date="2021-01-01",
        )

    def test_no_signprot_complex_returns_false(self):
        import tempfile
        from contactnetwork.models import InteractingResiduePair

        doc = {
            "schema_version": "engine2/1.0",
            "metadata": {"ligand_type": "protein", "receptor_chain": "R",
                         "ligand_chain": "A", "partner_category": "g_protein_alpha"},
            "interface_interactions": [{
                "type": "hydrogen_bond", "hbond_class": "hb_ss",
                "donor": {"chain_id": "R", "resname": "ASP", "resid": 100,
                          "atom_name": "N", "side": "selection1"},
                "acceptor": {"chain_id": "A", "resname": "GLU", "resid": 200,
                             "atom_name": "O", "side": "selection2"},
            }],
            "residue_pair_summaries": [],
        }
        with tempfile.TemporaryDirectory() as tmp:
            inst = os.path.join(tmp, self.PDB, f"{self.PDB}_R_A")
            os.makedirs(inst)
            with open(os.path.join(inst, f"{self.PDB}_R_A.yaml"), "w") as fh:
                yaml.safe_dump(doc, fh)
            with self.assertLogs("interaction.schrodinger_processor", level="WARNING"):
                ok = process_schrodinger_gprotein_interactions(
                    current_structure_obj=self.structure, ligand_chain="A",
                    pdb_code_str=self.PDB, schrodinger_interactions_dir_override=tmp,
                )
        self.assertFalse(ok)
        self.assertEqual(InteractingResiduePair.objects.count(), 0)
