"""Tests for the Schrödinger small-molecule interaction processor (Plan A).

Run with::

    docker exec gpcrdb-app python manage.py test interaction

protwis ships no pytest-django, and ``interaction.schrodinger_processor``
imports Django models at module load, so these tests run under Django's own
test runner. Pure-logic tests use ``SimpleTestCase`` (no test database; it
raises if a test touches the DB), which faithfully implements the Plan A rule
that parsing/mapping logic must be testable "纯逻辑、不碰 DB". DB-integration
tests (Steps A6 / A10) use ``TestCase`` with
``fixtures = ['interaction_types.json']`` (delivered by Plan B).

Fixture data: ``interaction/tests/schrodinger_data/`` holds a byte-for-byte
copy of the canonical Engine 1 output for 6LN2 / ligand 97Y, in its real
nested layout ``{PDB}/{PDB}/{HET}_{chain}_{resnum}/{HET}_{chain}_{resnum}.yaml``.
"""

import json
import os

import yaml
from django.test import SimpleTestCase, TestCase

from interaction.schrodinger_processor import (
    apply_backbone_override,
    apply_priority_dedup,
    get_receptor_pdb_block,
    is_standard_residue,
    locate_interaction_yamls,
    parse_receptor_residue,
    passes_chain_filter,
    process_schrodinger_sm_interactions,
    resolve_slug,
)

HERE = os.path.dirname(os.path.abspath(__file__))
SCHRODINGER_DATA_DIR = os.path.join(HERE, "schrodinger_data")
CANONICAL_YAML = os.path.join(
    SCHRODINGER_DATA_DIR, "6LN2_97Y", "6LN2_97Y", "97Y_A_503", "97Y_A_503.yaml"
)


FIXTURE_SLUGS_PATH = os.path.join(HERE, "..", "fixtures", "interaction_types.json")


def load_canonical_interactions():
    """Return the ``result.interactions`` list from the frozen 6LN2_97Y YAML."""
    with open(CANONICAL_YAML) as fh:
        return yaml.safe_load(fh)["result"]["interactions"]


def load_fixture_slugs():
    """The canonical slug set from Plan B's loaddata fixture (no DB needed)."""
    with open(FIXTURE_SLUGS_PATH) as fh:
        return {e["fields"]["slug"] for e in json.load(fh)}


class ReceptorResidueParsingTests(SimpleTestCase):
    """Step A1 — read the receptor residue number from a real YAML entry.

    Before the fix the processor read ``receptor_res_info["pdb_number"]`` while
    Engine 1 emits ``pdb_residue_number`` → ``KeyError: 'pdb_number'``.
    """

    def test_reads_receptor_residue_number(self):
        # interactions[0] in the frozen fixture: Hydroxyl / Donor, residue N406 chain A.
        entry = load_canonical_interactions()[0]
        parsed = parse_receptor_residue(entry["receptor_residue"])
        self.assertEqual(parsed["sequence_number"], 406)
        self.assertEqual(parsed["amino_acid"], "N")
        self.assertEqual(parsed["chain_id"], "A")
        self.assertEqual(parsed["insertion_code"], "")

    def test_parses_every_receptor_residue_in_fixture(self):
        # All six interactions must parse without KeyError and yield ints/strs.
        for entry in load_canonical_interactions():
            parsed = parse_receptor_residue(entry["receptor_residue"])
            self.assertIsInstance(parsed["sequence_number"], int)
            self.assertTrue(parsed["amino_acid"])
            self.assertEqual(parsed["chain_id"], "A")


class YamlLocationTests(SimpleTestCase):
    """Step A3 — locate the YAML in Engine 1's real nested layout.

    The original processor built a flat ``{PDB}_{HET}_interactions.yaml`` path,
    but Engine 1 emits
    ``{PDB}_{HET}/{PDB}_{HET}/{HET}_{chain}_{resnum}/{HET}_{chain}_{resnum}.yaml``.
    """

    def test_locate_yaml_for_ligand(self):
        found = locate_interaction_yamls(SCHRODINGER_DATA_DIR, "6LN2", "97Y")
        self.assertEqual(len(found), 1, f"expected exactly one instance, got {found}")
        self.assertTrue(found[0].endswith("97Y_A_503/97Y_A_503.yaml"))
        self.assertTrue(os.path.exists(found[0]))

    def test_lowercase_codes_are_normalised(self):
        # Callers may pass lowercase; layout dirs are uppercase.
        found = locate_interaction_yamls(SCHRODINGER_DATA_DIR, "6ln2", "97y")
        self.assertEqual(len(found), 1)

    def test_missing_ligand_returns_empty(self):
        self.assertEqual(
            locate_interaction_yamls(SCHRODINGER_DATA_DIR, "6LN2", "ZZZ"), []
        )


class ShallowLayoutLocationTests(SimpleTestCase):
    """Step 2c-1/S1 — locate the YAML in 2Y02's *shallow* layout.

    6LN2 nests the instance dirs one level deeper
    (``{PDB}_{HET}/{PDB}_{HET}/{HET}_*/...``), but 2Y02 places them directly
    under the top dir (``{PDB}_{HET}/{HET}_*/...``). ``locate_interaction_yamls``
    must handle both. 2Y02_WHJ has two ligand instances (chain A + chain B); the
    non-preferred chain A is dropped later by the chain filter (preferred=B).
    """

    def test_locates_both_2y02_instances(self):
        found = locate_interaction_yamls(SCHRODINGER_DATA_DIR, "2Y02", "WHJ")
        self.assertEqual(len(found), 2, f"expected chain A + B instances, got {found}")
        names = {os.path.basename(p) for p in found}
        self.assertEqual(names, {"WHJ_A_601.yaml", "WHJ_B_601.yaml"})
        for p in found:
            self.assertTrue(os.path.exists(p))

    def test_nested_6ln2_still_found(self):
        # Regression: adapting to the shallow layout must not break the nested one.
        found = locate_interaction_yamls(SCHRODINGER_DATA_DIR, "6LN2", "97Y")
        self.assertEqual(len(found), 1)
        self.assertTrue(found[0].endswith("97Y_A_503/97Y_A_503.yaml"))


class StandardResidueTests(SimpleTestCase):
    """Step A5 — non-standard ('X') receptor residues are recognised for skipping.

    Engine 1 writes ``name_1_letter: "X"`` for waters / ions / glycans / covalent
    modifications fused into the receptor (trap H5). They have no protwis Residue
    row, so the processor skips them with a warning instead of querying the DB.
    """

    def test_x_is_not_standard(self):
        self.assertFalse(is_standard_residue("X"))

    def test_twenty_standard_aas(self):
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            self.assertTrue(is_standard_residue(aa), aa)

    def test_lowercase_is_normalised(self):
        self.assertTrue(is_standard_residue("n"))

    def test_synthetic_x_entry_flagged(self):
        # A hand-built entry (the real 6LN2_97Y YAML has no X residue).
        x_entry = {
            "name_1_letter": "X",
            "pdb_residue_number": 9999,
            "chain_id": "A",
            "insertion_code": "",
        }
        parsed = parse_receptor_residue(x_entry)
        self.assertFalse(is_standard_residue(parsed["amino_acid"]))

    def test_canonical_fixture_all_standard(self):
        for entry in load_canonical_interactions():
            aa = parse_receptor_residue(entry["receptor_residue"])["amino_acid"]
            self.assertTrue(is_standard_residue(aa), aa)


class ChainFilterTests(SimpleTestCase):
    """Step A7 — only preferred_chain interactions are kept.

    protwis stores no per-residue chain; ``build_structures.py:1410`` drops any
    residue whose chain != ``Structure.preferred_chain`` (reduced to the first
    chain for multi-chain receptors at :248-249). We mirror that exactly.
    """

    def test_keeps_preferred_chain(self):
        self.assertTrue(passes_chain_filter("A", "A"))

    def test_drops_non_preferred_chain(self):
        self.assertFalse(passes_chain_filter("B", "A"))

    def test_multichain_preference_uses_first(self):
        # build_structures reduces "A,B" -> "A" before filtering.
        self.assertTrue(passes_chain_filter("A", "A,B"))
        self.assertFalse(passes_chain_filter("B", "A,B"))

    def test_empty_preference_keeps_all(self):
        self.assertTrue(passes_chain_filter("A", ""))

    def test_6ln2_all_chain_A(self):
        # The 6LN2_97Y fixture is entirely chain A; nothing is dropped (pref=A).
        for entry in load_canonical_interactions():
            chain = parse_receptor_residue(entry["receptor_residue"])["chain_id"]
            self.assertTrue(passes_chain_filter(chain, "A"))


class SlugResolutionTests(SimpleTestCase):
    """Steps A6a + A4 — (feature_family, direction) -> canonical slug via the
    delivered Plan C map (``interaction/interaction_type_map.yaml``).

    Routing ignores the chemical ``feature`` name, so comma-joined feature
    strings (trap H2 / §0.5 A5) can never corrupt a slug. Unknown combos raise
    (fail-loud) so non-default contacts surface instead of silently dropping.
    """

    EXPECTED_6LN2 = {
        ("Aromatic", "edge-to-face"): "aro_ef_protein",
        ("Donor", "ligand-donor"): "polar_acceptor_protein",
        ("Acceptor", "ligand-acceptor"): "polar_donor_protein",
    }

    def test_resolves_6ln2_combos(self):
        for (family, direction), slug in self.EXPECTED_6LN2.items():
            self.assertEqual(resolve_slug(family, direction), slug)

    def test_unknown_combo_raises(self):
        with self.assertRaises(ValueError):
            resolve_slug("Halogen", "halogen-bond")

    def test_feature_name_irrelevant_to_routing(self):
        # A4: resolve_slug takes no `feature` argument, so a comma-joined feature
        # is structurally incapable of producing a garbage slug.
        for entry in load_canonical_interactions():
            family, direction = entry["feature_family"], entry["direction"]
            self.assertEqual(
                resolve_slug(family, direction),
                self.EXPECTED_6LN2[(family, direction)],
            )

    def test_resolved_slugs_exist_in_plan_b_fixture(self):
        fixture_slugs = load_fixture_slugs()
        for entry in load_canonical_interactions():
            slug = resolve_slug(entry["feature_family"], entry["direction"])
            self.assertIn(slug, fixture_slugs)


class ReceptorPdbBlockTests(SimpleTestCase):
    """Step A8 — Fragment.pdbdata comes from the YAML receptor_pdb_block.

    The original code used ``sli.pdb_file`` (always None in this flow), gated by
    ``if not sli.pdb_file: continue`` — so Fragments were never created. Engine 1
    ships the receptor residue's atoms per interaction; we store that as opaque
    text (trap H3: non-standard PDB columns — never feed it to BioPython).
    """

    def test_extracts_block_text(self):
        entry = load_canonical_interactions()[0]  # ASN 406
        block = get_receptor_pdb_block(entry)
        self.assertIn("ASN A 406", block)
        self.assertTrue(block.lstrip().startswith("ATOM"))

    def test_every_entry_has_a_block(self):
        for entry in load_canonical_interactions():
            self.assertTrue(get_receptor_pdb_block(entry).strip())

    def test_missing_block_returns_empty(self):
        self.assertEqual(get_receptor_pdb_block({}), "")

    def test_fragment_not_gated_on_sli_pdb_file(self):
        # The B4 bug: `if not sli.pdb_file: continue` skipped every Fragment.
        import inspect
        from interaction import schrodinger_processor

        self.assertNotIn("if not sli.pdb_file", inspect.getsource(schrodinger_processor))


class TransactionAndIdempotencyTests(SimpleTestCase):
    """Step A9 — the per-PDB write is atomic and idempotent (delete-before-write).

    Behavioural proof (rollback leaves no dirty rows) is exercised by the A10
    characterization command, which runs the processor against real 6LN2 data
    inside a rolled-back transaction. These guards lock the implementation shape.
    """

    def _source(self):
        import inspect
        from interaction import schrodinger_processor

        return inspect.getsource(schrodinger_processor)

    def test_processor_is_atomic(self):
        self.assertIn("@transaction.atomic", self._source())

    def test_deletes_existing_rows_before_write(self):
        self.assertIn(
            "ResidueFragmentInteraction.objects.filter("
            "structure_ligand_pair=sli).delete()",
            self._source(),
        )


class LogicalEndToEndTests(SimpleTestCase):
    """Step A10 (logical) — the full parse -> skip-X -> chain-filter -> resolve ->
    backbone-override -> dedup chain over the real 6LN2_97Y YAML yields the
    expected (seq, slug) set, with no DB. The DB-writing e2e is the
    ``schrodinger_compare`` management command (run against the populated GPCRdb
    DB as a rolled-back dry run).

    Post Phase 1a-main (ADR-005 + A6b backbone override): the 6 YAML interactions
    yield 6 distinct rows. L401's H-bond is on a main-chain O and N407 has one
    entry on a main-chain N — both are reclassified to ``polar_backbone`` (views.py
    priority chain :1008-1019, backbone first). N407's other entry (on side-chain
    ND2) keeps ``polar_donor_protein``, so 407 ends up with two distinct rows:
    one backbone, one donor. Nothing dedups (all six (seq, slug) pairs distinct).
    """

    EXPECTED_NEW_SET = {
        (347, "aro_ef_protein"),
        (355, "polar_acceptor_protein"),       # T355 OG1 sidechain
        (401, "polar_backbone"),               # L401 O (backbone) — A6b override
        (406, "polar_acceptor_protein"),       # N406 OD1 sidechain
        (407, "polar_backbone"),               # N407 N (backbone) — A6b override
        (407, "polar_donor_protein"),          # N407 ND2 sidechain — untouched
    }

    def _run_logical_pipeline(self, preferred_chain="A"):
        result = set()
        for entry in load_canonical_interactions():
            parsed = parse_receptor_residue(entry["receptor_residue"])
            if not is_standard_residue(parsed["amino_acid"]):
                continue
            if not passes_chain_filter(parsed["chain_id"], preferred_chain):
                continue
            slug = resolve_slug(entry["feature_family"], entry["direction"])
            slug = apply_backbone_override(slug, entry.get("receptor_atom_name"))
            result.add((parsed["sequence_number"], slug))  # set membership = dedup
        return result

    def test_logical_pipeline_matches_expected(self):
        self.assertEqual(self._run_logical_pipeline(), self.EXPECTED_NEW_SET)

    def test_407_splits_into_backbone_and_donor(self):
        # Post-A6b: the two Acceptor entries on residue 407 carry different slugs
        # (one is on main-chain N → polar_backbone; the other on side-chain ND2 →
        # polar_donor_protein), so dedup does NOT collapse them.
        rows_407 = {slug for (seq, slug) in self._run_logical_pipeline() if seq == 407}
        self.assertEqual(rows_407, {"polar_backbone", "polar_donor_protein"})

    def test_401_promoted_to_backbone(self):
        # L401 has no polar sidechain (Leu); its H-bond is on the main-chain
        # carbonyl O → polar_backbone (Step A10 prediction validated).
        rows_401 = {slug for (seq, slug) in self._run_logical_pipeline() if seq == 401}
        self.assertEqual(rows_401, {"polar_backbone"})

    def test_sidechain_hbond_not_promoted(self):
        # N406 OD1 is a sidechain H-bond; A6b must not touch it.
        rows_406 = {slug for (seq, slug) in self._run_logical_pipeline() if seq == 406}
        self.assertEqual(rows_406, {"polar_acceptor_protein"})

    def test_all_new_slugs_are_canonical(self):
        canonical = load_fixture_slugs()
        for _, slug in self.EXPECTED_NEW_SET:
            self.assertIn(slug, canonical)


class BackboneOverrideTests(SimpleTestCase):
    """ADR-005 / A6b — main-chain N or O ⇒ ``polar_backbone`` overrides the
    donor/acceptor slug. This mirrors protwis views.py:1008-1019 where the
    elif-chain checks backbone first; we get the equivalent ordering by
    applying it as a post-routing override.

    Phase 1a-main: the deterministic A6b implementation depends on the new
    ``receptor_atom_name`` field emitted by the wrapper (ADR-005). Before that
    field landed, A6b had to be deferred — the old code couldn't tell a
    main-chain H-bond from a side-chain one without re-parsing the
    ``receptor_pdb_block`` and matching atoms by index (the ADR-005 Option 2
    geometry-reconstruction path, rejected because of the RDKit/Schrödinger
    dual-index hazard).
    """

    def test_main_chain_n_promotes_donor(self):
        # Acceptor entry whose protein partner is the backbone N becomes
        # polar_backbone (protein donates H via main-chain N).
        self.assertEqual(
            apply_backbone_override("polar_donor_protein", "N"), "polar_backbone"
        )

    def test_main_chain_o_promotes_acceptor(self):
        # Donor entry whose protein partner is the backbone carbonyl O.
        self.assertEqual(
            apply_backbone_override("polar_acceptor_protein", "O"), "polar_backbone"
        )

    def test_sidechain_o_not_promoted(self):
        # ASP OD1/OD2, GLU OE1/OE2, SER/THR OG/OG1, etc. — single-letter O
        # detection must NOT trigger; we match the exact atom name.
        for atom in ("OD1", "OD2", "OE1", "OG", "OG1", "OH"):
            self.assertEqual(
                apply_backbone_override("polar_acceptor_protein", atom),
                "polar_acceptor_protein",
                f"atom={atom} should not promote",
            )

    def test_sidechain_n_not_promoted(self):
        for atom in ("ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2", "NZ"):
            self.assertEqual(
                apply_backbone_override("polar_donor_protein", atom),
                "polar_donor_protein",
                f"atom={atom} should not promote",
            )

    def test_whitespace_stripped(self):
        # Some Schrödinger atom name fields ship with padding; .strip() guard
        # is essential (and is also done at the wrapper side via .pdbname.strip()).
        self.assertEqual(
            apply_backbone_override("polar_donor_protein", " N "), "polar_backbone"
        )

    def test_missing_atom_name_no_op(self):
        # Old fixtures without the ADR-005 field must not crash.
        self.assertEqual(
            apply_backbone_override("polar_donor_protein", None),
            "polar_donor_protein",
        )
        self.assertEqual(
            apply_backbone_override("polar_donor_protein", ""),
            "polar_donor_protein",
        )

    def test_carbon_atoms_not_promoted(self):
        # Aromatic/PiCat entries carry carbon atoms (CG / CD1 / CD2 / CZ). The
        # override must only fire for N/O — never CA / C / S.
        for atom in ("CA", "C", "CB", "CG", "CD2", "CZ", "SD", "SG"):
            self.assertEqual(
                apply_backbone_override("aro_ef_protein", atom),
                "aro_ef_protein",
                f"atom={atom} should not promote",
            )

    def test_only_polar_hbond_slugs_affected(self):
        # Salt/aromatic/picat slugs must be left alone even if atom name is N/O
        # (defensive — in practice they don't appear with backbone N/O, but the
        # override semantics are scoped to donor/acceptor by views.py:1008-1019).
        for slug in (
            "polar_double_neg_protein",
            "polar_double_pos_protein",
            "aro_ef_protein",
            "aro_ff_protein",
            "aro_ion_protein",
        ):
            self.assertEqual(
                apply_backbone_override(slug, "N"), slug, f"slug={slug} should pass through"
            )

    def test_polar_backbone_is_in_canonical_fixture(self):
        # Plan B fixture must seed the slug or A6 will fail-loud at DB insert.
        canonical = load_fixture_slugs()
        self.assertIn("polar_backbone", canonical)


class BuildStructuresIntegrationTests(TestCase):
    """Step A2 — build_structures wires in the processor under its correct name.

    The stale ``interaction`` branch called ``self.process_schrodinger_interactions(...)``
    (wrong name — missing ``_sm_`` — invoked as a method, never imported) which
    raised NameError/AttributeError and was a key reason that branch could not run.
    This branch rebuilds the seam from clean dev_build: a module-level import of the
    correctly named ``process_schrodinger_sm_interactions``.

    NB: ``TestCase`` (not ``SimpleTestCase``) because importing the build command
    module executes ``test_model_updates(initialize=True)`` at class-definition
    time, which queries the DB — so a real test DB must exist.
    """

    def test_build_structures_imports_processor(self):
        from build.management.commands import build_structures as bs

        self.assertTrue(
            hasattr(bs, "process_schrodinger_sm_interactions"),
            "build_structures must import process_schrodinger_sm_interactions",
        )
        self.assertIs(
            bs.process_schrodinger_sm_interactions, process_schrodinger_sm_interactions
        )

    def test_no_stale_misnamed_processor_call(self):
        # Guard against regressing to the §0.5 A2 typo (name without `_sm_`).
        import inspect
        from build.management.commands import build_structures as bs

        source = inspect.getsource(bs)
        self.assertNotIn("process_schrodinger_interactions(", source)


class PriorityDedupTests(SimpleTestCase):
    """A6c — identical-(seq, slug) dedup only (ADR-008 / ADR-009).

    Reverted 2026-05-28: the earlier 2c-1 implementation mirrored protwis
    views.py:998/1002 by dropping plain h-bond slugs on charged residues that
    also carried a charge slug. ADR-008 reverses that — charge-assisted h-bonds
    carry two independent physical forces (Coulomb + directional dipole) and
    Schrödinger's two detectors fire independently, so both rows are kept as a
    faithful record. ADR-009 generalises: DB layer preserves information; UI
    compactness is a front-end concern. Only identical (seq, slug) duplicates
    collapse — no information loss.
    """

    def _keys(self, records):
        return {(r["sequence_number"], r["slug"]) for r in apply_priority_dedup(records)}

    def test_charge_does_not_suppress_hbond_on_charged_residue(self):
        # ADR-008: D121 with both a PosCharge salt bridge and a Donor h-bond
        # keeps BOTH rows; the two identical Donor entries collapse to one.
        records = [
            {"sequence_number": 121, "amino_acid": "D", "slug": "polar_acceptor_protein"},
            {"sequence_number": 121, "amino_acid": "D", "slug": "polar_acceptor_protein"},
            {"sequence_number": 121, "amino_acid": "D", "slug": "polar_double_neg_protein"},
        ]
        self.assertEqual(
            self._keys(records),
            {(121, "polar_acceptor_protein"), (121, "polar_double_neg_protein")},
        )

    def test_noncharged_residue_keeps_both_hbond_directions(self):
        # ASN 329 both donates and accepts -> both rows survive.
        records = [
            {"sequence_number": 329, "amino_acid": "N", "slug": "polar_acceptor_protein"},
            {"sequence_number": 329, "amino_acid": "N", "slug": "polar_donor_protein"},
        ]
        self.assertEqual(
            self._keys(records),
            {(329, "polar_acceptor_protein"), (329, "polar_donor_protein")},
        )

    def test_charged_residue_without_charge_slug_keeps_hbond(self):
        # A lone h-bond on a charged residue (no salt entry) passes through —
        # consistent with ADR-008/009: no information loss either way.
        records = [
            {"sequence_number": 200, "amino_acid": "D", "slug": "polar_acceptor_protein"},
        ]
        self.assertEqual(self._keys(records), {(200, "polar_acceptor_protein")})

    def test_backbone_and_charge_both_survive(self):
        # ADR-008/009: distinct slugs all survive. No priority cascade is
        # applied in the DB layer; if a future polar_backbone slug ever lands
        # on the same residue as a charge slug, both are recorded as faithful
        # readings of independent detectors. (A6b backbone slug is still
        # deferred — see Inbox-5 — so this is forward-looking.)
        records = [
            {"sequence_number": 121, "amino_acid": "D", "slug": "polar_backbone"},
            {"sequence_number": 121, "amino_acid": "D", "slug": "polar_double_neg_protein"},
        ]
        self.assertEqual(
            self._keys(records),
            {(121, "polar_backbone"), (121, "polar_double_neg_protein")},
        )

    def test_dedup_identical_seq_slug(self):
        records = [
            {"sequence_number": 307, "amino_acid": "F", "slug": "aro_ef_protein"},
            {"sequence_number": 307, "amino_acid": "F", "slug": "aro_ef_protein"},
        ]
        self.assertEqual(len(apply_priority_dedup(records)), 1)

    def test_aromatic_channel_untouched_on_charged_residue(self):
        records = [
            {"sequence_number": 121, "amino_acid": "D", "slug": "polar_double_neg_protein"},
            {"sequence_number": 121, "amino_acid": "D", "slug": "aro_ion_protein"},
        ]
        self.assertEqual(
            self._keys(records),
            {(121, "polar_double_neg_protein"), (121, "aro_ion_protein")},
        )


def _load_2y02_records(preferred_chain="B"):
    """Parse -> skip-X -> chain-filter -> resolve -> backbone-override, aggregating
    both 2Y02_WHJ instances (chain A + chain B) the way the orchestrator does.
    No DB, no dedup. Mirrors processor:349 + 448 + backbone override (A6b)."""
    records = []
    for path in locate_interaction_yamls(SCHRODINGER_DATA_DIR, "2Y02", "WHJ"):
        with open(path) as fh:
            interactions = yaml.safe_load(fh)["result"]["interactions"]
        for entry in interactions:
            parsed = parse_receptor_residue(entry["receptor_residue"])
            if not is_standard_residue(parsed["amino_acid"]):
                continue
            if not passes_chain_filter(parsed["chain_id"], preferred_chain):
                continue
            slug = resolve_slug(entry["feature_family"], entry["direction"])
            slug = apply_backbone_override(slug, entry.get("receptor_atom_name"))
            records.append(
                {
                    "sequence_number": parsed["sequence_number"],
                    "amino_acid": parsed["amino_acid"],
                    "slug": slug,
                }
            )
    return records


class LogicalEndToEnd2Y02Tests(SimpleTestCase):
    """Step 2c-1 (logical e2e) — the full chain over 2Y02_WHJ, no DB, no golden.

    2Y02 (turkey β1-adrenergic receptor + carmoterol/WHJ) has RFI_rows=0 in the
    DB, so there is nothing to diff against; instead this verifies the new
    pipeline against known β-AR pharmacology. It exercises the two things 6LN2
    could not: ADR-002 salt-bridge mapping (D121/D3.32 -> polar_double_neg_protein)
    and A6c dedup. Per ADR-008/009 (2026-05-28), the previous "charge suppresses
    h-bond" branch was reverted: D121 now keeps BOTH the salt bridge and the
    h-bond row (total 10 rows, was 9).
    """

    EXPECTED_2Y02_SET = {
        (121, "polar_double_neg_protein"),  # Asp3.32 salt bridge (ADR-002)
        (121, "polar_acceptor_protein"),    # Asp3.32 charge-assisted H-bond (ADR-008: kept)
        (211, "polar_acceptor_protein"),    # Ser5.42 catechol H-bond
        (215, "polar_donor_protein"),       # Ser5.46 catechol H-bond
        (310, "polar_donor_protein"),       # Asn6.55
        (329, "polar_acceptor_protein"),    # Asn7.39 donates
        (329, "polar_donor_protein"),       # Asn7.39 accepts
        (117, "aro_ion_protein"),           # Trp3.28 pi-cation
        (201, "aro_ef_protein"),            # Phe45.52 (ECL2) edge-to-face
        (307, "aro_ef_protein"),            # Phe6.52 edge-to-face
    }

    def _survivor_keys(self, preferred_chain="B"):
        survivors = apply_priority_dedup(_load_2y02_records(preferred_chain))
        return {(r["sequence_number"], r["slug"]) for r in survivors}

    def test_pipeline_matches_expected(self):
        self.assertEqual(self._survivor_keys(), self.EXPECTED_2Y02_SET)

    def test_chain_A_instance_fully_dropped(self):
        # locate returns both WHJ_A_601 and WHJ_B_601; preferred=B keeps only B.
        self.assertTrue(self._survivor_keys("B"))
        self.assertEqual(self._survivor_keys("Z"), set())  # no chain Z -> nothing

    def test_d121_keeps_both_salt_and_hbond(self):
        # ADR-008: D121 carries BOTH the salt bridge and the directional h-bond
        # — two faithful readings of independent Schrödinger detectors.
        rows_121 = {slug for (seq, slug) in self._survivor_keys() if seq == 121}
        self.assertEqual(
            rows_121, {"polar_double_neg_protein", "polar_acceptor_protein"}
        )

    def test_a6c_only_collapses_identical_dupes(self):
        # Raw load has 11 (seq, slug) rows including N2's two identical
        # Donor->D121 entries (both -> polar_acceptor_protein). A6c collapses
        # only that exact duplicate; nothing else is removed -> 10 rows.
        records = _load_2y02_records("B")
        before_list = [(r["sequence_number"], r["slug"]) for r in records]
        after = {
            (r["sequence_number"], r["slug"]) for r in apply_priority_dedup(records)
        }
        # the (121, polar_acceptor_protein) pair appears twice pre-dedup
        self.assertEqual(before_list.count((121, "polar_acceptor_protein")), 2)
        # exactly that duplicate is collapsed -> 10 unique survivors
        self.assertEqual(len(after), 10)
        # set of unique pre-dedup rows == set of post-dedup rows (no info lost)
        self.assertEqual(set(before_list), after)

    def test_w117_picat_maps_aro_ion(self):
        rows_117 = {slug for (seq, slug) in self._survivor_keys() if seq == 117}
        self.assertEqual(rows_117, {"aro_ion_protein"})

    def test_all_slugs_canonical(self):
        canonical = load_fixture_slugs()
        for _, slug in self._survivor_keys():
            self.assertIn(slug, canonical)

    def test_a6b_does_not_touch_2y02(self):
        # Sanity check that the backbone override is correctly *scoped*. 2Y02's
        # H-bonds are all on side-chain oxygens (D121 OD1/OD2 carboxyl, S211/S215
        # OG hydroxyl, N310/N329 OD1/ND2 amide) or aromatic carbons (W117/F201/F307
        # CG/CD2). None are main-chain N or O. So A6b must promote zero rows
        # here — the EXPECTED_2Y02_SET above carries no ``polar_backbone`` entry,
        # which validates the ADR-005 / ADR-008 decoupling: charge-assisted h-bond
        # on D121 stays as polar_acceptor_protein + polar_double_neg_protein
        # rather than being collapsed into polar_backbone.
        backbones = [slug for (_, slug) in self._survivor_keys() if slug == "polar_backbone"]
        self.assertEqual(backbones, [])
