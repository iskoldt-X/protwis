from django.test import TestCase

from interaction.models import ResidueFragmentInteractionType


# The 21 canonical interaction-type slugs (18 original + 3 added 2026-05-28
# via Phase 1d: ADR-011 halogen_protein, ADR-012 water_bridge_protein,
# ADR-013 metal_coordination_protein). The old `parsecalculation` created the
# original 18 lazily via get_or_create during a structure build. The
# Schrodinger pipeline writes ResidueFragmentInteraction rows directly and
# therefore needs these foreign-key targets to already exist. The 3 new slugs
# are pre-seeded ahead of Plan D wrapper opt-in detection (halogen / water /
# metal contacts in Engine 1) per ADR-009 information preservation. See
# interaction/fixtures/interaction_types.json.
EXPECTED_SLUGS = {
    'acc',
    'hyd',
    'aro_ff',
    'aro_fe_protein',
    'aro_ef_protein',
    'aro_ion_protein',
    'polar_backbone',
    'polar_donor_protein',
    'polar_acceptor_protein',
    'polar_double_pos_protein',
    'polar_double_neg_protein',
    'polar_pos_ligand',
    'polar_neg_ligand',
    'polar_pos_protein',
    'polar_neg_protein',
    'polar_unknown_protein',
    'polar_unspecified',
    'Van der Waals',
    # ADR-011/012/013 (Phase 1d, 2026-05-28): protwis schema 18 -> 21.
    'halogen_protein',
    'water_bridge_protein',
    'metal_coordination_protein',
}

# The 3 slugs added by Phase 1d (ADR-011/012/013). Kept as a separate set so
# regressions on their fields are caught by name, not by the catch-all
# EXPECTED_SLUGS.
PHASE_1D_NEW_SLUGS = {
    'halogen_protein',
    'water_bridge_protein',
    'metal_coordination_protein',
}


class InteractionTypeFixtureTest(TestCase):
    """Lock the interaction_types.json seed fixture (Plan B, steps B2 & B3;
    extended 2026-05-28 by Phase 1d to 21 slugs).

    The ``fixtures`` attribute makes Django ``loaddata`` the file into a fresh
    test database before every test method, which is exactly the B2 scenario:
    an empty DB that must end up with the 21 canonical rows.
    """

    fixtures = ['interaction_types.json']

    # --- B2: fixture is loadable, count == 21 -----------------------------

    def test_fixture_loads(self):
        self.assertEqual(ResidueFragmentInteractionType.objects.count(), 21)

    def test_fixture_has_exact_slug_set(self):
        slugs = set(
            ResidueFragmentInteractionType.objects.values_list('slug', flat=True)
        )
        self.assertEqual(slugs, EXPECTED_SLUGS)

    # --- B3: key slug exactness -------------------------------------------

    def test_van_der_waals_has_space(self):
        # The slug literally contains spaces ("Van der Waals"). If this drifts
        # to e.g. "van-der-waals" the unique constraint lets old and new data
        # coexist, silently splitting the interaction table.
        self.assertTrue(
            ResidueFragmentInteractionType.objects.filter(
                slug='Van der Waals'
            ).exists()
        )

    def test_acc_is_hidden(self):
        # The UI explicitly excludes slug='acc'; its type must stay 'hidden'.
        acc = ResidueFragmentInteractionType.objects.get(slug='acc')
        self.assertEqual(acc.type, 'hidden')

    # --- Phase 1d (ADR-011/012/013): 3 new canonical slugs ----------------

    def test_phase_1d_slugs_exist(self):
        # All three Phase 1d slugs must exist after fixture load. Regression
        # of any one means the ADR-009 information-preservation guarantee has
        # been silently weakened (upstream Schrodinger detection has no slug
        # to land on -> info would be dropped).
        for slug in PHASE_1D_NEW_SLUGS:
            self.assertTrue(
                ResidueFragmentInteractionType.objects.filter(slug=slug).exists(),
                f"Phase 1d slug {slug!r} missing from fixture",
            )

    def test_phase_1d_slugs_have_polar_type_and_empty_direction(self):
        # ADR-011/012/013 all chose type='polar' (engineering pragmatic
        # categorization aligned with existing protwis 'type' enum) and
        # direction='' (chemistry is direction-less for these contact types,
        # analogous to polar_double_*_protein). Drift here would corrupt
        # downstream queries that filter by type or direction.
        for slug in PHASE_1D_NEW_SLUGS:
            row = ResidueFragmentInteractionType.objects.get(slug=slug)
            self.assertEqual(row.type, 'polar', f"{slug} type drift")
            self.assertEqual(row.direction, '', f"{slug} direction drift")
