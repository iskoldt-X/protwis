from django.test import TestCase

from interaction.models import ResidueFragmentInteractionType


# The 18 canonical interaction-type slugs that the old `parsecalculation`
# created lazily via get_or_create during a structure build. The Schrodinger
# pipeline writes ResidueFragmentInteraction rows directly and therefore needs
# these foreign-key targets to already exist. See
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
}


class InteractionTypeFixtureTest(TestCase):
    """Lock the interaction_types.json seed fixture (Plan B, steps B2 & B3).

    The ``fixtures`` attribute makes Django ``loaddata`` the file into a fresh
    test database before every test method, which is exactly the B2 scenario:
    an empty DB that must end up with the 18 canonical rows.
    """

    fixtures = ['interaction_types.json']

    # --- B2: fixture is loadable, count == 18 -----------------------------

    def test_fixture_loads(self):
        self.assertEqual(ResidueFragmentInteractionType.objects.count(), 18)

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
