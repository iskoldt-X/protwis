"""Tests for the InteractionCalculator strategy factory.

Phase 1c axis 2E. These tests cover:
 - Factory dispatch to known strategies.
 - fail-loud on unknown strategy names (ADR-009).
 - Subclass / instance shape of returned objects.
"""

from django.test import TestCase

from interaction.calculators import (
    InteractionCalculator,
    get_interaction_calculator,
)
from interaction.calculators.rdkit import RDKitInteractionCalculator
from interaction.calculators.schrodinger import SchrodingerInteractionCalculator


class CalculatorFactoryTests(TestCase):

    def test_rdkit_dispatch(self):
        calc = get_interaction_calculator('rdkit')
        self.assertIsInstance(calc, RDKitInteractionCalculator)
        self.assertIsInstance(calc, InteractionCalculator)

    def test_schrodinger_dispatch(self):
        calc = get_interaction_calculator('schrodinger')
        self.assertIsInstance(calc, SchrodingerInteractionCalculator)
        self.assertIsInstance(calc, InteractionCalculator)

    def test_unknown_name_fails_loud(self):
        # ADR-009: silent fallback would hide config typos.
        with self.assertRaises(ValueError) as ctx:
            get_interaction_calculator('rdket')
        msg = str(ctx.exception)
        self.assertIn("'rdket'", msg)
        self.assertIn('rdkit', msg)
        self.assertIn('schrodinger', msg)

    def test_rdkit_requires_legacy_kwargs(self):
        # The legacy strategy preserves the build_structures contract — it must
        # refuse to run without sd/command rather than silently no-op.
        calc = get_interaction_calculator('rdkit')
        with self.assertRaises(RuntimeError):
            calc.compute_interactions(structure=None)
