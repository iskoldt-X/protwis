"""Pluggable interaction calculators (Strategy Pattern).

See Phase 1c 生产导入路径设计 axis 2E and ADR-007 (北极星 — RFI 表是协议层).
"""

from .base import InteractionCalculator


def get_interaction_calculator(name):
    """Resolve a calculator strategy by name.

    Valid names: 'rdkit' (legacy real-time), 'schrodinger' (Phase 1+ YAML).

    fail-loud (ADR-009): unknown name raises ValueError instead of silently
    falling back to a default — silent fallbacks hide config bugs.
    """
    from .rdkit import RDKitInteractionCalculator
    from .schrodinger import SchrodingerInteractionCalculator

    calculators = {
        'rdkit': RDKitInteractionCalculator,
        'schrodinger': SchrodingerInteractionCalculator,
    }
    if name not in calculators:
        raise ValueError(
            "Unknown INTERACTION_CALCULATOR={!r}, expected one of {}".format(
                name, sorted(calculators)))
    return calculators[name]()
