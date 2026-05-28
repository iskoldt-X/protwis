"""Abstract InteractionCalculator (Strategy Pattern interface)."""

import abc


class InteractionCalculator(abc.ABC):
    """Strategy for computing residue-ligand interactions for a structure.

    Implementations:
      - RDKitInteractionCalculator: legacy real-time RDKit (preserves the
        existing build_structures behavior when settings.INTERACTION_CALCULATOR
        defaults to 'rdkit').
      - SchrodingerInteractionCalculator: Phase 1+ Schrodinger YAML consumer
        (entry point: interaction.schrodinger_processor).

    See ADR-007 (RFI table is the protocol layer / 北极星),
    ADR-009 (information preservation, fail-loud),
    Phase 1c 生产导入路径设计 axis 2E.
    """

    @abc.abstractmethod
    def compute_interactions(self, structure, **kwargs):
        """Compute and persist interactions for one Structure.

        Implementations MUST be idempotent within their own pipeline:
        repeated calls with the same inputs produce the same DB state
        (delete-then-insert is acceptable). See ADR-009.

        Strategy-specific kwargs are advertised by each implementation
        (e.g. RDKit needs sd + command for legacy build_structures contract;
        Schrodinger needs at most schrodinger_data_dir override).
        """
        raise NotImplementedError
