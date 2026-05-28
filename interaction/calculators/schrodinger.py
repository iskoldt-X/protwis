"""Schrödinger YAML-based interaction calculator.

Thin wrapper over interaction.schrodinger_processor.process_schrodinger_sm_interactions.
Iterates each StructureLigandInteraction belonging to the structure and dispatches
to the processor, which is already @transaction.atomic internally.

Phase 1c 生产导入路径设计 axis 2E.
"""

from .base import InteractionCalculator


class SchrodingerInteractionCalculator(InteractionCalculator):
    """Schrödinger calculator.

    Accepts an optional `schrodinger_data_dir` kwarg to override the YAML
    base directory per-call (used by import_schrodinger_interactions
    --schrodinger-data-dir). When None, the processor falls back to
    settings.SCHRODINGER_INTERACTIONS_DIR.
    """

    def compute_interactions(self, structure, schrodinger_data_dir=None, **kwargs):
        from interaction.models import (
            ResidueFragmentInteraction,
            StructureLigandInteraction,
        )
        from interaction.schrodinger_processor import process_schrodinger_sm_interactions

        pdb_code_str = structure.pdb_code.index
        slis = StructureLigandInteraction.objects.filter(structure=structure)

        results = []
        for sli in slis:
            het_code = sli.pdb_reference
            if not het_code:
                results.append((sli.id, None, False, 'missing_pdb_reference'))
                continue

            # Phase 1c axis 5A: explicit per-SLI delete BEFORE the processor.
            # The processor's own delete-then-insert short-circuits on an empty
            # YAML (defensive: don't wipe legacy rows on a bad file), which means
            # a covalent ligand like 1F88 RET would leave 29 legacy RDKit rows
            # in place. For the production import path we want the new pipeline
            # to be the source of truth — so the import command's calculator
            # makes the delete explicit and unconditional.
            ResidueFragmentInteraction.objects.filter(structure_ligand_pair=sli).delete()

            ok = process_schrodinger_sm_interactions(
                current_structure_obj=structure,
                current_ligand_db_obj=sli.ligand,
                ligand_pdb_het_code=het_code,
                pdb_code_str=pdb_code_str,
                schrodinger_interactions_dir_override=schrodinger_data_dir,
            )
            results.append((sli.id, het_code, ok, None))
        return results
