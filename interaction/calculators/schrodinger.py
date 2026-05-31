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
        from interaction.schrodinger_processor import (
            locate_engine2_yamls,
            process_schrodinger_gprotein_interactions,
            process_schrodinger_peptide_interactions,
            process_schrodinger_sm_interactions,
        )

        pdb_code_str = structure.pdb_code.index

        results = []

        # --- Engine 1: small-molecule SLIs (RFI rows) -----------------------
        slis = StructureLigandInteraction.objects.filter(structure=structure)
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

        # --- Engine 2: peptide / G-protein interfaces -----------------------
        # engine2/1.0 schema. Two target models, routed off the YAML's
        # partner_category (the truth, not ligand_type -- 7RYC carries BOTH a
        # G-alpha instance (chain D) and a legacy oxytocin peptide (chain L), and
        # they land in different models):
        #   * '' (legacy peptide)  -> InteractingPeptideResiduePair (LPS-keyed)
        #   * 'g_protein_alpha'    -> InteractingResiduePair (SignprotComplex-keyed)
        #   * 'arrestin'           -> deferred (separate StructureExtraProteins +
        #                             "_arrestin" bridge, not yet consumed)
        # Each processor does its own per-instance idempotency. Missing Engine 2
        # YAMLs are normal for SM-only PDBs.
        import yaml as _yaml
        from django.conf import settings as _settings
        from ligand.models import LigandPeptideStructure

        base_dir = schrodinger_data_dir or getattr(
            _settings, 'SCHRODINGER_INTERACTIONS_DIR', None)

        # Map each interface chain to its partner_category by peeking at the YAMLs.
        chain_category = {}
        if base_dir:
            for path in locate_engine2_yamls(base_dir, pdb_code_str):
                try:
                    with open(path) as fh:
                        meta = (_yaml.safe_load(fh) or {}).get('metadata') or {}
                except (OSError, _yaml.YAMLError):
                    continue
                ch = (meta.get('ligand_chain') or '').strip()
                if ch:
                    chain_category[ch] = (meta.get('partner_category') or '').strip()

        # Legacy peptide path (LPS-keyed). Skip chains the YAML tags as a
        # G-alpha interface (routed below) or arrestin (deferred -- recorded,
        # not run).
        lpss = LigandPeptideStructure.objects.filter(structure=structure)
        for lps in lpss:
            chain = (lps.chain or '').strip()
            if not chain:
                results.append((lps.id, None, False, 'missing_peptide_chain'))
                continue
            category = chain_category.get(chain, '')
            if category == 'g_protein_alpha':
                continue  # handled by the G-protein loop below
            if category == 'arrestin':
                results.append(
                    (lps.id, 'arrestin:{}'.format(chain), False, 'arrestin_phase2b'))
                continue
            ok = process_schrodinger_peptide_interactions(
                current_structure_obj=structure,
                ligand_chain=chain,
                pdb_code_str=pdb_code_str,
                schrodinger_interactions_dir_override=schrodinger_data_dir,
            )
            results.append((lps.id, 'peptide:{}'.format(chain), ok, None))

        # G-protein interface path (SignprotComplex-keyed -> InteractingResiduePair).
        # Driven off the YAMLs, NOT LigandPeptideStructure: the G-alpha chain may
        # have no LPS row (cube.py's do_complexes likewise keys off SignprotComplex).
        for chain in sorted(chain_category):
            if chain_category[chain] != 'g_protein_alpha':
                continue
            ok = process_schrodinger_gprotein_interactions(
                current_structure_obj=structure,
                ligand_chain=chain,
                pdb_code_str=pdb_code_str,
                schrodinger_interactions_dir_override=schrodinger_data_dir,
            )
            results.append(
                ('gp:{}'.format(chain), 'g_protein_alpha:{}'.format(chain), ok, None))

        return results
