"""Legacy RDKit-based real-time interaction calculator.

Wraps the existing per-ligand loop from
build.management.commands.build_structures.Command.main_func (the block
that calls runcalculation_2022 + parsecalculation). Preserves byte-for-byte
behavior under settings.INTERACTION_CALCULATOR='rdkit' (the default).

Phase 1c 生产导入路径设计 axis 2E: this is the legacy strategy. It exists so
build_structures keeps working as before while a parallel Schrödinger path
is wired up via the import_schrodinger_interactions management command.
"""

import time

from .base import InteractionCalculator


class RDKitInteractionCalculator(InteractionCalculator):
    """Legacy RDKit calculator. Required kwargs: sd, command.

    sd:      parsed structure dict from ParseStructureCSV (must contain
             'pdb' and the ligand list reachable as 'ligand').
    command: the build_structures.Command instance (provides logger,
             interaction_errors collector, and the static parsecalculation
             method that writes RFI rows).
    """

    def compute_interactions(self, structure, sd=None, command=None, **kwargs):
        if sd is None or command is None:
            raise RuntimeError(
                "RDKitInteractionCalculator requires sd= and command= kwargs "
                "(legacy build_structures contract). See Phase 1c axis 2E.")

        from interaction.views import runcalculation_2022

        ligand_list = sd.get('ligand') or []
        if not isinstance(ligand_list, list):
            ligand_list = [ligand_list]

        pdb_code_str = sd['pdb']
        entry_name = structure.protein_conformation.protein.entry_name

        for ligand in ligand_list:
            if ligand['type'].strip() not in ['small molecule', 'protein', 'peptide']:
                continue
            if not ligand.get('in_structure'):
                continue
            try:
                current = time.time()
                peptide_chain = ""
                if ligand['chain'] != '':
                    peptide_chain = ligand['chain']
                data_results = runcalculation_2022(pdb_code_str, peptide_chain)
                if 'NAG' in data_results:
                    del data_results['NAG']
                command.parsecalculation(pdb_code_str, data_results, ligand['name'], False)
                end = time.time()
                diff = round(end - current, 1)
                print('Interaction calculations done for {}. {} seconds.'.format(entry_name, diff))
                command.logger.info(
                    'Interaction calculations done for {}. {} seconds.'.format(entry_name, diff))
            except Exception as msg:
                print(msg)
                print('ERROR WITH INTERACTIONS {}'.format(pdb_code_str))
                command.logger.error(
                    'Error parsing interactions output for {}'.format(pdb_code_str))
                command.interaction_errors.append(structure)
