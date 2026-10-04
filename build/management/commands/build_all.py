from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command
import os
from django.conf import settings

import datetime

# Where the Engine 1 products and their per-PDB chain maps are delivered,
# relative to DATA_DIR.
ENGINE1_DIR = os.sep.join(['structure_data', 'schrodinger', 'engine1'])
# Where the Engine 2 products and their per-PDB peptide maps are delivered.
ENGINE2_DIR = os.sep.join(['structure_data', 'schrodinger', 'engine2'])

# Where each import run leaves its accounting, relative to BASE_DIR. Not in
# DATA_DIR: that is a git checkout of shared, versioned input, and a per-run,
# per-machine output does not belong in it -- `git clean` there would take the
# delivery with it. logs/ is this project's own runtime output directory and is
# ignored by git in full. One directory per run, so a later build cannot
# overwrite the record of an earlier one.
ENGINE1_RUN_DIR = os.sep.join(['logs', 'engine1_import'])
ENGINE2_RUN_DIR = os.sep.join(['logs', 'engine2_peptide_import'])


class Command(BaseCommand):
    help = 'Runs all build functions'

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
                            type=int,
                            action='store',
                            dest='proc',
                            default=1,
                            help='Number of processes to run')
        parser.add_argument('-t', '--test',
                            action='store_true',
                            dest='test',
                            default=False,
                            help='Include only a subset of data for testing')
        # parser.add_argument('--hommod',
        #                     action='store_true',
        #                     dest='hommod',
        #                     default=False,
        #                     help='Include build of homology models')
        parser.add_argument('--no_reload',
                            action='store_true',
                            dest='no_reload',
                            default=False,
                            help='Skip ligand dump reload scripts')
        parser.add_argument('--phase',
                            type=int,
                            action='store',
                            dest='phase',
                            default=None,
                            help='Specify build phase to run (1 or 2, default: None)')
        parser.add_argument('--engine1_data_dir',
                            action='store',
                            dest='engine1_data_dir',
                            default=None,
                            help='Engine 1 product tree; default DATA_DIR/' + ENGINE1_DIR)
        parser.add_argument('--engine1_report_dir',
                            action='store',
                            dest='engine1_report_dir',
                            default=None,
                            help='Where this run leaves its Engine 1 import accounting; '
                                 'default BASE_DIR/' + ENGINE1_RUN_DIR + '/<timestamp>')
        parser.add_argument('--engine2_data_dir',
                            action='store',
                            dest='engine2_data_dir',
                            default=None,
                            help='Engine 2 product tree; default DATA_DIR/' + ENGINE2_DIR)
        parser.add_argument('--engine2_report_dir',
                            action='store',
                            dest='engine2_report_dir',
                            default=None,
                            help='Where this run leaves its Engine 2 peptide import accounting; '
                                 'default BASE_DIR/' + ENGINE2_RUN_DIR + '/<timestamp>')
        parser.add_argument('--skip_ligand_import',
                            action='store_true',
                            dest='skip_ligand_import',
                            default=False,
                            help='Do not import the Schrodinger ligand interactions (Engine 1 '
                                 'small molecules, Engine 2 "pep" chains). build_structures '
                                 'does not compute them, so the build then has none')

    def engine1_dir(self, options):
        return options['engine1_data_dir'] or os.sep.join([settings.DATA_DIR, ENGINE1_DIR])

    def engine2_dir(self, options):
        return options['engine2_data_dir'] or os.sep.join([settings.DATA_DIR, ENGINE2_DIR])

    def ligand_import_steps(self, options):
        """The ligand interactions: both imports dry-run first, then both for real.

        build_structures computes no ligand interaction; these two imports are
        where the ligand tables come from. Engine 1 serves anchors named by a
        HET code, Engine 2's peptide lane the "pep" chains. Each dry run rolls
        every structure back, and both run before either import writes, so a
        structure that would fail stops the build with nothing imported rather
        than half imported: recovering is a re-run, not an excavation. Each
        run writes its accounting next to the others of the same build.
        """
        if options['skip_ligand_import']:
            print('{} SKIPPING the ligand imports: this build has no ligand '
                  'interactions'.format(datetime.datetime.strftime(
                      datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')))
            return []
        stamp = datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
        lanes = []
        for command, data_dir, report_dir, default_dir in (
                ('import_schrodinger_interactions', self.engine1_dir(options),
                 options['engine1_report_dir'], ENGINE1_RUN_DIR),
                ('import_schrodinger_peptides', self.engine2_dir(options),
                 options['engine2_report_dir'], ENGINE2_RUN_DIR)):
            run_dir = report_dir or os.sep.join([settings.BASE_DIR, default_dir, stamp])
            print('{} {} accounting goes to {}'.format(
                datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'),
                command, run_dir))
            lanes.append((command, data_dir, run_dir))
        dry = [[command,
                {'data_dir': data_dir, 'dry_run': True,
                 'anomaly_csv': os.path.join(run_dir, 'anomalies.dryrun.csv'),
                 'report_json': os.path.join(run_dir, 'report.dryrun.json')}]
               for command, data_dir, run_dir in lanes]
        real = [[command,
                 {'data_dir': data_dir,
                  'anomaly_csv': os.path.join(run_dir, 'anomalies.csv'),
                  'report_json': os.path.join(run_dir, 'report.json')}]
                for command, data_dir, run_dir in lanes]
        return dry + real

    def handle(self, *args, **options):
        if options['test']:
            print('Running in test mode')

        phase1 = [
            ['clear_cache'],
            ['build_common'],
            ['build_citations'],
            ['build_human_proteins'],
            ['build_blast_database'],
            ['build_other_proteins', {'constructs_only': options['test'] ,'proc': options['proc']}], # build only constructs in test mode
            ['build_classification_annotations'],
            ['build_annotation', {'proc': options['proc']}],
            ['build_blast_database'],
            ['build_links'],
            ['build_construct_proteins'],
            ['build_experimental_data_light', {'test_run': options['test'], 'no_reload': options['no_reload']}],
            # ['build_all_gtp_ligands', {'test_run': options['test']}],
            # ['build_endogenous_data_from_gtp_source', {'test_run': options['test']}],
            ['build_bias_preprocess_data', {'test_run': options['test']}],
            #['build_balanced_ligands', {'test_run': options['test']}],
            # ['build_chembl_data', {'test_run': options['test']}],
            ['build_mutant_data', {'test_run': options['test']}],
            ['build_structures', {'proc': options['proc'], 'skip_cn': options['test']}],
            # The ligand interactions come from the Schrodinger deliveries,
            # imported before anything reads them.
            *(self.ligand_import_steps(options) if options['phase'] in (None, 1) else []),
            ['build_consensus_sequences', {'proc': options['proc']}],
            ['build_g_proteins'],
            ['build_consensus_sequences', {'proc': options['proc'], 'signprot': 'Alpha'}],
            ['build_arrestins'],
            ['build_coupling_data'],
            ['build_consensus_sequences', {'proc': options['proc'], 'signprot': 'Arrestin'}],
            ['build_signprot_complex'],
            ['build_g_protein_structures', {'proc': options['proc']}],
            ['build_arrestin_structures'],
            ['build_structure_extra_proteins'],
            ['build_structure_model_rmsd'],
            ['build_gain_domain'],
            ['build_blast_database']
        ]
        phase2 = [
            ['build_structure_angles', {'proc': options['proc']}],
            ['build_construct_data', {'proc': options['proc']}],
            ['update_construct_mutations', {'proc': options['proc']}],
            ['build_protein_sets'],
            ['build_drugs_updated'],
            ['build_mutational_landscape'],
            ['build_residue_sets'],
            ['build_dynamine_annotation', {'proc': options['proc']}],
            ['build_complex_interactions'],
            ['assign_structure_states'],
            ['build_contact_representative'],
            ['build_mammalian_representative'],
            ['upload_excel_bias_pathways'],
            ['build_receptor_similarity'],
            ['build_treenetwork'],
            ['build_structure_similarity'],
            ['build_clustercoord'],
            ['build_ligand_search'],
            ['build_text'],
            ['build_frontend_table_datasources'],
        ]
        phase3 = [
            ['build_complex_models', {'proc': options['proc'], 'parser' : 'alphafoldcomplex', 'model_set_name' : 'AlphaFold_multimer_non_phys', 'cleaned_seq_csv' : os.sep.join([settings.DATA_DIR, 'structure_data', 'AlphaFold_multimer_non_phys', 'cleaned_seqs.csv']) }],
            ['build_complex_models', {'proc': options['proc'], 'parser' : 'alphafoldcomplex', 'model_set_name' : 'AlphaFold_multimer_phys' }],
            ['build_complex_models', {'proc': options['proc'], 'parser' : 'alphafoldcomplex', 'model_set_name' : 'AlphaFold_multimer_G_protein' }],
            ['build_complex_models', {'proc': options['proc'], 'parser' : 'alphafoldcomplex', 'model_set_name' : 'Arrestins_AF_models', "deposition_date": '2024-10-31'}],
            ['build_complex_models', {'proc': options['proc'], 'parser' : 'boltztwocomplex', 'model_set_name' : 'boltz2_complex', "deposition_date": '2026-03-01'}],
            ['build_rfaa_models'],
            ### build_homology_models --alphafold -r {active pdbs} -p ### build refined structures for new G prot coupled structures
            ['build_homology_models_zip', {'proc': options['proc']}],
            ['build_homology_models_zip', {'proc': options['proc'], 'c': True}],
            ['foldseek_db', {'proc': options['proc']}],
            ['build_release_notes'],
        ]

        if options['phase']:
            if options['phase']==1:
                commands = phase1
            elif options['phase']==2:
                commands = phase2
            elif options['phase']==3:
                commands = phase3
        else:
            commands = phase1+phase2+phase3

        for command, label, data_dir in (
                ('import_schrodinger_interactions', 'Engine 1', self.engine1_dir(options)),
                ('import_schrodinger_peptides', 'Engine 2', self.engine2_dir(options))):
            if any(c[0] == command for c in commands) and not os.path.isdir(data_dir):
                raise CommandError(
                    '{} products are missing: {} is not a directory. Deliver them, or pass '
                    '--skip_ligand_import to build without ligand interactions.'.format(
                        label, data_dir))

        for c in commands:
            print('{} Running {}'.format(
                datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), c[0]))
            if len(c) == 2:
                call_command(c[0], **c[1])
            elif len(c) == 3:
                call_command(c[0], *c[1], **c[2])
            else:
                call_command(c[0])

        print('{} Build completed'.format(datetime.datetime.strftime(
            datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')))
