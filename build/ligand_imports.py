"""The Schrodinger ligand imports, as build_all and build_all_interactions run them.

Neither command computes ligand interactions with the legacy calculation; they
come from two imports of the Schrodinger deliveries: Engine 1 serves the
anchors named by a HET code (import_schrodinger_interactions), Engine 2's
peptide lane the "pep" chains (import_schrodinger_peptides). Both dry-run
first, then both run for real, so a structure that would fail stops the
caller with nothing imported rather than half imported: recovering is a
re-run, not an excavation.
"""

import datetime
import os

from django.conf import settings
from django.core.management import call_command
from django.core.management.base import CommandError

# Where the Engine 1 products and their per-PDB chain maps are delivered,
# relative to DATA_DIR.
ENGINE1_DIR = os.sep.join(['structure_data', 'schrodinger', 'engine1'])
# Where the Engine 2 products and their per-PDB peptide maps are delivered.
ENGINE2_DIR = os.sep.join(['structure_data', 'schrodinger', 'engine2'])

# Where each import run leaves its accounting, relative to BASE_DIR. Not in
# DATA_DIR: that is a git checkout of shared, versioned input, and a per-run,
# per-machine output does not belong in it -- `git clean` there would take the
# delivery with it. logs/ is this project's own runtime output directory and is
# ignored by git in full. One directory per run, so a later run cannot
# overwrite the record of an earlier one.
ENGINE1_RUN_DIR = os.sep.join(['logs', 'engine1_import'])
ENGINE2_RUN_DIR = os.sep.join(['logs', 'engine2_peptide_import'])

COMMANDS = ('import_schrodinger_interactions', 'import_schrodinger_peptides')


def _now():
    return datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')


def add_arguments(parser):
    """The options both callers take."""
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
                             'small molecules, Engine 2 "pep" chains) in this run; nothing '
                             'else computes them, so the ligand tables keep what they hold')


def engine1_dir(options):
    return options['engine1_data_dir'] or os.sep.join([settings.DATA_DIR, ENGINE1_DIR])


def engine2_dir(options):
    return options['engine2_data_dir'] or os.sep.join([settings.DATA_DIR, ENGINE2_DIR])


def steps(options):
    """[[command, options]]: both imports dry-run first, then both for real."""
    if options['skip_ligand_import']:
        print('{} SKIPPING the ligand imports: no ligand interactions are '
              'written'.format(_now()))
        return []
    stamp = datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')
    lanes = []
    for command, data_dir, report_dir, default_dir in (
            (COMMANDS[0], engine1_dir(options), options['engine1_report_dir'], ENGINE1_RUN_DIR),
            (COMMANDS[1], engine2_dir(options), options['engine2_report_dir'], ENGINE2_RUN_DIR)):
        run_dir = report_dir or os.sep.join([settings.BASE_DIR, default_dir, stamp])
        print('{} {} accounting goes to {}'.format(_now(), command, run_dir))
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


def check_deliveries(options, command_names):
    """Refuse to start when a delivery an import in ``command_names`` reads is missing."""
    for command, label, data_dir in ((COMMANDS[0], 'Engine 1', engine1_dir(options)),
                                     (COMMANDS[1], 'Engine 2', engine2_dir(options))):
        if command in command_names and not os.path.isdir(data_dir):
            raise CommandError(
                '{} products are missing: {} is not a directory. Deliver them, or pass '
                '--skip_ligand_import to go on without ligand interactions.'.format(
                    label, data_dir))


def run(options):
    """Run the import steps in order; a failing import raises and stops the caller."""
    planned = steps(options)
    check_deliveries(options, [c for c, _o in planned])
    for command, kwargs in planned:
        print('{} Running {}{}'.format(_now(), command,
                                       ' (dry run)' if kwargs.get('dry_run') else ''))
        call_command(command, **kwargs)
