from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.core.management import call_command
from django.conf import settings
from django.db import connection
from structure.models import Structure

from contactnetwork.cube import *
from build import ligand_imports

import logging, json, os

class Command(BaseBuild):

    help = ("Recompute the interactions of all experimental GPCR structures: the "
            "intra-receptor contact network (legacy calculation; no Schrodinger lane yet), "
            "then the ligand interactions by importing the Schrodinger deliveries, as "
            "build_all does (Engine 1 for HET anchors, Engine 2 for \"pep\" chains).")

    logger = logging.getLogger(__name__)
    pdbs = Structure.objects.filter(structure_type__origin='experiment').values_list('pdb_code__index', flat=True)


    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
            type=int,
            action='store',
            dest='proc',
            default=1,
            help='Number of processes to run')
        ligand_imports.add_arguments(parser)

    def handle(self, *args, **options):
        # Refuse before the long contact-network pass, not after it.
        if not options['skip_ligand_import']:
            ligand_imports.check_deliveries(options, ligand_imports.COMMANDS)
        try:
            self.logger.info('CREATING ALL INTERACTIONS')
            self.prepare_input(options['proc'], self.pdbs)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)
        # Ligand interactions: the same imports build_all runs; a failing
        # import raises, so the command exits non-zero.
        ligand_imports.run(options)
        self.logger.info('COMPLETED ALL INTERACTIONS')

    def main_func(self, positions, iteration,count,lock):
        pdbs = self.pdbs
        while count.value<len(pdbs):
            with lock:
                pdb = pdbs[count.value]
                count.value +=1
            try:
                # compute_interactions(pdb, True)
                # The contact network only: the receptor x peptide pairs are imported
                # by ligand_imports.run (import_schrodinger_peptides) in handle().
                compute_interactions(pdb, do_interactions=True, do_peptide_ligand=False, save_to_db=True)
            except:
                print('Issue making interactions for',pdb)
