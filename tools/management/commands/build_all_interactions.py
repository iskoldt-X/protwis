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
            "and the ligand interactions by importing the Schrodinger deliveries, as "
            "build_all does (Engine 1 for HET anchors, Engine 2 for \"pep\" chains). "
            "The maps both imports read are built and both imports dry-run before the "
            "contact network, and they import after it; do not change the deliveries "
            "while it runs.")

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
        # Ligand interactions: the same imports build_all runs; a failing
        # import raises, so the command exits non-zero. Refuse before the long
        # contact-network pass, not after it: the deliveries must exist, the
        # maps be built and both dry runs pass first. The contact network writes only
        # interacting_residue_pair and interaction, which neither import reads
        # or writes, so what the dry runs found in the database still holds
        # after it. The deliveries are read again by the imports, so they must
        # not change while this command runs.
        before, imports = ligand_imports.split(ligand_imports.steps(options))
        ligand_imports.check_deliveries(options, [c for c, _o in before + imports])
        ligand_imports.run(before)
        try:
            self.logger.info('CREATING ALL INTERACTIONS')
            self.prepare_input(options['proc'], self.pdbs)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)
        ligand_imports.run(imports)
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
                # in handle() (import_schrodinger_peptides, via ligand_imports.run).
                compute_interactions(pdb, do_interactions=True, do_peptide_ligand=False, save_to_db=True)
            except:
                print('Issue making interactions for',pdb)
