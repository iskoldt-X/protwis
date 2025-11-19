from build.management.commands.base_build import Command as BaseBuild

from django.db.models import F,Q


from protein.models import Protein, ProteinSegment, ProteinFamily, Species
from alignment.models import ClassSimilarity, ClassSimilarityTie, ClassSimilarityType, ClassRepresentativeSpecies

from common.alignment import Alignment

from collections import OrderedDict
import os
import sys
import logging

from datetime import datetime, date
import time



starttime = datetime.now()
logger = logging.getLogger('class_similarity')
hdlr = logging.FileHandler('./logs/class_similarity.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)

from protein.models import CLASSLESS_PARENT_GPCR_SLUGS

matrix_header_max_length = 20

build_date = date.today()

import warnings
warnings.filterwarnings("ignore")

class Command(BaseBuild):
    help = 'Build cross-class similarity and identity matrices and stores them into GPCRdb database purging first the existing similarity and identity entries. Also stores the most similar and the highest identity GPCR pairs.'
    

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--verbose', help='Prints progress in stdout.', default=False, action='store_true')
        parser.add_argument('--print_results', help='Prints results in stdout.', default=False, action='store_true')
        parser.add_argument('--limit',type=int, help='For testing purposes. Use only any indicated number of GPCRs per class and per batch. Default 10.', default=False, action='store')
       
    def get_parent_gpcr_families(self,exclude_classless_artificial_class=True,include_classless_natural_classes=True):
        parent_family = ProteinFamily.objects.get(slug='000') 
        parent_gpcr_families = ProteinFamily.objects.filter(parent_id=parent_family.pk, slug__startswith='0').exclude(pk=parent_family.pk)
        if exclude_classless_artificial_class:
            for slug in CLASSLESS_PARENT_GPCR_SLUGS:
                parent_gpcr_families = parent_gpcr_families.exclude(slug__startswith=slug)
        if include_classless_natural_classes:
            classless_protein_families = self.get_classless_bottom_protein_families()
        parent_gpcr_families = list(parent_gpcr_families)+classless_protein_families
        return sorted(parent_gpcr_families,key=lambda f: (int(f.slug.split('_')[0])))
    
    def get_human_species(self):
        return Species.objects.get(common_name__iexact='Human')
    
    def get_yeast_species(self):
        return Protein.objects.filter(entry_name__iendswith='_yeast')[0].species
    
    def __slug_tree_branch(self, slug_parts,slug_tree_dict):
        #print('hola0:',slug_tree_dict,slug_parts)

        if len(slug_parts) == 1:
            slug_tree_dict[slug_parts[0]] = None
            #print('hola1:',slug_tree_dict)
            return slug_tree_dict
        else:
            if slug_parts[0] in slug_tree_dict:
                slug_subtree_dict = slug_tree_dict[slug_parts[0]]
                if slug_subtree_dict is None:
                    slug_subtree_dict = {}
            else:
                slug_subtree_dict = {}
            slug_tree_dict[slug_parts[0]] = self.__slug_tree_branch(slug_parts[1:],slug_subtree_dict)
            #print('hola3:',slug_tree_dict)
        return slug_tree_dict

    def __sort_slug_tree_branch(self, slug_tree_dict):
        slug_tree_ordered_dict = OrderedDict()
        for slug in sorted(sorted(slug_tree_dict.keys(),key = lambda x: int(x[1:])),key = lambda x: x[0]):
            slug_subtree_dict = slug_tree_dict[slug]
            if slug_subtree_dict is not None:
                slug_tree_ordered_dict[slug] = self.__sort_slug_tree_branch(slug_subtree_dict)
            else:
                slug_tree_ordered_dict[slug] = None
        return slug_tree_ordered_dict
    
    def __parse_slug_tree_(self, slug_tree_dict,slug_list_list,slug_list):
        for slug,subtree in slug_tree_dict.items():
            slug_list.append(slug)
            if subtree is not None:
                self.__parse_slug_tree_(subtree,slug_list_list,slug_list)
            else:
                slug_list_list.append(slug_list.copy())
            slug_list.pop()
    
    def get_classless_bottom_protein_families(self):
        classless_parent_gpcrs_slugs_list = sorted(sorted(CLASSLESS_PARENT_GPCR_SLUGS,key = lambda x: int(x[1:])),key = lambda x: x[0])
        parent_gpcr_families = ProteinFamily.objects.filter(slug__startswith=classless_parent_gpcrs_slugs_list[0])
        
        for slug in classless_parent_gpcrs_slugs_list[1:]:
            parent_gpcr_families = parent_gpcr_families.filter(slug__startswith=slug)
        
        slug_2_family_dict = {}
        for f in parent_gpcr_families:
            slug_2_family_dict[f.slug] = f
        family_slug_tree_dict = {}
        for slug in slug_2_family_dict.keys():
            slug_parts = slug.split('_')
            self.__slug_tree_branch(slug_parts,family_slug_tree_dict)
        family_slug_tree_ordered_dict = self.__sort_slug_tree_branch(family_slug_tree_dict)
        slug_list_list = []
        slug_list = []
        self.__parse_slug_tree_(family_slug_tree_ordered_dict,slug_list_list,slug_list)
        classless_bottom_slugs_list = ['_'.join(slug_list) for slug_list in slug_list_list]
        return [slug_2_family_dict[slug] for slug in classless_bottom_slugs_list]
    
    def filter_out_non_species_parent_gpcr_families(self,parent_gpcr_families,species):
        """ Filters out parent GPCR families as a list of ProteinFamily objects that belong to a species.
            parent_gpcr_families: a list of protein.ProteinFamily objects
            species: protein.Species object
        """
        new_parent_gpcr_families_slugs_set = set()
        species2 = species
        try:
            species_iterator = iter(species)
        except TypeError as te:
             species2 = [species]
        parent_gpcr_families_slugs = []
        for family in parent_gpcr_families:
            parent_gpcr_families_slugs.append(family.slug)

        for slug in parent_gpcr_families_slugs:
            q = Protein.objects.annotate(family_slug=F('family__slug')).filter(family_slug__startswith=slug,species__in=species2)
            if q.exists():
                new_parent_gpcr_families_slugs_set.add(slug)
        return [f for f in parent_gpcr_families if f.slug in new_parent_gpcr_families_slugs_set]
        
    
    def filter_out_non_human_parent_gpcr_families(self,parent_gpcr_families):
        return self.filter_out_non_species_parent_gpcr_families(parent_gpcr_families,self.get_human_species())

    def filter_out_non_yeast_parent_gpcr_families(self,parent_gpcr_families):
        return self.filter_out_non_species_parent_gpcr_families(parent_gpcr_families,self.get_yeast_species())
    
    def handle(self, *args, **options):
        initial_step1 = 380 #If alignment fails, please, set this to a lower value
        initial_step2 = 380 #If alignment fails, please, set this to a lower value
        start_time = time.time()

        parent_families = self.get_parent_gpcr_families(exclude_classless_artificial_class=True,include_classless_natural_classes=True)
        human_parent_gpcr_families = self.filter_out_non_human_parent_gpcr_families(parent_families)
        human_species = self.get_human_species()

        ClassRepresentativeSpecies.custom_objects.truncate_table() #custom_objects is a custom manager, see ClassRepresentativeSpecies model definition
        class_representative_species_list = []
        for gpcr_family in human_parent_gpcr_families:
            class_representative_species_list.append(ClassRepresentativeSpecies(protein_family = gpcr_family, species = human_species))
        ClassRepresentativeSpecies.objects.bulk_create(class_representative_species_list)

        yeast_parent_gpcr_families = self.filter_out_non_yeast_parent_gpcr_families(parent_families)
        yeast_non_human_parent_gpcr_families_set = set(yeast_parent_gpcr_families) - set(human_parent_gpcr_families)
        yeast_non_human_parent_gpcr_families = [gpcr_family for gpcr_family in yeast_parent_gpcr_families if gpcr_family in yeast_non_human_parent_gpcr_families_set]
        yeast_species = self.get_yeast_species()

        class_representative_species_list = []
        for gpcr_family in yeast_non_human_parent_gpcr_families:
            class_representative_species_list.append(ClassRepresentativeSpecies(protein_family = gpcr_family, species = yeast_species))
        ClassRepresentativeSpecies.objects.bulk_create(class_representative_species_list)


        gpcr_segments = ProteinSegment.objects.filter(Q(proteinfamily='GPCR') & (Q(slug__regex='TM[1-7]') | Q(slug='H8')))

        cross_class_similarities = OrderedDict()
        proteins = None
        i = 1

        step1=int(initial_step1) #If alignment fails, please, set this to a lower value
        step2=int(initial_step2) #If alignment fails, please, set this to a lower value
        step_halved = False 


        human_parent_gpcr_families_protein = {}
        human_parent_gpcr_families_protein_num = {}
        for gpcr_class in human_parent_gpcr_families:
            gpcr_class_proteins = Protein.objects.all().annotate(family_slug=F('family__slug')).filter(species=human_species,family_slug__startswith=gpcr_class.slug)
            gpcr_class_proteins = gpcr_class_proteins.exclude(accession=None).order_by('family_slug','entry_name')
            human_parent_gpcr_families_protein[gpcr_class] = list(gpcr_class_proteins)
            human_parent_gpcr_families_protein_num[gpcr_class] = len(gpcr_class_proteins)

        yeast_non_human_parent_gpcr_families_protein = {}
        yeast_non_human_parent_gpcr_families_protein_num = {} 
        for gpcr_class in yeast_non_human_parent_gpcr_families:
            gpcr_class_proteins = Protein.objects.all().annotate(family_slug=F('family__slug')).filter(species=yeast_species,family_slug__startswith=gpcr_class.slug)
            gpcr_class_proteins = gpcr_class_proteins.exclude(accession=None).order_by('family_slug','entry_name')
            yeast_non_human_parent_gpcr_families_protein[gpcr_class] = list(gpcr_class_proteins)
            yeast_non_human_parent_gpcr_families_protein_num[gpcr_class] = len(gpcr_class_proteins)

        selected_parent_gpcr_families = human_parent_gpcr_families + yeast_non_human_parent_gpcr_families
        selected_parent_gpcr_families_protein = {}
        selected_parent_gpcr_families_protein_num = {}
        for gpcr_class in selected_parent_gpcr_families:
            if gpcr_class in human_parent_gpcr_families_protein and gpcr_class in yeast_non_human_parent_gpcr_families_protein:
                selected_parent_gpcr_families_protein[gpcr_class] = human_parent_gpcr_families_protein[gpcr_class] + yeast_non_human_parent_gpcr_families_protein[gpcr_class]
                selected_parent_gpcr_families_protein_num[gpcr_class] = human_parent_gpcr_families_protein_num[gpcr_class] + yeast_non_human_parent_gpcr_families_protein_num[gpcr_class]
            elif gpcr_class in human_parent_gpcr_families_protein:
                selected_parent_gpcr_families_protein[gpcr_class] = human_parent_gpcr_families_protein[gpcr_class]
                selected_parent_gpcr_families_protein_num[gpcr_class] = human_parent_gpcr_families_protein_num[gpcr_class]
            elif gpcr_class in yeast_non_human_parent_gpcr_families_protein:
                selected_parent_gpcr_families_protein[gpcr_class] = yeast_non_human_parent_gpcr_families_protein[gpcr_class]
                selected_parent_gpcr_families_protein_num[gpcr_class] = yeast_non_human_parent_gpcr_families_protein_num[gpcr_class]


        for gpcr_class in selected_parent_gpcr_families[:-1]:
            while_loop_continue = False
            while True:
                cross_class_similarities2 = OrderedDict()
                for clim in range(0,selected_parent_gpcr_families_protein_num[gpcr_class],step1):
                    
                    gpcr_class_proteins = selected_parent_gpcr_families_protein[gpcr_class]
                    gpcr_class_proteins = gpcr_class_proteins[clim:clim+step1]
                    if options['limit']:
                        gpcr_class_proteins = gpcr_class_proteins[:options['limit']]
                    
                    for gpcr_class2 in selected_parent_gpcr_families[i:]:
                        
                        for clim2 in range(0,selected_parent_gpcr_families_protein_num[gpcr_class2],step2):
                            if options['verbose']:
                                if options['limit']:
                                    print(gpcr_class,"from:"+str(clim+1),"to:"+str(clim+options['limit']),"(of:{})".format(selected_parent_gpcr_families_protein_num[gpcr_class]),\
                                    'vs',gpcr_class2,"from:"+str(clim2+1),"to:"+str(clim2+options['limit']),"(of:{})".format(selected_parent_gpcr_families_protein_num[gpcr_class2]))
                                else:
                                    print(gpcr_class,"from:"+str(clim+1),"to:"+str(clim+step1),"(of:{})".format(selected_parent_gpcr_families_protein_num[gpcr_class]),\
                                    'vs',gpcr_class2,"from:"+str(clim2+1),"to:"+str(clim2+step2),"(of:{})".format(selected_parent_gpcr_families_protein_num[gpcr_class2]))
                            gpcr_class2_proteins = selected_parent_gpcr_families_protein[gpcr_class2]
                            gpcr_class2_proteins = gpcr_class2_proteins[clim2:clim2+step2]

                            if options['limit']:
                                gpcr_class2_proteins = gpcr_class2_proteins[:options['limit']]
                            proteins = gpcr_class_proteins + gpcr_class2_proteins

                            cs_alignment = Alignment()
                            cs_alignment.load_proteins(proteins)
                            cs_alignment.load_segments(gpcr_segments)
                            build_alignment_return_value = cs_alignment.build_alignment()
                            if build_alignment_return_value == "Too large":
                                print('Alignment too large. Retrying...', file=sys.stderr)
                                while_loop_continue = True
                                break
                            cs_alignment.remove_non_generic_numbers_from_alignment() 
                            cs_alignment.calculate_similarity_matrix()

                            entry_name_2_position = {}
                            pos = 0
                            for protein in cs_alignment.proteins:
                                entry_name_2_position[protein.protein.entry_name] = pos
                                pos += 1


                            max_similarity = float('-inf')
                            i_similarity_value = float('-inf')
                            max_identity = 0
                            s_identity_value = 0
                            for protein1 in gpcr_class_proteins:
                                for protein2 in gpcr_class2_proteins:
                                    pos_s = entry_name_2_position[protein1.entry_name]
                                    similarity_value = int(cs_alignment.similarity_matrix[protein2.entry_name]['values'][pos_s][0])
                                    pos_i = entry_name_2_position[protein2.entry_name]
                                    identity_value = int(cs_alignment.similarity_matrix[protein1.entry_name]['values'][pos_i][0])
                                    if identity_value > max_identity or (identity_value == max_identity and similarity_value > i_similarity_value):
                                        max_identity = identity_value
                                        i_protein_tuple = (protein1,protein2)
                                        i_similarity_value = similarity_value
                                        i_protein_tuple_w = []
                                    elif identity_value == max_identity and similarity_value == i_similarity_value:
                                        i_protein_tuple_w.append((protein1,protein2))
                                        # text = protein1.entry_name+' vs '+protein2.entry_name+' with highest identity tied with '+' vs '.join([p.entry_name for p in i_protein_tuple])+'.'
                                        # print('WARNING: '+text, file=sys.stderr)
                                        # logger.warning(text)
                                            
                                    if  similarity_value > max_similarity or (similarity_value == max_similarity and identity_value > s_identity_value):
                                        max_similarity = similarity_value
                                        s_protein_tuple = (protein1,protein2)
                                        s_identity_value = identity_value
                                        s_protein_tuple_w = []
                                    elif similarity_value == max_similarity and identity_value == s_identity_value:
                                        s_protein_tuple_w.append((protein1,protein2))
                                        # text = protein1.entry_name+' vs '+protein2.entry_name+' with highest similarity tied with '+' vs '.join([p.entry_name for p in s_protein_tuple])+'.'
                                        # print('WARNING: '+text, file=sys.stderr)
                                        # logger.warning(text)
                                        
                            if gpcr_class2 in cross_class_similarities2:
                                if cross_class_similarities2[gpcr_class2]['identity'] < max_identity or (cross_class_similarities2[gpcr_class2]['identity'] == max_identity 
                                                                                                and cross_class_similarities2[gpcr_class2]['i_similarity_value'] < i_similarity_value):
                                    cross_class_similarities2[gpcr_class2]['identity'] = max_identity
                                    cross_class_similarities2[gpcr_class2]['identity_gpcr_pair'] = i_protein_tuple
                                    cross_class_similarities2[gpcr_class2]['i_similarity_value']  = i_similarity_value
                                    cross_class_similarities2[gpcr_class2]['identity_gpcr_pair_w'] = []
                                elif (cross_class_similarities2[gpcr_class2]['identity'] == max_identity
                                                                                                and cross_class_similarities2[gpcr_class2]['i_similarity_value'] == i_similarity_value):
                                    cross_class_similarities2[gpcr_class2]['identity_gpcr_pair_w'] += i_protein_tuple_w
                                    # text = 'vs_'.join([p.entry_name for p in cross_class_similarities2[gpcr_class2]['identity_gpcr_pair']])+' with highest identity tied with '+' vs '.join([p.entry_name for p in i_protein_tuple])+'.'
                                    # print('WARNING: '+text, file=sys.stderr)
                                    # logger.warning(text)
                                if cross_class_similarities2[gpcr_class2]['similarity'] < max_similarity or (cross_class_similarities2[gpcr_class2]['similarity'] == max_similarity 
                                                                                                and cross_class_similarities2[gpcr_class2]['s_identity_value'] < s_identity_value):
                                    cross_class_similarities2[gpcr_class2]['similarity'] = max_similarity
                                    cross_class_similarities2[gpcr_class2]['similarity_gpcr_pair'] = s_protein_tuple
                                    cross_class_similarities2[gpcr_class2]['s_identity_value']  = s_identity_value
                                    cross_class_similarities2[gpcr_class2]['similarity_gpcr_pair_w'] = []
                                elif (cross_class_similarities2[gpcr_class2]['similarity'] == max_similarity
                                                                                                and cross_class_similarities2[gpcr_class2]['s_identity_value'] == s_identity_value):
                                    cross_class_similarities2[gpcr_class2]['similarity_gpcr_pair_w'] += s_protein_tuple_w

                                    
                            else:
                                # cross_class_similarities2[gpcr_class2] = {'identity':None,'similarity':None,
                                #              'identity_gpcr_pair': None,'similarity_gpcr_pair': None}
                                cross_class_similarities2[gpcr_class2] = {}
                                cross_class_similarities2[gpcr_class2]['identity'] = max_identity
                                cross_class_similarities2[gpcr_class2]['identity_gpcr_pair'] = i_protein_tuple
                                cross_class_similarities2[gpcr_class2]['i_similarity_value']  = i_similarity_value
                                cross_class_similarities2[gpcr_class2]['identity_gpcr_pair_w'] = i_protein_tuple_w
                                cross_class_similarities2[gpcr_class2]['similarity'] = max_similarity
                                cross_class_similarities2[gpcr_class2]['similarity_gpcr_pair'] = s_protein_tuple
                                cross_class_similarities2[gpcr_class2]['s_identity_value']  = s_identity_value
                                cross_class_similarities2[gpcr_class2]['similarity_gpcr_pair_w'] = s_protein_tuple_w

                        if while_loop_continue:
                            break
                    del proteins
                    del gpcr_class_proteins
                    if while_loop_continue:
                        break
                    selected_parent_gpcr_families_protein[gpcr_class][clim:clim+step1] = [None for i1 in range(clim,clim+step1)]
                if while_loop_continue:
                    if options['verbose']: print("Halving step1 and step2...") 
                    step1 = step1 // 2
                    step2 = step2 // 2
                    step_halved = True
                    while_loop_continue = False
                    if options['verbose']: print("Retrying last class similarity computation with the new steps...")
                    continue
                elif step_halved:
                    print("Restoring initial step1 and step2...")
                    step1=initial_step1
                    step2=initial_step2
                    step_halved = False
                cross_class_similarities[gpcr_class] = cross_class_similarities2
                del selected_parent_gpcr_families_protein[gpcr_class]
                for gpcr_class2 in cross_class_similarities2:
                    for type in ('identity','similarity'):
                        if len(cross_class_similarities2[gpcr_class2][type+'_gpcr_pair_w']) > 0:
                            text = ' vs '.join([p.entry_name for p in cross_class_similarities2[gpcr_class2][type+'_gpcr_pair']])+' with highest '+type+' tied with '+\
                            ', '.join([' vs '.join([p.entry_name for p in protein_tuple]) for protein_tuple in cross_class_similarities2[gpcr_class2][type+'_gpcr_pair_w']])+'.'
                            print('WARNING: '+text, file=sys.stderr)
                            logger.warning(text)
                i += 1
                break    

        # for key in  cross_class_similarities:
        #   print(key,cross_class_similarities[key])
        if options['print_results']:
            print("[Similarity_protein_pairs]")
            for gpcr_class in cross_class_similarities:      
                for gpcr_class2 in cross_class_similarities[gpcr_class]:
                    protein1,protein2 = cross_class_similarities[gpcr_class][gpcr_class2]['similarity_gpcr_pair']
                    print(';'.join((gpcr_class.name,gpcr_class2.name,protein1.entry_name,protein2.entry_name)))
            print("[End Similarity_protein_pairs]")
            print("[Identity_protein_pairs]")
            for gpcr_class in cross_class_similarities:      
                for gpcr_class2 in cross_class_similarities[gpcr_class]:
                    protein1,protein2 = cross_class_similarities[gpcr_class][gpcr_class2]['identity_gpcr_pair']
                    print(';'.join((gpcr_class.name,gpcr_class2.name,protein1.entry_name,protein2.entry_name)))
            print("[End Identity_protein_pairs]")

            cross_class_similarity_matrix = OrderedDict()
            i = 0
            for gpcr_class in selected_parent_gpcr_families:
                j = 0
                row = []
                for gpcr_class2 in selected_parent_gpcr_families:
                    if gpcr_class == gpcr_class2:
                        row.append('-')
                        j += 1
                        continue
                    if i > j:
                        value_type = 'similarity'
                    else:
                        value_type = 'identity'
                    if gpcr_class in cross_class_similarities:
                        if gpcr_class2 in cross_class_similarities[gpcr_class]:
                            row.append(str(cross_class_similarities[gpcr_class][gpcr_class2][value_type]))
                        else:
                            row.append(str(cross_class_similarities[gpcr_class2][gpcr_class][value_type]))
                    else:
                        row.append(str(cross_class_similarities[gpcr_class2][gpcr_class][value_type]))
                    j += 1
                cross_class_similarity_matrix[gpcr_class.name] = row
                i += 1
            print("[Matrix]")
            headers = cross_class_similarity_matrix.keys()
            print(';'+';'.join(list(headers)))
            for key in  cross_class_similarity_matrix:
                print(("{};"+';'.join(["{}" for i in range(0,len(cross_class_similarity_matrix[key]))])).format(key,*cross_class_similarity_matrix[key]))
            print("[End Matrix]")
            if options['verbose']: 
                print("[Human readable Similarity_protein_pairs]")
                for gpcr_class in cross_class_similarities:      
                    for gpcr_class2 in cross_class_similarities[gpcr_class]:
                        protein1,protein2 = cross_class_similarities[gpcr_class][gpcr_class2]['similarity_gpcr_pair']
                        print('{:{w}s} vs {:{w}s}'.format(gpcr_class.name[:matrix_header_max_length],gpcr_class2.name[:matrix_header_max_length],\
                                                        w=matrix_header_max_length),'{:>{w}s} vs {:>{w}s}'.format(protein1.entry_name,protein2.entry_name,\
                                                                                                                w=25))
                print("[End Human readable Similarity__pairs]")
                print("[Human readable Identity_protein_pairs]")
                for gpcr_class in cross_class_similarities:      
                    for gpcr_class2 in cross_class_similarities[gpcr_class]:
                        protein1,protein2 = cross_class_similarities[gpcr_class][gpcr_class2]['identity_gpcr_pair']
                        print('{:{w}s} vs {:{w}s}'.format(gpcr_class.name[:matrix_header_max_length],gpcr_class2.name[:matrix_header_max_length],\
                                                        w=matrix_header_max_length),'{:>{w}s} vs {:>{w}s}'.format(protein1.entry_name,protein2.entry_name,\
                                                                                                                w=25))
                print("[End Human readable Identity__pairs]")

                print("[Human readable Matrix]")
                headers = cross_class_similarity_matrix.keys()
                print('; '+'; '.join(list(headers)))
                for key in  cross_class_similarity_matrix:
                    nkey = key[:matrix_header_max_length]
                    print(("{:{header_max_length}s} "+' '.join(["{: >3}" for i in range(0,len(cross_class_similarity_matrix[key]))]))\
                        .format(nkey,*cross_class_similarity_matrix[key],header_max_length=matrix_header_max_length))
                print("[End Human readable Matrix]")

        if options['verbose']: print('Truncating old table...')
        ClassSimilarity.custom_objects.truncate_table() #custom_objects is a custom manager, see ClassSimilarity model definition
        if options['verbose']: print('Saving new table...')
        class_similarity_list = []
        for gpcr_class in cross_class_similarities:      
            for gpcr_class2 in cross_class_similarities[gpcr_class]:
                s_protein1,s_protein2 = cross_class_similarities[gpcr_class][gpcr_class2]['similarity_gpcr_pair']
                i_protein1,i_protein2 = cross_class_similarities[gpcr_class][gpcr_class2]['identity_gpcr_pair']
                class_similarity_list.append(ClassSimilarity(protein_family1=gpcr_class,protein_family2=gpcr_class2,similar_protein1=s_protein1,similar_protein2=s_protein2,
                                                             ident_protein1=i_protein1,ident_protein2=i_protein2,
                                                             similarity=cross_class_similarities[gpcr_class][gpcr_class2]['similarity'],
                                                               identity=cross_class_similarities[gpcr_class][gpcr_class2]['identity']))
        ClassSimilarity.objects.bulk_create(class_similarity_list)

        class_similarity_ties_list = []
        
        for gpcr_class in cross_class_similarities:      
            for gpcr_class2 in cross_class_similarities[gpcr_class]:
                class_similarity = ClassSimilarity.objects.get(protein_family1=gpcr_class,protein_family2=gpcr_class2)
                for type in ('identity','similarity'):
                    for protein_pair in cross_class_similarities[gpcr_class][gpcr_class2][type+'_gpcr_pair_w']:
                        protein1,protein2 = protein_pair
                        class_similarity_ties_list.append(ClassSimilarityTie(class_similarity=class_similarity,protein1=protein1,protein2=protein2,type=ClassSimilarityType.__dict__[type.upper()]))
        ClassSimilarityTie.objects.bulk_create(class_similarity_ties_list)                
        end_time = time.time()
        elapsed_time = end_time - start_time
        if options['verbose']:
            print('Execution time:', elapsed_time, 'seconds')
            print('Done.')

     





             



