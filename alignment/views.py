from django.shortcuts import render, redirect
from django.conf import settings
from django.http import HttpResponse
from django.views.generic import TemplateView
from django.db.models import Q, F, Case, When, IntegerField, Prefetch, Max
from django.db.models.functions import Least, Greatest
from django.core.cache import cache
from django.core.cache import caches
from django.utils.html import escape
from django.core.files.base import File
from django.http import JsonResponse
from django.views import View

try:
    cache_alignment = caches['alignments']
except:
    cache_alignment = cache

from alignment.functions import get_proteins_from_selection
from common import definitions
from common.selection import Selection, SelectionItem
from common.views import AbsTargetSelection, AbsTargetSelectionTable, AbsSegmentSelection, AbsMiscSelection
from structure.functions import BlastSearch
from protwis.context_processors import site_title

# from common.alignment_SITE_NAME import Alignment
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')
from protein.models import Protein, ProteinSegment, ProteinFamily, ProteinSet, Gene
from residue.models import ResidueNumberingScheme, ResiduePositionSet
from alignment.models import ClassSimilarity, ClassSimilarityTie, ClassSimilarityType, ClassRepresentativeSpecies
from seqsign.sequence_signature import SequenceSignature, signature_score_excel
from protein.models import CLASSLESS_PARENT_GPCR_SLUGS
from mapper.views import DataMapperHome
from alignment.models import ReceptorSimilarity
from common.models import WebLink
from ligand.models import Endogenous_GTP
from structure.models import Structure, StructureModel

from collections import OrderedDict, defaultdict
from copy import deepcopy
import hashlib
import inspect
from io import BytesIO
import itertools
import json
import numpy as np
import os
import xlsxwriter
from xlsxwriter.utility import xl_range_abs
import xlrd
import re
import pandas as pd
from string import Template
from sklearn.manifold import TSNE

strain_re = re.compile(r'\bstrain\b', flags=re.I)
class_fungal_re = re.compile(r'(\(Ste2-like)(\s+)(fungal)(\s+)(pheromone\))', flags=re.I)
class_fullname_re = re.compile(r'^(Class\s+\w+)(\s+)(\(.*\))', flags=re.I)
class_prefix_re = re.compile(r'^(Class)\s+', flags=re.I)

# class TargetSelection(AbsTargetSelection):
#     step = 1
#     number_of_steps = 2
#     filter_tableselect = False
#     docs = 'sequences.html#structure-based-alignments'
#     selection_boxes = OrderedDict([
#         ('reference', False),
#         ('targets', True),
#         ('segments', False),
#     ])
#     buttons = {
#         'continue': {
#             'label': 'Continue to next step',
#             'url': '/alignment/segmentselection',
#             'color': 'success',
#         },
#     }

class TargetSelection(AbsTargetSelectionTable):
    step = 1
    number_of_steps = 2
    filter_tableselect = False
    docs = 'sequences.html#structure-based-alignments'
    title = "SELECT RECEPTORS"
    description = 'Select receptors in the table (below) or browse the classification tree (right). You can select entire' \
        + ' families or individual receptors.\n\nOnce you have selected all your receptors, click the green button.'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Next',
            'onclick': "submitSelection('/alignment/segmentselection');",
            'color': 'success',
        },
    }

class TargetSelectionGprotein(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    psets = False
    filters = True
    filter_gprotein = True

    docs = 'sequences.html#structure-based-alignments'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/alignment/segmentselectiongprot',
            'color': 'success',
        },
    }
    try:
        if ProteinFamily.objects.filter(slug="100_001").exists():
            ppf = ProteinFamily.objects.get(slug="100_001")
            pfs = ProteinFamily.objects.filter(parent=ppf.id)
            ps = Protein.objects.filter(family=ppf)

            tree_indent_level = []
            action = 'expand'
            # remove the parent family (for all other families than the root of the tree, the parent should be shown)
            del ppf
    except Exception as e:
        pass

class TargetSelectionArrestin(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    psets = False
    filters = True
    filter_gprotein = True

    docs = 'sequences.html#structure-based-alignments'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/alignment/segmentselectionarrestin',
            'color': 'success',
        },
    }
    try:
        if ProteinFamily.objects.filter(slug="200_000").exists():
            ppf = ProteinFamily.objects.get(slug="200_000")
            pfs = ProteinFamily.objects.filter(parent=ppf.id)
            ps = Protein.objects.filter(family=ppf)

            tree_indent_level = []
            action = 'expand'
            # remove the parent family (for all other families than the root of the tree, the parent should be shown)
            del ppf
    except Exception as e:
        pass

class SegmentSelection(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show alignment',
            'url': '/alignment/render',
            'color': 'success',
        },
    }

class SegmentSelectionGprotein(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
    description = 'Select sequence segments in the middle column for G proteins. You can expand every structural element and select individual' \
        + ' residues by clicking on the down arrows next to each helix, sheet or loop.\n\n You can select the full sequence or show all structured regions at the same time.\n\nSelected segments will appear in the' \
        + ' right column, where you can edit the list.\n\nOnce you have selected all your segments, click the green' \
        + ' button.'

    template_name = 'common/segmentselection.html'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show alignment',
            'url': '/alignment/render',
            'color': 'success',
        },
    }

    position_type = 'gprotein'
    rsets = ResiduePositionSet.objects.filter(name__in=['Gprotein Barcode', 'YM binding site']).prefetch_related('residue_position')

    ss = ProteinSegment.objects.filter(partial=False, proteinfamily='Alpha').prefetch_related('generic_numbers')
    ss_cats = ss.values_list('category').order_by('category').distinct('category')


class SegmentSelectionArrestin(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
    description = 'Select sequence segments in the middle column for beta and visual arrestins. You can expand every structural element and select individual' \
        + ' residues by clicking on the down arrows next to each helix, sheet or loop.\n\n You can select the full sequence or show all structured regions at the same time.\n\nSelected segments will appear in the' \
        + ' right column, where you can edit the list.\n\nOnce you have selected all your segments, click the green' \
        + ' button.'

    template_name = 'common/segmentselection.html'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show alignment',
            'url': '/alignment/render',
            'color': 'success',
        },
    }

    position_type = 'arrestin'

    ## Add some Arrestin specific positions
    rsets = ResiduePositionSet.objects.filter(name__in=['Arrestin interface']).prefetch_related('residue_position')

    ## ProteinSegment for different proteins
    ss = ProteinSegment.objects.filter(partial=False, proteinfamily='Arrestin').prefetch_related('generic_numbers')
    ss_cats = ss.values_list('category').order_by('category').distinct('category')


class BlastSearchInput(AbsMiscSelection):
    step = 1
    number_of_steps = 1
    docs = 'sequences.html#similarity-search-blast'
    title = 'BLAST search'
    description = 'Enter a sequence into the text box and press the green button.'
    buttons = {
        'continue': {
            'label': 'BLAST',
            'onclick': 'document.getElementById(\'form\').submit()',
            'color': 'success',
        },
    }
    selection_boxes = {}
    blast_input = True


class BlastSearchResults(TemplateView):
    """
    An interface for blast similarity search of the input sequence.
    """
    template_name="blast/blast_search_results.html"

    def post(self, request, *args, **kwargs):

        if 'human' in request.POST.keys():
            blast = BlastSearch(blastdb=os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'protwis_human_blastdb']), top_results=50)
            blast_out = blast.run(request.POST['input_seq'])
        else:
            blast = BlastSearch(top_results=50)
            blast_out = blast.run(request.POST['input_seq'])

        context = {}
        context['results'] = [(Protein.objects.get(pk=x[0]), x[1]) for x in blast_out]
        context["input"] = request.POST['input_seq']

        return render(request, self.template_name, context)


def render_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)
    if simple_selection == False or not simple_selection.targets:
        return redirect("/alignment/targetselection")

    # create an alignment object
    a = Alignment()

    # load data from selection into the alignment
    # only show wildtype protein entries if selection type is family
    if len([t for t in simple_selection.targets if t.type=='family'])>0:
        a.load_proteins_from_selection(simple_selection, only_wildtype=True)
    else:
        a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    key = "ALIGNMENT_" + a.get_hash()
    return_html = cache_alignment.get(key)

    if return_html==None:
        # build the alignment data matrix
        check = a.build_alignment()
        if check == 'Too large':
            return render(request, 'alignment/error.html', {'proteins': len(a.proteins), 'residues':a.number_of_residues_total})

        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()
        num_of_sequences = len(a.proteins)
        num_residue_columns = len(a.positions) + len(a.segments)
        # segment_headers_to_hide = [seg for seg, posis in a.segments.items() if len(posis)==0]

        return_html = render(request, 'alignment/alignment.html', {'a': a, 'num_of_sequences': num_of_sequences,
            'num_residue_columns': num_residue_columns})

    cache_alignment.set(key, return_html, 60*60*24*7) #set alignment cache one week

    return return_html

def render_family_alignment(request, slug):
    # create an alignment object
    a = Alignment()

    # fetch proteins and segments
    proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt')

    if len(proteins)>50 and len(slug.split("_"))<4:
        # If alignment is going to be too big, only pick human.
        proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt', species__latin_name='Homo sapiens')

    if slug.startswith('100'):

        gsegments = definitions.G_PROTEIN_SEGMENTS

        preserved = Case(*[When(slug=pk, then=pos) for pos, pk in enumerate(gsegments['Full'])])
        segments = ProteinSegment.objects.filter(slug__in=gsegments['Full'], partial=False).order_by(preserved)

    elif slug.startswith('200'):
        arrsegments = definitions.ARRESTIN_SEGMENTS

        preserved = Case(*[When(slug=pk, then=pos) for pos, pk in enumerate(arrsegments['Full'])])
        segments = ProteinSegment.objects.filter(slug__in=arrsegments['Full'], partial=False).order_by(preserved)

    else:
        segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR')
        if len(proteins)>50:
            # if a lot of proteins, exclude some segments
            segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(slug__in=['N-term','C-term'])
        if len(proteins)>200:
            # if many more proteins exluclude more segments
            segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(slug__in=['N-term','C-term']).exclude(category='loop')

    protein_ids = []
    for p in proteins:
        protein_ids.append(p.pk)
    protein_list = ','.join(str(x) for x in sorted(protein_ids))

    #create unique proteins_id
    segments_ids = []
    for s in segments:
        segments_ids.append(s.slug)
    segments_list = ','.join(str(x) for x in sorted(segments_ids))

    # Store proteins and segments as selection to enable Fasta/Excel/CSV downloads
    selection = Selection()
    for prot in proteins:
        selection.add('targets', 'protein', SelectionItem('protein', prot))
    for segment in segments:
        selection.add('segments', 'protein_segment', SelectionItem('protein_segment', segment))
    request.session['selection'] = selection.exporter()

    s = str(protein_list+"_"+segments_list)
    key = "ALIGNMENT_"+hashlib.md5(s.encode('utf-8')).hexdigest()
    return_html = cache_alignment.get(key)
    if return_html==None:
        # load data into the alignment
        a.load_proteins(proteins)
        a.load_segments(segments)

        # build the alignment data matrix
        a.build_alignment()

        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()

        num_of_sequences = len(a.proteins)
        num_residue_columns = len(a.positions) + len(a.segments)

        return_html = render(request, 'alignment/alignment.html', {'a': a, 'num_of_sequences': num_of_sequences,
        'num_residue_columns': num_residue_columns})

    #update it if used
    cache_alignment.set(key,return_html, 60*60*24*7) #set alignment cache one week

    return return_html

def render_fasta_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)

    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # load data from selection into the alignment
    a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    # build the alignment data matrix
    a.build_alignment()

    response = render(request, 'alignment/alignment_fasta.html', context={'a': a}, content_type='text/fasta')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.fasta"
    return response

def render_fasta_family_alignment(request, slug):
    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # fetch proteins and segments
    proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt')
    segments = ProteinSegment.objects.filter(partial=False)

    # load data into the alignment
    a.load_proteins(proteins)
    a.load_segments(segments)

    # build the alignment data matrix
    a.build_alignment()

    response = render(request, 'alignment/alignment_fasta.html', context={'a': a}, content_type='text/fasta')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.fasta"
    return response

def render_csv_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)

    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # load data from selection into the alignment
    a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    # build the alignment data matrix
    a.build_alignment()

    # calculate consensus sequence + amino acid and feature frequency
    a.calculate_statistics()

    response = render(request, 'alignment/alignment_csv.html', context={'a': a}, content_type='text/csv')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.csv"
    return response

# Excel download based on seq. signature tool
def render_alignment_excel(request):

    # Grab all targets
    targets = request.session.get('selection', False)

    # create placeholder seq signature
    signature = SequenceSignature()
    signature.setup_alignments_from_selection(targets, targets)

    # calculate the signture
    signature.calculate_signature()
    signature.calculate_zscales_signature()

    outstream = BytesIO()
    wb = xlsxwriter.Workbook(outstream, {'in_memory': True})

    # Sequence alignment of targets
    signature.prepare_excel_worksheet(
        wb,
        'Alignment',
        'positive',
        'alignment'
    )

    # Residue properties stats
    signature.prepare_excel_worksheet(
        wb,
        'Property_conservation',
        'positive',
        'features'
    )
    # Z-scales
    signature.zscales_excel(
        wb,
        "Z-scales",
        'positive'
    )
    wb.close()
    outstream.seek(0)
    response = HttpResponse(
        outstream.read(),
        content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.xlsx"

    return response

render_class_similarity_csv_matrix_urls = ['human_only_without_classless', 'human_only_with_classless','all_without_classless','all_with_classless'] # look this up to see where it is used.

def retrieve_class_similarity_matrix(output_type='html',classless=False,human_only=False):
    cross_class_similarities = OrderedDict()
    class_representative_species_list = ClassRepresentativeSpecies.objects.all().prefetch_related()
    if human_only:
        class_representative_species_list = class_representative_species_list.filter(species__common_name__iexact='human')
        human_classes_set = set()
        for class_representative_species in class_representative_species_list:
            human_classes_set.add(class_representative_species.protein_family)
    if output_type == 'html':
        gpcr_family_name_2_class_representative_species = {}
        for class_representative_species in class_representative_species_list:
            gpcr_family_name = class_representative_species.protein_family.name
            species_name = class_representative_species.species.common_name
            if species_name.lower() != 'human':
                if strain_re.search(species_name):
                    species_name = '<i>'+escape(class_representative_species.species.latin_name)+'</i> <br>'+escape(species_name)
                    gpcr_family_name_2_class_representative_species[gpcr_family_name] = species_name
                else:
                    species_name = class_representative_species.species.latin_name
                    gpcr_family_name_2_class_representative_species[gpcr_family_name] = '<i>'+escape(species_name)+'</i>'
            else:
                gpcr_family_name_2_class_representative_species[gpcr_family_name] = escape(species_name)

    ties_list = ClassSimilarityTie.objects.all().select_related('protein1','protein2').values('class_similarity_id',
                                                                    'protein1__entry_name','protein2__entry_name',
                                                                    'protein1__name','protein2__name','type').order_by('protein1__name','protein2__name')

    class_similarity_id_2_ties = {}
    for tie in ties_list:
        class_similarity_id = tie['class_similarity_id']
        if class_similarity_id not in class_similarity_id_2_ties:
            class_similarity_id_2_ties[class_similarity_id] = []
        class_similarity_id_2_ties[class_similarity_id].append(tie)

    class_sim = ClassSimilarity.objects.all().select_related('protein_family1','protein_family2')
    if human_only:
        class_sim = class_sim.filter(protein_family1__in=human_classes_set,protein_family2__in=human_classes_set)


    if not classless:
        for slug in list(CLASSLESS_PARENT_GPCR_SLUGS):
            class_sim = class_sim.exclude(protein_family1__slug__startswith=slug)
            class_sim = class_sim.exclude(protein_family2__slug__startswith=slug)
    else:
        # NOT USED: code for getting ClassSimilarity.id for pairs containing a classless GPCR
        # classless_sim1 = classless_sim2 = class_sim
        # for i,slug in enumerate(list(CLASSLESS_PARENT_GPCR_SLUGS)):
        #     classless_sim1 = classless_sim1.filter(protein_family1__slug__startswith=slug)
        #     classless_sim2 = classless_sim2.filter(protein_family2__slug__startswith=slug)
        #     if i == 0:
        #         classless_sim = classless_sim1 | classless_sim2
        #     else:
        #         classless_sim = classless_sim | classless_sim1 | classless_sim2
            
        # classless_sim_ids = set(classless_sim.values_list('id',flat=True))
        ProteinFamily.objects.filter()
        for i,slug in enumerate(list(CLASSLESS_PARENT_GPCR_SLUGS)):
            protein_family_0 = ProteinFamily.objects.filter(slug__startswith=slug)
            if i == 0:
                protein_family = protein_family_0
            else:
                protein_family = protein_family | protein_family_0
        classless_protein_family_names = set(protein_family .values_list('name',flat=True))

    class_pairs = class_sim.values('id','protein_family1__slug','protein_family1__name',
                                                                   'protein_family2__slug','protein_family2__name',
                                                                    'similar_protein1__entry_name','similar_protein2__entry_name',
                                                                    'similar_protein1__name','similar_protein2__name',
                                                                    'ident_protein1__entry_name','ident_protein2__entry_name',
                                                                    'ident_protein1__name','ident_protein2__name','identity','similarity')\
                                                                    .order_by('protein_family1__name','protein_family2__name')

    for class_pair in class_pairs:
        gpcr_class1_name = class_pair['protein_family1__name']
        gpcr_class2_name = class_pair['protein_family2__name']
        if gpcr_class1_name not in cross_class_similarities:
            cross_class_similarities[gpcr_class1_name] = OrderedDict()
        if gpcr_class2_name not in cross_class_similarities[gpcr_class1_name]:
            cross_class_similarities[gpcr_class1_name][gpcr_class2_name] = {}
        cross_class_sim = cross_class_similarities[gpcr_class1_name][gpcr_class2_name]
        cross_class_sim['identity'] = class_pair['identity']
        cross_class_sim['similarity'] = class_pair['similarity']
        if output_type == 'html':
            name_type = 'name'
            ident_protein1_name = Protein(name=class_pair['ident_protein1__'+name_type]).short().replace('<i>','').replace('</i>','')
            ident_protein2_name = Protein(name=class_pair['ident_protein2__'+name_type]).short().replace('<i>','').replace('</i>','')
            similar_protein1_name = Protein(name=class_pair['similar_protein1__'+name_type]).short().replace('<i>','').replace('</i>','')
            similar_protein2_name = Protein(name=class_pair['similar_protein2__'+name_type]).short().replace('<i>','').replace('</i>','')
            name_type = 'entry_name'
            ident_protein1_entry_name = class_pair['ident_protein1__'+name_type]
            ident_protein2_entry_name = class_pair['ident_protein2__'+name_type]
            similar_protein1_entry_name = class_pair['similar_protein1__'+name_type]
            similar_protein2_entry_name = class_pair['similar_protein2__'+name_type]
            cross_class_sim['identity_gpcr_pair_entry_name'] = (ident_protein1_entry_name,ident_protein2_entry_name)
            cross_class_sim['identity_gpcr_pair_w_entry_name'] = []
            cross_class_sim['similarity_gpcr_pair'] = (similar_protein1_name,similar_protein2_name)
            cross_class_sim['similarity_gpcr_pair_entry_name'] = (similar_protein1_entry_name,similar_protein2_entry_name)
            cross_class_sim['similarity_gpcr_pair_w_entry_name'] = []
        else:
            name_type = 'entry_name'
            ident_protein1_name = class_pair['ident_protein1__'+name_type]
            ident_protein2_name = class_pair['ident_protein2__'+name_type]
            similar_protein1_name = class_pair['similar_protein1__'+name_type]
            similar_protein2_name = class_pair['similar_protein2__'+name_type]
        cross_class_sim['identity_gpcr_pair'] = (ident_protein1_name,ident_protein2_name)
        cross_class_sim['identity_gpcr_pair_w'] = []
        cross_class_sim['similarity_gpcr_pair'] = (similar_protein1_name,similar_protein2_name)
        cross_class_sim['similarity_gpcr_pair_w'] = []


        if class_pair['id'] in class_similarity_id_2_ties:
            for tie in class_similarity_id_2_ties[class_pair['id']]:
                if tie['type'] == ClassSimilarityType.IDENTITY:
                    type = 'identity'
                elif tie['type'] == ClassSimilarityType.SIMILARITY:
                    type = 'similarity'
                
                if output_type == 'html':
                    name_type = 'name'
                    protein1 = Protein(name=tie['protein1__'+name_type]).short().replace('<i>','').replace('</i>','')
                    protein2 = Protein(name=tie['protein2__'+name_type]).short().replace('<i>','').replace('</i>','')
                    name_type = 'entry_name'
                    protein1_entry_name= tie['protein1__'+name_type]
                    protein2_entry_name = tie['protein2__'+name_type]
                    cross_class_sim[type+'_gpcr_pair_w_entry_name'].append((protein1_entry_name,protein2_entry_name))
                else:
                    name_type = 'entry_name'
                    protein1 = tie['protein1__'+name_type]
                    protein2 = tie['protein2__'+name_type]

                cross_class_sim[type+'_gpcr_pair_w'].append((protein1,protein2))

    selected_parent_gpcr_families_names1 = OrderedDict()
    selected_parent_gpcr_families_names2 = OrderedDict()
    for key in cross_class_similarities:
        selected_parent_gpcr_families_names1[key] = None
        for key2 in cross_class_similarities[key]:
            selected_parent_gpcr_families_names2[key2] = None
    for key in selected_parent_gpcr_families_names2:
        selected_parent_gpcr_families_names1[key] = None

    selected_parent_gpcr_families_names = list(selected_parent_gpcr_families_names1.keys())
    if classless:
        selected_p_gpcr_families_names1 = [name for name in selected_parent_gpcr_families_names if name not in classless_protein_family_names]
        selected_p_gpcr_families_names1.sort()
        selected_p_gpcr_families_names2 = [name for name in selected_parent_gpcr_families_names if name in classless_protein_family_names]
        selected_p_gpcr_families_names2.sort()
        selected_parent_gpcr_families_names = selected_p_gpcr_families_names1 + selected_p_gpcr_families_names2
        del selected_p_gpcr_families_names1
        del selected_p_gpcr_families_names2

    else:
        selected_parent_gpcr_families_names.sort()
    
    cross_class_similarity_matrix = OrderedDict()
    if output_type=='xls' or output_type=='xlsx':
        cross_class_similarities_xls = OrderedDict()
    i = 0
    for key in selected_parent_gpcr_families_names:
        row = []
        j = 0
        if output_type=='xls' or output_type=='xlsx':
            cross_class_similarities_xls[key] = OrderedDict()
        for key2 in selected_parent_gpcr_families_names:
            if key == key2:
                row.append(['-','-',''])
                continue
            if i > j:
                type = 'similarity'
            else:
                type = 'identity'
            
            if key in cross_class_similarities:
                if key2 in cross_class_similarities[key]:
                    sim = cross_class_similarities[key][key2]
                    invert_pairs = False
                else:
                    sim = cross_class_similarities[key2][key]
                    invert_pairs = True
            else:
                sim = cross_class_similarities[key2][key]
                invert_pairs = True
            if invert_pairs:
                pairs = [(sim[type+'_gpcr_pair'][1],sim[type+'_gpcr_pair'][0])]
                pairsw = [(p[1],p[0]) for p in sim[type+'_gpcr_pair_w']]
            else:
                pairs = [sim[type+'_gpcr_pair']]
                pairsw = sim[type+'_gpcr_pair_w']
            if len(pairsw) > 0:
                pairs += pairsw
            if output_type == 'html':
                if invert_pairs:
                    pairs_entry_name = [(sim[type+'_gpcr_pair_entry_name'][1],sim[type+'_gpcr_pair_entry_name'][0])]
                    pairsw_entry_name = [(p[1],p[0]) for p in sim[type+'_gpcr_pair_w_entry_name']]
                else:
                    pairs_entry_name = [sim[type+'_gpcr_pair_entry_name']]
                    pairsw_entry_name = sim[type+'_gpcr_pair_w_entry_name']    
                if len(pairsw_entry_name) > 0:
                    pairs_entry_name += pairsw_entry_name
                pairs_entry_name_order_list = [(i,p) for i,p in enumerate(pairs_entry_name)]
                pairs_entry_name_order_list.sort(key=lambda p : p[1][1])
                pairs_entry_name_order_list.sort(key=lambda p : p[1][0])
                pairs_entry_name = [e[1] for e in pairs_entry_name_order_list]
                pairs = [pairs[e[0]] for e in pairs_entry_name_order_list]
                row.append([str(sim[type]),str(sim[type]//10),pairs,pairs_entry_name])
                species_name = gpcr_family_name_2_class_representative_species[key] 
            else:
                pairs.sort(key=lambda p : p[1])
                pairs.sort(key=lambda p : p[0])
                row.append([str(sim[type]),str(sim[type]//10),pairs])
            if output_type=='xls' or output_type=='xlsx':
                sim2 = deepcopy(sim)
                if type == 'identity':
                    rec12 = pairs
                elif type == 'similarity':
                    rec12 = [sim2[type+'_gpcr_pair']] + sim2[type+'_gpcr_pair_w']
                    rec12.sort(key=lambda p : p[1])
                    rec12.sort(key=lambda p : p[0])
                sim2[type+'_gpcr_pair'] = rec12[0]
                if len(rec12) > 1:
                    sim2[type+'_gpcr_pair_w'] = rec12[1:]
                cross_class_similarities_xls[key][key2] = sim2

            j += 1
        
        name = key.replace('<i>','').replace('</i>','')
        if output_type == 'html':
            name = class_fullname_re.sub(r'\1<br>\3', name)
            name = class_fungal_re.sub(r'(Ste2 \3<br>\5', name)
            cross_class_similarity_matrix[key] = {'name':name,'values':row, 'species':species_name}
        else:
            cross_class_similarity_matrix[key] = {'name':name,'values':row}
        i += 1
    if output_type=='xls' or output_type=='xlsx':
        return (cross_class_similarity_matrix,selected_parent_gpcr_families_names,cross_class_similarities_xls)
    else:
        return (cross_class_similarity_matrix,selected_parent_gpcr_families_names)

def render_class_similarity_matrix(request):
    r_human_only_without_classless = retrieve_class_similarity_matrix(classless=False,human_only=True)
    r_human_only_with_classless = retrieve_class_similarity_matrix(classless=True,human_only=True)
    r_all_without_classless = retrieve_class_similarity_matrix(classless=False,human_only=False)
    r_all_with_classless = retrieve_class_similarity_matrix(classless=True,human_only=False)

    list_tmp = [r_human_only_without_classless,r_human_only_with_classless,r_all_without_classless,r_all_with_classless]

    h_list = ['Human only', 'Human only (including Classless)','All (including non-human)','All (including non-human & Classless)']
    csv_list = [{'classless': '0','human_only':'1'},{'classless': '1','human_only':'1'},{'classless': '0','human_only':'0'},{'classless': '1','human_only':'0'}]

    return render(request, 'class_similarity/matrix.html', {'csv_list':csv_list,'h_list': h_list, 'p_list': [r[1] for r in list_tmp],'m_list': [r[0] for r in list_tmp]})



def render_class_similarity_csv_matrix(request):
    classless = True
    human_only = False
    try:
        classless_g = int(request.GET.get('classless', 1))
        if classless_g == 0:
            classless = False
    except ValueError as e:
        pass
    try:
        human_only_g = int(request.GET.get('human_only', 0))
        if abs(human_only_g) > 0:
            human_only = True
    except ValueError as e:
        pass

    r = retrieve_class_similarity_matrix(output_type='csv',classless=classless,human_only=human_only)

    response = render(request, 'class_similarity/matrix_csv.html', {'p': r[1],'m': r[0]}, content_type='text/csv')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_similaritymatrix.csv"
    return response

def render_class_similarity_xlsx_matrix(request):
    classless = True
    human_only = False
    try:
        classless_g = int(request.GET.get('classless', 1))
        if classless_g == 0:
            classless = False
    except ValueError as e:
        pass
    try:
        human_only_g = int(request.GET.get('human_only', 0))
        if abs(human_only_g) > 0:
            human_only = True
    except ValueError as e:
        pass

    r = retrieve_class_similarity_matrix(output_type='xlsx',classless=classless,human_only=human_only)
    m, selected_parent_gpcr_families_names, cross_class_similarities = r

    xlsx_output = BytesIO()
    
    workbook = xlsxwriter.Workbook(xlsx_output, {'in_memory': True})

    axis = True

    v_axis_row = 1    #starting row index of vertical axis
    v_axis_col = 0    #column index of identity axis
    h_axis_row = 0    #row index of horitzontal axis
    h_axis_col = 1    #starting column index of indentity axis
    row_header_row = v_axis_row + 1    #starting row index of the row header
    row_header_col = v_axis_col + 1    #column index of the row header
    col_header_row = h_axis_row + 1    #starting column index of the column header
    col_header_col = h_axis_col + 1    #row index of the column header
    if not axis:
        row_header_row = 1    #starting row index of the row header
        row_header_col = 0    #column index of the row header
        col_header_row = 0    #starting column index of the column header
        col_header_col = 1    #row index of the column header
    
    
    first_data_row = row_header_row    #first data cell row index
    first_data_col = col_header_col    #first data cell column index

    n_header = len(m.keys())

    last_data_row = n_header - 1 + first_data_row    #last data cell row index
    last_data_col = n_header - 1 + first_data_col    #last data cell column index

    #format

    table_headers = [v['name'] for p, v in m.items()]
    worksheet_info = workbook.add_worksheet('Info')
    worksheet_matrix = workbook.add_worksheet('Matrix')
    worksheet_matrix_p = workbook.add_worksheet('Matrix_pairs')
    worksheet_pairs_id = workbook.add_worksheet('Pairs_Id')
    worksheet_pairs_sim = workbook.add_worksheet('Pairs_Sim')

    table_headers_format = workbook.add_format()            #row and column headers format
    table_h_axis_format = workbook.add_format()             #horitzontal axis format
    table_v_axis_format = workbook.add_format()             #vertical axis format
    table_headers_numbers_format = workbook.add_format()    #row and column headers format for cells containing numeric data
    numbers_format = workbook.add_format()                  #data cell format for cells containing numeric data                                        
    cell_format = workbook.add_format()                     #data cell format for cells containing non-numeric data
    pairs_header_format = workbook.add_format()             #receptor pairs column headers format
    info_header_format = workbook.add_format()              #info column headers format

    table_headers_format.set_bold(True)
    table_headers_numbers_format.set_bold(True)
    table_h_axis_format.set_bold(True)
    table_v_axis_format.set_bold(True)
    pairs_header_format.set_bold(True)
    info_header_format.set_bold(True)

    table_headers_format.set_align('vcenter')
    table_headers_numbers_format.set_align('vcenter')
    table_h_axis_format.set_align('vcenter')
    table_v_axis_format.set_align('vcenter')
    numbers_format.set_align('vcenter')
    cell_format.set_align('vcenter')

    table_headers_numbers_format.set_align('center')
    table_h_axis_format.set_align('center')
    table_v_axis_format.set_align('center')
    numbers_format.set_align('center')

    numbers_format.set_border(1)
    numbers_format.set_border_color('#DDDDDD')
    table_v_axis_format.set_rotation(90)
    table_headers_format.set_text_wrap()
    table_headers_numbers_format.set_text_wrap()

    worksheet_pairs_id.freeze_panes(1, 0)
    worksheet_pairs_sim.freeze_panes(1, 0)

    #Info
    worksheet_info.write_row(0 , 0, ['Tabs', 'Description'], info_header_format)
    worksheet_info.write_row(1 , 0, ['Matrix', 'Cross-class GPCR similarity matrix, where cross-class similarity/identity is defined as the highest similarity/identity among each pair of receptors.'])
    worksheet_info.write_row(2 , 0, ['Matrix_Pairs', 'Highest cross-class GPCR similarity/identity pairs of receptors matrix.'])
    worksheet_info.write_row(3 , 0, ['Pairs_Id', 'List of cross-class identity of each pair of receptors with the highest identity.'])
    worksheet_info.write_row(4 , 0, ['Pairs_Sim', 'List of cross-class similarity of each pair of receptors with the highest similarity.'])


    try:
        # Requires xlsxwriter > 3.0.8
        worksheet_info.autofit()
    except Exception as e:
        worksheet_info.set_column(0, 0, 10)
        worksheet_info.set_column(1, 1, 110)
        pass

    #axis titles

    if axis:
            worksheet_matrix.merge_range(h_axis_row , h_axis_col, h_axis_row, last_data_col, 'Identity (%)',table_h_axis_format)
            worksheet_matrix_p.merge_range(h_axis_row , h_axis_col, h_axis_row, last_data_col, 'Identity',table_h_axis_format)
    if axis:
            worksheet_matrix.merge_range(v_axis_row , v_axis_col, last_data_row, v_axis_col, 'Similarity (%)',table_v_axis_format)
            worksheet_matrix_p.merge_range(v_axis_row , v_axis_col, last_data_row, v_axis_col, 'Similarity',table_v_axis_format)

    #column header
    worksheet_matrix.write_row(col_header_row , col_header_col, table_headers, table_headers_numbers_format)
    worksheet_matrix_p.write_row(col_header_row, col_header_col, table_headers, table_headers_format)

    row = first_data_row
    for p, v in m.items():
        worksheet_matrix.write(row,row_header_col,v['name'],table_headers_numbers_format)
        col = first_data_col 
        for v_col in v['values']:
            cell_value = v_col[0]
            if cell_value  != '-':
                cell_value = int(cell_value )
            worksheet_matrix.write(row,col,cell_value,numbers_format)
            col += 1
        row += 1


    row = first_data_row
    for p, v in m.items():
        worksheet_matrix_p.write(row,row_header_col,v['name'],table_headers_format)
        col = first_data_col
        for v_col in v['values']:
            cell_value = ''
            for counter, pair in enumerate(v_col[2]):
                cell_value += pair[0]+' vs '+pair[1]
                if counter != len(v_col[2]) - 1:
                    cell_value += ':'
            if row == col:
                cell_value = '-'
            worksheet_matrix_p.write(row,col,cell_value,cell_format)
            col += 1
        row += 1
    try:
        # Requires xlsxwriter > 3.0.8

        worksheet_matrix.autofit()
        worksheet_matrix_p.autofit()
    except Exception as e:
        if axis:
            first_col = v_axis_col
        else:
            first_col = row_header_col
        worksheet_matrix_p.set_column(first_col, last_data_col, 12)
        worksheet_matrix_p.set_column(first_data_col, last_data_col, 81)
        pass
    worksheet_matrix.set_column(0, last_data_col, 12)
    for row in range(0,last_data_col + 1):
        worksheet_matrix.set_row(row, 50)
    worksheet_matrix.freeze_panes(2, 2)
    worksheet_matrix_p.freeze_panes(2, 2)
    multi_range_list = []
    # compute cell multi-range for color scale

    # identity
    multi_range_list = []
    for i in range(1,last_data_row - first_data_row + 1):
        first_row = first_data_row
        last_row = first_data_row + i - 1 # do not color scale the diagonal
        first_col = last_col = first_data_col + i
        multi_range_list.append(xl_range_abs(first_row,first_col,last_row,last_col))
    worksheet_matrix.conditional_format(multi_range_list[0],
                                {'type': '2_color_scale','min_color': '#FFFFFF',
                                'max_color': '#999999','multi_range': ' '.join(multi_range_list)})

    # similarity
    multi_range_list = []
    for i in range(0,last_data_row - first_data_row + 1):
        first_row = first_data_row + 1 + i # do not color scale the diagonal
        last_row = last_data_row
        first_col = last_col = first_data_col + i
        multi_range_list.append(xl_range_abs(first_row,first_col,last_row,last_col))
    worksheet_matrix.conditional_format(multi_range_list[0],
                                  {'type': '2_color_scale','min_color': '#FFFFFF',
                                    'max_color': '#999999','multi_range': ' '.join(multi_range_list)})
    worksheet_matrix_p.write_row(col_header_row, col_header_col, table_headers, table_headers_format)


    #Pairs
    pairs_worksheet_dict_list = [
        {'worksheet': worksheet_pairs_id,'values_field_header':'Id (%)','type':'identity'},
        {'worksheet': worksheet_pairs_sim,'values_field_header':'Sim (%)','type':'similarity'}
    ]
    worksheet_pairs_headers = ['Class 1','Class 2','Rec 1','Rec 2']
    for pairs_worksheet_dict in pairs_worksheet_dict_list:
        worksheet = pairs_worksheet_dict['worksheet']
        worksheet_headers = worksheet_pairs_headers + [pairs_worksheet_dict['values_field_header']]
        type = pairs_worksheet_dict['type']
        worksheet.write_row(0 , 0, worksheet_headers,pairs_header_format)
        unique_keys_set = set()
        row = 1
        for gpcr_class1_name in selected_parent_gpcr_families_names:
            # NOT USED: retrieve_class_similarity_matrix() returns a cross_class_similarities with reverse combinations of classes
            # if gpcr_class1_name not in cross_class_similarities:
            #     continue
            gpcr_class1_name_d = class_prefix_re.sub(r'',gpcr_class1_name.replace('<i>','').replace('</i>',''))
            for gpcr_class2_name in selected_parent_gpcr_families_names:
                if gpcr_class2_name == gpcr_class1_name:
                    continue
                # skip duplicated combinations of classes
                unique_list = [gpcr_class1_name,gpcr_class2_name]
                unique_list.sort()
                unique_key = '@'.join(unique_list)
                if unique_key in unique_keys_set:
                    continue
                # NOT USED: retrieve_class_similarity_matrix() returns a cross_class_similarities with reverse combinations of classes
                # if gpcr_class2_name not in cross_class_similarities[gpcr_class1_name]:
                #     continue

                gpcr_class2_name_d = class_prefix_re.sub(r'',gpcr_class2_name.replace('<i>','').replace('</i>',''))
                cross_class_sim = cross_class_similarities[gpcr_class1_name][gpcr_class2_name]
                
                value = cross_class_sim[type]
                rec12 = [cross_class_sim[type+'_gpcr_pair']] + cross_class_sim[type+'_gpcr_pair_w']
                # NOT USED: already sorted in retrieve_class_similarity_matrix()
                # rec12.sort(key=lambda p : p[1])
                # rec12.sort(key=lambda p : p[0])
                rec1 = ':'.join([p[0] for p in rec12])
                rec2 = ':'.join([p[1] for p in rec12])
                worksheet.write_row(row , 0, [gpcr_class1_name_d,gpcr_class2_name_d,rec1,rec2,value])
                row += 1
                unique_keys_set.add(unique_key)
        try:
            # Requires xlsxwriter > 3.0.8

            worksheet.autofit()
        except Exception as e:
            worksheet.set_column(0, 1, 25)
            worksheet.set_column(2, 3, 40)
            worksheet.set_column(4, 4, 6)
            pass
        
        

    workbook.close()
    xlsx_file = File(xlsx_output)
    xlsx_file_size = xlsx_file.size
    
    response = HttpResponse(xlsx_file ,'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
    response['Content-Length'] = xlsx_file_size
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_similaritymatrix.xlsx"
    return response

class ClassesAndCounts(TemplateView):
    template_name = "class_similarity/ClassesAndCounts.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # ---- file load (local to this view) ----
        data_folder = "protein_data"
        file_name = "Ligand type update plus sense column.xlsx"
        file_path = os.path.join(settings.DATA_DIR, data_folder, file_name)

        try:
            df = pd.read_excel(file_path)
        except FileNotFoundError:
            context["error"] = f"File not found: {file_path}"
            return context

        df.columns = df.columns.str.strip().str.replace("\n", " ")

        # ---- mapping (local to this view) ----
        mapping = {
            "A":  (["Class A (Rhodopsin)"], "Rhodopsin", "Rhodopsin", None),
            "B1": (["Class B1 (Secretin)"], "Secretin", "Secretin", None),
            "B2": (["Class B2 (Adhesion)"], "Adhesion", "Adhesion", None),
            "C":  (["Class C (Glutamate)"], "Glutamate", "Glutamate", None),
            "D":  (["Class D (Fungal pheromone)"], "- (fungal)", "Fungal pheromone", None),
            "E":  (["Class E (Yeast cAMP)"], "- (yeast)", "Yeast cAMP", "2"),
            "F":  (["Class F (Frizzled)"], "Frizzled/Taste2", "Frizzled", None),
            "T2": (["Class T2 (Taste 2)"], "Frizzled/Taste2", "Taste 2", None),
            "OR": (["Class O1 (fish-like odorant)", "Class O2 (tetrapod specific odorant)"],
                   "Rhodopsin", "Odorant (not olfactory)", None),
            "V?": (["Class V? (Vomeronasal/pheromone?)"],
                   "- (non-functional in human)", "Vomeronasal or pheromone?", "5"),
            "Cl": (["Other GPCRs"], "-", "Classless", None),
        }
        nonhuman_symbols = {"D", "E", "V?"}

        # ---- build class/counts table ----
        table_data = []
        for symbol, (class_names, grafs_family, display_name, fixed_sensory) in mapping.items():
            class_df = df[df["Class"].isin(class_names)]

            receptor_families_count = class_df["Receptor family"].nunique()
            members_total_count = class_df["GPCRs (UniProt)"].nunique()
            non_sensory_count = (class_df["Sense"] == "Non-sensory").sum()

            if fixed_sensory is not None:
                sensory_str = fixed_sensory
            else:
                sensory_counts = []
                for sense_type in ["Vision", "Taste", "Odorant"]:
                    c = (class_df["Sense"] == sense_type).sum()
                    if c > 0:
                        sensory_counts.append(f"{c} {sense_type}")
                sensory_str = ", ".join(sensory_counts) or 0

            orphan_count = (class_df["Sense"] == "Unknown").sum()

            entries = {
                "receptor_families": sorted(class_df["Receptor family"].dropna().unique().tolist()),
                "members_total": class_df[
                    ["GPCRs (Gene name)", "Receptor family", "Ligand type", "Sense"]
                ].drop_duplicates().to_dict(orient="records"),
                "members_non_sensory": class_df[class_df["Sense"] == "Non-sensory"][
                    ["GPCRs (Gene name)", "Receptor family", "Ligand type", "Sense"]
                ].drop_duplicates().to_dict(orient="records"),
                "members_sensory": class_df[class_df["Sense"].isin(["Vision", "Taste", "Odorant"])][
                    ["GPCRs (Gene name)", "Receptor family", "Ligand type", "Sense"]
                ].drop_duplicates().to_dict(orient="records"),
                "orphan": class_df[class_df["Sense"] == "Unknown"][
                    ["GPCRs (Gene name)", "Receptor family", "Ligand type", "Sense"]
                ].drop_duplicates().to_dict(orient="records"),
            }

            table_data.append({
                "symbol": symbol,
                "class_name": display_name,
                "grafs_family": grafs_family,
                "receptor_families_count": int(receptor_families_count),
                "members_total_count": int(members_total_count),
                "members_non_sensory_count": int(non_sensory_count),
                "members_sensory_count": sensory_str,
                "orphan_count": int(orphan_count),
                "entries": entries,
            })

        context["human_table_data"] = [r for r in table_data if r["symbol"] not in nonhuman_symbols]
        context["nonhuman_table_data"] = [r for r in table_data if r["symbol"] in nonhuman_symbols]
        return context


class Classification(TemplateView):
    template_name = "class_similarity/Classification.html"

    # mapping from “Excel Class” → (symbol, display_name)
    CLASS_MAPPING = {
        "Class A (Rhodopsin)":                      ("A",  "Rhodopsin"),
        "Class B1 (Secretin)":                      ("B1", "Secretin"),
        "Class B2 (Adhesion)":                      ("B2", "Adhesion"),
        "Class C (Glutamate)":                      ("C",  "Glutamate"),
        "Class F (Frizzled)":                       ("F",  "Frizzled"),
        "Class T2 (Taste 2)":                       ("T2", "Taste 2"),
        "Class O1 (fish-like odorant)":             ("O1", "Fish-like olfactory receptors"),
        "Class O2 (tetrapod specific odorant)":     ("O2", "Tetrapod-specific olfactory receptors"),
        "Other GPCRs":                              ("Cl", "Classless"),
        # non-human:
        "Class D (Fungal pheromone)":               ("D1", "Fungal pheromone"),
        "Class E (Yeast cAMP)":                     ("E",  "Yeast cAMP"),
        "Class V? (Vomeronasal/pheromone?)":        ("V?", "Vomeronasal or pheromone?"),
        # you can add V1/V2 etc later if they exist in backbone
    }

    # mapping “Ligand type” → “Ligand type group”
    LIGAND_GROUP_MAP = {
        "Ion receptors":            "Ion receptors",
        "Peptide receptors":        "Polypeptide receptors",
        "Protein receptors":        "Polypeptide receptors",
        "Alicarboxylic acid receptors": "Small molecule receptors",
        "Aminergic receptors":      "Small molecule receptors",
        "Amino acid receptors":     "Small molecule receptors",
        "Lipid receptors":          "Small molecule receptors",
        "Melatonin receptors":      "Small molecule receptors",
        "Nucleotide receptors":     "Small molecule receptors",
        "Odorant receptors":        "Small molecule receptors",
        "Opsin receptors":          "Small molecule receptors",
        "Pheromone receptors":      "Pheromone receptors",
        "Steroid receptors":        "Small molecule receptors",
        "Tastant receptors":        "Small molecule receptors",
        "Mechano-activated receptors": "Tethered peptide receptors",
        "Protease-activated receptors": "Tethered peptide receptors",
        "Orphan receptors":         "Orphan receptors",
        # fall-back for anything not listed:
        #   → we’ll label as "Other / unknown" in code
    }

    SENSORY_SENSES = {"vision", "taste", "odorant"}

    def _load_df(self):
        folder = "protein_data"
        fname = "Ligand type update plus sense column.xlsx"
        path = os.path.join(settings.DATA_DIR, folder, fname)
        df = pd.read_excel(path)
        df.columns = df.columns.str.strip().str.replace("\n", " ")
        return df

    def get_context_data(self, **kwargs):
        ctx = super().get_context_data(**kwargs)

        try:
            df = self._load_df()
        except FileNotFoundError as e:
            ctx["error"] = f"File not found: {e}"
            return ctx

        # only keep backbone columns we care about
        cols_needed = [
            "GPCRs (Gene name)",
            "GPCRs (UniProt)",
            "Class",
            "Receptor family",
            "Ligand type",
            "Sense",
        ]
        df = df[cols_needed].copy()

        # map Excel class → symbol (A, B1, O1, …)
        class_to_symbol = {}
        for excel_cls, (symbol, _name) in self.CLASS_MAPPING.items():
            class_to_symbol[excel_cls] = symbol

        df["Class_symbol"] = df["Class"].map(class_to_symbol)

        # drop rows that don’t map to a symbol (just to be safe)
        df = df.dropna(subset=["Class_symbol"])

        # ---------- 1) Ligand type table ----------
        lt_agg = {}  # ligand_type -> {"group":..., "classes": set([...])}
        for _, row in df.iterrows():
            lt = str(row["Ligand type"]).strip()
            if not lt or lt.lower() == "nan":
                continue
            symbol = row["Class_symbol"]
            group = self.LIGAND_GROUP_MAP.get(lt, "Other / unknown")
            entry = lt_agg.setdefault(lt, {"group": group, "classes": set()})
            # keep the first group if mapping disagrees; or you can assert
            if entry["group"] != group:
                # you *could* log or harmonise here; for now we keep the first
                pass
            entry["classes"].add(symbol)

        # nice ordered list of classes
        class_order = ["A", "B1", "B2", "C", "D1", "D2", "E", "F",
                       "T2", "O1", "O2", "V1", "V2", "Cl"]
        def sort_classes(s):
            return sorted(s, key=lambda x: (class_order.index(x)
                                            if x in class_order else 999, x))

        ligand_type_rows = []
        for lt in sorted(lt_agg.keys(), key=str.lower):
            entry = lt_agg[lt]
            cls_list = sort_classes(entry["classes"])
            ligand_type_rows.append({
                "ligand_type": lt,
                "ligand_group": entry["group"],
                "classes": ", ".join(cls_list),
            })

        # ---------- 2) Receptor families: non-sensory vs sensory ----------

        non_sens_triples = set()       # (class_symbol, family, ligand_type)
        sensory_map = {}              # (family, ligand_type) -> set(classes)

        for _, row in df.iterrows():
            fam = str(row["Receptor family"]).strip()
            if not fam or fam.lower() == "nan":
                continue
            lt = str(row["Ligand type"]).strip()
            symbol = row["Class_symbol"]
            sense = str(row["Sense"]).strip().lower()

            if sense in self.SENSORY_SENSES:
                key = (fam, lt)
                sensory_map.setdefault(key, set()).add(symbol)
            elif sense == "non-sensory":
                non_sens_triples.add((symbol, fam, lt))
            else:
                # "unknown" or anything else – you can decide where to put these;
                # for now we ignore them for the receptor-family tables
                pass

        # non-sensory: Class / Receptor family / Ligand type
        rf_non_rows = [
            {
                "class_symbol": cs,
                "receptor_family": fam,
                "ligand_type": lt,
            }
            for (cs, fam, lt) in sorted(
                non_sens_triples,
                key=lambda t: (class_order.index(t[0])
                               if t[0] in class_order else 999,
                               t[0].lower(), t[1].lower())
            )
        ]

        # sensory: Receptor family / Ligand type / Found in classes
        rf_sens_rows = []
        for (fam, lt), classes in sensory_map.items():
            cls_list = sort_classes(classes)
            rf_sens_rows.append({
                "receptor_family": fam,
                "ligand_type": lt,
                "classes": ", ".join(cls_list),
            })

        rf_sens_rows.sort(key=lambda r: (r["receptor_family"].lower(),
                                         r["ligand_type"].lower()))

        ctx["ligand_type_rows"] = ligand_type_rows
        ctx["rf_non_rows"] = rf_non_rows
        ctx["rf_sens_rows"] = rf_sens_rows

        return ctx


class GPCRBrowser(TemplateView):
    template_name = "class_similarity/GPCRBrowser.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # ---- file load (local to this view) ----
        data_folder = "protein_data"
        file_name = "Ligand type update plus sense column.xlsx"
        file_path = os.path.join(settings.DATA_DIR, data_folder, file_name)

        try:
            df = pd.read_excel(file_path)
        except FileNotFoundError:
            context["error"] = f"File not found: {file_path}"
            return context

        df.columns = df.columns.str.strip().str.replace("\n", " ")

        # ---- mapping (local to this view) ----
        mapping = {
            "A":  (["Class A (Rhodopsin)"], "Rhodopsin", "Rhodopsin", None),
            "B1": (["Class B1 (Secretin)"], "Secretin", "Secretin", None),
            "B2": (["Class B2 (Adhesion)"], "Adhesion", "Adhesion", None),
            "C":  (["Class C (Glutamate)"], "Glutamate", "Glutamate", None),
            "D":  (["Class D (Fungal pheromone)"], "- (fungal)", "Fungal pheromone", None),
            "E":  (["Class E (Yeast cAMP)"], "- (yeast)", "Yeast cAMP", "2"),
            "F":  (["Class F (Frizzled)"], "Frizzled/Taste2", "Frizzled", None),
            "T2": (["Class T2 (Taste 2)"], "Frizzled/Taste2", "Taste 2", None),
            "OR": (["Class O1 (fish-like odorant)", "Class O2 (tetrapod specific odorant)"],
                   "Rhodopsin", "Odorant (not olfactory)", None),
            "V?": (["Class V? (Vomeronasal/pheromone?)"],
                   "- (non-functional in human)", "Vomeronasal or pheromone?", "5"),
            "Cl": (["Other GPCRs"], "-", "Classless", None),
        }
        nonhuman_symbols = {"D", "E", "V?"}

        # ---- placeholder ligand data ----
        # when ready, build ligand-table here from df + mapping
        context["ligand_table_data"] = []
        return context



class ClassificationWheel(TemplateView):
    template_name = 'class_similarity/ClassificationWheel.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # --- Step 1: Load Excel metadata ---

        data_folder = 'protein_data'
        file_name = 'Ligand type update plus sense column.xlsx'
        file_path = os.path.join(settings.DATA_DIR, data_folder, file_name)
        df = pd.read_excel(file_path)

        # Clean column names
        df.columns = df.columns.str.strip().str.replace("\n", " ")

        # Build a lookup dict by UniProt entry name
        meta_lookup = {}
        for _, row in df.iterrows():
            entry = str(row["GPCRs (UniProt)"]).strip()
            if entry:
                meta_lookup[entry] = {
                    "Class": row.get("Class", ""),
                    "Ligand type": row.get("Ligand type", ""),
                    "Receptor family": row.get("Receptor family", ""),
                    "Sense": row.get("Sense", "")
                }

        # --- Step 2: Helper to inject metadata into wheel structure ---
        def enrich_wheel_with_metadata(wheelstructure):
            def recurse(node, current_class=None):
                if isinstance(node, dict):
                    for k, v in node.items():
                        if isinstance(v, dict):
                            # If we're inside a Circle_X, the keys here are actual classes (A, B1, etc.)
                            if k.startswith("Circle_"):
                                recurse(v, current_class=None)  # reset class at start of a circle
                            elif current_class is None and not "EntryName" in v:
                                # This k is the class code (A, B1, etc.)
                                recurse(v, current_class=k)
                            elif "EntryName" in v:
                                entry_code = v["EntryName"]
                                meta = meta_lookup.get(entry_code, {})
                                v.update(meta)

                                # Add the class from one level above (A, B1, etc.)
                                v["Class"] = current_class

                                if "Color" not in v:
                                    v["Color"] = "#FFFFFF"
                                if "Data" not in v:
                                    v["Data"] = ""
                            else:
                                recurse(v, current_class=current_class)
                elif isinstance(node, list):
                    for item in node:
                        recurse(item, current_class=current_class)

            recurse(wheelstructure.get("Data", {}))
            return wheelstructure


        # --- Step 3: Build the wheels ---
        odorant_wheel = DataMapperHome.GenerateGPCRomeDataStructure(data_type="Odorant")
        classic_wheel = DataMapperHome.GenerateGPCRomeDataStructure(data_type="Classic")

        # Inject metadata into both
        updated_odorant = enrich_wheel_with_metadata(odorant_wheel)
        updated_classic = enrich_wheel_with_metadata(classic_wheel)

        # --- Step 4: Pass to template ---
        context['GPCRomeData'] = json.dumps(updated_classic['Data'])
        context['GPCRomeOdorantData'] = json.dumps(updated_odorant['Data'])

        return context

class CrossClassSimilarity(TemplateView):
    template_name = 'class_similarity/CrossClassSimilarity.html'

    # Core fixed classes (original order)
    CLASS_ORDER = [
        "Class A (Rhodopsin)",
        "Class B1 (Secretin)",
        "Class B2 (Adhesion)",
        "Class C (Glutamate)",
        "Class F (Frizzled)",
        "Class O1 (fish-like)",
        "Class O2 (tetrapod specific)",
        "Class T2 (Taste 2)",
    ]

    # Map display names -> top-level family slug codes
    CLASS_CODE_BY_NAME = {
        "Class A (Rhodopsin)": "001",
        "Class B1 (Secretin)": "002",
        "Class B2 (Adhesion)": "003",
        "Class C (Glutamate)": "004",
        "Class F (Frizzled)":  "006",
        "Class O1 (fish-like)": "007",
        "Class O2 (tetrapod specific)": "008",
        "Class T2 (Taste 2)": "009",
    }

    # Five single-protein “Classless” items as separate groups (display order)
    SINGLE_PROTEIN_LABELS = ["GPR107", "GPR137", "TPRA1", "GPR143", "GPR157"]

    # Exact entry_name per label (case-insensitive)
    SINGLE_PROTEIN_ENTRYNAMES = {
        "GPR107": "gp107_human",
        "GPR137": "g137a_human",
        "TPRA1":  "tpra1_human",
        "GPR143": "gp143_human",
        "GPR157": "gp157_human",
    }

    # ---------- helpers ----------
    @staticmethod
    def clean_name(nm: str) -> str:
        if not nm:
            return "-"
        return (nm.replace("receptor", "")
                  .replace("-adrenoceptor", "")
                  .replace("<i>", "").replace("</i>", "")
                  .strip()) or "-"

    @staticmethod
    def primary_gene_of(p: 'Protein') -> str:
        if getattr(p, "primary_genes_self", None):
            return p.primary_genes_self[0].name
        if p.entry_name:
            return p.entry_name.split("_")[0].upper()
        return "-"

    def build_gtop_url(self, wl):
        try:
            return Template(wl.web_resource.url).substitute(index=wl.index)
        except Exception:
            return None

    def pack_hover(self, p: 'Protein') -> dict:
        """Return per-protein metadata for the tooltip."""
        wl = p.gtop_links_self[0] if getattr(p, "gtop_links_self", None) else None
        return {
            "display_name": self.clean_name(p.name),
            "gtopdb_link": self.build_gtop_url(wl) or "",
            "uniprot": p.entry_name or "",
            "gene": self.primary_gene_of(p),
            "gpcrdb_link": f"/protein/{p.entry_name}" if p.entry_name else "",
            "uniprot_link": (f"https://www.uniprot.org/uniprot/{getattr(p, 'accession', '')}"
                             if getattr(p, "accession", None) else ""),
        }

    # ---- resolvers
    def _resolve_family_ids(self, slug_codes):
        qs = ProteinFamily.objects.filter(slug__in=slug_codes).only('id', 'slug', 'name')
        return {f.slug: f.id for f in qs}

    def _fetch_by_entry_names(self, entry_names_lower):
        if not entry_names_lower:
            return {}
        qs = (Protein.objects
              .filter(entry_name__in=entry_names_lower)
              .only('id', 'entry_name', 'name', 'accession'))
        return {(p.entry_name or "").lower(): p for p in qs}

    def _resolve_single_proteins(self):
        """Resolve SINGLE_PROTEIN_LABELS by entry_name first, then primary gene."""
        wanted_lc = {
            lbl: (self.SINGLE_PROTEIN_ENTRYNAMES.get(lbl) or "").lower()
            for lbl in self.SINGLE_PROTEIN_LABELS
        }
        entry_to_label = {en: lbl for lbl, en in wanted_lc.items() if en}

        found = {}
        if entry_to_label:
            by_en = self._fetch_by_entry_names(list(entry_to_label.keys()))
            for en, prot in by_en.items():
                lbl = entry_to_label.get(en)
                if lbl:
                    found[lbl] = prot

        missing = [lbl for lbl in self.SINGLE_PROTEIN_LABELS if lbl not in found]
        if missing:
            wanted_genes = set(missing)
            gene_qs = (
                Protein.objects
                .prefetch_related(
                    Prefetch('genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self')
                )
                .only('id', 'entry_name', 'name', 'accession')
            )
            for p in gene_qs:
                g = self.primary_gene_of(p)
                if g in wanted_genes and g not in found:
                    found[g] = p
        return found

    # ---- main
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # 1) Build display list (no extra/non-human groups)
        base_names   = list(self.CLASS_ORDER)
        single_names = [f"{lab} (Classless)" for lab in self.SINGLE_PROTEIN_LABELS]
        display_names = base_names + single_names

        # 2) Resolve base class families for fast class↔class aggregation
        code_to_famid = self._resolve_family_ids(list(self.CLASS_CODE_BY_NAME.values()))
        name_to_famid = {
            name: code_to_famid[self.CLASS_CODE_BY_NAME[name]]
            for name in self.CLASS_ORDER
            if self.CLASS_CODE_BY_NAME[name] in code_to_famid
        }
        allowed_class_ids = list(name_to_famid.values())

        # 3) Resolve the 5 classless singles
        resolved_singles = self._resolve_single_proteins()
        classless_name_to_protein = {}
        for lab in self.SINGLE_PROTEIN_LABELS:
            key = f"{lab} (Classless)"
            p = resolved_singles.get(lab)
            if p:
                classless_name_to_protein[key] = p

        # 4) Final groups (skip unresolved safely)
        groups = []
        for name in display_names:
            if name in name_to_famid:
                groups.append({"display": name, "kind": "class", "id": name_to_famid[name]})
            elif name in classless_name_to_protein:
                groups.append({"display": name, "kind": "protein", "id": classless_name_to_protein[name].id})
        n = len(groups)

        # ------------------------------ OPTIMIZED AGGREGATION ------------------------------
        # A) class↔class maxima (2 queries)
        base_pairs = (
            ReceptorSimilarity.objects
            .filter(ref_class_id__in=allowed_class_ids, target_class_id__in=allowed_class_ids)
            .annotate(
                pair_a=Least('ref_class_id', 'target_class_id'),
                pair_b=Greatest('ref_class_id', 'target_class_id'),
            )
            .values('pair_a', 'pair_b')
        )
        cc_max = {
            (row['pair_a'], row['pair_b']): (row['max_id'], row['max_sim'])
            for row in base_pairs.annotate(
                max_id=Max('identity'),
                max_sim=Max('similarity')
            )
        }

        # B) class↔protein maxima (1 query)
        protein_ids = [g['id'] for g in groups if g['kind'] == 'protein']
        cp_max = {}
        if allowed_class_ids and protein_ids:
            qs_cp = (
                ReceptorSimilarity.objects
                .filter(
                    (Q(ref_class_id__in=allowed_class_ids, protein_target_id__in=protein_ids)) |
                    (Q(target_class_id__in=allowed_class_ids, protein_ref_id__in=protein_ids))
                )
                .annotate(
                    canon_class_id=Case(
                        When(ref_class_id__in=allowed_class_ids, then='ref_class_id'),
                        default='target_class_id',
                        output_field=IntegerField()
                    ),
                    canon_protein_id=Case(
                        When(protein_target_id__in=protein_ids, then='protein_target_id'),
                        default='protein_ref_id',
                        output_field=IntegerField()
                    ),
                )
                .values('canon_class_id', 'canon_protein_id')
                .annotate(
                    max_id=Max('identity'),
                    max_sim=Max('similarity')
                )
            )
            cp_max = {
                (row['canon_class_id'], row['canon_protein_id']): (row['max_id'], row['max_sim'])
                for row in qs_cp
            }

        # C) protein↔protein maxima (1 query)
        pp_max = {}
        if len(protein_ids) >= 2:
            qs_pp = (
                ReceptorSimilarity.objects
                .filter(protein_ref_id__in=protein_ids, protein_target_id__in=protein_ids)
                .annotate(
                    pair_a=Least('protein_ref_id', 'protein_target_id'),
                    pair_b=Greatest('protein_ref_id', 'protein_target_id'),
                )
                .values('pair_a', 'pair_b')
                .annotate(
                    max_id=Max('identity'),
                    max_sim=Max('similarity')
                )
            )
            pp_max = {
                (row['pair_a'], row['pair_b']): (row['max_id'], row['max_sim'])
                for row in qs_pp
            }

        # 5) Build value-only matrix
        matrix = [[None for _ in range(n)] for _ in range(n)]

        def best_for(a, b, metric):
            if a['kind'] == 'class' and b['kind'] == 'class':
                key = (min(a['id'], b['id']), max(a['id'], b['id']))
                tup = cc_max.get(key)
            elif a['kind'] == 'class' and b['kind'] == 'protein':
                tup = cp_max.get((a['id'], b['id']))
            elif a['kind'] == 'protein' and b['kind'] == 'class':
                tup = cp_max.get((b['id'], a['id']))
            else:
                key = (min(a['id'], b['id']), max(a['id'], b['id']))
                tup = pp_max.get(key)
            if not tup:
                return None
            return tup[0] if metric == 'identity' else tup[1]

        # Collect tie fetch specs so we can pull them in big batches later
        tie_specs = []
        for i in range(n):
            for j in range(n):
                if i == j:
                    matrix[i][j] = None
                    continue
                a, b = groups[i], groups[j]
                metric = 'identity' if i < j else 'similarity'
                best = best_for(a, b, metric)
                matrix[i][j] = {"value": int(best) if best is not None else None,
                                "type": metric,
                                "items": []}
                if best is not None:
                    if a['kind'] == 'class' and b['kind'] == 'class':
                        tie_specs.append(('cc', (min(a['id'], b['id']), max(a['id'], b['id'])), metric, int(best)))
                    elif a['kind'] == 'class' and b['kind'] == 'protein':
                        tie_specs.append(('cp', (a['id'], b['id']), metric, int(best)))
                    elif a['kind'] == 'protein' and b['kind'] == 'class':
                        tie_specs.append(('cp', (b['id'], a['id']), metric, int(best)))
                    else:
                        tie_specs.append(('pp', (min(a['id'], b['id']), max(a['id'], b['id'])), metric, int(best)))

        # 6) Batch-fetch tie rows (kept separate per metric to avoid mixing)
        def bucket_specs(specs):
            buckets = defaultdict(list)
            for kind, ids, metric, best in specs:
                buckets[(kind, metric)].append((ids, best))
            return buckets

        buckets = bucket_specs(tie_specs)

        # Collectors: dict-of-dicts keyed by metric
        cc_ties = {'identity': defaultdict(list), 'similarity': defaultdict(list)}
        cp_ties = {'identity': defaultdict(list), 'similarity': defaultdict(list)}
        pp_ties = {'identity': defaultdict(list), 'similarity': defaultdict(list)}

        # weblink prefetch (GtoPdb)
        gtop_links_qs = WebLink.objects.select_related('web_resource').filter(web_resource__slug='gtop')

        def fetch_cc(metric):
            pairs = buckets.get(('cc', metric), [])
            if not pairs:
                return
            q = Q()
            for (a_id, b_id), best in pairs:
                cond = (Q(ref_class_id=a_id, target_class_id=b_id) |
                        Q(ref_class_id=b_id, target_class_id=a_id))
                cond &= Q(**{metric: best})
                q |= cond
            if not q.children:
                return
            rows = (
                ReceptorSimilarity.objects
                .filter(q)
                .select_related('protein_ref', 'protein_target')
                .only(
                    'identity', 'similarity',
                    'protein_ref__id', 'protein_ref__entry_name', 'protein_ref__name', 'protein_ref__accession',
                    'protein_target__id', 'protein_target__entry_name', 'protein_target__name', 'protein_target__accession',
                    'ref_class_id', 'target_class_id'
                )
                .prefetch_related(
                    Prefetch('protein_ref__genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self'),
                    Prefetch('protein_target__genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self'),
                    Prefetch('protein_ref__web_links',
                             queryset=gtop_links_qs,
                             to_attr='gtop_links_self'),
                    Prefetch('protein_target__web_links',
                             queryset=gtop_links_qs,
                             to_attr='gtop_links_self'),
                )
            )
            for r in rows:
                a = min(r.ref_class_id, r.target_class_id)
                b = max(r.ref_class_id, r.target_class_id)
                cc_ties[metric][(a, b)].append(r)

        def fetch_cp(metric):
            pairs = buckets.get(('cp', metric), [])
            if not pairs:
                return
            q = Q()
            for (cls_id, prot_id), best in pairs:
                cond = (
                    Q(ref_class_id=cls_id, protein_target_id=prot_id) |
                    Q(target_class_id=cls_id, protein_ref_id=prot_id)
                )
                cond &= Q(**{metric: best})
                q |= cond
            if not q.children:
                return
            rows = (
                ReceptorSimilarity.objects
                .filter(q)
                .select_related('protein_ref', 'protein_target')
                .only(
                    'identity', 'similarity',
                    'protein_ref__id', 'protein_ref__entry_name', 'protein_ref__name', 'protein_ref__accession',
                    'protein_target__id', 'protein_target__entry_name', 'protein_target__name', 'protein_target__accession',
                    'ref_class_id', 'target_class_id', 'protein_ref_id', 'protein_target_id'
                )
                .prefetch_related(
                    Prefetch('protein_ref__genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self'),
                    Prefetch('protein_target__genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self'),
                    Prefetch('protein_ref__web_links',
                             queryset=gtop_links_qs,
                             to_attr='gtop_links_self'),
                    Prefetch('protein_target__web_links',
                             queryset=gtop_links_qs,
                             to_attr='gtop_links_self'),
                )
            )
            for r in rows:
                if r.ref_class_id is not None and r.protein_target_id is not None:
                    key = (r.ref_class_id, r.protein_target_id)
                else:
                    key = (r.target_class_id, r.protein_ref_id)
                cp_ties[metric][key].append(r)

        def fetch_pp(metric):
            pairs = buckets.get(('pp', metric), [])
            if not pairs:
                return
            q = Q()
            for (a_id, b_id), best in pairs:
                cond = (
                    Q(protein_ref_id=a_id, protein_target_id=b_id) |
                    Q(protein_ref_id=b_id, protein_target_id=a_id)
                )
                cond &= Q(**{metric: best})
                q |= cond
            if not q.children:
                return
            rows = (
                ReceptorSimilarity.objects
                .filter(q)
                .select_related('protein_ref', 'protein_target')
                .only(
                    'identity', 'similarity',
                    'protein_ref__id', 'protein_ref__entry_name', 'protein_ref__name', 'protein_ref__accession',
                    'protein_target__id', 'protein_target__entry_name', 'protein_target__name', 'protein_target__accession',
                    'protein_ref_id', 'protein_target_id'
                )
                .prefetch_related(
                    Prefetch('protein_ref__genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self'),
                    Prefetch('protein_target__genes',
                             queryset=Gene.objects.filter(position=0),
                             to_attr='primary_genes_self'),
                    Prefetch('protein_ref__web_links',
                             queryset=gtop_links_qs,
                             to_attr='gtop_links_self'),
                    Prefetch('protein_target__web_links',
                             queryset=gtop_links_qs,
                             to_attr='gtop_links_self'),
                )
            )
            for r in rows:
                a = min(r.protein_ref_id, r.protein_target_id)
                b = max(r.protein_ref_id, r.protein_target_id)
                pp_ties[metric][(a, b)].append(r)

        # Execute the 6 batched tie fetches
        for m in ('identity', 'similarity'):
            fetch_cc(m)
            fetch_cp(m)
            fetch_pp(m)

        # 7) Fill items for tooltips (include identity & similarity per row)
        def pack_rows(rows):
            return [{
                "ref":     self.pack_hover(r.protein_ref),
                "target":  self.pack_hover(r.protein_target),
                "identity": r.identity,
                "similarity": r.similarity,
            } for r in rows]

        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                a, b = groups[i], groups[j]
                metric = matrix[i][j]["type"]            # 'identity' for upper, 'similarity' for lower
                val = matrix[i][j]["value"]
                if val is None:
                    continue

                if a['kind'] == 'class' and b['kind'] == 'class':
                    key = (min(a['id'], b['id']), max(a['id'], b['id']))
                    rows = cc_ties[metric].get(key, [])
                elif a['kind'] == 'class' and b['kind'] == 'protein':
                    rows = cp_ties[metric].get((a['id'], b['id']), [])
                elif a['kind'] == 'protein' and b['kind'] == 'class':
                    rows = cp_ties[metric].get((b['id'], a['id']), [])
                else:
                    key = (min(a['id'], b['id']), max(a['id'], b['id']))
                    rows = pp_ties[metric].get(key, [])

                matrix[i][j]["items"] = pack_rows(rows)

        # 8) Context
        context["classes"] = [g["display"] for g in groups]
        context["matrix_json"] = json.dumps(matrix)
        return context


class OrphanSimilarity(TemplateView):
    template_name = 'class_similarity/OrphanSimilarity.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        def clean_gtop_name(nm):
            if not nm:
                return "-"
            s = (nm
                 .replace("receptor", "")
                 .replace("-adrenoceptor", "")
                 .replace("<i>", "").replace("</i>", "")
                 .strip())
            return s or "-"

        orphans_qs = (
            Protein.objects
            .filter(
                parent_id__isnull=True,
                species_id=1,
                family__parent__parent__name__iexact='Orphan receptors',
            )
            .prefetch_related(
                Prefetch('genes',
                         queryset=Gene.objects.filter(position=0),
                         to_attr='primary_genes_self')
            )
            .order_by('entry_name')
        )

        data = []
        for p in orphans_qs:
            gene = (p.primary_genes_self[0].name
                    if getattr(p, 'primary_genes_self', None)
                    else (p.entry_name.split('_')[0].upper() if p.entry_name else "-"))
            data.append({
                "id": p.id,
                "text": clean_gtop_name(p.name),
                "name": clean_gtop_name(p.name),
                "entry_name": p.entry_name,
                "gene": gene,
            })

        context['orphans_select2'] = json.dumps(data)
        return context
    
class SimilarityTopAPI(View):
    """
    GET /class_similarity/api/similarity?ref=<protein_id>

    Rule:
      1) Sort ALL neighbors by similarity (desc) – do NOT use identity.
      2) Compute the cutoff from the 10th **liganded** row (or the last liganded row if <10 exist).
         Keep **all** liganded rows with similarity >= cutoff (ties included).
      3) Also include **orphan** rows with similarity >= the same cutoff.
      4) If there are NO liganded rows, fall back to the 10th overall row (or last) as cutoff.

    Returns each row with:
      Gene (primary M2M gene, position=0), entry_name, links (GPCRdb/UniProt/IUPHAR),
      family/ligand_type/class, similarity & identity, endogenous ligands (id/name) + types.
    """
    def get(self, request):
        from string import Template
        ref_raw = request.GET.get('ref')
        try:
            ref_id = int(ref_raw)
        except (TypeError, ValueError):
            return JsonResponse({"error": "Missing or invalid 'ref' parameter"}, status=400)

        # 1) All neighbors sorted by similarity only (no identity in ordering)
        pairs_qs = (
            ReceptorSimilarity.objects
            .filter(Q(protein_ref_id=ref_id) | Q(protein_target_id=ref_id))
            .annotate(
                other_id=Case(
                    When(protein_ref_id=ref_id, then=F('protein_target_id')),
                    default=F('protein_ref_id'),
                    output_field=IntegerField(),
                )
            )
            .order_by('-similarity')   # <— identity NOT used
            .values('other_id', 'similarity', 'identity')
        )
        pairs = list(pairs_qs)
        if not pairs:
            return JsonResponse({"results": []})

        other_ids_all = [p['other_id'] for p in pairs]

        # 2) For cutoff we only need to know who is orphan. Fetch ligand type name cheaply.
        LT_ORPHAN = 'Orphan receptors'
        lt_map = dict(
            Protein.objects
                   .filter(id__in=other_ids_all)
                   .values_list('id', 'family__parent__parent__name')   # ligand type name
        )
        def is_orphan(pid):
            return (lt_map.get(pid) or '').strip().lower() == LT_ORPHAN.lower()

        liganded_rows = [p for p in pairs if not is_orphan(p['other_id'])]

        # 3) Compute the similarity cutoff from liganded rows; if none, fall back to overall.
        if liganded_rows:
            if len(liganded_rows) >= 10:
                cutoff = liganded_rows[9]['similarity']  # 10th liganded row
            else:
                cutoff = liganded_rows[-1]['similarity'] # last liganded row
        else:
            # No liganded hits: use overall 10th (or last) similarity as cutoff
            cutoff = pairs[min(9, len(pairs) - 1)]['similarity']

        # 4) Keep everyone (liganded + orphan) with similarity >= cutoff
        kept = [p for p in pairs if p['similarity'] >= cutoff]
        kept_ids = [p['other_id'] for p in kept]

        # 5) Fetch full protein meta only for kept ids
        gtop_links_qs = WebLink.objects.select_related('web_resource').filter(web_resource__slug='gtop')
        proteins_qs = (
            Protein.objects
            .filter(id__in=kept_ids)
            .select_related('family__parent__parent__parent')
            .prefetch_related(
                Prefetch('genes',
                         queryset=Gene.objects.filter(position=0),
                         to_attr='primary_genes_self'),
                Prefetch('web_links',
                         queryset=gtop_links_qs,
                         to_attr='gtop_links_self'),
            )
        )
        proteins = {p.id: p for p in proteins_qs}

        # 6) Batch endogenous ligands for kept ids
        endo_qs = (
            Endogenous_GTP.objects
            .filter(receptor_id__in=kept_ids)
            .select_related('ligand', 'ligand__ligand_type')
        )
        endo_by_receptor = {}
        for e in endo_qs:
            if e.ligand:
                endo_by_receptor.setdefault(e.receptor_id, []).append(e)

        # helpers
        def clean_iuphar_name(nm):
            if not nm:
                return "-"
            s = nm.replace("receptor", "").replace("-adrenoceptor", "").replace("<i>", "").replace("</i>", "").strip()
            return s or "-"

        def build_gtop_url(wl):
            try:
                return Template(wl.web_resource.url).substitute(index=wl.index)
            except Exception:
                return None

        sim_map = {p['other_id']: p['similarity'] for p in kept}
        idn_map = {p['other_id']: p['identity']   for p in kept}

        # 7) Build results in similarity-desc order
        results = []
        for row in kept:
            pid = row['other_id']
            p = proteins.get(pid)
            if not p:
                continue

            # Gene: primary (position=0) or fallback to entry code
            gene_name = p.primary_genes_self[0].name if getattr(p, 'primary_genes_self', None) else (
                (p.entry_name.split('_')[0].upper()) if p.entry_name else "-"
            )

            entry_name  = p.entry_name or None
            gpcrdb_link = f"/protein/{entry_name}" if entry_name else "-"
            uniprot_link = f"https://www.uniprot.org/uniprot/{p.accession}" if p.accession else None

            wl_self = p.gtop_links_self[0] if getattr(p, 'gtop_links_self', None) else None
            iuphar_link = build_gtop_url(wl_self) if wl_self else None
            iuphar_name = clean_iuphar_name(p.name)

            family_name = getattr(getattr(p.family, "parent", None), "name", None)
            ligand_type = getattr(getattr(getattr(p.family, "parent", None), "parent", None), "name", None)
            clazz       = getattr(getattr(getattr(getattr(p.family, "parent", None), "parent", None), "parent", None), "name", None)

            # endogenous ligands (dedup by ligand.id)
            lig_items = endo_by_receptor.get(pid, [])
            seen, endo_ligands, lig_types = set(), [], set()
            for e in lig_items:
                lig = e.ligand
                if not lig:
                    continue
                if lig.id not in seen:
                    seen.add(lig.id)
                    endo_ligands.append({"id": lig.id, "name": lig.name})
                if lig.ligand_type:
                    lig_types.add(lig.ligand_type.name)
            endo_type = "<br>".join(sorted(lig_types)) if lig_types else "-"

            results.append({
                "other_id": pid,
                "Gene": gene_name,
                "entry_name": entry_name,
                "gpcrdb_link": gpcrdb_link,
                "uniprot_link": uniprot_link,
                "iuphar_name": iuphar_name,
                "iuphar_link": iuphar_link,
                "family": family_name,
                "ligand_type": ligand_type,
                "class": clazz,
                "similarity": sim_map.get(pid, 0),
                "identity": idn_map.get(pid, 0),
                "endo_ligands": endo_ligands,
                "endo_type": endo_type,
            })

        # already in similarity desc because `kept` was built from `pairs` order
        return JsonResponse({"results": results})

class OrhanSimilarityClustering(TemplateView):
    template_name = 'class_similarity/OrhanSimilarityClustering.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        def clean_gtop_name(nm):
            if not nm:
                return "-"
            s = (nm
                 .replace("receptor", "")
                 .replace("-adrenoceptor", "")
                 .replace("<i>", "").replace("</i>", "")
                 .strip())
            return s or "-"

        orphans_qs = (
            Protein.objects
            .filter(
                parent_id__isnull=True,
                species_id=1,
                family__parent__parent__name__iexact='Orphan receptors',
            )
            .prefetch_related(
                Prefetch(
                    'genes',
                    queryset=Gene.objects.filter(position=0),
                    to_attr='primary_genes_self'
                )
            )
            .order_by('entry_name')
        )

        data = []
        for p in orphans_qs:
            gene = (
                p.primary_genes_self[0].name
                if getattr(p, 'primary_genes_self', None)
                else (p.entry_name.split('_')[0].upper() if p.entry_name else "-")
            )
            data.append({
                "id": p.id,
                "text": clean_gtop_name(p.name),
                "name": clean_gtop_name(p.name),
                "entry_name": p.entry_name,
                "gene": gene,
            })

        context['orphans_select2'] = json.dumps(data)
        return context

class SimilarityEmbeddingAPI(View):
    ORPHAN_LT = "Orphan receptors"

    @staticmethod
    def _truthy(v):
        return str(v).lower() not in ("", "0", "false", "no", "off", "none")

    @staticmethod
    def _lab(en):
        return (en or "").replace("_human", "")

    def get_orphan_ids(self):
        return set(
            Protein.objects
                   .filter(species_id=1,
                           family__parent__parent__name__iexact=self.ORPHAN_LT)
                   .values_list("id", flat=True)
        )

    def get(self, request):
        # ---- params ----
        try:
            ref_id = int(request.GET.get("ref"))
        except (TypeError, ValueError):
            return JsonResponse({"error": "Missing or invalid 'ref' parameter"}, status=400)

        top_n = request.GET.get("top_n")
        try:
            top_n = int(top_n) if top_n is not None else 50
        except Exception:
            top_n = 50
        top_n = max(10, min(200, top_n))

        method = (request.GET.get("method") or "umap").lower()
        metric = (request.GET.get("metric") or "identity").lower()
        if metric not in ("identity", "similarity"):
            metric = "identity"
        exclude_orphans = self._truthy(request.GET.get("exclude_orphans", "true"))

        # ---- 1) neighbors of ref (similarity desc) ----
        pairs = list(
            ReceptorSimilarity.objects
            .filter(Q(protein_ref_id=ref_id) | Q(protein_target_id=ref_id))
            .annotate(
                other_id=Case(
                    When(protein_ref_id=ref_id, then=F("protein_target_id")),
                    default=F("protein_ref_id"),
                    output_field=IntegerField(),
                )
            )
            .order_by("-similarity")
            .values("other_id", "similarity", "identity")
        )
        if not pairs:
            return JsonResponse({"points": [], "meta": {"note": "No neighbors for ref"}})

        # ---- 2) optionally drop other orphans ----
        if exclude_orphans:
            orphan_ids = self.get_orphan_ids()
            pairs = [p for p in pairs if p["other_id"] not in orphan_ids]

        ranked = pairs[:top_n]
        if not ranked:
            return JsonResponse({"points": [], "meta": {"note": "No non-orphan neighbors found"}})

        kept_ids = [ref_id] + [p["other_id"] for p in ranked]

        # ---- 3) protein metadata ----
        proteins = (
            Protein.objects
            .filter(id__in=kept_ids)
            .select_related("family__parent__parent__parent")
        )
        pmap = {p.id: p for p in proteins}
        if ref_id not in pmap:
            return JsonResponse({"error": "Reference protein not found"}, status=404)

        labels, ordered_ids = [], []
        for pid in kept_ids:
            p = pmap.get(pid)
            if not p:
                continue
            lab = self._lab(p.entry_name or "")
            if lab and lab not in labels:
                labels.append(lab)
                ordered_ids.append(pid)

        N = len(labels)
        if N < 3:
            return JsonResponse({"points": [], "meta": {"note": "Too few points", "n_points": N}})

        id_to_idx = {pid: i for i, pid in enumerate(ordered_ids)}

        # ---- 4) all pairs among kept_ids (one small query) ----
        subpairs = list(
            ReceptorSimilarity.objects
            .filter(protein_ref_id__in=kept_ids, protein_target_id__in=kept_ids)
            .values("protein_ref_id", "protein_target_id", "similarity", "identity")
        )
        # undirected map
        pv = {}
        for r in subpairs:
            a, b = r["protein_ref_id"], r["protein_target_id"]
            if a == b:
                continue
            key = (a, b) if a < b else (b, a)
            if key not in pv:
                pv[key] = r

        D = np.full((N, N), np.nan, dtype=float)
        np.fill_diagonal(D, 0.0)

        def to_dist(rec):
            val = rec["identity"] if metric == "identity" else rec["similarity"]
            try:
                v = float(val)
            except Exception:
                return np.nan
            return max(0.0, min(1.0, 1.0 - v / 100.0))

        for (a, b), rec in pv.items():
            if a in id_to_idx and b in id_to_idx:
                i, j = id_to_idx[a], id_to_idx[b]
                d = to_dist(rec)
                D[i, j] = d
                D[j, i] = d

        if np.isnan(D).any():
            col_med = np.nanmedian(D, axis=0)
            inds = np.where(np.isnan(D))
            D[inds] = np.take(col_med, inds[1])
            D = 0.5 * (D + D.T)
            np.fill_diagonal(D, 0.0)

        # ---- 5) embed (t-SNE only) ----

        # perplexity must be strictly < N; keep it reasonable for tiny N
        perplexity = max(
            1.0,
            min(40.0, (N - 1) / 3.0, N - 1 - 1e-9),
        )

        tsne = TSNE(
            n_components=2,
            metric="precomputed",
            perplexity=perplexity,
            random_state=42,
            init="random",        # required with precomputed distances
            learning_rate="auto", # silences future warning / modern default
            square_distances=True # modern behavior; no deprecation warning
        )
        coords = tsne.fit_transform(D)
        used_method = "tsne"


        # ---- 6) scalar vs ref for Gradient color ----
        ref_field = "identity" if metric == "identity" else "similarity"
        val_vs_ref = {}
        for r in pairs:  # original ref↔other list
            if r["other_id"] in id_to_idx:
                val = r.get(ref_field)
                if val is not None:
                    lab = self._lab(pmap[r["other_id"]].entry_name or "")
                    val_vs_ref[lab] = float(val)

        # ---- 7) response points ----
        points = []
        for i, pid in enumerate(ordered_ids):
            p = pmap[pid]
            fam = getattr(p.family, "parent", None)
            lig = getattr(fam, "parent", None) if fam else None
            cls = getattr(lig, "parent", None) if lig else None

            clazz = getattr(cls, "name", "") if cls else ""
            lig_t = getattr(lig, "name", "") if lig else ""
            fam_n = getattr(fam, "name", "") if fam else ""
            lab = labels[i]

            fill = 100.0 if pid == ref_id else val_vs_ref.get(lab)

            points.append({
                "id": pid,
                "label": lab,
                "x": float(coords[i, 0]),
                "y": float(coords[i, 1]),
                "Class": clazz,
                "Ligand type": lig_t,
                "Receptor family": fam_n,
                "fill": float(fill) if fill is not None else None,
                "is_ref": (pid == ref_id),
            })

        return JsonResponse({
            "points": points,
            "ref": {"id": ref_id, "label": labels[0]},
            "meta": {"method": used_method, "metric": metric, "top_n": len(points), "n_points": N}
        })

# ----------------------------- Shared helpers ------------------------------

class OrphanSelect2Mixin:
    ORPHAN_LT = 'Orphan receptors'

    @staticmethod
    def _clean_gtop_name(nm):
        if not nm:
            return "-"
        s = (nm
             .replace("receptor", "")
             .replace("-adrenoceptor", "")
             .replace("<i>", "").replace("</i>", "")
             .strip())
        return s or "-"

    def get_orphans_select2(self):
        # pull the whole lineage to read "Class ..."
        orphans_qs = (
            Protein.objects
            .filter(
                parent_id__isnull=True,
                species_id=1,
                family__parent__parent__name__iexact=self.ORPHAN_LT,
            )
            .select_related('family__parent__parent__parent')  # <-- add this
            .prefetch_related(
                Prefetch('genes',
                         queryset=Gene.objects.filter(position=0),
                         to_attr='primary_genes_self')
            )
            .order_by('entry_name')
        )

        data = []
        for p in orphans_qs:
            gene = (
                p.primary_genes_self[0].name
                if getattr(p, 'primary_genes_self', None)
                else (p.entry_name.split('_')[0].upper() if p.entry_name else "-")
            )
            nm = self._clean_gtop_name(p.name)

            # read raw class name from lineage (e.g. "Class A orphans")
            fam = getattr(p.family, 'parent', None)
            lig = getattr(fam, 'parent', None) if fam else None
            cls = getattr(lig, 'parent', None) if lig else None
            raw_class = getattr(cls, 'name', '') or ''

            # normalize directly in Python
            if raw_class.lower().startswith("other gpcr"):
                raw_class = "Classless"
            elif raw_class.lower().startswith("class "):
                # keep only the letter, e.g. "Class A orphans" → "Class A"
                m = re.search(r"class\s*([a-z])", raw_class, re.I)
                raw_class = f"Class {m.group(1).upper()}" if m else raw_class

            data.append({
                "id": p.id,
                "text": nm,
                "name": nm,
                "entry_name": p.entry_name,
                "gene": gene,
                "class": raw_class,
            })
        return data


def _get_ref_class_info(ref_id):
    """
    Return (class_id, class_name) for the reference protein, inferred directly
    from ReceptorSimilarity (works regardless of row direction).
    """
    row = (
        ReceptorSimilarity.objects
        .filter(Q(protein_ref_id=ref_id) | Q(protein_target_id=ref_id))
        .values('protein_ref_id', 'protein_target_id', 'ref_class_id', 'target_class_id')
        .first()
    )
    if not row:
        return (None, None)

    if row['protein_ref_id'] == ref_id:
        class_id = row['ref_class_id']
    else:
        class_id = row['target_class_id']

    class_name = None
    if class_id:
        try:
            cls = ProteinFamily.objects.only('id', 'name').get(id=class_id)
            # Normalize like your Select2 (Classless / "Class X")
            name = (cls.name or '').strip()
            if name.lower().startswith('other gpcr'):
                name = 'Classless'
            else:
                m = re.search(r'class\s*([a-z])', name, re.I)
                if m:
                    name = f'Class {m.group(1).upper()}'
            class_name = name
        except ProteinFamily.DoesNotExist:
            pass

    return (class_id, class_name)

def _build_similarity_rows(ref_id):
    """
    Reuses the exact logic from SimilarityTopAPI to produce the table rows.
    Returns: list[dict] (the 'results' list you already send today)
    """
    from string import Template

    pairs_qs = (
        ReceptorSimilarity.objects
        .filter(Q(protein_ref_id=ref_id) | Q(protein_target_id=ref_id))
        .annotate(
            other_id=Case(
                When(protein_ref_id=ref_id, then=F('protein_target_id')),
                default=F('protein_ref_id'),
                output_field=IntegerField(),
            )
        )
        .order_by('-similarity')
        .values('other_id', 'similarity', 'identity')
    )
    pairs = list(pairs_qs)
    if not pairs:
        return []

    other_ids_all = [p['other_id'] for p in pairs]

    LT_ORPHAN = 'Orphan receptors'
    lt_map = dict(
        Protein.objects
               .filter(id__in=other_ids_all)
               .values_list('id', 'family__parent__parent__name')
    )
    def is_orphan(pid):
        return (lt_map.get(pid) or '').strip().lower() == LT_ORPHAN.lower()

    liganded_rows = [p for p in pairs if not is_orphan(p['other_id'])]

    if liganded_rows:
        cutoff = liganded_rows[9]['similarity'] if len(liganded_rows) >= 10 else liganded_rows[-1]['similarity']
    else:
        cutoff = pairs[min(9, len(pairs) - 1)]['similarity']

    kept = [p for p in pairs if p['similarity'] >= cutoff]
    kept_ids = [p['other_id'] for p in kept]

    gtop_links_qs = WebLink.objects.select_related('web_resource').filter(web_resource__slug='gtop')
    proteins_qs = (
        Protein.objects
        .filter(id__in=kept_ids)
        .select_related('family__parent__parent__parent')
        .prefetch_related(
            Prefetch('genes',
                     queryset=Gene.objects.filter(position=0),
                     to_attr='primary_genes_self'),
            Prefetch('web_links',
                     queryset=gtop_links_qs,
                     to_attr='gtop_links_self'),
        )
    )
    proteins = {p.id: p for p in proteins_qs}

    endo_qs = (
        Endogenous_GTP.objects
        .filter(receptor_id__in=kept_ids)
        .select_related('ligand', 'ligand__ligand_type')
    )
    endo_by_receptor = {}
    for e in endo_qs:
        if e.ligand:
            endo_by_receptor.setdefault(e.receptor_id, []).append(e)

    def clean_iuphar_name(nm):
        if not nm:
            return "-"
        s = nm.replace("receptor", "").replace("-adrenoceptor", "").replace("<i>", "").replace("</i>", "").strip()
        return s or "-"

    def build_gtop_url(wl):
        try:
            return Template(wl.web_resource.url).substitute(index=wl.index)
        except Exception:
            return None

    sim_map = {p['other_id']: p['similarity'] for p in kept}
    idn_map = {p['other_id']: p['identity']   for p in kept}

    results = []
    for row in kept:
        pid = row['other_id']
        p = proteins.get(pid)
        if not p:
            continue

        gene_name = p.primary_genes_self[0].name if getattr(p, 'primary_genes_self', None) else (
            (p.entry_name.split('_')[0].upper()) if p.entry_name else "-")

        entry_name  = p.entry_name or None
        gpcrdb_link = f"/protein/{entry_name}" if entry_name else "-"
        uniprot_link = f"https://www.uniprot.org/uniprot/{p.accession}" if p.accession else None

        wl_self = p.gtop_links_self[0] if getattr(p, 'gtop_links_self', None) else None
        iuphar_link = build_gtop_url(wl_self) if wl_self else None
        iuphar_name = clean_iuphar_name(p.name)

        family_name = getattr(getattr(p.family, "parent", None), "name", None)
        ligand_type = getattr(getattr(getattr(p.family, "parent", None), "parent", None), "name", None)
        clazz       = getattr(getattr(getattr(getattr(p.family, "parent", None), "parent", None), "parent", None), "name", None)

        lig_items = endo_by_receptor.get(pid, [])
        seen, endo_ligands, lig_types = set(), [], set()
        for e in lig_items:
            lig = e.ligand
            if not lig:
                continue
            if lig.id not in seen:
                seen.add(lig.id)
                endo_ligands.append({"id": lig.id, "name": lig.name})
            if lig.ligand_type:
                lig_types.add(lig.ligand_type.name)
        endo_type = "<br>".join(sorted(lig_types)) if lig_types else "-"

        results.append({
            "other_id": pid,
            "Gene": gene_name,
            "entry_name": entry_name,
            "gpcrdb_link": gpcrdb_link,
            "uniprot_link": uniprot_link,
            "iuphar_name": iuphar_name,
            "iuphar_link": iuphar_link,
            "family": family_name,
            "ligand_type": ligand_type,
            "class": clazz,
            "similarity": sim_map.get(pid, 0),
            "identity": idn_map.get(pid, 0),
            "endo_ligands": endo_ligands,
            "endo_type": endo_type,
        })

    return results

def _build_embedding_payload(ref_id, *, top_n=50, metric='identity', exclude_orphans=True):
    """
    Embedding selection logic (per your new rules), then the same t-SNE pipeline.
    Rules:
      - Classless  -> top N across all classes
      - Class C    -> ALL Class C
      - Other      -> top N within the same class as ref
    """
    # 0) Figure out the reference class once (using the RS table)
    ref_class_id, ref_class_name = _get_ref_class_info(ref_id)

    # If we couldn't deduce the class, fall back to the broad top N across all
    if not ref_class_id:
        base_qs = (
            ReceptorSimilarity.objects
            .filter(Q(protein_ref_id=ref_id) | Q(protein_target_id=ref_id))
            .order_by('-similarity')[:max(10, min(200, int(top_n)))]
        )
    else:
        # 1) Build class-based neighbor selection
        #    Use indexed filters directly on RS:
        #    - same-class rows regardless of direction
        same_class_q = (
            Q(protein_ref_id=ref_id, target_class_id=ref_class_id) |
            Q(protein_target_id=ref_id, ref_class_id=ref_class_id)
        )

        if ref_class_name == 'Classless':
            # Top N across all classes
            base_qs = (
                ReceptorSimilarity.objects
                .filter(Q(protein_ref_id=ref_id) | Q(protein_target_id=ref_id))
                .order_by('-similarity')[:max(10, min(200, int(top_n)))]
            )
        elif ref_class_name == 'Class C':
            # ALL Class C (no slice)
            base_qs = (
                ReceptorSimilarity.objects
                .filter(same_class_q)
                .order_by('-similarity')
            )
        else:
            # Same class only, Top N
            base_qs = (
                ReceptorSimilarity.objects
                .filter(same_class_q)
                .order_by('-similarity')[:max(10, min(200, int(top_n)))]
            )

    # 2) Convert to a uniform shape with other_id + values (works both directions)
    pairs = list(
        base_qs.annotate(
            other_id=Case(
                When(protein_ref_id=ref_id, then=F('protein_target_id')),
                default=F('protein_ref_id'),
                output_field=IntegerField(),
            )
        ).values('other_id', 'similarity', 'identity')
    )

    if not pairs:
        return {"points": [], "ref": {"id": ref_id, "label": ""}, "meta": {"note": "No neighbors for ref"}}

    # Note: the old 'exclude_orphans' flag is redundant here because the
    # class-based selection already governs inclusion. Kept for compatibility.

    # 3) Cap the "ranked" set only if we have a slice-less case above (Class C keeps all already)
    ranked = pairs  # already sliced where needed
    kept_ids = [ref_id] + [p["other_id"] for p in ranked]

    # 4) Fetch protein annotation for labels & legend
    proteins = (
        Protein.objects
        .filter(id__in=kept_ids)
        .select_related("family__parent__parent__parent")
    )
    pmap = {p.id: p for p in proteins}
    if ref_id not in pmap:
        return {"error": "Reference protein not found"}

    def _lab(en):
        return (en or "").replace("_human", "")

    labels, ordered_ids = [], []
    for pid in kept_ids:
        p = pmap.get(pid)
        if not p:
            continue
        lab = _lab(p.entry_name or "")
        if lab and lab not in labels:
            labels.append(lab)
            ordered_ids.append(pid)

    N = len(labels)
    if N < 3:
        return {"points": [], "ref": {"id": ref_id, "label": labels[0] if labels else ""}, "meta": {"note": "Too few points", "n_points": N}}

    id_to_idx = {pid: i for i, pid in enumerate(ordered_ids)}

    # 5) Build dense pair set among kept ids for distances
    subpairs = list(
        ReceptorSimilarity.objects
        .filter(protein_ref_id__in=kept_ids, protein_target_id__in=kept_ids)
        .values("protein_ref_id", "protein_target_id", "similarity", "identity")
    )
    pv = {}
    for r in subpairs:
        a, b = r["protein_ref_id"], r["protein_target_id"]
        if a == b:
            continue
        key = (a, b) if a < b else (b, a)
        if key not in pv:
            pv[key] = r

    D = np.full((N, N), np.nan, dtype=float)
    np.fill_diagonal(D, 0.0)

    metric = (metric or "identity").lower()
    if metric not in ("identity", "similarity"):
        metric = "identity"

    def to_dist(rec):
        val = rec["identity"] if metric == "identity" else rec["similarity"]
        try:
            v = float(val)
        except Exception:
            return np.nan
        return max(0.0, min(1.0, 1.0 - v / 100.0))

    for (a, b), rec in pv.items():
        if a in id_to_idx and b in id_to_idx:
            i, j = id_to_idx[a], id_to_idx[b]
            d = to_dist(rec)
            D[i, j] = d
            D[j, i] = d

    if np.isnan(D).any():
        col_med = np.nanmedian(D, axis=0)
        inds = np.where(np.isnan(D))
        D[inds] = np.take(col_med, inds[1])
        D = 0.5 * (D + D.T)
        np.fill_diagonal(D, 0.0)

    perplexity = max(1.0, min(40.0, (N - 1) / 3.0, N - 1 - 1e-9))
    tsne = TSNE(
        n_components=2,
        metric="precomputed",
        perplexity=perplexity,
        random_state=42,
        init="random",
        learning_rate="auto",
        square_distances=True
    )
    coords = tsne.fit_transform(D)
    used_method = "tsne"

    # 6) Values vs ref (for hover/gradient)
    ref_field = "identity" if metric == "identity" else "similarity"
    val_vs_ref = {}
    for r in ranked:
        oid = r["other_id"]
        if oid in id_to_idx:
            p = pmap.get(oid)
            if p:
                lab = _lab(p.entry_name or "")
                val = r.get(ref_field)
                if val is not None:
                    val_vs_ref[lab] = float(val)

    # 7) Build points w/ annotation
    points = []
    for i, pid in enumerate(ordered_ids):
        p = pmap[pid]
        fam = getattr(p.family, "parent", None)
        lig = getattr(fam, "parent", None) if fam else None
        cls = getattr(lig, "parent", None) if lig else None

        clazz = getattr(cls, "name", "") if cls else ""
        lig_t = getattr(lig, "name", "") if lig else ""
        fam_n = getattr(fam, "name", "") if fam else ""
        lab = labels[i]

        # Normalize class label like before
        if clazz.lower().startswith('other gpcr'):
            clazz = 'Classless'
        else:
            m = re.search(r'class\s*([a-z])', clazz, re.I)
            if m:
                clazz = f'Class {m.group(1).upper()}'

        fill = 100.0 if pid == ref_id else val_vs_ref.get(lab)

        points.append({
            "id": pid,
            "label": lab,
            "x": float(coords[i, 0]),
            "y": float(coords[i, 1]),
            "Class": clazz,
            "Ligand type": lig_t,
            "Receptor family": fam_n,
            "fill": float(fill) if fill is not None else None,
            "is_ref": (pid == ref_id),
        })

    return {
        "points": points,
        "ref": {"id": ref_id, "label": labels[0]},
        "meta": {"method": used_method, "metric": metric, "top_n": len(points), "n_points": N, "ref_class": ref_class_name}
    }

# ------------------------------ New merged page -----------------------------

class OrphanSimilarityExplorer(OrphanSelect2Mixin, TemplateView):
    """
    Single page that will host tabs: (1) Neighbor Table, (2) Cluster Embedding.
    The template (to be added) will keep a hidden wrapper, and after a selection
    it will call the bundle API once, show BusyLoad, then reveal the tabs.
    """
    template_name = 'class_similarity/OrphanSimilarityExplorer.html'

    def get_context_data(self, **kwargs):
        ctx = super().get_context_data(**kwargs)
        ctx['orphans_select2'] = json.dumps(self.get_orphans_select2())
        return ctx


# ------------------------------- New merged API -----------------------------

class SimilarityBundleAPI(View):
    """
    GET /class_similarity/api/bundle?ref=<protein_id>&top_n=50&metric=identity&exclude_orphans=true
    Returns BOTH:
      - table.results[]   (same structure as SimilarityTopAPI)
      - embedding.{points,ref,meta} (same structure as SimilarityEmbeddingAPI)
    Optional query flags:
      - only=table   -> return only table part
      - only=embed   -> return only embedding part
    """

    @staticmethod
    def _truthy(v):
        return str(v).lower() not in ("", "0", "false", "no", "off", "none")

    def get(self, request):
        # ---- params ----
        ref_raw = request.GET.get('ref')
        try:
            ref_id = int(ref_raw)
        except (TypeError, ValueError):
            return JsonResponse({"error": "Missing or invalid 'ref' parameter"}, status=400)

        top_n = request.GET.get("top_n")
        try:
            top_n = int(top_n) if top_n is not None else 50
        except Exception:
            top_n = 50
        top_n = max(10, min(200, top_n))

        metric = (request.GET.get("metric") or "identity").lower()
        if metric not in ("identity", "similarity"):
            metric = "identity"

        exclude_orphans = self._truthy(request.GET.get("exclude_orphans", "true"))

        only = (request.GET.get("only") or "").strip().lower()

        payload = {}

        # Build parts according to 'only'
        if only in ("", "table"):
            table_rows = _build_similarity_rows(ref_id)
            payload["table"] = {"results": table_rows}

            # If user asked only table, short-circuit
            if only == "table":
                return JsonResponse(payload)

        if only in ("", "embed"):
            embed = _build_embedding_payload(
                ref_id,
                top_n=top_n,
                metric=metric,
                exclude_orphans=exclude_orphans,
            )
            payload["embedding"] = embed

            if only == "embed":
                return JsonResponse(payload)

        return JsonResponse(payload)



class ReceptorSimilarityExportExcel(View):
    """
    GET /class_similarity/ReceptorSimilarityExportExcel
    -> returns an .xlsx file with:
       ref_id, ref_entry, target_id, target_entry,
       ref_class, target_class, identity, similarity
    """

    def get(self, request, *args, **kwargs):
        # Allow override of filename via ?filename=...
        filename = request.GET.get('filename', 'receptor_similarity.xlsx')

        # Query the data we need
        qs = (
            ReceptorSimilarity.objects
            .select_related('protein_ref', 'protein_target', 'ref_class', 'target_class')
            .values(
                'protein_ref_id',
                'protein_ref__entry_name',
                'protein_target_id',
                'protein_target__entry_name',
                'ref_class__name',
                'target_class__name',
                'identity',
                'similarity',
            )
        )

        # Build Excel in memory
        output = BytesIO()
        workbook = xlsxwriter.Workbook(output, {'in_memory': True})
        worksheet = workbook.add_worksheet('ReceptorSimilarity')

        # Headers as requested
        headers = [
            'ref_id',
            'ref_entry',
            'target_id',
            'target_entry',
            'ref_class',
            'target_class',
            'identity',
            'similarity',
        ]

        # Write header row
        for col, header in enumerate(headers):
            worksheet.write(0, col, header)

        # Write data
        row_idx = 1
        for row in qs:
            # ref
            worksheet.write(row_idx, 0, row['protein_ref_id'])
            worksheet.write(row_idx, 1, row['protein_ref__entry_name'])

            # target
            worksheet.write(row_idx, 2, row['protein_target_id'])
            worksheet.write(row_idx, 3, row['protein_target__entry_name'])

            # classes (nullable)
            ref_class_name = row['ref_class__name'] or ''
            target_class_name = row['target_class__name'] or ''

            worksheet.write(row_idx, 4, ref_class_name)
            worksheet.write(row_idx, 5, target_class_name)

            # numbers
            worksheet.write(row_idx, 6, row['identity'])
            worksheet.write(row_idx, 7, row['similarity'])

            row_idx += 1

        workbook.close()
        output.seek(0)

        # HTTP response
        response = HttpResponse(
            output.read(),
            content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
        )
        response['Content-Disposition'] = 'attachment; filename="%s"' % filename
        return response


# ------------------------------ Structure similarity -----------------------------

class StructureSim(TemplateView):
    template_name = 'class_similarity/StructureSim.html'

    # --- config ---
    DATA_FOLDER = 'structure_data'
    FILES = {
        'inactive': 'GPCR_structure_clustering_inactiveRep.csv',
        'active':   'GPCR_structure_clustering_activeStructuresRep.csv',
    }
    DEFAULT_STATE = 'inactive'

    # NEW: ligand/sense annotation Excel
    LIGAND_META_FOLDER = 'protein_data'
    LIGAND_META_FILE = 'Ligand type update plus sense column.xlsx'

    CACHE_TIMEOUT = 60 * 15  # 15 minutes
    CACHE_NS = 'structuresim:data'

    # ---------- helpers ----------
    def _file_path(self, state: str):
        """Get full path to the CSV for the given state ('inactive' or 'active')."""
        try:
            file_name = self.FILES[state]
        except KeyError:
            raise ValueError(
                f"Unknown state '{state}'. Expected one of: {', '.join(self.FILES.keys())}"
            )
        return os.path.join(settings.DATA_DIR, self.DATA_FOLDER, file_name)

    def _ligand_meta_path(self):
        """Full path to the Excel with ligand-type + sense annotations."""
        return os.path.join(
            settings.DATA_DIR,
            self.LIGAND_META_FOLDER,
            self.LIGAND_META_FILE,
        )

    def _cache_key(self, state: str, extra: str = ''):
        """Auto-bust on CSV mtime; allow manual ?v=...; separated per state."""
        p = self._file_path(state)
        try:
            mtime = int(os.path.getmtime(p))
        except Exception:
            mtime = 0
        manual = self.request.GET.get('v', '')
        # bump version for new physio logic
        return f"{self.CACHE_NS}:v9:{state}:{mtime}:{manual}:{extra}"

    def _family_lineage_names(self, protein):
        """
        Returns (class_name, ligand_type, receptor_family) by walking up family parents:
          parent^3 = Class, parent^2 = Ligand type, parent^1 = Receptor family
        """
        f = getattr(protein, 'family', None)
        if not f:
            return None, None, None
        p1 = getattr(f, 'parent', None)
        p2 = getattr(p1, 'parent', None) if p1 else None
        p3 = getattr(p2, 'parent', None) if p2 else None
        receptor_family = getattr(p1, 'name', None)
        ligand_type     = getattr(p2, 'name', None)
        class_name      = getattr(p3, 'name', None)
        return class_name, ligand_type, receptor_family

    def _canonical_protein_from_structure(self, s, models_by_template):
        """
        Prefer canonical protein rows for names/labels:
          1) structure's protein if it has accession
          2) else its parent if it has accession
          3) else StructureModel(main_template=s) protein with accession
          4) else fallback to structure's protein
        """
        p = s.protein_conformation.protein
        if getattr(p, 'accession', None):
            return p
        parent = getattr(p, 'parent', None)
        if parent and getattr(parent, 'accession', None):
            return parent
        for sm in models_by_template.get(s.id, ()):
            mp = sm.protein
            if getattr(mp, 'accession', None):
                return mp
        return p

    def _compute_tsne(self, D):
        """
        D: numpy (n x n) distance matrix, symmetric, zeros on diagonal.
        Returns coords (n x 2), with robust defaults for perplexity.
        """
        n = D.shape[0]
        perplexity = max(5.0, min(40.0, (n - 1) / 3.0, n - 2.0))
        tsne = TSNE(
            n_components=2,
            metric="precomputed",
            perplexity=perplexity,
            random_state=42,
            init="random",
            learning_rate="auto",
            square_distances=True,
        )
        coords = tsne.fit_transform(D)
        return coords

    @staticmethod
    def _norm_lt(name: str) -> str:
        """Normalize a ligand type string to coarse buckets: 'peptide', 'small', or 'other'."""
        s = (name or '').strip().lower()
        if 'peptide' in s or 'protein' in s:
            return 'peptide'
        if ('small' in s and 'molecule' in s) or s == 'small-molecule' or s == 'small molecule':
            return 'small'
        return 'other'

    @staticmethod
    def _uniprot_from_entry(entry_name):
        """
        Convert GPCRdb entry_name like 'adra1a_human' → 'ADRA1A'.
        """
        if not entry_name:
            return None
        s = str(entry_name).strip()
        if not s:
            return None
        s = s.split('_')[0]
        return s.upper()

    def _load_ligand_meta(self):
        """
        Load Excel with ligand type / sense info, return:
          { UNIPROT_CODE (upper) : { 'Ligand type': ..., 'Sense': ..., 'Receptor family': ..., 'Class': ... } }
        Cached via Django cache.
        """
        cache_key = f"{self.CACHE_NS}:ligandmeta:v1"
        cached = cache.get(cache_key)
        if cached is not None:
            return cached

        path = self._ligand_meta_path()
        if not os.path.exists(path):
            cache.set(cache_key, {}, self.CACHE_TIMEOUT)
            return {}

        df = pd.read_excel(path)

        # Find column names, being tolerant to line breaks / spacing
        def pick(*cands):
            cands = {c.strip() for c in cands}
            for name in df.columns:
                s = str(name).strip()
                if s in cands:
                    return name
            return None

        col_uni    = pick('GPCRs (UniProt)', 'GPCRs\n(UniProt)')
        col_lt     = pick('Ligand type')
        col_sense  = pick('Sense')
        col_family = pick('Receptor family')
        col_class  = pick('Class')

        if not col_uni:
            cache.set(cache_key, {}, self.CACHE_TIMEOUT)
            return {}

        mapping = {}
        for _, row in df.iterrows():
            uni = str(row[col_uni]).strip()
            if not uni or uni.lower() == 'nan':
                continue
            key = uni.upper()
            rec = {}

            if col_lt is not None:
                v = row[col_lt]
                if pd.notna(v):
                    rec['Ligand type'] = str(v).strip()

            if col_sense is not None:
                v = row[col_sense]
                if pd.notna(v):
                    rec['Sense'] = str(v).strip()

            if col_family is not None:
                v = row[col_family]
                if pd.notna(v):
                    rec['Receptor family'] = str(v).strip()

            if col_class is not None:
                v = row[col_class]
                if pd.notna(v):
                    rec['Class'] = str(v).strip()

            if rec:
                mapping[key] = rec

        cache.set(cache_key, mapping, self.CACHE_TIMEOUT)
        return mapping

    def _build_payload(self, state: str, want_embed: bool = False):
        """
        Read CSV for given state, validate, enrich with DB metadata (canonical protein),
        compute physiological ligand-type consensus per receptor from ALL Endogenous_GTP,
        and optionally compute t-SNE embedding on the distance matrix.
        """
        cache_key = self._cache_key(state, extra=('embed' if want_embed else 'noembed'))
        cached = cache.get(cache_key)
        if cached:
            return cached

        # --- load CSV ---
        file_path = self._file_path(state)
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        df = pd.read_csv(file_path, index_col=0)

        # --- basic validation ---
        if df.shape[0] != df.shape[1]:
            raise ValueError("CSV must be a square distance matrix (rows == columns).")

        labels = df.index.astype(str).tolist()
        cols = df.columns.astype(str).tolist()

        if cols != labels:
            try:
                df = df.loc[labels, labels]
            except Exception:
                raise ValueError("CSV columns don't match index (PDB codes).")

        # force numeric, impute, symmetrize, zero diag
        df = df.apply(pd.to_numeric, errors='coerce')
        arr = df.values.astype(float)
        if np.isnan(arr).any():
            col_med = np.nanmedian(arr, axis=0)
            inds = np.where(np.isnan(arr))
            arr[inds] = np.take(col_med, inds[1])
        arr = 0.5 * (arr + arr.T)
        np.fill_diagonal(arr, 0.0)
        matrix = arr.tolist()

        # --- DB fetch: map PDB -> Structure ---
        codes_upper = [c.upper() for c in labels]
        structures_qs = (
            Structure.objects
            .filter(pdb_code__index__in=codes_upper)
            .select_related(
                'pdb_code',
                'structure_type',
                'state',
                'protein_conformation__protein',
                'protein_conformation__protein__family',
                'protein_conformation__protein__parent',
                'protein_conformation__protein__parent__family',
            )
        )
        structures = list(structures_qs)
        by_pdb = {s.pdb_code.index.upper(): s for s in structures}

        # preload StructureModel by main_template (fallback canonicalization)
        template_ids = [s.id for s in structures]
        models_by_template = {}
        if template_ids:
            for sm in (
                StructureModel.objects
                .filter(main_template_id__in=template_ids)
                .select_related('protein')
            ):
                models_by_template.setdefault(sm.main_template_id, []).append(sm)

        # --- build proteins dict + collect receptor ids for consensus lookup ---
        proteins = {}
        unmatched = []
        label_to_receptor_id = {}

        for original in labels:
            s = by_pdb.get(original.upper())
            if not s:
                unmatched.append(original)
                continue

            p = self._canonical_protein_from_structure(s, models_by_template)

            stype = getattr(s.structure_type, "type_short", None)
            stype = stype() if callable(stype) else getattr(s.structure_type, "name", None)

            class_name, ligand_type, receptor_family = self._family_lineage_names(p)

            proteins[original] = {
                "pdb_code": s.pdb_code.index.upper(),
                "protein_entry": getattr(p, "entry_name", None),
                "protein_name": getattr(p, "name", None),
                "protein_class": class_name,
                "ligand_type": ligand_type,
                "protein_family": receptor_family,
                "receptor_family": receptor_family,
                "state": getattr(s.state, "slug", None),
                "structure_type": stype,
            }

            if getattr(p, "id", None):
                label_to_receptor_id[original] = p.id

        # --- consensus physiological ligand type per receptor (ALL Endogenous_GTP) ---
        consensus_map = {}
        raw_types_map = {}
        receptor_ids = list(set(label_to_receptor_id.values()))
        if receptor_ids:
            lt_rows = (
                Endogenous_GTP.objects
                .filter(receptor_id__in=receptor_ids)
                .values('receptor_id', 'ligand__ligand_type__name')
                .distinct()
            )
            agg = {}
            raw = {}
            for row in lt_rows:
                rid = row['receptor_id']
                lt = (row['ligand__ligand_type__name'] or '').strip()
                if rid not in agg:
                    agg[rid] = set()
                    raw[rid] = set()
                agg[rid].add(self._norm_lt(lt))
                if lt:
                    raw[rid].add(lt)

            for rid, kinds in agg.items():
                has_pep = ('peptide' in kinds)
                has_sml = ('small' in kinds)
                if has_pep and has_sml:
                    cat = 'Peptide/protein & small-molecule'
                elif has_pep:
                    cat = 'Peptide/protein'
                elif has_sml:
                    cat = 'Small-molecule'
                else:
                    cat = 'Other'
                consensus_map[rid] = cat
                raw_types_map[rid] = sorted(raw.get(rid, []))

        for lab, rid in label_to_receptor_id.items():
            if rid in consensus_map:
                proteins[lab]["physio_ligand_consensus"] = consensus_map[rid]
                proteins[lab]["physio_ligand_types_raw"] = raw_types_map.get(rid, [])

        # --- apply Excel ligand-type + sense overrides ---
        ligand_meta = self._load_ligand_meta()
        if ligand_meta:
            for lab, meta in proteins.items():
                entry = meta.get("protein_entry")
                uni = self._uniprot_from_entry(entry)
                if not uni:
                    continue
                ann = ligand_meta.get(uni)
                if not ann:
                    continue

                lt = ann.get("Ligand type")
                if lt:
                    meta["ligand_type"] = lt

                fam = ann.get("Receptor family")
                if fam:
                    meta["protein_family"] = fam
                    meta["receptor_family"] = fam

                cls = ann.get("Class")
                if cls:
                    meta["protein_class"] = cls

                sense = ann.get("Sense")
                if sense is not None:
                    meta["sense"] = sense
            
        # --- final physio-ligand cleanup based on UPDATED ligand_type ---
        for lab, meta in proteins.items():
            lt = (meta.get("ligand_type") or "").strip()
            lt_low = lt.lower()
            physio = (meta.get("physio_ligand_consensus") or "").strip()

            # 1) Special cases that should always override

            #    Adhesion receptors → Cleaved-endterm | PPI
            if lt_low == "adhesion receptors":
                meta["physio_ligand_consensus"] = "Tethered ligand | PPI"

            #    Ion receptors → Ion
            elif lt_low == "ion receptors":
                meta["physio_ligand_consensus"] = "Ion"

            #    Peptide / amino acid / protein receptors → Peptide/protein
            elif lt_low in {
                "peptide receptors",
                "amino acid receptors",
                "protein receptors",
            }:
                meta["physio_ligand_consensus"] = "Peptide/protein"

            #    Light / odorant / tastant receptors → Small-molecule
            elif lt_low in {
                "light receptors",
                "odorant receptors",
                "tastant receptors",
            }:
                meta["physio_ligand_consensus"] = "Small-molecule"

            #    Unknown receptors → Orphan
            elif lt_low == "unknown receptors":
                meta["physio_ligand_consensus"] = "Orphan"

            # 2) Anything still missing after all the above → Orphan
            if not meta.get("physio_ligand_consensus"):
                meta["physio_ligand_consensus"] = "Orphan"


        payload = {
            "state": state,
            "labels": labels,
            "matrix": matrix,
            "proteins": proteins,
            "unmatched": unmatched,
        }

        # --- optional embedding (t-SNE on the distance matrix) ---
        if want_embed:
            try:
                D = np.array(matrix, dtype=float)
                coords = self._compute_tsne(D)
                points = []
                for i, lab in enumerate(labels):
                    meta = proteins.get(lab, {})
                    points.append({
                        "label": lab,
                        "x": float(coords[i, 0]),
                        "y": float(coords[i, 1]),
                        "Class": meta.get("protein_class") or "",
                        "Ligand type": meta.get("ligand_type") or "",
                        "Receptor family": meta.get("receptor_family") or meta.get("protein_family") or "",
                        "Sense": meta.get("sense") or "",
                    })
                payload["tsne"] = {
                    "method": "tsne",
                    "points": points,
                    "n": len(points),
                }
            except Exception as e:
                payload["tsne"] = {"error": str(e)}

        cache.set(cache_key, payload, self.CACHE_TIMEOUT)
        return payload

    # ---------- TemplateView overrides ----------
    def get(self, request, *args, **kwargs):
        """
        Serve HTML by default; JSON when requested; t-SNE via ?embed=1.
        Dataset is selected via ?state=inactive|active (default: inactive).
        """
        state_param = request.GET.get('state')
        if state_param in self.FILES:
            state = state_param
        else:
            state = self.DEFAULT_STATE

        want_json = (
            request.GET.get('format') == 'json'
            or request.GET.get('data') == '1'
            or 'application/json' in request.headers.get('Accept', '')
        )
        if want_json:
            want_embed = (request.GET.get('embed') in ('1', 'true', 'yes'))
            try:
                payload = self._build_payload(state=state, want_embed=want_embed)
            except Exception as e:
                return JsonResponse({"error": str(e)}, status=400)
            return JsonResponse(payload, safe=True)

        return super().get(request, *args, **kwargs)

    def get_context_data(self, **kwargs):
        """
        Expose separate URLs for inactive & active datasets so JS can toggle between them.
        """
        ctx = super().get_context_data(**kwargs)
        base = self.request.build_absolute_uri(self.request.path)

        ctx['inactive_embed_url'] = f"{base}?format=json&state=inactive&embed=1"
        ctx['active_embed_url']   = f"{base}?format=json&state=active&embed=1"

        return ctx
