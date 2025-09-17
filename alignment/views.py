from django.shortcuts import render, redirect
from django.conf import settings
from django.http import HttpResponse
from django.views.generic import TemplateView
from django.db.models import Q, F, Case, When, IntegerField, Prefetch
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

from collections import OrderedDict
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
    template_name = 'class_similarity/ClassesAndCounts.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        data_folder = 'protein_data'
        file_name = 'Ligand type update plus sense column.xlsx'
        file_path = os.path.join(settings.DATA_DIR, data_folder, file_name)

        try:
            df = pd.read_excel(file_path)
        except FileNotFoundError:
            context['error'] = f"File not found: {file_path}"
            return context

        df.columns = df.columns.str.strip().str.replace("\n", " ")

        # symbol -> ([Excel class names], GRAFS family, display name, fixed sensory value or None)
        mapping = {
            "A": (["Class A (Rhodopsin)"], "Rhodopsin", "Rhodopsin", None),
            "B1": (["Class B1 (Secretin)"], "Secretin", "Secretin", None),
            "B2": (["Class B2 (Adhesion)"], "Adhesion", "Adhesion", None),
            "C": (["Class C (Glutamate)"], "Glutamate", "Glutamate", None),
            "D": (["Class D (Fungal pheromone)"], "- (fungal)", "Fungal pheromone", None),
            "E": (["Class E (Yeast cAMP)"], "- (yeast)", "Yeast cAMP", "2"),  # fixed value
            "F": (["Class F (Frizzled)"], "Frizzled/Taste2", "Frizzled", None),
            "T2": (["Class T2 (Taste 2)"], "Frizzled/Taste2", "Taste 2", None),
            "OR": (["Class O1 (fish-like odorant)", "Class O2 (tetrapod specific odorant)"], "Rhodopsin", "Odorant (not olfactory)", None),
            "V?": (["Class V? (Vomeronasal/pheromone?)"], "- (non-functional in human)", "Vomeronasal or pheromone?", "5"),  # fixed value
            "Cl?": (["Other GPCRs"], "-", "Classless", None),
        }

        table_data = []

        for symbol, (class_names, grafs_family, display_name, fixed_sensory) in mapping.items():
            class_df = df[df["Class"].isin(class_names)]

            # counts
            receptor_families_count = class_df["Receptor family"].nunique()
            members_total_count = class_df["GPCRs (UniProt)"].nunique()
            non_sensory_count = (class_df["Sense"] == "non-sensory").sum()

            if fixed_sensory is not None:
                sensory_str = fixed_sensory
            else:
                sensory_counts = []
                for sense_type in ["vision", "taste", "odorant"]:
                    count = (class_df["Sense"] == sense_type).sum()
                    if count > 0:
                        sensory_counts.append(f"{count} {sense_type}")
                sensory_str = ", ".join(sensory_counts)

            orphan_count = (class_df["Sense"] == "unknown").sum()

            # entries for popups
            entries = {
                "receptor_families": sorted(class_df["Receptor family"].dropna().unique().tolist()),
                "members_total": class_df[["GPCRs (Gene name)", "Receptor family", "Ligand type", "Sense"]]
                                    .drop_duplicates().to_dict(orient="records"),
                "members_non_sensory": class_df[class_df["Sense"] == "non-sensory"]
                                        [["GPCRs (Gene name)", "Receptor family", "Ligand type", "Sense"]]
                                        .drop_duplicates().to_dict(orient="records"),
                "members_sensory": class_df[class_df["Sense"].isin(["vision", "taste", "odorant"])]
                                    [["GPCRs (Gene name)", "Receptor family", "Ligand type", "Sense"]]
                                    .drop_duplicates().to_dict(orient="records"),
                "orphan": class_df[class_df["Sense"] == "unknown"]
                           [["GPCRs (Gene name)", "Receptor family", "Ligand type", "Sense"]]
                           .drop_duplicates().to_dict(orient="records"),
            }

            row = {
                "symbol": symbol,
                "class_name": display_name,
                "grafs_family": grafs_family,
                "receptor_families_count": int(receptor_families_count),
                "members_total_count": int(members_total_count),
                "members_non_sensory_count": int(non_sensory_count),
                "members_sensory_count": sensory_str if sensory_str else 0,
                "orphan_count": int(orphan_count),
                "entries": entries,
            }

            table_data.append(row)

        # split into human vs non-human
        nonhuman_symbols = ["D", "E", "V?"]
        context["human_table_data"] = [row for row in table_data if row["symbol"] not in nonhuman_symbols]
        context["nonhuman_table_data"] = [row for row in table_data if row["symbol"] in nonhuman_symbols]

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
    
class OrphanSimilarity(TemplateView):
    template_name = 'class_similarity/OrphanSimilarity.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Build the orphan list once and embed it (93 items ≈ perfect for local Select2)
        orphans = (
            Protein.objects
            .filter(
                parent_id__isnull=True,
                species_id=1,
                family__parent__parent__name__iexact='Orphan receptors',
            )
            .select_related('family__parent__parent__parent')
            .annotate(
                family_name=F('family__parent__name'),
                ligand_type=F('family__parent__parent__name'),
                clazz=F('family__parent__parent__parent__name'),  # see note below
            )
            .values('id', 'entry_name', 'name', 'family_name', 'ligand_type', 'clazz')
            .order_by('entry_name')
        )

        # Shape for Select2: use only "name" as the visible text.
        context['orphans_select2'] = json.dumps([
            {
                "id": p["id"],
                "text": p["name"],                 # <- Select2 label and search text
                "entry_name": p["entry_name"],     # extra fields for your UI
                "name": p["name"],
                "family": p["family_name"],
                "ligand_type": p["ligand_type"],
                "class": p["clazz"],               # map clazz -> "class" for the client
            }
            for p in orphans
        ])
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
