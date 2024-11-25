from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings
from django.db.models import Count, Max, Q, F, Value, CharField, Case, When, IntegerField
from django.db.models import Count, Max
from django.core.cache import cache
from django.db import connection, reset_queries
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView
from common.views import AbsReferenceSelectionTable, getReferenceTable, getLigandTable, getLigandCountTable, AbsTargetSelection
from drugs.models import Drugs, Drugs2024, Indication
from protein.models import Protein, ProteinFamily, TissueExpression
from structure.models import Structure
from drugs.models import Drugs, Drugs2024, Indication, ATCCodes
from protein.models import Protein, ProteinFamily, Tissues, TissueExpression
from mutational_landscape.models import NHSPrescribings
from mapper.views import LandingPage

import re
import json
import numpy as np
from collections import OrderedDict
from copy import deepcopy
import pandas as pd
import os



def get_spaced_colors(n):
    max_value = 16581375 #255**3
    interval = int(max_value / n)
    colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

    return ['#'+i for i in colors] # HEX colors
    # return [(int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors] # RGB colors

def striphtml(data):
    p = re.compile(r'<.*?>')
    return p.sub('', data)

@cache_page(60 * 60 * 24 * 28)
def drugstatistics(request):

    # ===== drugtargets =====
    drugtargets_raw_approved = Protein.objects.filter(drugs__status='approved').values('entry_name').annotate(value=Count('drugs__name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugtargets_raw_approved))
    drugtargets_approved = []
    for i, drugtarget in enumerate(drugtargets_raw_approved):
        drugtarget['label'] = drugtarget['entry_name'].replace("_human","").upper()
        # drugtarget['color'] = str(list_of_hec_colors[i])
        del drugtarget['entry_name']
        drugtargets_approved.append(drugtarget)

    drugtargets_raw_trials = Protein.objects.filter(drugs__status__in=['in trial'], drugs__clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']).values('entry_name').annotate(value=Count('drugs__name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugtargets_raw_trials))
    drugtargets_trials = []
    for i, drugtarget in enumerate(drugtargets_raw_trials):
        drugtarget['label'] = drugtarget['entry_name'].replace("_human","").upper()
        # drugtarget['color'] = str(list_of_hec_colors[i])
        del drugtarget['entry_name']
        drugtargets_trials.append(drugtarget)


    all_human_GPCRs = Protein.objects.filter(species_id=1, sequence_type_id=1, family__slug__startswith='0').distinct()

    in_trial = Protein.objects.filter(drugs__status__in=['in trial'] ).exclude(drugs__status='approved').distinct() #drugs__clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']

    not_targeted = len(all_human_GPCRs) - len(drugtargets_approved) - len(in_trial)

    # ===== drugfamilies =====
    drugfamilies_raw_approved = Protein.objects.filter(drugs__status='approved').values('family_id__parent__name').annotate(value=Count('drugs__name', distinct = True))

    list_of_hec_colors = get_spaced_colors(len(drugfamilies_raw_approved))
    drugfamilies_approved = []
    for i, drugfamily in enumerate(drugfamilies_raw_approved):
        drugfamily['label'] = striphtml(drugfamily['family_id__parent__name']).replace(" receptors","")
        drugfamily['color'] = str(list_of_hec_colors[i])
        del drugfamily['family_id__parent__name']
        drugfamilies_approved.append(drugfamily)

    drugfamilies_raw_trials = Protein.objects.exclude(drugs__status='approved').values('family_id__parent__name').annotate(value=Count('drugs__name', distinct = True))

    list_of_hec_colors = get_spaced_colors(len(drugfamilies_raw_trials))
    drugfamilies_trials = []
    for i, drugfamily in enumerate(drugfamilies_raw_trials):
        drugfamily['label'] = striphtml(drugfamily['family_id__parent__name']).replace(" receptors","")
        drugfamily['color'] = str(list_of_hec_colors[i])
        del drugfamily['family_id__parent__name']
        drugfamilies_trials.append(drugfamily)

    # ===== drugclas =====
    drugclasses_raw_approved = Protein.objects.filter(drugs__status='approved').values('family_id__parent__parent__parent__name').annotate(value=Count('drugs__name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugclasses_raw_approved)+1)
    drugClasses_approved = []
    for i, drugclas in enumerate(drugclasses_raw_approved):
        drugclas['label'] = drugclas['family_id__parent__parent__parent__name']
        drugclas['color'] = str(list_of_hec_colors[i+1])
        del drugclas['family_id__parent__parent__parent__name']
        drugClasses_approved.append(drugclas)

    drugclasses_raw_trials = Protein.objects.filter(drugs__status__in=['in trial'], drugs__clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']).values('family_id__parent__parent__parent__name').annotate(value=Count('drugs__name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugclasses_raw_trials)+1)
    drugClasses_trials = []
    for i, drugclas in enumerate(drugclasses_raw_trials):
        drugclas['label'] = drugclas['family_id__parent__parent__parent__name']
        drugclas['color'] = str(list_of_hec_colors[i+1])
        del drugclas['family_id__parent__parent__parent__name']
        drugClasses_trials.append(drugclas)

    # ===== drugtypes =====
    drugtypes_raw_approved = Drugs.objects.values('drugtype').filter(status='approved').annotate(value=Count('name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugtypes_raw_approved)+20)
    drugtypes_approved = []
    for i, drugtype in enumerate(drugtypes_raw_approved):
        drugtype['label'] = drugtype['drugtype']
        drugtype['color'] = str(list_of_hec_colors[i])
        del drugtype['drugtype']
        drugtypes_approved.append(drugtype)

    drugtypes_raw_trials = Drugs.objects.values('drugtype').filter(status='in trial', clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']).annotate(value=Count('name', distinct = True)).order_by('-value')

    # list_of_hec_colors = get_spaced_colors(len(drugtypes_raw_trials)+5)
    drugtypes_trials = []
    for i, drugtype in enumerate(drugtypes_raw_trials):
        drugtype['label'] = drugtype['drugtype']
        drugtype['color'] = str(list_of_hec_colors[i])
        del drugtype['drugtype']
        drugtypes_trials.append(drugtype)

    drugtypes_raw_not_estab = Drugs.objects.values('drugtype').filter(novelty='not established').annotate(value=Count('name', distinct = True)).order_by('-value')

    # list_of_hec_colors = get_spaced_colors(len(drugtypes_raw_not_estab)+5)
    drugtypes_not_estab = []
    for i, drugtype in enumerate(drugtypes_raw_not_estab):
        drugtype['label'] = drugtype['drugtype']
        drugtype['color'] = str(list_of_hec_colors[i])
        del drugtype['drugtype']
        drugtypes_not_estab.append(drugtype)

    drugtypes_raw_estab = Drugs.objects.values('drugtype').filter(novelty='established').annotate(value=Count('name', distinct = True)).order_by('-value')

    # list_of_hec_colors = get_spaced_colors(len(drugtypes_raw_estab)+5)
    drugtypes_estab = []
    for i, drugtype in enumerate(drugtypes_raw_estab):
        drugtype['label'] = drugtype['drugtype']
        drugtype['color'] = str(list_of_hec_colors[i])
        del drugtype['drugtype']
        drugtypes_estab.append(drugtype)

    # ===== modes of action =====
    moas_raw_approved = Drugs.objects.values('moa').filter(status='approved').annotate(value=Count('name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(moas_raw_approved)+5)
    moas_approved = []
    for i, moa in enumerate(moas_raw_approved):
        moa['label'] = moa['moa']
        moa['color'] = str(list_of_hec_colors[i])
        del moa['moa']
        moas_approved.append(moa)

    moa_raw_trials = Drugs.objects.values('moa').filter(status='in trial', clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']).annotate(value=Count('name', distinct = True)).order_by('-value')

    # list_of_hec_colors = get_spaced_colors(len(moa_raw_trials)+5)
    moas_trials = []
    for i, moa in enumerate(moa_raw_trials):
        moa['label'] = moa['moa']
        moa['color'] = str(list_of_hec_colors[i])
        del moa['moa']
        moas_trials.append(moa)

    # ===== Phase distributions =====
    # Distinguish between different Clinical Status
    phases_raw_active = Drugs.objects.values('phase').filter(status='in trial', clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']).annotate(value=Count('name', distinct = True)).order_by('-value')

    phase_trials = []
    list_of_hec_colors = ["#88df8c", "#43A047", "#b0f2b2"]
    for i, phase in enumerate(phases_raw_active):
        phase['label'] = 'Phase ' + phase['phase']
        phase['color'] = str(list_of_hec_colors[i])
        del phase['phase']
        phase_trials.append(phase)

    phases_raw_inactive = Drugs.objects.values('phase').filter(status='in trial', clinicalstatus__in=['terminated','discontinued','unknown','withdrawn']).annotate(value=Count('name', distinct = True)).order_by('-value')

    phase_trials_inactive = []
    list_of_hec_colors = ["#88df8c", "#43A047", "#b0f2b2"]
    for i, phase in enumerate(phases_raw_inactive):
        phase['label'] = 'Phase ' + phase['phase']
        phase['color'] = str(list_of_hec_colors[i])
        del phase['phase']
        phase_trials_inactive.append(phase)

    # ===== drugindications =====
    drugindications_raw_approved = Drugs.objects.values('indication').filter(status='approved').annotate(value=Count('name', distinct = True)).order_by('-value')

    list_of_hec_colors = get_spaced_colors(len(drugindications_raw_approved)+10)
    drugindications_approved = []
    for i, drugindication in enumerate(drugindications_raw_approved):
        drugindication['label'] = drugindication['indication']
        drugindication['color'] = str(list_of_hec_colors[i])
        del drugindication['indication']
        drugindications_approved.append(drugindication)

    drugindications_raw_trials = Drugs.objects.values('indication').filter(status='in trial', clinicalstatus__in=['completed','not open yet','ongoing','recruiting','suspended']).annotate(value=Count('name', distinct = True)).order_by('-value')

    # list_of_hec_colors = get_spaced_colors(len(drugindications_raw_trials))
    drugindications_trials = []
    for i, drugindication in enumerate(drugindications_raw_trials):
        drugindication['label'] = drugindication['indication']
        drugindication['color'] = str(list_of_hec_colors[i])
        del drugindication['indication']
        drugindications_trials.append(drugindication)

    # ===== drugtimes =====
    drugtime_raw = Drugs.objects.values('approval').filter(status='approved').annotate(y=Count('name', distinct = True)).order_by('approval')

    drugtimes = []
    running_total = 0

    for i, time in enumerate(range(1942,2017,1)):
        if str(time) in [i['approval'] for i in drugtime_raw]:
            y = [i['y'] for i in drugtime_raw if i['approval']==str(time)][0] + running_total
            x = time
            running_total = y
        else:
            x = time
            y = running_total
        if time % 2 == 0:
            drugtimes.append({'x':x,'y':y})

    drugs_over_time = [{"values": drugtimes, "yAxis": "1", "key": "GPCRs"}, {'values': [{'y': 2, 'x': '1942'}, {'x': '1944', 'y': 2}, {'y': 6, 'x': '1946'}, {'y': 9, 'x': '1948'}, {'y': 18, 'x': '1950'}, {'y': 30, 'x': '1952'}, {'y': 55, 'x': '1954'}, {'y': 72, 'x': '1956'}, {'y': 98, 'x': '1958'}, {'y': 131, 'x': '1960'}, {'y': 153, 'x': '1962'}, {'y': 171, 'x': '1964'}, {'y': 188, 'x': '1966'}, {'y': 205, 'x': '1968'}, {'y': 224, 'x': '1970'}, {'y': 242, 'x': '1972'}, {'y': 265, 'x': '1974'}, {'y': 300, 'x': '1976'}, {'y': 340, 'x': '1978'}, {'y': 361, 'x': '1980'}, {'y': 410, 'x': '1982'}, {'y': 442, 'x': '1984'}, {'y': 499, 'x': '1986'}, {'y': 542, 'x': '1988'}, {'y': 583, 'x': '1990'}, {'y': 639, 'x': '1992'}, {'y': 686, 'x': '1994'}, {'y': 779, 'x': '1996'}, {'y': 847, 'x': '1998'}, {'y': 909, 'x': '2000'}, {'y': 948, 'x': '2002'}, {'y': 1003, 'x': '2004'}, {'y': 1041, 'x': '2006'}, {'y': 1078, 'x': '2008'}, {'y': 1115, 'x': '2010'}, {'y': 1177, 'x': '2012'}, {'y': 1239, 'x': '2014'}, {'y': 1286, 'x': '2016'}], 'key': 'All FDA drugs', 'yAxis': '1'}]

    # ===== drugtimes =====


    return render(request, 'drugstatistics.html', {'drugtypes_approved':drugtypes_approved, 'drugtypes_trials':drugtypes_trials,  'drugtypes_estab':drugtypes_estab,  'drugtypes_not_estab':drugtypes_not_estab, 'drugindications_approved':drugindications_approved, 'drugindications_trials':drugindications_trials, 'drugtargets_approved':drugtargets_approved, 'drugtargets_trials':drugtargets_trials, 'phase_trials':phase_trials, 'phase_trials_inactive': phase_trials_inactive, 'moas_trials':moas_trials, 'moas_approved':moas_approved, 'drugfamilies_approved':drugfamilies_approved, 'drugfamilies_trials':drugfamilies_trials, 'drugClasses_approved':drugClasses_approved, 'drugClasses_trials':drugClasses_trials, 'drugs_over_time':drugs_over_time, 'in_trial':len(in_trial), 'not_targeted':not_targeted})

class DrugBrowser(TemplateView):
    template_name = 'drugbrowser.html'

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        drugs = Drugs.objects.all().prefetch_related('target', 'target__family__parent__parent__parent', 'publication', 'publication__web_link', 'publication__web_link__web_resource', 'publication__journal')
        drugs_NHS_names = list(NHSPrescribings.objects.values_list('drugname__name', flat=True).distinct())
        context_data = list()

        def get_pmid(publication):
            try:
                pmid = publication.web_link.index if publication.web_link.web_resource.slug == "pubmed" else None
            except AttributeError:
                pmid = None
            return pmid

        # Create the string for the References column under APA rules
        def format_title(title):
            return re.sub(r'[\'"]', '-', title) if title else ''

        def format_authors(authors):
            return f"<b>{authors},</b>" if authors else ''

        def format_year(year):
            return f"<b>({year})</b><br/>" if year else ''

        def format_journal(journal):
            return f"<i>{journal.name}</i>" if journal else ''

        def format_volume_and_pages(reference):
            if reference and ':' in reference:
                volume, pages = reference.split(':')
                return f"<b>{volume}:</b><b>{pages}</b>"
            else:
                return ''

        def format_pmid_link(pmid):
            return f"[PMID: <a href='https://pubmed.ncbi.nlm.nih.gov/{pmid}' target='_blank'>{pmid}</a>]<br/>" if pmid else ''

        def format_publication(pub):
            pmid = get_pmid(pub)
            return f"{format_authors(pub.authors)} {format_year(pub.year)} {format_title(pub.title)}<br/> {format_journal(pub.journal)}, {format_volume_and_pages(pub.reference)} {format_pmid_link(pmid)}"

        for drug in drugs:
            drugname = drug.name
            NHS = 'yes' if drugname in drugs_NHS_names else 'no'
            target_list = drug.target.all()

            publications = drug.publication.all()
            publication_info = [format_publication(pub) for pub in publications]

            publication_info_string = '<br>'.join(publication_info)

            regex = r'\bClass\b' # Regular expression to extract text after the word 'Class'

            for protein in target_list:

                # Get the hierarchical class name from the protein family
                hierarchical_class_name = str(protein.family.parent.parent.parent.name)

                # Remove the word "Class" from the hierarchical class name using the provided regular expression
                class_name = re.sub(regex, '', hierarchical_class_name)

                family = str(protein.family.parent.name)

                jsondata = {
                    'name': drugname,
                    'target': str(protein),
                    'phase': drug.phase,
                    'approval': drug.approval,
                    'class': class_name,
                    'family': family,
                    'indication': drug.indication,
                    'status': drug.status,
                    'drugtype': drug.drugtype,
                    'moa': drug.moa,
                    'novelty': drug.novelty,
                    'targetlevel': drug.targetlevel,
                    'clinicalstatus': drug.clinicalstatus,
                    'NHS': NHS,
                    'publications': publication_info_string
                }

                context_data.append(jsondata)

        context['drugdata'] = context_data

        return context

def DrugsVenn(request):
    return Venn(request, "drugs")

def TargetVenn(request):
    return Venn(request, "targets")


def Venn(request, origin="both"):
    name_of_cache = 'venn_' + origin
    context = cache.get(name_of_cache)
    # NOTE cache disabled for development only!
    # context = None
    if context == None:
        context = OrderedDict()
    # Here we need to generate fata for three different Venn diagrams plus associated tables:
    # if origin is drugs, we need a Venn diagram showing drugs across different clinical phases
    # plus one comparing drugs that are in phase 1-3 and those in phase 4, and potential overlap
    # if origin is targets, we need a single Venn diagram showing targets across different clinical phases
    # the Venn design has to be the same of the SignProt pages
        phases_dict = {}
        key_mapping = {
            1: "Phase I",
            2: "Phase II",
            3: "Phase III",
            4: "Phase IV"
        }
        if origin == "drugs":
            # Call to get drugs in each maximum phase
            drug_phases = Drugs2024.objects.all().values_list('indication_max_phase','ligand_id__name').distinct()
            for item in drug_phases:
                if item[0] not in phases_dict.keys():
                    phases_dict[item[0]] = []
                phases_dict[item[0]].append(item[1])
            phases_dict = {key_mapping[k]: phases_dict[k] for k in key_mapping if k in phases_dict}
            for key in phases_dict.keys():
                phases_dict[key] = '\n'.join(phases_dict[key])
        else:
            # Call to get receptors in each maximum phase
            receptor_phases = Drugs2024.objects.all().values_list('indication_max_phase','target_id__entry_name').distinct()
            for item in receptor_phases:
                if item[0] not in phases_dict.keys():
                    phases_dict[item[0]] = []
                phases_dict[item[0]].append(item[1])
            phases_dict = {key_mapping[k]: phases_dict[k] for k in key_mapping if k in phases_dict}
            for key in phases_dict.keys():
                phases_dict[key] = '\n'.join(phases_dict[key])

        context["phases_dict"] = phases_dict
        context["phases_dict_keys"] = list(phases_dict.keys())

        #  Fetch table data with all related information
        table_data = Drugs2024.objects.select_related(
            'target__family__parent__parent__parent',  # All target info
            'ligand__ligand_type',
            'indication__code',
            'moa',
            'disease_association'
        ).values(
            'target', # target ID
            'target__entry_name',  # Gene name
            'target__name',  # Protein name
            'target__family__parent__name',  # Receptor family
            'target__family__parent__parent__name',  # Ligand type
            'target__family__parent__parent__parent__name',  # Class
            'ligand__id', #Ligand ID
            'ligand__name',  # Agent/Drug
            'ligand__ligand_type__name',  # Modality
            'moa__name',  # Mode of action
            'indication__title',  # Disease name
            'indication__code',  # Disease ICD11 code
            'indication_max_phase',  # Max phase
            'disease_association__association_score', # Disease association score
            'drug_status',  # Approval
        )

        # Convert the table_data queryset to a list of dictionaries
        table_data_list = list(table_data)

        # Convert the list of dictionaries to a pandas DataFrame
        df = pd.DataFrame(table_data_list)

        # Rename the columns to your desired format
        df.rename(columns={
            'target': 'TargetID',
            'target__entry_name': "Gene name",
            'target__name': "Protein name",
            'target__family__parent__name': 'Receptor family',
            'target__family__parent__parent__name': 'Ligand type',
            'target__family__parent__parent__parent__name': 'Class',
            'ligand__id': 'LigandID',
            'ligand__name': 'Ligand name',
            'ligand__ligand_type__name': 'Drug type',
            'moa__name': 'Modality',
            'indication__title': 'Indication name',
            'indication__code': 'ICD11',
            'indication_max_phase': 'Phase',
            'disease_association__association_score' : 'Association score',
            'drug_status': 'Approved'
        }, inplace=True)
            

        # Convert 'Approved' from integer to 'Yes'/'No'
        df['Approved'] = df['Approved'].apply(lambda x: 'Yes' if x == "Approved" else 'No')

        # Drop rows with missing LigandID
        df.dropna(subset=['LigandID'], inplace=True)

        # Fill missing values for other columns
        df.fillna({
            'Ligand name': 'Unknown',  # Fill missing Ligand name with 'Unknown'
            'Phase': 'N/A',  # Replace missing phases with 'N/A'
            'Approved': 'No',  # Default approval status to 'No'
            # Add additional columns as needed...
        }, inplace=True)

        # Convert DataFrame to JSON
        json_records = df.to_json(orient='records')

        # Pass the JSON data to the template context
        context['Full_data'] = json_records
    context["layout"] = origin

    return render(request,'venn_diagrams.html',context)

##########################################
#### Drugs, Indications, and Targets #####
##########################################

class DrugSectionSelection(TemplateView):
    # variable setup #
    template_name = 'common/selection_drugs.html'

    page = 'Drugs'
    title = "Drug search"
    description = 'Search by drug name'  # Default description

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Ensure that both title and description are initialized
        title = self.title
        description = self.description  # Set default description
        page = self.page

         # Modify the description and fetch data based on the page type
        if page == 'Drugs':
            description = 'Search by agent / drug name'
            # Fetch distinct ligands and create a dictionary of {ligand.name: ligand.id}
            search_data = Drugs2024.objects.all().prefetch_related('ligand').distinct('ligand')
            search_dict = {drug.ligand.name: drug.ligand.id for drug in search_data}

            # Fetch ATC codes for the table
            ATC_data = ATCCodes.objects.select_related(
                'code'
            ).values(
                'ligand', # Agent/Drug id
                'code__index' # ATC codes
            )

            # Convert ATC_data queryset to a list of dictionaries
            ATC_data_list = list(ATC_data)

            # Create a DataFrame from the ATC data
            atc_df = pd.DataFrame(ATC_data_list)

            # Group ATC codes by 'Ligand ID' and concatenate them
            atc_df_grouped = atc_df.groupby('ligand')['code__index'].apply(lambda x: ', '.join(x)).reset_index()

            # Rename columns for the ATC DataFrame
            atc_df_grouped.rename(columns={'ligand': 'LigandID', 'code__index': 'ATC'}, inplace=True)

            # Fetch table data with all related information
            table_data = Drugs2024.objects.select_related(
                'ligand', 
                'target__family__parent__parent__parent', # All target info
                'indication__code'
                'disease_association'
            ).values(
                'ligand', # Agent/Drug id
                'ligand__name', # Agent/Drug name
                'target__entry_name', # Gene name
                'target__name', # Protein name
                'target__family__parent__name',  # Receptor family
                'target__family__parent__parent__name',  # Ligand type
                'target__family__parent__parent__parent__name',  # Class
                'indication__title', # Disease name
                'indication__code', # Disease ICD11 code
                'indication_max_phase', # Max phase
                'drug_status', # Approval
                'disease_association__association_score' # Disease association score
            )

            # Convert the table_data queryset to a list of dictionaries
            table_data_list = list(table_data)

            # Convert the list of dictionaries to a pandas DataFrame
            df = pd.DataFrame(table_data_list)

            # Rename the columns to your desired format
            df.rename(columns={
                'ligand': 'LigandID',
                'ligand__name': 'Ligand name',
                'target__entry_name': 'Gene name',
                'target__name': 'Protein name',
                'target__family__parent__name': 'Receptor family',
                'target__family__parent__parent__name': 'Ligand type',
                'target__family__parent__parent__parent__name': 'Class',
                'indication__title': 'Indication name',
                'indication__code': 'ICD11',
                'indication_max_phase': 'Phase',
                'disease_association__association_score' : 'Association score',
                'drug_status': 'Status'}, inplace=True)

            # Merge the ATC data into the main DataFrame (df) on 'Ligand ID'
            df = df.merge(atc_df_grouped, on='LigandID', how='left')
            # Fill NaN values in the 'ATC' column with None
            df['ATC'] = df['ATC'].fillna("")
            
            # Precompute Phase and Status Information
            df['Is_Approved'] = (df['Status'] == 'Approved').astype(int)

            # Perform GroupBy and Aggregate
            Modified_df = df.groupby(
                ['LigandID', 'Gene name', 'Indication name', 'Ligand name', 'Protein name', 'Receptor family', 'Ligand type', 'Class', 'ICD11', 'ATC', 'Association score']
            ).agg(
                Highest_phase=('Phase', 'max'),  # Get the highest phase for each group
                Approved=('Status', lambda x: 1 if 'Approved' in x.values else 0)  # Mark as 1 if any row has 'Approved'
            ).reset_index()

            # Convert 'Approved' from integer to 'Yes'/'No'
            Modified_df['Approved'] = Modified_df['Approved'].apply(lambda x: 'Yes' if x == 1 else 'No')

            # Convert 'Approved' from integer to 'Yes'/'No'
            Modified_df['Approved'] = Modified_df['Approved'].apply(lambda x: 'Yes' if x == 1 else 'No')

            # Convert DataFrame to JSON
            json_records = Modified_df.to_json(orient='records')

            # Pass the JSON data to the template context
            context['Full_data'] = json_records

        elif page == 'Targets':
            description = 'Search by target name'
            # Fetch distinct targets and create a dictionary of {target.name: target.id}
            search_data = Drugs2024.objects.all().prefetch_related('target').distinct('target')
            search_dict = {drug.target.name: drug.target.id for drug in search_data}

            # Fetch ATC codes for the table
            ATC_data = ATCCodes.objects.select_related(
                'code'
            ).values(
                'ligand', # Agent/Drug id
                'code__index' # ATC codes
            )

            # Convert ATC_data queryset to a list of dictionaries
            ATC_data_list = list(ATC_data)

            # Create a DataFrame from the ATC data
            atc_df = pd.DataFrame(ATC_data_list)

            # Group ATC codes by 'Ligand ID' and concatenate them
            atc_df_grouped = atc_df.groupby('ligand')['code__index'].apply(lambda x: ', '.join(x)).reset_index()

            # Rename columns for the ATC DataFrame
            atc_df_grouped.rename(columns={'ligand': 'LigandID', 'code__index': 'ATC'}, inplace=True)


            # Fetch table data with all related information
            table_data = Drugs2024.objects.select_related(
                'target',
                'ligand__ligand_type',
                'indication__code',
                'moa',
                'disease_association'
            ).values(
                'target',  # Target ID
                'target__entry_name', # Gene Name
                'target__name',  # Target name
                'ligand__name',  # Agent/Drug
                'ligand__ligand_type__name',  # Modality
                'moa__name',  # Mode of action
                'indication__title',  # Disease name
                'indication__code',  # Disease ICD11 code
                'indication_max_phase',  # Max phase
                'ligand', # ligand ID
                'disease_association__association_score', # Disease association score
                'drug_status',  # Approval
            )

            # Convert the table_data queryset to a list of dictionaries
            table_data_list = list(table_data)

            # Convert the list of dictionaries to a pandas DataFrame
            df = pd.DataFrame(table_data_list)

            # Rename the columns to your desired format
            df.rename(columns={
                'target': 'Target ID',
                'target__entry_name': 'Gene name',
                'target__name': 'Target name',
                'ligand__name': 'Ligand name',
                'ligand__ligand_type__name': 'Modality',
                'moa__name': 'Mode of action',
                'indication__title': 'Indication name',
                'indication__code': 'ICD11',
                'indication_max_phase': 'Phase',
                'ligand': 'LigandID',
                'disease_association__association_score' : 'Association score',
                'drug_status': 'Status'
            }, inplace=True)
            
            # Merge the ATC data into the main DataFrame (df) on 'Ligand ID'
            df = df.merge(atc_df_grouped, on='LigandID', how='left')
            # Fill NaN values in the 'ATC' column with None
            df['ATC'] = df['ATC'].fillna("")

            # Precompute flags for phase and approval status
            df['Is_Approved'] = (df['Status'] == 'Approved').astype(int)

            # Group by necessary columns and perform aggregation in one go
            grouped = df.groupby(
                ['Target ID', 'Target name','Gene name', 'LigandID', 'Ligand name', 'Indication name', 'Modality', 'Mode of action', 'ICD11', 'ATC', 'Association score']
            )

            # Perform aggregation
            Modified_df = grouped.agg(
                Highest_phase=('Phase', 'max'),  # Get the highest phase for each group
                Approved=('Is_Approved', 'max')  # Check if any row has 'Approved' status (max of binary flag)
            ).reset_index()

            # Convert 'Approved' from integer to 'Yes'/'No'
            Modified_df['Approved'] = Modified_df['Approved'].apply(lambda x: 'Yes' if x == 1 else 'No')

            # Convert DataFrame to JSON
            json_records = Modified_df.to_json(orient='records')

            # Pass the JSON data to the template context
            context['Full_data'] = json_records


        elif page == 'Indications':
            description = 'Search by indication name'
            # Fetch distinct indications and create a dictionary of {indication.name: indication.id}
            search_data = Drugs2024.objects.all().prefetch_related('indication').distinct('indication')
            search_dict = {drug.indication.title: drug.indication.id for drug in search_data}

            # ###########################
            # Single Data Query
            # ###########################
            # Fetch all data in a single query
            table_data = Drugs2024.objects.select_related(
                'indication__code',
                'ligand__ligand_type',
                'target__family__parent__parent__parent',  # All target info
                'moa',
                'disease_association'
            ).values(
                'indication',  # Indication ID
                'indication__title',  # Indication name
                'indication__code',  # Disease ICD11 code
                'ligand__id',  # Ligand ID
                'ligand__name',  # Ligand name
                'indication_max_phase',  # Max phase
                'drug_status',  # Approval
                'ligand__ligand_type__name',  # Molecule type
                'moa__name',  # Mode of action
                'target__entry_name',  # Gene name
                'target__name',  # Protein name
                'target__family__parent__name',  # Receptor family
                'target__family__parent__parent__name',  # Ligand type
                'target__family__parent__parent__parent__name',  # Class
                'disease_association__association_score', # Disease association score
                "disease_association__ot_genetics_portal", # OT Genetics
                "disease_association__gene_burden", # Gene Burden
                "disease_association__eva", # ClinVar
                "disease_association__genomics_england", # GEL PanelApp
                "disease_association__gene2phenotype", # Gene2phenotype
                "disease_association__uniprot_literature", # UniProt literature
                "disease_association__uniprot_variants", # UniProt curated variants
                "disease_association__orphanet", # Orphanet
                "disease_association__clingen", # ClinGen
                "disease_association__cancer_gene_census", # Cancer Gene Census
                "disease_association__intogen", # IntOGen
                "disease_association__eva_somatic", # Clinvar (somatic)
                "disease_association__cancer_biomarkers", # Cancer Biomarkers
                "disease_association__chembl", # ChEMBL
                "disease_association__crispr_screen", # CRISPR Screens
                "disease_association__crispr", # Project Score
                "disease_association__slapenrich", # SLAPenrich
                "disease_association__reactome", # Reactome
                "disease_association__sysbio", # Gene signatures
                "disease_association__europepmc", # Europe PMC
                "disease_association__expression_atlas", # Expression Atlas
                "disease_association__impc" # IMPC
            )

            # Convert the table_data queryset to a list of dictionaries
            table_data_list = list(table_data)

            # Convert the list of dictionaries to a pandas DataFrame
            df = pd.DataFrame(table_data_list)

            # Rename the columns to your desired format
            df.rename(columns={
                'indication': 'Indication ID',
                'indication__title': 'Indication name',
                'indication__code': 'ICD11',
                'ligand__id': 'LigandID',
                'ligand__name': 'Drug name',
                'indication_max_phase': 'Phase',
                'drug_status': 'Status',
                'ligand__ligand_type__name': 'Molecule_type',
                'moa__name': 'Mode of action',
                'target__entry_name': 'Gene name',
                'target__name': 'Protein name',
                'target__family__parent__name': 'Receptor family',
                'target__family__parent__parent__name': 'Ligand type',
                'target__family__parent__parent__parent__name': 'Class',
                'disease_association__association_score' : 'Association score',
                "disease_association__ot_genetics_portal": "OT Genetics",
                "disease_association__gene_burden": "Gene Burden",
                "disease_association__eva": "ClinVar",
                "disease_association__genomics_england": "GEL PanelApp",
                "disease_association__gene2phenotype": "Gene2phenotype",
                "disease_association__uniprot_literature": "UniProt literature",
                "disease_association__uniprot_variants": "UniProt curated variants",
                "disease_association__orphanet": "Orphanet",
                "disease_association__clingen": "ClinGen",
                "disease_association__cancer_gene_census": "Cancer Gene Census",
                "disease_association__intogen": "IntOGen",
                "disease_association__eva_somatic": "Clinvar (somatic)",
                "disease_association__cancer_biomarkers": "Cancer Biomarkers",
                "disease_association__chembl":"ChEMBL",
                "disease_association__crispr_screen": "CRISPR Screens",
                "disease_association__crispr": "Project Score",
                "disease_association__slapenrich": "SLAPenrich",
                "disease_association__reactome": "Reactome",
                "disease_association__sysbio": "Gene signatures",
                "disease_association__europepmc": "Europe PMC",
                "disease_association__expression_atlas": "Expression Atlas",
                "disease_association__impc": "IMPC"
            }, inplace=True)

            # Define MOA categories
            stim_moa = ['Partial agonist', 'Agonist', 'PAM']
            inhib_moa = ['Antagonist', 'Inverse agonist', 'NAM']

            # Split the DataFrame into two: one for targets and one for drugs
            df_targets = df.copy()
            df_drugs = df.copy()

            # ###########################
            # Data Aggregation for Targets
            # ###########################
            # Update group_cols to include 'Indication name'
            group_cols = ['Indication ID', 'ICD11', 'Gene name', 'Indication name', 'Protein name', 'Receptor family', 'Ligand type', 'Class']

            # Define disease association columns to be added to the grouping
            disease_cols = [
                'Association score', 'OT Genetics', 'Gene Burden', 'ClinVar', 'GEL PanelApp', 
                'Gene2phenotype', 'UniProt literature', 'UniProt curated variants', 'Orphanet', 
                'ClinGen', 'Cancer Gene Census', 'IntOGen', 'Clinvar (somatic)', 'Cancer Biomarkers', 
                'ChEMBL', 'CRISPR Screens', 'Project Score', 'SLAPenrich', 'Reactome', 'Gene signatures', 
                'Europe PMC', 'Expression Atlas', 'IMPC'
            ]

            # Add these to group_cols
            group_cols += disease_cols

            df_targets[disease_cols] = df[disease_cols].fillna("")

            # Precompute classifications
            df_targets['Classification'] = df_targets['Status'].apply(lambda x: 'Drug' if x == 'Approved' else 'Agent')

            # Pre-filter the DataFrame for stimulatory and inhibitory MOAs
            stim_moa_df = df_targets[df_targets['Mode of action'].isin(stim_moa)]
            inhib_moa_df = df_targets[df_targets['Mode of action'].isin(inhib_moa)]


            # Create a function to get the max phase safely
            def get_max_phase(df, group_cols):
                if df.empty:
                    return pd.Series(dtype='float64')  # Return an empty Series if DataFrame is empty
                return df.groupby(group_cols)['Phase'].max()

            # Group the targets DataFrame by the updated columns
            grouped_targets = df_targets.groupby(group_cols)

            # Aggregate data using efficient vectorized operations
            agg_data_targets = grouped_targets.agg(
                All_max_phase=('Phase', 'max'),
                All_Drugs=('Classification', lambda x: (x == 'Drug').sum()),
                All_Agents=('Classification', lambda x: (x == 'Agent').sum())
            ).reset_index()

            # Compute the stimulatory and inhibitory max phase and counts efficiently
            stim_max_phase = get_max_phase(stim_moa_df, group_cols)
            inhib_max_phase = get_max_phase(inhib_moa_df, group_cols)

            # Add the precomputed data into the DataFrame
            agg_data_targets = agg_data_targets.merge(stim_max_phase.rename('Stimulatory_max_phase'), on=group_cols, how='left')
            agg_data_targets = agg_data_targets.merge(inhib_max_phase.rename('Inhibitory_max_phase'), on=group_cols, how='left')

           # Compute Stimulatory and Inhibitory Drugs and Agents counts
            # Reindex to ensure both 'Drug' and 'Agent' columns are present, even if they don't exist in the data
            stim_counts = stim_moa_df.groupby(group_cols)['Classification'].value_counts().unstack(fill_value=0).reindex(
                columns=['Drug', 'Agent'], fill_value=0).reset_index()

            # Rename columns to make sure they have meaningful names
            stim_counts.rename(columns={'Drug': 'Stimulatory_Drugs', 'Agent': 'Stimulatory_Agents'}, inplace=True)

            # Reindex to ensure both 'Drug' and 'Agent' columns are present for inhibitory counts
            inhib_counts = inhib_moa_df.groupby(group_cols)['Classification'].value_counts().unstack(fill_value=0).reindex(
                columns=['Drug', 'Agent'], fill_value=0).reset_index()

            # Rename columns for inhibitory MOAs
            inhib_counts.rename(columns={'Drug': 'Inhibitory_Drugs', 'Agent': 'Inhibitory_Agents'}, inplace=True)

            # Merge the counts back into the aggregated data
            agg_data_targets = agg_data_targets.merge(stim_counts[group_cols + ['Stimulatory_Drugs', 'Stimulatory_Agents']], on=group_cols, how='left').fillna(0)
            agg_data_targets = agg_data_targets.merge(inhib_counts[group_cols + ['Inhibitory_Drugs', 'Inhibitory_Agents']], on=group_cols, how='left').fillna(0)

            # Convert DataFrame to JSON
            json_records_targets = agg_data_targets.to_json(orient='records')
            context['Full_data_targets'] = json_records_targets

            # ###########################
            # Data Aggregation for Drugs
            # ###########################
            group_cols_drugs = ['Indication ID', 'ICD11', 'Gene name', 'Drug name', 'LigandID', 'Indication name', 'Protein name', 'Receptor family', 'Ligand type', 'Class', 'Molecule_type','Mode of action']

            # Precompute relevant columns
            df_drugs['Is_Approved'] = df_drugs['Status'].apply(lambda x: 1 if x == 'Approved' else 0)
            df_drugs['Is_Phase_I'] = (df_drugs['Phase'] == 1).astype(int)
            df_drugs['Is_Phase_II'] = (df_drugs['Phase'] == 2).astype(int)
            df_drugs['Is_Phase_III'] = (df_drugs['Phase'] == 3).astype(int)

            # Group the DataFrame only once
            grouped_drugs = df_drugs.groupby(group_cols_drugs)

            # Aggregate data using vectorized operations
            agg_data_drugs = grouped_drugs.agg(
                Highest_phase=('Phase', 'max'),
                Phase_I_trials=('Is_Phase_I', 'sum'),
                Phase_II_trials=('Is_Phase_II', 'sum'),
                Phase_III_trials=('Is_Phase_III', 'sum'),
                Approved=('Is_Approved', 'max')
            ).reset_index()

            # Convert 'Approved' from integer to 'Yes'/'No'
            agg_data_drugs['Approved'] = agg_data_drugs['Approved'].apply(lambda x: 'Yes' if x == 1 else 'No')

            # Convert DataFrame to JSON
            json_records_drugs = agg_data_drugs.to_json(orient='records')
            context['Full_data_drugs'] = json_records_drugs


        # Convert to JSON string and pass to context
        context['search_data'] = json.dumps(search_dict)
        context['page'] = page
        context['title'] = title
        context['description'] = description
        # context['table_data'] = table_data

        return context

class DruggedGPCRome(TemplateView):
    """
    Per class statistics of known ligands.
    """

    template_name = 'drugged_gpcrome.html'

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)

        drug_count_receptor_dict = {}

        # Use iterator to process the queryset in chunks
        drug_data = Drugs2024.objects.all().values_list('target__entry_name', 'indication_max_phase').distinct()

        all_proteins = Protein.objects.filter(species_id=1, parent_id__isnull=True, accession__isnull=False, family_id__slug__startswith='0').exclude(
                                            family_id__slug__startswith='007'
                                        ).exclude(
                                            family_id__slug__startswith='008'
                                        )
        for prot in all_proteins:
            drug_count_receptor_dict[prot.entry_name] = 0

        for pair in drug_data:
            if pair[0] in drug_count_receptor_dict.keys():
                if int(pair[1]) > int(drug_count_receptor_dict[pair[0]]):
                    drug_count_receptor_dict[pair[0]] = int(pair[1])

        proteins = list(Protein.objects.filter(entry_name__in=drug_count_receptor_dict.keys()
        ).values('entry_name', 'name').order_by('entry_name'))

        names_conversion_dict = {item['entry_name']: item['name'] for item in proteins}

        names = list(names_conversion_dict.values())
        IUPHAR_to_uniprot_dict = {item['name']: item['entry_name'] for item in proteins}

        families = ProteinFamily.objects.all()
        datatree = {}
        conversion = {}

        for item in families:
            if len(item.slug) == 3 and item.slug not in datatree.keys():
                datatree[item.slug] = {}
                conversion[item.slug] = item.name
            if len(item.slug) == 7 and item.slug not in datatree[item.slug[:3]].keys():
                datatree[item.slug[:3]][item.slug[:7]] = {}
                conversion[item.slug] = item.name
            if len(item.slug) == 11 and item.slug not in datatree[item.slug[:3]][item.slug[:7]].keys():
                datatree[item.slug[:3]][item.slug[:7]][item.slug[:11]] = []
                conversion[item.slug] = item.name
            if len(item.slug) == 15 and item.slug not in datatree[item.slug[:3]][item.slug[:7]][item.slug[:11]]:
                datatree[item.slug[:3]][item.slug[:7]][item.slug[:11]].append(item.name)

        datatree2 = LandingPage.convert_keys(datatree, conversion)
        datatree2.pop('Parent family', None)
        datatree3 = LandingPage.filter_dict(datatree2, names)
        data_converted = {names_conversion_dict[key]: {'Value1':value} for key, value in drug_count_receptor_dict.items()}
        Data_full = {"NameList": datatree3, "DataPoints": data_converted, "LabelConversionDict":IUPHAR_to_uniprot_dict}
        context['GPCRome_data'] = json.dumps(Data_full["NameList"])
        context['GPCRome_data_variables'] = json.dumps(Data_full['DataPoints'])
        context['GPCRome_Label_Conversion'] = json.dumps(Data_full['LabelConversionDict'])

        #TREE SECTION
        drug_data = Drugs2024.objects.all().values_list('target__entry_name', 'ligand__name','indication_max_phase')
        
        drug_dict = {}
        
        # Populate the dictionary
        for drug in drug_data:
            target, ligand, phase = drug
            if target not in drug_dict:
                # Initialize with empty lists
                drug_dict[target] = {'Outer1': [], 'Outer2': [], 'Outer3': [], 'Outer4': [], 'Inner': []}
            
            if phase in [1, 2, 3, 4]:
                outer_key = f"Outer{phase}"
                drug_dict[target][outer_key].append(ligand)  # Add ligand to the corresponding Outer list
            else:
                drug_dict[target]['Inner'].append(ligand)  # Add ligand to the Inner list

        # Convert lists to unique counts
        for target in drug_dict:
            for key in drug_dict[target]:
                unique_entries = set(drug_dict[target][key])  # Get unique ligands
                if target == 'drd2_human' and key == 'Outer4':
                    print(len(unique_entries),sorted(unique_entries))
                drug_dict[target][key] = len(unique_entries)  # Replace the list with the count

        tree, tree_options, circles, receptors = LandingPage.generate_tree_plot(drug_dict)
        #Remove 0 circles
        for key, outer_dict in circles.items():
            circles[key] = {k: v for k, v in outer_dict.items() if v != 0}

        context['tree'] = json.dumps(tree)
        context['tree_options'] = tree_options
        context['circles'] = json.dumps(circles)
        context['whole_dict'] = json.dumps(receptors)

        return context

class TargetSelectionTool(TemplateView):
    # Template using this class #
    template_name = 'TargetSelectionTool.html'
    # Get context for hmtl usage #
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        
        # Fetch all data in a single query
        table_data = Drugs2024.objects.select_related(
                'target__family__parent__parent__parent',  # All target info
                'moa'
            ).values(
                'target', # Target ID
                'target__entry_name',  # Gene name
                'target__name',  # Protein name
                'target__family__parent__name',  # Receptor family
                'target__family__parent__parent__name',  # Ligand type
                'target__family__parent__parent__parent__name',  # Class
                'ligand',
                'novelty_score', # Novelty score
                'drug_status',  # Approval
                'moa__name',  # Mode of action
                'indication_max_phase' # Max pahse
            )
        

        # Convert the table_data queryset to a list of dictionaries
        table_data_list = list(table_data)

        # Convert the list of dictionaries to a pandas DataFrame
        df = pd.DataFrame(table_data_list)

        # Drop duplicates based on all columns
        # df.drop_duplicates(inplace=True)

        # Rename the columns to your desired format
        df.rename(columns={
            'target': 'Target ID', # Target ID
            'target__entry_name': 'Gene name', # Gene name
            'target__name': 'Protein name', # Protein name
            'target__family__parent__name': 'Receptor family', # Receptor family
            'target__family__parent__parent__name': 'Ligand type', # Ligand type
            'target__family__parent__parent__parent__name': 'Class', # Class
            'ligand': "Ligand ID",
            'novelty_score': 'Novelty (Pharos)', #Only novelty we have
            'drug_status': 'Status', #Approval
            'moa__name': 'Mode of action', #Modality
            'indication_max_phase': 'Phase'
        }, inplace=True)

        # Assuming `target_ids` is a list of IDs from the Drugs data
        target_ids = df['Target ID'].unique()

        # Fetch only relevant `Structure` data for matching `target_ids`
        structure_data = Structure.objects.filter(
            protein_conformation__protein__parent__id__in=target_ids
        ).exclude(structure_type__slug__startswith='af-').values(
            'protein_conformation__protein__parent__id'
        ).annotate(
            Inactive=Count(Case(When(state_id=1, then=1), output_field=IntegerField())),
            Active=Count(Case(When(state_id=2, then=1), output_field=IntegerField())),
            Intermediate=Count(Case(When(state_id=3, then=1), output_field=IntegerField())),
            Total=Count('id')
        ).order_by('protein_conformation__protein__parent__id')

        # Convert structure data to a DataFrame
        structure_df = pd.DataFrame(list(structure_data))

        # Rename columns for better clarity
        structure_df.rename(columns={
            'protein_conformation__protein__parent__id': 'Target ID'
        }, inplace=True)

        # Merge structure counts into the drugs DataFrame on `Target ID`
        df = df.merge(structure_df, on='Target ID', how='left')

        # Fill any missing values for states with 0 (in case a target ID has no structure data)
        df[['Active', 'Inactive', 'Intermediate', 'Total']] = df[['Active', 'Inactive', 'Intermediate', 'Total']].fillna(0)

        # Define stimulatory and inhibitory MOA categories
        stim_moa = ['Partial agonist', 'Agonist', 'PAM']
        inhib_moa = ['Antagonist', 'Inverse agonist', 'NAM']

        # Precompute classifications: Label as 'Drug' if status is 'Approved', otherwise as 'Agent'
        df['Classification'] = df['Status'].apply(lambda x: 'Drug' if x == 'Approved' else 'Agent')

       # Precompute classifications: Label as 'Drug' if status is 'Approved', otherwise as 'Agent'
        df['Classification'] = df['Status'].apply(lambda x: 'Drug' if x == 'Approved' else 'Agent')

        # Filter data for stimulatory and inhibitory based on Mode of Action
        stim_moa = ['Partial agonist', 'Agonist', 'PAM']
        inhib_moa = ['Antagonist', 'Inverse agonist', 'NAM']

        stim_moa_df = df[df['Mode of action'].isin(stim_moa)]
        inhib_moa_df = df[df['Mode of action'].isin(inhib_moa)]

        # Define a helper function to compute the max phase safely
        def get_max_phase(df, group_cols):
            if df.empty:
                return pd.Series(dtype='float64')  # Return an empty Series if DataFrame is empty
            return df.groupby(group_cols)['Phase'].max()

        # Define a helper function to compute unique counts for drugs/agents
        def compute_unique_counts(df, group_cols, value_col, classification_col):
            """
            Aggregates the unique count of `value_col` grouped by `group_cols` 
            and splits them by their classification.
            """
            grouped = df.groupby(group_cols)[[value_col, classification_col]].apply(
                lambda x: pd.Series({
                    'Drug': len(set(x.loc[x[classification_col] == 'Drug', value_col])),
                    'Agent': len(set(x.loc[x[classification_col] == 'Agent', value_col]))
                })
            )
            return grouped.reset_index()

        # Group columns for aggregation
        group_cols = ['Target ID']

        # Compute All_Max_Phase for each Target ID as the maximum phase across all drugs and agents
        all_max_phase = df.groupby(group_cols)['Phase'].max().reset_index().rename(columns={'Phase': 'All_Max_Phase'})

        # Compute stimulatory max phases
        stim_max_phase = get_max_phase(stim_moa_df, group_cols).rename('Stimulatory_max_phase')

        # Compute inhibitory max phases
        inhib_max_phase = get_max_phase(inhib_moa_df, group_cols).rename('Inhibitory_max_phase')

        # Compute unique stimulatory counts
        stim_counts_unique = compute_unique_counts(
            stim_moa_df, group_cols, 'Ligand ID', 'Classification'
        ).rename(columns={'Drug': 'Stimulatory_Drugs', 'Agent': 'Stimulatory_Agents'})

        # Compute unique inhibitory counts
        inhib_counts_unique = compute_unique_counts(
            inhib_moa_df, group_cols, 'Ligand ID', 'Classification'
        ).rename(columns={'Drug': 'Inhibitory_Drugs', 'Agent': 'Inhibitory_Agents'})

        # Merge these into the aggregated data
        agg_data_targets = all_max_phase.merge(stim_max_phase, on=group_cols, how='left')
        agg_data_targets = agg_data_targets.merge(inhib_max_phase, on=group_cols, how='left')
        agg_data_targets = agg_data_targets.merge(stim_counts_unique, on=group_cols, how='left').fillna(0)
        agg_data_targets = agg_data_targets.merge(inhib_counts_unique, on=group_cols, how='left').fillna(0)

        # Compute total unique drugs and agents
        agg_data_targets['All_Drugs'] = agg_data_targets['Stimulatory_Drugs'] + agg_data_targets['Inhibitory_Drugs']
        agg_data_targets['All_Agents'] = agg_data_targets['Stimulatory_Agents'] + agg_data_targets['Inhibitory_Agents']

        # Ensure NaNs are replaced with 0 for all aggregated fields
        agg_data_targets.fillna({
            'Stimulatory_max_phase': 0, 'Inhibitory_max_phase': 0,
            'Stimulatory_Drugs': 0, 'Stimulatory_Agents': 0,
            'Inhibitory_Drugs': 0, 'Inhibitory_Agents': 0,
            'All_Drugs': 0, 'All_Agents': 0,
        }, inplace=True)

        # Merge `agg_data_targets` back into the original DataFrame on `Target ID`
        df_final = df.merge(agg_data_targets, on=group_cols, how='left')
        
        # Define the columns to add with "Coming soon" as the default value
        coming_soon_columns = [
            "No. Publications", "IDG target level", "Novelty (OpenTargets)",
            "A No. Ligands", "B No. Ligands", "C No. Ligands"
        ]

        # Add each column to df_final with "Coming soon" as the value for all rows
        for column in coming_soon_columns:
            df_final[column] = "Coming soon"

        # Keep only the specified columns in df_final
        keep_col_names = [
            'Target ID', 'Gene name', 'Protein name', 'Receptor family', 'Ligand type',
            'Class', 'No. Publications', 'Novelty (Pharos)', 'IDG target level', 'Novelty (OpenTargets)',
            'Total', 'Active', 'Inactive', 'All_Max_Phase', 'All_Drugs', 'All_Agents', 'A No. Ligands',
            'Stimulatory_max_phase', 'Stimulatory_Drugs', 'Stimulatory_Agents', 'B No. Ligands',
            'Inhibitory_max_phase', 'Inhibitory_Drugs', 'Inhibitory_Agents', 'C No. Ligands'
        ]

        # Keep only the specified columns in df_final
        df_final = df_final[keep_col_names]

        # Drop duplicates based on the reduced set of columns
        df_final.drop_duplicates(inplace=True)

        # Convert the final DataFrame to JSON
        json_records_targets = df_final.to_json(orient='records')
        context['Full_data'] = json_records_targets
        return context

        # # Get data - server side - Queries #
        # TissueExp = TissueExpression.objects.all().prefetch_related('protein','tissue')
        # Target_drug_data = Drugs2024.objects.all().prefetch_related('target','ligand','indication','indication__code')
        # target_ids = [drug.target.id for drug in Target_drug_data if drug.target]
        # #Structures_data = Structure.objects.all().prefetch_related('state', 'protein_conformation__protein')
        # Structure_data = Structure.objects.all().exclude(structure_type__slug__startswith='af-').values('protein_conformation__protein__parent__id').annotate(
        #             Inactive=Count(Case(When(state_id=1, then=1), output_field=IntegerField())),
        #             Active=Count(Case(When(state_id=2, then=1), output_field=IntegerField())),
        #             Intermediate=Count(Case(When(state_id=3, then=1), output_field=IntegerField())),
        #             Other=Count(Case(When(state_id=4, then=1), output_field=IntegerField())),
        #         ).order_by('protein_conformation__protein__parent__id') # Fetches only human structures --> make sure we dont want other species


        # #'target__family__parent__parent','target__family__parent__parent__parent','indication','indication__code','ligand','moa'
        # # Context lists for pushing data #
        # context_target_selection = list()
        # context_data_tissue = list()
        # # Dicts to modulate the data from the server side #
        # structure_dict = {}
        # target_selection_dict = {}
        # target_indication_dict = {}
        # Tissue_expression_dict = {}
        # index_dict = {}
        # Drugs_target_ids = []
        # # Create dict for structure total and state (active, inactive and intermediate)
        # for entry in Structure_data:
        #     protein_id = str(entry['protein_conformation__protein__parent__id'])
        #     if protein_id in target_ids:
        #         if protein_id not in structure_dict:
        #             structure_dict[protein_id] = {}
        #             structure_dict[protein_id]['Active'] = entry['Active']
        #             structure_dict[protein_id]['Inactive'] = entry['Inactive']
        #             structure_dict[protein_id]['Intermediate'] = entry['Intermediate']
        #             structure_dict[protein_id]['Total'] = int(entry['Active'])+entry['Inactive']+entry['Intermediate']
        #         else:
        #             print("Something should be terribly wrong if there is two identical ids")
        #             break
        #     else:
        #         structure_dict[protein_id] = {}
        #         structure_dict[protein_id]['Active'] = 0
        #         structure_dict[protein_id]['Inactive'] = 0
        #         structure_dict[protein_id]['Intermediate'] = 0
        #         structure_dict[protein_id]['Total'] = 0
        # # Create Target selection browser #
        # for entry in Target_drug_data:
        #     # Ids and keys
        #     Protein_id = str(entry.target.id)
        #     Indication_id = str(entry.indication.code.index)
        #     Target_indication_pair = "{}___{}".format(Protein_id,Indication_id)
        #     # Static values #
        #     Protein_uniprot = str(entry.target.accession)
        #     Protein_name = str(entry.target.entry_name).split("_")[0].upper()
        #     Indication = str(entry.indication.name)
        #     Drug_name  = str(entry.ligand.name)
        #     Drug_status = str(entry.drug_status)
        #     Novelty_score = float(entry.novelty_score)
        #     # Dicts #
        #     if Target_indication_pair not in target_indication_dict:
        #         target_indication_dict[Target_indication_pair] = {}
        #         target_indication_dict[Target_indication_pair]['information'] = [Protein_id,Indication_id,Indication]
        #         target_indication_dict[Target_indication_pair]['Novelty_score'] = Novelty_score
        #         #target_indication_dict[Target_indication_pair]['Drugs_total'] = 1
        #         target_indication_dict[Target_indication_pair]['Drugs'] = []
        #         target_indication_dict[Target_indication_pair]['Drugs'].append(Drug_name)
        #         # Handle drugs
        #         target_indication_dict[Target_indication_pair]['Drug__status'] = {}
        #         target_indication_dict[Target_indication_pair]['Drug__status']['Active'] = 0
        #         target_indication_dict[Target_indication_pair]['Drug__status']['Approved'] = 0
        #         target_indication_dict[Target_indication_pair]['Drug__status']['Discontinued'] = 0
        #         target_indication_dict[Target_indication_pair]['Drug__status'][Drug_status] += 1
        #         # Handle structures #
        #         # if Protein_id in structure_dict:
        #         #     target_indication_dict[Target_indication_pair]['Structures'] = {}
        #         #     target_indication_dict[Target_indication_pair]['Structures']['Total'] = structure_dict[Protein_id]['Total']
        #         #     target_indication_dict[Target_indication_pair]['Structures']['Active'] = structure_dict[Protein_id]['Active']
        #         #     target_indication_dict[Target_indication_pair]['Structures']['Inactive'] = structure_dict[Protein_id]['Inactive']
        #         #     target_indication_dict[Target_indication_pair]['Structures']['Intermediate'] = structure_dict[Protein_id]['Intermediate']
        #         # else:
        #         #     print("Something fishy")
        #         #     target_indication_dict[Target_indication_pair]['Structures'] = {}
        #         #     target_indication_dict[Target_indication_pair]['Structures']['Total'] = 0
        #         #     target_indication_dict[Target_indication_pair]['Structures']['Active'] = 0
        #         #     target_indication_dict[Target_indication_pair]['Structures']['Inactive'] = 0
        #         #     target_indication_dict[Target_indication_pair]['Structures']['Intermediate'] = 0
        #     else:
        #         # Drugs #
        #         target_indication_dict[Target_indication_pair]['Drugs'].append(Drug_name)
        #         target_indication_dict[Target_indication_pair]['Drug__status'][Drug_status] += 1

        #     if Protein_id not in target_selection_dict:
        #         target_selection_dict[Protein_id] = [Protein_uniprot,Protein_name]
        # for key in target_indication_dict:
        #     key_id = str(target_indication_dict[key]['information'][0])
        #     Drugs_target_ids.append(key_id)
        #     jsondata_TargetSelectionTool = {
        #             'Index_number': key_id,
        #             'Target_name': target_selection_dict[key_id][1],
        #             'Target_uniprot': target_selection_dict[key_id][0],
        #             'Indication_name': target_indication_dict[key]['information'][2],
        #             'Indication_id': target_indication_dict[key]['information'][1],
        #             'Novelty_score': str(target_indication_dict[key]['Novelty_score']),
        #             'IDG': "Coming soon",
        #             'Drugs_approved_names': target_indication_dict[key]['Drugs'],
        #             'Drugs_total': int(len(target_indication_dict[key]['Drugs'])),
        #             'Drugs_approved': int(target_indication_dict[key]['Drug__status']['Approved']),
        #             'Drugs_in_trial': int(target_indication_dict[key]['Drug__status']['Active']),
        #             'Drugs_discontinued': int(target_indication_dict[key]['Drug__status']['Discontinued'])
        #             # 'Structures_total': int(target_indication_dict[key]['Structures']['Total']),
        #             # 'Structures_active' : int(target_indication_dict[key]['Structures']['Active']),
        #             # 'Structures_inactive' : int(target_indication_dict[key]['Structures']['Inactive']),
        #             # 'Structures_intermediate' : int(target_indication_dict[key]['Structures']['Intermediate'])
        #     }
        #     context_target_selection.append(jsondata_TargetSelectionTool)
        # context['Drug_Targets'] = context_target_selection
        # # Go through server side data and modulate into a dict #
        # for entry in TissueExp:
        #     # string values for Tissue expression table #
        #     protein_id = entry.protein.entry_name
        #     value = entry.value
        #     ### Get rid of the nan values ### <--- should be removed through the build
        #     if value != value:
        #         value = None
        #     Tissue_id = entry.tissue.name
        #     # Index key #
        #     index_key = str(entry.protein.id)
        #     if index_key not in Drugs_target_ids:
        #         pass
        #     else:
        #         if protein_id not in index_dict:
        #             index_dict[str(protein_id)] = index_key
        #         # Expression value linked to protein / target #
        #         if protein_id not in Tissue_expression_dict:
        #             Tissue_expression_dict[str(protein_id)] = {}
        #             Tissue_expression_dict[str(protein_id)][str(Tissue_id)] = float(value) if value else '-'
        #         else:
        #             Tissue_expression_dict[str(protein_id)][str(Tissue_id)] = float(value) if value else '-'
        # # Run through dict and assign the correct values into the context data #
        # for key in Tissue_expression_dict:
        #     jsondata_tissue = {
        #             'Index_number': index_dict[key],
        #             'ProteinID': key,
        #             'ProteinName': str(key).split("_")[0].upper(),
        #             'adipose_tissue': Tissue_expression_dict[key]['adipose_tissue'],
        #             'adrenal_gland': Tissue_expression_dict[key]['adrenal_gland'],
        #             'amygdala': Tissue_expression_dict[key]['amygdala'],
        #             'appendix': Tissue_expression_dict[key]['appendix'],
        #             'basal_ganglia': Tissue_expression_dict[key]['basal_ganglia'],
        #             'bone_marrow': Tissue_expression_dict[key]['bone_marrow'],
        #             'breast': Tissue_expression_dict[key]['breast'],
        #             'cerebellum': Tissue_expression_dict[key]['cerebellum'],
        #             'cerebral_cortex': Tissue_expression_dict[key]['cerebral_cortex'],
        #             'cervix': Tissue_expression_dict[key]['cervix'],
        #             'choroid_plexus': Tissue_expression_dict[key]['choroid_plexus'],
        #             'colon': Tissue_expression_dict[key]['colon'],
        #             'duodenum': Tissue_expression_dict[key]['duodenum'],
        #             'endometrium': Tissue_expression_dict[key]['endometrium_1'], #Should be updated
        #             'epididymis': Tissue_expression_dict[key]['epididymis'],
        #             'esophagus': Tissue_expression_dict[key]['esophagus'],
        #             'fallopian_tube': Tissue_expression_dict[key]['fallopian_tube'],
        #             'gallbladder': Tissue_expression_dict[key]['gallbladder'],
        #             'heart_muscle': Tissue_expression_dict[key]['heart_muscle'],
        #             'hippocampal_formation': Tissue_expression_dict[key]['hippocampal_formation'],
        #             'hypothalamus': Tissue_expression_dict[key]['hypothalamus'],
        #             'kidney': Tissue_expression_dict[key]['kidney'],
        #             'liver': Tissue_expression_dict[key]['liver'],
        #             'lung': Tissue_expression_dict[key]['lung'],
        #             'lymph_node': Tissue_expression_dict[key]['lymph_node'],
        #             'midbrain': Tissue_expression_dict[key]['midbrain'],
        #             'ovary': Tissue_expression_dict[key]['ovary'],
        #             'pancreas': Tissue_expression_dict[key]['pancreas'],
        #             'parathyroid_gland': Tissue_expression_dict[key]['parathyroid_gland'],
        #             'pituitary_gland': Tissue_expression_dict[key]['pituitary_gland'],
        #             'placenta': Tissue_expression_dict[key]['placenta'],
        #             'prostate': Tissue_expression_dict[key]['prostate'],
        #             'rectum': Tissue_expression_dict[key]['rectum'],
        #             'retina': Tissue_expression_dict[key]['retina'],
        #             'salivary_gland': Tissue_expression_dict[key]['salivary_gland'],
        #             'seminal_vesicle': Tissue_expression_dict[key]['seminal_vesicle'],
        #             'skeletal_muscle': Tissue_expression_dict[key]['skeletal_muscle'],
        #             'skin': Tissue_expression_dict[key]['skin_1'],
        #             'small_intestine': Tissue_expression_dict[key]['small_intestine'],
        #             'smooth_muscle': Tissue_expression_dict[key]['smooth_muscle'],
        #             'spinal_cord': Tissue_expression_dict[key]['spinal_cord'],
        #             'spleen': Tissue_expression_dict[key]['spleen'],
        #             'stomach': Tissue_expression_dict[key]['stomach_1'],
        #             'testis': Tissue_expression_dict[key]['testis'],
        #             'thymus': Tissue_expression_dict[key]['thymus'],
        #             'thyroid_gland': Tissue_expression_dict[key]['thyroid_gland'],
        #             'tongue': Tissue_expression_dict[key]['tongue'],
        #             'tonsil': Tissue_expression_dict[key]['tonsil'],
        #             'urinary_bladder': Tissue_expression_dict[key]['urinary_bladder'],
        #             'vagina': Tissue_expression_dict[key]['vagina']
        #         }
        #     # Append context data into list #
        #     context_data_tissue.append(jsondata_tissue)
        # # Create context data for tissue expression data #
        # context['Tissue_data'] = context_data_tissue
        # context['Tissue_header_list'] = ['Adipose tissue',
        #                                'Adrenal gland',
        #                                'Amygdala',
        #                                'Appendix',
        #                                'Basal ganglia',
        #                                'Bone marrow',
        #                                'Breast',
        #                                'Cerebellum',
        #                                'Cerebral cortex',
        #                                'Cervix',
        #                                'Choroid plexus',
        #                                'Colon',
        #                                'Duodenum',
        #                                'Endometrium',
        #                                'Epididymis',
        #                                'Esophagus',
        #                                'Fallopian tube',
        #                                'Gallbladder',
        #                                'Heart muscle',
        #                                'Hippocampal formation',
        #                                'Hypothalamus',
        #                                'Kidney',
        #                                'Liver',
        #                                'Lung',
        #                                'Lymph node',
        #                                'Midbrain',
        #                                'Ovary',
        #                                'Pancreas',
        #                                'Parathyroid gland',
        #                                'Pituitary gland',
        #                                'Placenta',
        #                                'Prostate',
        #                                'Rectum',
        #                                'Retina',
        #                                'Salivary gland',
        #                                'Seminal vesicle',
        #                                'Skeletal muscle',
        #                                'Skin',
        #                                'Small intestine',
        #                                'Smooth muscle',
        #                                'Spinal cord',
        #                                'Spleen',
        #                                'Stomach',
        #                                'Testis',
        #                                'Thymus',
        #                                'Thyroid gland',
        #                                'Tongue',
        #                                'Tonsil',
        #                                'Urinary bladder',
        #                                'Vagina']
        # context['Tissue_header_list'] = [x.replace(' ', '<br>') for x in context['Tissue_header_list']]
        # context['Structure_data'] = Structure_data
        # # Lastly return context for html usage #
        # return context








@cache_page(60 * 60 * 24 * 28)
def drugmapping(request):
    context = dict()

    families = ProteinFamily.objects.all()
    lookup = {}
    for f in families:
        lookup[f.slug] = f.name.replace("receptors","").replace(" receptor","").replace(" hormone","").replace("/neuropeptide","/").replace(" (G protein-coupled)","").replace(" factor","").replace(" (LPA)","").replace(" (S1P)","").replace("GPR18, GPR55 and GPR119","GPR18/55/119").replace("-releasing","").replace(" peptide","").replace(" and oxytocin","/Oxytocin").replace("Adhesion class orphans","Adhesion orphans").replace("muscarinic","musc.").replace("-concentrating","-conc.")

    class_proteins = Protein.objects.filter(family__slug__startswith="0",source__name='SWISSPROT', species_id=1).prefetch_related('family').order_by('family__slug')

    temp = OrderedDict([
                    ('name',''),
                    ('trials', 0),
                    ('maxphase', 0),
                    ('approved', 0),
                    ('family_sum_approved', 0),
                    ('family_sum_trials' , 0),
                    ('establishment', 2),
                    ('children', OrderedDict())
                    ])

    coverage = OrderedDict()

    # Make the scaffold
    for p in class_proteins:
        #print(p,p.family.slug)
        fid = p.family.slug.split("_")
        if fid[0] not in coverage:
            coverage[fid[0]] = deepcopy(temp)
            coverage[fid[0]]['name'] = lookup[fid[0]]
        if fid[1] not in coverage[fid[0]]['children']:
            coverage[fid[0]]['children'][fid[1]] = deepcopy(temp)
            coverage[fid[0]]['children'][fid[1]]['name'] = lookup[fid[0]+"_"+fid[1]]
        if fid[2] not in coverage[fid[0]]['children'][fid[1]]['children']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]] = deepcopy(temp)
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['name'] = lookup[fid[0]+"_"+fid[1]+"_"+fid[2]][:28]
        if fid[3] not in coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]] = deepcopy(temp)
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['name'] = p.entry_name.split("_")[0] #[:10]

    # # POULATE WITH DATA
    total_approved = 0
    drugtargets_approved_class = Protein.objects.filter(drugs__status='approved').values('family_id__parent__parent__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_approved_class:
        fid = i['family_id__parent__parent__parent__slug'].split("_")
        coverage[fid[0]]['family_sum_approved'] += i['value']
        total_approved += i['value']

    drugtargets_approved_type = Protein.objects.filter(drugs__status='approved').values('family_id__parent__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_approved_type:
        fid = i['family_id__parent__parent__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['family_sum_approved'] += i['value']

    drugtargets_approved_family = Protein.objects.filter(drugs__status='approved').values('family_id__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_approved_family:
        fid = i['family_id__parent__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['family_sum_approved'] += i['value']
    drugtargets_approved_target = Protein.objects.filter(drugs__status='approved').values('family_id__slug').annotate(value=Count('drugs__name', distinct = True)).annotate(maxphase=Max('drugs__phase'))
    for i in drugtargets_approved_target:
        fid = i['family_id__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['approved'] += i['value']
        if int(i['maxphase']) > coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['maxphase']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['maxphase'] = int(i['maxphase'])
        if i['value'] > 0:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['establishment'] = 4

    total_trials = 0
    drugtargets_trials_class = Protein.objects.filter(drugs__status__in=['in trial']).values('family_id__parent__parent__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_trials_class:
        fid = i['family_id__parent__parent__parent__slug'].split("_")
        coverage[fid[0]]['family_sum_trials'] += i['value']
        total_trials += i['value']

    drugtargets_trials_type = Protein.objects.filter(drugs__status__in=['in trial']).values('family_id__parent__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_trials_type:
        fid = i['family_id__parent__parent__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['family_sum_trials'] += i['value']

    drugtargets_trials_family = Protein.objects.filter(drugs__status__in=['in trial']).values('family_id__parent__slug').annotate(value=Count('drugs__name', distinct = True))
    for i in drugtargets_trials_family:
        fid = i['family_id__parent__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['family_sum_trials'] += i['value']

    drugtargets_trials_target = Protein.objects.filter(drugs__status__in=['in trial']).values('family_id__slug').annotate(value=Count('drugs__name', distinct = True)).annotate(maxphase=Max('drugs__phase'))
    # add highest reached trial here
    for i in drugtargets_trials_target:
        fid = i['family_id__slug'].split("_")
        coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['trials'] += i['value']
        if int(i['maxphase']) > coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['maxphase']:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['maxphase'] = int(i['maxphase'])
        if i['value'] > 0 and coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['establishment'] == 2:
            coverage[fid[0]]['children'][fid[1]]['children'][fid[2]]['children'][fid[3]]['establishment'] = 7

    # MAKE THE TREE
    tree = OrderedDict({'name':'GPCRome', 'family_sum_approved': total_approved, 'family_sum_trials': total_trials,'children':[]})
    i = 0
    n = 0
    for c,c_v in coverage.items():
        c_v['name'] = c_v['name'].split("(")[0]
        if c_v['name'].strip() == 'Other GPCRs':
            # i += 1
            continue
            # pass
        children = []
        for lt,lt_v in c_v['children'].items():
            children_rf = []
            for rf,rf_v in lt_v['children'].items():
                rf_v['name'] = rf_v['name'].split("<")[0]
                # if rf_v['name'].strip() == 'Taste 2':
                    # continue
                children_r = []
                for r,r_v in rf_v['children'].items():
                    r_v['sort'] = n
                    children_r.append(r_v)
                    n += 1
                rf_v['children'] = children_r
                rf_v['sort'] = n
                children_rf.append(rf_v)
            lt_v['children'] = children_rf
            lt_v['sort'] = n
            children.append(lt_v)
        c_v['children'] = children
        c_v['sort'] = n
        tree['children'].append(c_v)
        #tree = c_v
        #break
        i += 1

    jsontree = json.dumps(tree)

    context["drugdata"] = jsontree

    return render(request, 'drugmapping.html', {'drugdata':context})

@cache_page(60 * 60 * 24 * 28)
def nhs_drug(request, slug):

    nhs_data = NHSPrescribings.objects.filter(drugname__name=slug.lower()).order_by('date')

    data_dic = {}
    sections = []
    query_translate = {}
    for i in nhs_data:
        prescription_name = i.op_name +' (' + i.drugCode + ')'
        queryname = i.drugname.name

        if not prescription_name in data_dic:
            data_dic[prescription_name] = []
            sections.append(i.bnf_section)
        dic = {}
        dic['x'] = str(i.date)
        dic['y'] = int(i.actual_cost)
        data_dic[prescription_name].append(dic)

        if not prescription_name in query_translate:
            query_translate[prescription_name] = queryname

    data = []
    for nhs_name in data_dic.keys():
        data.append({'values': data_dic[nhs_name], 'query_key':str(query_translate[nhs_name]), 'key':nhs_name})

    return render(request, 'nhs.html', {'data':data, 'drug':slug, 'section':list(set(sections))})

# @cache_page(60 * 60 * 24 * 28)
def indication_detail(request, code):

    code = code.upper()
    context = dict()
    #code = 'EFO_0003843'
    indication_data = Drugs2024.objects.filter(indication__code=code).prefetch_related('ligand',
                                                                                        'target',
                                                                                        'indication',
                                                                                        'indication__uri')

    indication_name = Indication.objects.filter(code=code).values_list('title', flat=True).distinct()[0]

    sankey = {"nodes": [],
              "links": []}
    caches = {'indication':[],
              'ligands': [],
              'targets': [],
              'entries': []}

    node_counter = 0
    for record in indication_data:
        #assess the values for indication/ligand/protein
        indication_code = record.indication.title.capitalize()
        indication_uri = record.indication.uri.index
        ligand_name = record.ligand.name.capitalize()
        ligand_id = record.ligand.id
        protein_name = record.target.name
        target_name = record.target.entry_name
        #check for each value if it exists and retrieve the source node value
        if indication_code not in caches['indication']:
            sankey['nodes'].append({"node": node_counter, "name": indication_code, "url":'https://icd.who.int/browse/2024-01/mms/en#'+indication_uri})
            node_counter += 1
            caches['indication'].append(indication_code)
        indi_node = next((item['node'] for item in sankey['nodes'] if item['name'] == indication_code), None)
        if [ligand_name, ligand_id] not in caches['ligands']:
            sankey['nodes'].append({"node": node_counter, "name": ligand_name, "url":'/ligand/'+str(ligand_id)+'/info'})
            node_counter += 1
            caches['ligands'].append([ligand_name, ligand_id])
        lig_node = next((item['node'] for item in sankey['nodes'] if item['name'] == ligand_name), None)
        if protein_name not in caches['targets']:
            sankey['nodes'].append({"node": node_counter, "name": protein_name, "url":'/protein/'+str(target_name)})
            node_counter += 1
            caches['targets'].append(protein_name)
            caches['entries'].append(target_name)
        prot_node = next((item['node'] for item in sankey['nodes'] if item['name'] == protein_name), None)
        #append connection between indication and ligand
        sankey['links'].append({"source":indi_node, "target":lig_node, "value":1, "ligtrace": ligand_name, "prottrace": None})
        #append connection between ligand and target
        sankey['links'].append({"source":lig_node, "target":prot_node, "value":1, "ligtrace": ligand_name, "prottrace": protein_name})

    #Fixing redundancy in sankey['links']
    unique_combinations = {}

    for d in sankey['links']:
        # Create a key based on source and target for identifying unique combinations
        key = (d['source'], d['target'])

        if key in unique_combinations:
            # If the combination exists, add the value to the existing entry
            unique_combinations[key]['value'] += d['value']
        else:
            # If it's a new combination, add it to the dictionary
            unique_combinations[key] = d

    # Convert the unique_combinations back to a list of dictionaries
    sankey['links'] = list(unique_combinations.values())
    total_points = len(caches['targets']) + len(caches['targets']) + 1;
    if len(caches['ligands']) > len(caches['targets']):
        context['nodes_nr'] = len(caches['ligands'])
    else:
        context['nodes_nr'] = len(caches['targets'])
    context['indication_code'] = code
    context['indication_uri'] = indication_uri
    context['indication'] = indication_name.capitalize()
    context['sankey'] = json.dumps(sankey)
    context['points'] = total_points
    context['targets'] = list(caches['entries'])
    context['ligands'] = list(caches['ligands'])
    return render(request, 'indication_detail.html', context)

@cache_page(60 * 60 * 24 * 28)
def nhs_section(request, slug):

    nhs_data = NHSPrescribings.objects.filter(bnf_section=slug).order_by('date')

    data_dic = {}
    sections = []
    query_translate = {}
    for i in nhs_data:
        prescription_name = i.op_name +' (' + i.drugCode + ')'
        queryname = i.drugname.name

        if not prescription_name in data_dic:
            data_dic[prescription_name] = []
            sections.append(i.bnf_section)

        dic = {}
        dic['x'] = str(i.date)
        dic['y'] = int(i.actual_cost)
        data_dic[prescription_name].append(dic)

        if not prescription_name in query_translate:
            query_translate[prescription_name] = queryname

    data = []
    for nhs_name in data_dic.keys():
        data.append({'values': data_dic[nhs_name], 'query_key':str(query_translate[nhs_name]), 'key':nhs_name})

    return render(request, 'nhs.html', {'data':data, 'drug':slug, 'section':list(set(sections))})
