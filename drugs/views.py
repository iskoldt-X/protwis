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
from drugs.models import Drugs, Drugs2024, Indication
from protein.models import Protein, ProteinFamily, Tissues, TissueExpression
from mutational_landscape.models import NHSPrescribings
from mapper.views import LandingPage

import re
import json
import numpy as np
from collections import OrderedDict
from copy import deepcopy



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


    all_human_GPCRs = Protein.objects.filter(species_id=1, sequence_type_id=1, family__slug__startswith='00').distinct()

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
        print('Check 1')
        phases_dict = {}
        key_mapping = {
            1: "Phase I",
            2: "Phase II",
            3: "Phase III",
            4: "Phase IV"
        }
        if origin == "drugs":
            print('Check 2')
            # Call to get drugs in each maximum phase
            drug_phases = Drugs2024.objects.all().values_list('indication_max_phase','ligand_id__name').distinct()
            for item in drug_phases:
                if item[0] not in phases_dict.keys():
                    phases_dict[item[0]] = []
                phases_dict[item[0]].append(item[1])
            phases_dict = {key_mapping[k]: phases_dict[k] for k in key_mapping if k in phases_dict}
            print('Check 3')
            for key in phases_dict.keys():
                phases_dict[key] = '\n'.join(phases_dict[key])
            print('Check 4')
        else:
            # Call to get receptors in each maximum phase
            receptor_phases = Drugs2024.objects.all().values_list('indication_max_phase','target_id__entry_name').distinct()
            for item in receptor_phases:
                if item[0] not in phases_dict.keys():
                    phases_dict[item[0]] = []
                phases_dict[item[0]].append(item[1].split("_")[0].upper())
            phases_dict = {key_mapping[k]: phases_dict[k] for k in key_mapping if k in phases_dict}
            for key in phases_dict.keys():
                phases_dict[key] = '\n'.join(phases_dict[key])

        print('Check 5')

        context["phases_dict"] = phases_dict
        context["phases_dict_keys"] = list(phases_dict.keys())

        # Collect drugs information
        drugs_panel = Drugs2024.objects.all().select_related(
                                                            "ligand",
                                                            "moa",
                                                            "target__family__parent",
                                                            "target",
                                                            "indication__code"
                                                            )
        print('Check 6')
        drug_dictionary = {}
        for p in drugs_panel:
            # Collect receptor data
            lig_name = p.ligand.name.lower().capitalize()
            lig_type = p.moa.name
            rec_family = p.target.family.parent.short()
            rec_uniprot = p.target.entry_short()
            rec_iuphar = p.target.family.name.replace("receptor", '').replace("<i>","").replace("</i>","").strip()
            clinical_phase = p.indication_max_phase
            indication_name = p.indication.name
            indication_code = p.indication.code.index
            # Create a tuple to store the values
            drug_entry = (lig_name, lig_type, rec_family, rec_uniprot, rec_iuphar, clinical_phase, indication_name, indication_code)
            # Initialize key if it doesn't exist
            if origin == 'drugs':
                key = lig_name
            else:
                key = rec_uniprot.lower().capitalize()
            if key not in drug_dictionary:
                drug_dictionary[key] = []
            # Check for duplicates before adding
            if drug_entry not in drug_dictionary[key]:
                drug_dictionary[key].append(drug_entry)

        print('Check 7')
        context["drug_dictionary"] = json.dumps(drug_dictionary)
        print(drug_dictionary)
    # cache.set(name_of_cache, context, 60 * 60 * 24 * 7)  # seven days timeout on cache
    context["layout"] = origin

    return render(request,'venn_diagrams.html',context)


class DrugSectionSelection(AbsTargetSelection):
    # Left panel
    page = 'drugs'
    step = 1
    number_of_steps = 1
    template_name = 'common/selection_drugs.html'
    filters = False
    import_export_box = False
    target_input = False
    psets = False
    family_tree = False
    type_of_selection = 'ligands'
    selection_only_receptors = False
    title = "Drug search"
    description = 'Search by drug name or database ID (GPCRdb, GtP, ChEMBL)'

    buttons = {
        'continue' : {
            'label' : 'Show ligand information',
            'url' : '',
            'color' : 'success',
            }
        }


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

        # tree = PhylogeneticTreeGenerator()
        # class_a_data = tree.get_tree_data(ProteinFamily.objects.get(name='Class A (Rhodopsin)'))
        # context['class_a_options'] = deepcopy(tree.d3_options)
        # context['class_a_options']['anchor'] = 'class_a'
        # context['class_a_options']['leaf_offset'] = 50
        # context['class_a_options']['label_free'] = []
        # # section to remove Orphan from Class A tree and apply to a different tree
        # whole_class_a = class_a_data.get_nodes_dict(self.page)
        # for item in whole_class_a['children']:
        #     if item['name'] == 'Orphan':
        #         orphan_data = OrderedDict(
        #             [('name', ''), ('value', 3000), ('color', ''), ('children', [item])])
        #         whole_class_a['children'].remove(item)
        #         break
        # context['class_a'] = json.dumps(whole_class_a)
        # class_b1_data = tree.get_tree_data(
        #     ProteinFamily.objects.get(name__startswith='Class B1 (Secretin)'))
        # context['class_b1_options'] = deepcopy(tree.d3_options)
        # context['class_b1_options']['anchor'] = 'class_b1'
        # context['class_b1_options']['branch_trunc'] = 60
        # context['class_b1_options']['label_free'] = [1, ]
        # context['class_b1'] = json.dumps(
        #     class_b1_data.get_nodes_dict(self.page))
        # class_b2_data = tree.get_tree_data(
        #     ProteinFamily.objects.get(name__startswith='Class B2 (Adhesion)'))
        # context['class_b2_options'] = deepcopy(tree.d3_options)
        # context['class_b2_options']['anchor'] = 'class_b2'
        # context['class_b2_options']['label_free'] = [1, ]
        # context['class_b2'] = json.dumps(
        #     class_b2_data.get_nodes_dict(self.page))
        # class_c_data = tree.get_tree_data(
        #     ProteinFamily.objects.get(name__startswith='Class C (Glutamate)'))
        # context['class_c_options'] = deepcopy(tree.d3_options)
        # context['class_c_options']['anchor'] = 'class_c'
        # context['class_c_options']['branch_trunc'] = 50
        # context['class_c_options']['label_free'] = [1, ]
        # context['class_c'] = json.dumps(class_c_data.get_nodes_dict(self.page))
        # class_f_data = tree.get_tree_data(
        #     ProteinFamily.objects.get(name__startswith='Class F (Frizzled)'))
        # context['class_f_options'] = deepcopy(tree.d3_options)
        # context['class_f_options']['anchor'] = 'class_f'
        # context['class_f_options']['label_free'] = [1, ]
        # context['class_f'] = json.dumps(class_f_data.get_nodes_dict(self.page))
        # class_t2_data = tree.get_tree_data(
        #     ProteinFamily.objects.get(name__startswith='Class T (Taste 2)'))
        # context['class_t2_options'] = deepcopy(tree.d3_options)
        # context['class_t2_options']['anchor'] = 'class_t2'
        # context['class_t2_options']['label_free'] = [1, ]
        # context['class_t2'] = json.dumps(
        #     class_t2_data.get_nodes_dict(self.page))
        # # definition of the class a orphan tree
        # context['orphan_options'] = deepcopy(tree.d3_options)
        # context['orphan_options']['anchor'] = 'orphan'
        # context['orphan_options']['label_free'] = [1, ]
        # context['orphan'] = json.dumps(orphan_data)
        #
        # whole_receptors = Protein.objects.prefetch_related(
        #     "family", "family__parent__parent__parent")
        # whole_rec_dict = {}
        # for rec in whole_receptors:
        #     rec_uniprot = rec.entry_short()
        #     rec_iuphar = rec.family.name.replace("receptor", '').replace(
        #         "<i>", "").replace("</i>", "").strip()
        #     if (rec_iuphar[0].isupper()) or (rec_iuphar[0].isdigit()):
        #         whole_rec_dict[rec_uniprot] = [rec_iuphar]
        #     else:
        #         whole_rec_dict[rec_uniprot] = [rec_iuphar.capitalize()]
        #
        # context["whole_receptors"] = json.dumps(whole_rec_dict)
        #
        #
        # circle_data = BiasedData.objects.filter(physiology_biased__isnull=False).values_list(
        #               "physiology_biased", "receptor_id__entry_name", "ligand_id").order_by(
        #               "physiology_biased", "receptor_id__entry_name", "ligand_id").distinct(
        #               "physiology_biased", "receptor_id__entry_name", "ligand_id")
        #
        # circles = {}
        # label_converter = {'Arrestin-2': "β-Arr",
        #                    'Arrestin-3': "β-Arr 2",
        #                    'Gaq/i-chimera': "Gq/i-chim",
        #                    'Minigi': "Mini-Gi"}
        # endpoint = 0
        # for data in circle_data:
        #     # if data[1].split('_')[1] == 'human':
        #     key = data[1].split('_')[0].upper()
        #     val = data[0].split(' (')[0].capitalize()
        #     if val in label_converter.keys():
        #         val = label_converter[val]
        #     if key not in circles.keys():
        #         circles[key] = {}
        #     if val not in circles[key].keys():
        #         circles[key][val] = 1
        #     circles[key][val] += 1
        #     if circles[key][val] > endpoint:
        #         endpoint = circles[key][val]
        #
        # context["circles_data"] = json.dumps(circles)
        # context["endpoint"] = endpoint

        return context



class NewDrugsBrowser(TemplateView):
    # Template using this class #
    template_name = 'NewDrugsBrowser.html'
    # Get context for hmtl usage #
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Get data - server side - Queries #
        Drug_browser_data = Drugs2024.objects.all().prefetch_related('target','target__family__parent__parent','target__family__parent__parent__parent','indication','indication__code','ligand','ligand__ligand_type','moa')
        # GtoPdb =
        # Possible addons: 'target__family__parent','target__family__parent__parent','target__family__parent__parent__parent','indication','ligand','moa'
        # initialize context list for pushing data to html #
        context_data_browser = list()
        i_breaker = 0
        #proteins = list(TissueExpression.objects.all().values_list('protein__entry_name').distinct())
        #drugs = Drugs.objects.all().prefetch_related('target', 'target__family__parent__parent__parent', 'publication', 'publication__web_link', 'publication__web_link__web_resource', 'publication__journal')
        #drugs_NHS_names = list(NHSPrescribings.objects.values_list('drugname__name', flat=True).distinct())
        for entry in Drug_browser_data:
            # i_breaker += 1
            # if i_breaker >= 10:
            #     break
            # else:
            #     pass
            ## For the drug browser table ##
            # Protein id, uniprot, and receptor name
            Protein_id = str(entry.target.id)
            Protein_uniprot = str(entry.target.accession)
            Protein_name = str(entry.target.entry_name).split("_")[0].upper()
            Protein_receptor_name = str(entry.target.name)
            Protein_family = str(entry.target.family.parent)
            Protein_class = str(entry.target.family.parent.parent.parent)
            Drug_name = str(entry.ligand.name)
            Indication = str(entry.indication.name)
            Indication_ID = str(entry.indication.code.index)
            Clinical_drug_status = str(entry.drug_status)
            Clinical_max_phase = int(entry.indication_max_phase)
            Clinical_approval_year = str(entry.approval_year)
            Clinical_drug_type = str(entry.ligand.ligand_type.name)
            Clinical_moa = str(entry.moa.name)

            # create and append to context data
            jsondata_browser = {
                'Index_number': Protein_id,
                'Drug': Drug_name,
                'Indication': Indication,
                'Indication_ID': Indication_ID,
                'Protein_uniprot': Protein_uniprot,
                'Protein_name': Protein_name,
                'Protein_receptor': Protein_receptor_name,
                'Protein_class': Protein_class,
                'Protein_family': Protein_family,
                'Drug_status': Clinical_drug_status,
                'Indication_max_phase': Clinical_max_phase,
                'Approval_year': Clinical_approval_year,
                'Drug_type': Clinical_drug_type,
                'Moa': Clinical_moa
            }
            context_data_browser.append(jsondata_browser)
        context['drug_data'] = context_data_browser
        return context

class TargetSelectionTool(TemplateView):
    # Template using this class #
    template_name = 'TargetSelectionTool.html'
    # Get context for hmtl usage #
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # Get data - server side - Queries #
        TissueExp = TissueExpression.objects.all().prefetch_related('protein','tissue')
        Target_drug_data = Drugs2024.objects.all().prefetch_related('target','ligand','indication','indication__code')
        target_ids = [drug.target.id for drug in Target_drug_data if drug.target]
        #Structures_data = Structure.objects.all().prefetch_related('state', 'protein_conformation__protein')
        Structure_data = Structure.objects.all().exclude(structure_type__slug__startswith='af-').values('protein_conformation__protein__parent__id').annotate(
                    Inactive=Count(Case(When(state_id=1, then=1), output_field=IntegerField())),
                    Active=Count(Case(When(state_id=2, then=1), output_field=IntegerField())),
                    Intermediate=Count(Case(When(state_id=3, then=1), output_field=IntegerField())),
                    Other=Count(Case(When(state_id=4, then=1), output_field=IntegerField())),
                ).order_by('protein_conformation__protein__parent__id') # Fetches only human structures --> make sure we dont want other species


        #'target__family__parent__parent','target__family__parent__parent__parent','indication','indication__code','ligand','moa'
        # Context lists for pushing data #
        context_target_selection = list()
        context_data_tissue = list()
        # Dicts to modulate the data from the server side #
        structure_dict = {}
        target_selection_dict = {}
        target_indication_dict = {}
        Tissue_expression_dict = {}
        index_dict = {}
        # Create dict for structure total and state (active, inactive and intermediate)
        for entry in Structure_data:
            protein_id = str(entry['protein_conformation__protein__parent__id'])
            if protein_id in target_ids:
                if protein_id not in structure_dict:
                    structure_dict[protein_id] = {}
                    structure_dict[protein_id]['Active'] = entry['Active']
                    structure_dict[protein_id]['Inactive'] = entry['Inactive']
                    structure_dict[protein_id]['Intermediate'] = entry['Intermediate']
                    structure_dict[protein_id]['Total'] = int(entry['Active'])+entry['Inactive']+entry['Intermediate']
                else:
                    print("Something should be terribly wrong if there is two identical ids")
                    break
            else:
                structure_dict[protein_id] = {}
                structure_dict[protein_id]['Active'] = 0
                structure_dict[protein_id]['Inactive'] = 0
                structure_dict[protein_id]['Intermediate'] = 0
                structure_dict[protein_id]['Total'] = 0
        # Create Target selection browser #
        for entry in Target_drug_data:
            # Ids and keys
            Protein_id = str(entry.target.id)
            Indication_id = str(entry.indication.code.index)
            Target_indication_pair = "{}___{}".format(Protein_id,Indication_id)
            # Static values #
            Protein_uniprot = str(entry.target.accession)
            Protein_name = str(entry.target.entry_name).split("_")[0].upper()
            Indication = str(entry.indication.name)
            Drug_name  = str(entry.ligand.name)
            Drug_status = str(entry.drug_status)
            Novelty_score = float(entry.novelty_score)
            # Dicts #
            if Target_indication_pair not in target_indication_dict:
                target_indication_dict[Target_indication_pair] = {}
                target_indication_dict[Target_indication_pair]['information'] = [Protein_id,Indication_id,Indication]
                target_indication_dict[Target_indication_pair]['Novelty_score'] = Novelty_score
                #target_indication_dict[Target_indication_pair]['Drugs_total'] = 1
                target_indication_dict[Target_indication_pair]['Drugs'] = []
                target_indication_dict[Target_indication_pair]['Drugs'].append(Drug_name)
                # Handle drugs
                target_indication_dict[Target_indication_pair]['Drug__status'] = {}
                target_indication_dict[Target_indication_pair]['Drug__status']['Active'] = 0
                target_indication_dict[Target_indication_pair]['Drug__status']['Approved'] = 0
                target_indication_dict[Target_indication_pair]['Drug__status']['Discontinued'] = 0
                target_indication_dict[Target_indication_pair]['Drug__status'][Drug_status] += 1
                # Handle structures #
                if Protein_id in structure_dict:
                    target_indication_dict[Target_indication_pair]['Structures'] = {}
                    target_indication_dict[Target_indication_pair]['Structures']['Total'] = structure_dict[Protein_id]['Total']
                    target_indication_dict[Target_indication_pair]['Structures']['Active'] = structure_dict[Protein_id]['Active']
                    target_indication_dict[Target_indication_pair]['Structures']['Inactive'] = structure_dict[Protein_id]['Inactive']
                    target_indication_dict[Target_indication_pair]['Structures']['Intermediate'] = structure_dict[Protein_id]['Intermediate']
                else:
                    print("Something fishy")
                    target_indication_dict[Target_indication_pair]['Structures'] = {}
                    target_indication_dict[Target_indication_pair]['Structures']['Total'] = 0
                    target_indication_dict[Target_indication_pair]['Structures']['Active'] = 0
                    target_indication_dict[Target_indication_pair]['Structures']['Inactive'] = 0
                    target_indication_dict[Target_indication_pair]['Structures']['Intermediate'] = 0
            else:
                # Drugs #
                target_indication_dict[Target_indication_pair]['Drugs'].append(Drug_name)
                target_indication_dict[Target_indication_pair]['Drug__status'][Drug_status] += 1

            if Protein_id not in target_selection_dict:
                target_selection_dict[Protein_id] = [Protein_uniprot,Protein_name]
        for key in target_indication_dict:
            key_id = str(target_indication_dict[key]['information'][0])
            jsondata_TargetSelectionTool = {
                    'Index_number': key_id,
                    'Target_name': target_selection_dict[key_id][1],
                    'Target_uniprot': target_selection_dict[key_id][0],
                    'Indication_name': target_indication_dict[key]['information'][2],
                    'Indication_id': target_indication_dict[key]['information'][1],
                    'Novelty_score': target_indication_dict[key]['Novelty_score'],
                    'IDG': "Coming soon",
                    'Drugs_approved_names': target_indication_dict[key]['Drugs'],
                    'Drugs_total': int(len(target_indication_dict[key]['Drugs'])),
                    'Drugs_approved': int(target_indication_dict[key]['Drug__status']['Approved']),
                    'Drugs_in_trial': int(target_indication_dict[key]['Drug__status']['Active']),
                    'Drugs_discontinued': int(target_indication_dict[key]['Drug__status']['Discontinued']),
                    'Structures_total': int(target_indication_dict[key]['Structures']['Total']),
                    'Structures_active' : int(target_indication_dict[key]['Structures']['Active']),
                    'Structures_inactive' : int(target_indication_dict[key]['Structures']['Inactive']),
                    'Structures_intermediate' : int(target_indication_dict[key]['Structures']['Intermediate'])
            }
            context_target_selection.append(jsondata_TargetSelectionTool)
        context['Target_data'] = context_target_selection
        # Go through server side data and modulate into a dict #
        for entry in TissueExp:
            # string values for Tissue expression table #
            protein_id = entry.protein.entry_name
            value = entry.value
            Tissue_id = entry.tissue.name
            # Index key #
            index_key = entry.protein.id
            if protein_id not in index_dict:
                index_dict[str(protein_id)] = index_key
            # Expression value linked to protein / target #
            if protein_id not in Tissue_expression_dict:
                Tissue_expression_dict[str(protein_id)] = {}
                Tissue_expression_dict[str(protein_id)][str(Tissue_id)] = float(value)
            else:
                Tissue_expression_dict[str(protein_id)][str(Tissue_id)] = float(value)
        # Run through dict and assign the correct values into the context data #
        for key in Tissue_expression_dict:
            jsondata_tissue = {
                    'Index_number': index_dict[key],
                    'ProteinID': key,
                    'adipose_tissue': Tissue_expression_dict[key]['adipose_tissue'],
                    'adrenal_gland': Tissue_expression_dict[key]['adrenal_gland'],
                    'amygdala': Tissue_expression_dict[key]['amygdala'],
                    'appendix': Tissue_expression_dict[key]['appendix'],
                    'basal_ganglia': Tissue_expression_dict[key]['basal_ganglia'],
                    'bone_marrow': Tissue_expression_dict[key]['bone_marrow'],
                    'breast': Tissue_expression_dict[key]['breast'],
                    'cerebellum': Tissue_expression_dict[key]['cerebellum'],
                    'cerebral_cortex': Tissue_expression_dict[key]['cerebral_cortex'],
                    'cervix': Tissue_expression_dict[key]['cervix'],
                    'choroid_plexus': Tissue_expression_dict[key]['choroid_plexus'],
                    'colon': Tissue_expression_dict[key]['colon'],
                    'duodenum': Tissue_expression_dict[key]['duodenum'],
                    'endometrium': Tissue_expression_dict[key]['endometrium_1'], #Should be updated
                    'epididymis': Tissue_expression_dict[key]['epididymis'],
                    'esophagus': Tissue_expression_dict[key]['esophagus'],
                    'fallopian_tube': Tissue_expression_dict[key]['fallopian_tube'],
                    'gallbladder': Tissue_expression_dict[key]['gallbladder'],
                    'heart_muscle': Tissue_expression_dict[key]['heart_muscle'],
                    'hippocampal_formation': Tissue_expression_dict[key]['hippocampal_formation'],
                    'hypothalamus': Tissue_expression_dict[key]['hypothalamus'],
                    'kidney': Tissue_expression_dict[key]['kidney'],
                    'liver': Tissue_expression_dict[key]['liver'],
                    'lung': Tissue_expression_dict[key]['lung'],
                    'lymph_node': Tissue_expression_dict[key]['lymph_node'],
                    'midbrain': Tissue_expression_dict[key]['midbrain'],
                    'ovary': Tissue_expression_dict[key]['ovary'],
                    'pancreas': Tissue_expression_dict[key]['pancreas'],
                    'parathyroid_gland': Tissue_expression_dict[key]['parathyroid_gland'],
                    'pituitary_gland': Tissue_expression_dict[key]['pituitary_gland'],
                    'placenta': Tissue_expression_dict[key]['placenta'],
                    'prostate': Tissue_expression_dict[key]['prostate'],
                    'rectum': Tissue_expression_dict[key]['rectum'],
                    'retina': Tissue_expression_dict[key]['retina'],
                    'salivary_gland': Tissue_expression_dict[key]['salivary_gland'],
                    'seminal_vesicle': Tissue_expression_dict[key]['seminal_vesicle'],
                    'skeletal_muscle': Tissue_expression_dict[key]['skeletal_muscle'],
                    'skin': Tissue_expression_dict[key]['skin_1'],
                    'small_intestine': Tissue_expression_dict[key]['small_intestine'],
                    'smooth_muscle': Tissue_expression_dict[key]['smooth_muscle'],
                    'spinal_cord': Tissue_expression_dict[key]['spinal_cord'],
                    'spleen': Tissue_expression_dict[key]['spleen'],
                    'stomach': Tissue_expression_dict[key]['stomach_1'],
                    'testis': Tissue_expression_dict[key]['testis'],
                    'thymus': Tissue_expression_dict[key]['thymus'],
                    'thyroid_gland': Tissue_expression_dict[key]['thyroid_gland'],
                    'tongue': Tissue_expression_dict[key]['tongue'],
                    'tonsil': Tissue_expression_dict[key]['tonsil'],
                    'urinary_bladder': Tissue_expression_dict[key]['urinary_bladder'],
                    'vagina': Tissue_expression_dict[key]['vagina']
                }
            # Append context data into list #
            context_data_tissue.append(jsondata_tissue)
        # Create context data for tissue expression data #
        context['Tissue_data'] = context_data_tissue
        # Lastly return context for html usage #
        return context


@cache_page(60 * 60 * 24 * 28)
def drugmapping(request):
    context = dict()

    families = ProteinFamily.objects.all()
    lookup = {}
    for f in families:
        lookup[f.slug] = f.name.replace("receptors","").replace(" receptor","").replace(" hormone","").replace("/neuropeptide","/").replace(" (G protein-coupled)","").replace(" factor","").replace(" (LPA)","").replace(" (S1P)","").replace("GPR18, GPR55 and GPR119","GPR18/55/119").replace("-releasing","").replace(" peptide","").replace(" and oxytocin","/Oxytocin").replace("Adhesion class orphans","Adhesion orphans").replace("muscarinic","musc.").replace("-concentrating","-conc.")

    class_proteins = Protein.objects.filter(family__slug__startswith="00",source__name='SWISSPROT', species_id=1).prefetch_related('family').order_by('family__slug')

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
    indication_data = Drugs2024.objects.filter(indication__code__index=code).prefetch_related('ligand',
                                                                                              'target',
                                                                                              'indication',
                                                                                              'indication__code')

    indication_name = Indication.objects.filter(code__index=code).values_list('name', flat=True).distinct()[0]

    sankey = {"nodes": [],
              "links": []}
    caches = {'indication':[],
              'ligands': [],
              'targets': [],
              'entries': []}

    node_counter = 0
    for record in indication_data:
        #assess the values for indication/ligand/protein
        indication_code = record.indication.name.capitalize()
        ligand_name = record.ligand.name.capitalize()
        ligand_id = record.ligand.id
        protein_name = record.target.name
        target_name = record.target.entry_name
        #check for each value if it exists and retrieve the source node value
        if indication_code not in caches['indication']:
            sankey['nodes'].append({"node": node_counter, "name": indication_code, "url":'https://www.ebi.ac.uk/ols4/ontologies/efo/classes?short_form='+code})
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
