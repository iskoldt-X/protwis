from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings
from django.db.models import Count, Max, Q, F, Value, CharField, Case, When, IntegerField
from django.db.models import Count, Max
from django.core.cache import cache
from django.db import connection, reset_queries
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView
from django.utils.safestring import mark_safe
from common.views import AbsReferenceSelectionTable, getReferenceTable, getLigandTable, getLigandCountTable, AbsTargetSelection
from protein.models import Protein, ProteinFamily, TissueExpression
from structure.models import Structure
from drugs.models import Drugs, Indication, ATCCodes
from protein.models import Protein, ProteinFamily, Tissues, TissueExpression
from mapper.views import LandingPage

import re
import json
import numpy as np
from collections import OrderedDict
from copy import deepcopy
import pandas as pd
import os

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
    # Here we need to generate data for three different Venn diagrams plus associated tables:
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
            drug_phases = Drugs.objects.all().values_list('indication_max_phase','ligand_id__name').distinct()
            for item in drug_phases:
                if item[0] not in phases_dict.keys():
                    phases_dict[item[0]] = []
                phases_dict[item[0]].append(item[1])
            phases_dict = {key_mapping[k]: phases_dict[k] for k in key_mapping if k in phases_dict}
            for key in phases_dict.keys():
                phases_dict[key] = '\n'.join(phases_dict[key])
        else:
            # Call to get receptors in each maximum phase
            receptor_phases = Drugs.objects.all().values_list('indication_max_phase','target_id__entry_name').distinct()
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
        table_data = Drugs.objects.select_related(
            'target__family__parent__parent__parent',  # All target info
            'ligand__ligand_type',
            'indication',
            'moa',
            'disease_association'
        ).values(
            'target__entry_name',  # Gene name
            'target__name',  # Protein name
            'target__family__parent__name',  # Receptor family
            'target__family__parent__parent__name',  # Ligand type
            'target__family__parent__parent__parent__name',  # Class
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
            'target__entry_name': "Gene name",
            'target__name': "Protein name",
            'target__family__parent__name': 'Receptor family',
            'target__family__parent__parent__name': 'Ligand type',
            'target__family__parent__parent__parent__name': 'Class',
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

        # Convert DataFrame to JSON
        json_records = df.to_json(orient='records')

        # Pass the JSON data to the template context
        context['Full_data'] = json_records

        # Collect drugs information
        # drugs_panel = Drugs.objects.all().select_related(
        #                                                     "ligand",
        #                                                     "moa",
        #                                                     "target__family__parent",
        #                                                     "target",
        #                                                     "indication__code"
        #                                                     )

        # drug_dictionary = {}
        # for p in drugs_panel:
        #     # Collect receptor data
        #     lig_name = p.ligand.name.lower().capitalize()
        #     lig_type = p.moa.name
        #     rec_family = p.target.family.parent.short()
        #     rec_uniprot = p.target.entry_short()
        #     rec_iuphar = p.target.family.name.replace("receptor", '').replace("<i>","").replace("</i>","").strip()
        #     clinical_phase = p.indication_max_phase
        #     indication_name = p.indication.name
        #     indication_code = p.indication.code.index
        #     # Create a tuple to store the values
        #     drug_entry = (lig_name, lig_type, rec_family, rec_uniprot, rec_iuphar, clinical_phase, indication_name, indication_code)
        #     # Initialize key if it doesn't exist
        #     if origin == 'drugs':
        #         key = lig_name
        #     else:
        #         key = rec_uniprot.lower().capitalize()
        #     if key not in drug_dictionary:
        #         drug_dictionary[key] = []
        #     # Check for duplicates before adding
        #     if drug_entry not in drug_dictionary[key]:
        #         drug_dictionary[key].append(drug_entry)

        # context["drug_dictionary"] = json.dumps(drug_dictionary)
    # cache.set(name_of_cache, context, 60 * 60 * 24 * 7)  # seven days timeout on cache
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
            search_data = Drugs.objects.all().prefetch_related('ligand').distinct('ligand')
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
            atc_df_grouped.rename(columns={'ligand': 'Ligand ID', 'code__index': 'ATC'}, inplace=True)

            # Fetch table data with all related information
            table_data = Drugs.objects.select_related(
                'ligand',
                'target__family__parent__parent__parent', # All target info
                'indication'
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
                'ligand': 'Ligand ID',
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
            df = df.merge(atc_df_grouped, on='Ligand ID', how='left')
            # Fill NaN values in the 'ATC' column with None
            df['ATC'] = df['ATC'].fillna("")

            # Precompute Phase and Status Information
            df['Is_Phase_I'] = (df['Phase'] == 1).astype(int)
            df['Is_Phase_II'] = (df['Phase'] == 2).astype(int)
            df['Is_Phase_III'] = (df['Phase'] == 3).astype(int)
            df['Is_Approved'] = (df['Status'] == 'Approved').astype(int)

            # Perform GroupBy and Aggregate
            Modified_df = df.groupby(
                ['Ligand ID', 'Gene name', 'Indication name', 'Ligand name', 'Protein name', 'Receptor family', 'Ligand type', 'Class', 'ICD11', 'ATC', 'Association score']
            ).agg(
                Highest_phase=('Phase', 'max'),  # Get the highest phase for each group
                Phase_I_trials=('Is_Phase_I', 'sum'),  # Count Phase I trials
                Phase_II_trials=('Is_Phase_II', 'sum'),  # Count Phase II trials
                Phase_III_trials=('Is_Phase_III', 'sum'),  # Count Phase III trials
                Approved=('Is_Approved', 'max')  # Check if any row has 'Approved' status
            ).reset_index()

            # Convert 'Approved' from integer to 'Yes'/'No'
            Modified_df['Approved'] = Modified_df['Approved'].apply(lambda x: 'Yes' if x == 1 else 'No')

            # Convert DataFrame to JSON
            json_records = Modified_df.to_json(orient='records')

            # Pass the JSON data to the template context
            context['Full_data'] = json_records

        elif page == 'Targets':
            description = 'Search by target name'
            # Fetch distinct targets and create a dictionary of {target.name: target.id}
            search_data = Drugs.objects.all().prefetch_related('target').distinct('target')
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
            atc_df_grouped.rename(columns={'ligand': 'Ligand ID', 'code__index': 'ATC'}, inplace=True)


            # Fetch table data with all related information
            table_data = Drugs.objects.select_related(
                'target',
                'ligand__ligand_type',
                'indication',
                'moa',
                'disease_association'
            ).values(
                'target',  # Target ID
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
                'target__name': 'Target name',
                'ligand__name': 'Ligand name',
                'ligand__ligand_type__name': 'Modality',
                'moa__name': 'Mode of action',
                'indication__title': 'Indication name',
                'indication__code': 'ICD11',
                'indication_max_phase': 'Phase',
                'ligand': 'Ligand ID',
                'disease_association__association_score' : 'Association score',
                'drug_status': 'Status'
            }, inplace=True)

            # Merge the ATC data into the main DataFrame (df) on 'Ligand ID'
            df = df.merge(atc_df_grouped, on='Ligand ID', how='left')
            # Fill NaN values in the 'ATC' column with None
            df['ATC'] = df['ATC'].fillna("Unavailable")

            # Precompute flags for phase and approval status
            df['Is_Phase_I'] = (df['Phase'] == 1).astype(int)
            df['Is_Phase_II'] = (df['Phase'] == 2).astype(int)
            df['Is_Phase_III'] = (df['Phase'] == 3).astype(int)
            df['Is_Approved'] = (df['Status'] == 'Approved').astype(int)

            # Group by necessary columns and perform aggregation in one go
            grouped = df.groupby(
                ['Target ID', 'Target name', 'Ligand name', 'Indication name', 'Modality', 'Mode of action', 'ICD11', 'ATC', 'Association score']
            )

            # Perform aggregation
            Modified_df = grouped.agg(
                Highest_phase=('Phase', 'max'),  # Get the highest phase for each group
                Phase_I_trials=('Is_Phase_I', 'sum'),  # Sum up the Phase I trials
                Phase_II_trials=('Is_Phase_II', 'sum'),  # Sum up the Phase II trials
                Phase_III_trials=('Is_Phase_III', 'sum'),  # Sum up the Phase III trials
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
            search_data = Drugs.objects.all().prefetch_related('indication').distinct('indication')
            search_dict = {drug.indication.title: drug.indication.id for drug in search_data}

            # ###########################
            # Single Data Query
            # ###########################
            # Fetch all data in a single query
            table_data = Drugs.objects.select_related(
                'indication',
                'ligand__ligand_type',
                'target__family__parent__parent__parent',  # All target info
                'moa',
                'disease_association'
            ).values(
                'indication',  # Indication ID
                'indication__title',  # Indication name
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
                'ligand__id': 'Ligand ID',
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
            group_cols = ['Indication ID', 'Gene name', 'Indication name', 'Protein name', 'Receptor family', 'Ligand type', 'Class']

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
            group_cols_drugs = ['Indication ID', 'Gene name', 'Drug name', 'Indication name', 'Protein name', 'Receptor family', 'Ligand type', 'Class', 'Molecule_type','Mode of action']

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

            #### GPCRome Indication Stuff START ####
            data_dir = os.sep.join([settings.DATA_DIR, 'drug_data'])
            filepath = os.sep.join([data_dir, 'short_titles_ICD.csv'])
            titles = pd.read_csv(filepath, sep=';', low_memory=False)
            title_conversion = {key: [value] for key, value in zip(titles['title'], titles['title_short'])}

            indication_levels_01 = Indication.objects.filter(level__in=[0,1])
            indication_tree = {}
            conversion = {}
            wheel_data = {}
            wheel_slugs = {}
            crunch = {}

            for item in indication_levels_01:
                if item.title == 'Symptoms, signs or clinical findings, not elsewhere classified':
                    item.title = 'Symptoms, signs or clinical findings'
                elif item.title == 'Certain conditions originating in the perinatal period':
                    item.title = 'Certain conditions originating in perinatal period'
                elif item.title == 'Injury, poisoning or certain other consequences of external causes':
                    item.title = 'Injury, poisoning or other external causes'
                elif item.title == 'Pregnancy, childbirth or the puerperium':
                    item.title = 'Pregnancy, childbirth or puerperium'
                elif item.title == 'Diseases of the blood or blood-forming organs':
                    item.title = 'Diseases of the blood or related organs'

                if (item.level == 0) and (item.title.split(' ')[0] not in ['Supplementary', 'Extension', 'External', 'Factors']):
                    indication_tree[item.slug] = []
                    conversion[item.slug] = item.title
                if (item.level == 1) and (item.parent.title.split(' ')[0] not in ['Supplementary', 'Extension', 'External', 'Factors']):
                    root = item.slug[:4]
                    if root not in indication_tree.keys():
                        indication_tree[root] = []
                    indication_tree[root].append(item.title)
                    conversion[item.slug] = item.title
                    wheel_data[item.title] = {'Value1': 0}
                    wheel_slugs[item.slug] = {'Value1': 0}
                    crunch[item.title] = {1: 0, 2: 0, 3: 0, 4: 0, 'unique': []}

            indication_tree2 = LandingPage.convert_keys(indication_tree, conversion)

            #Now get the drug data
            indication_drug_data = Drugs.objects.all().prefetch_related('indication')

            for item in indication_drug_data:
                try:
                    title = item.indication.get_level_1().title
                    slug = item.indication.get_level_1().slug
                    phase = item.indication_max_phase
                    wheel_data[title]['Value1'] +=1
                    wheel_slugs[slug]['Value1'] +=1
                    crunch[title][phase] += 1
                    crunch[title]['unique'].append(item.ligand_id)
                except:
                    continue

            for key in title_conversion.keys():
                title_conversion[key].append(wheel_data[key]['Value1'])
                title_conversion[key].append(crunch[key][1])
                title_conversion[key].append(crunch[key][2])
                title_conversion[key].append(crunch[key][3])
                title_conversion[key].append(crunch[key][4])
                uniq = len(set(crunch[key]['unique']))
                title_conversion[key].append(uniq)

            indication_full = {"NameList": indication_tree2, "DataPoints": wheel_data}
            context['GPCRome_data'] = json.dumps(indication_full["NameList"])
            context['GPCRome_data_variables'] = json.dumps(indication_full['DataPoints'])
            context['Title_conversion'] = json.dumps(title_conversion)
            #### GPCRome Indication Stuff END   ####

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
        drug_data = Drugs.objects.all().values_list('target__entry_name', 'indication_max_phase').distinct()

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
        drug_data = Drugs.objects.all().values_list('target__entry_name', 'indication_max_phase')
        drug_dict = {}
        for drug in drug_data:
            if drug[0] not in drug_dict.keys():
                drug_dict[drug[0]] = {'Outer1': 0, 'Outer2': 0, 'Outer3': 0, 'Outer4': 0, 'Inner': 0}
            if drug[1] in [1, 2, 3, 4]:
                outer_key = f"Outer{drug[1]}"
                drug_dict[drug[0]][outer_key] += 1
            if drug[1] > drug_dict[drug[0]]['Inner']:
                drug_dict[drug[0]]['Inner'] = drug[1]

        tree, tree_options, circles, receptors = LandingPage.generate_tree_plot(drug_dict)
        #Remove 0 circles
        for key, outer_dict in circles.items():
            circles[key] = {k: v for k, v in outer_dict.items() if v != 0}

        context['tree'] = json.dumps(tree)
        context['tree_options'] = tree_options
        context['circles'] = json.dumps(circles)
        context['whole_dict'] = json.dumps(receptors)

        #REPURPOSED TREE SECTION
        # Red: IndicationMaxPhase < 4 and MaxPhase (compound) < 4 [New agents in trial (lack approval)]
        # Purple: IndicationMaxPhase < 4 and MaxPhase = 4         [Drugs being repurposed in trials (have approval)]
        # Blue: IndicationMaxPhase and MaxPhase = 4               [All drugs (both those being repurposed in trials and not)]

        drug_data = Drugs.objects.all().values_list('ligand__name', 'target__entry_name', 'indication_max_phase').distinct()
        phase4 = {}
        for drug in drug_data:
            if drug[0] not in phase4.keys():
                phase4[drug[0]] = []
            if drug[2] == 4:
                phase4[drug[0]].append(drug[1])
        phase4 = {k: v for k, v in phase4.items() if v != []}
        # Red: IndicationMaxPhase < 4 and MaxPhase (compound) < 4 [New agents in trial (lack approval)]
        # Purple: IndicationMaxPhase < 4 and MaxPhase = 4         [Drugs being repurposed in trials (have approval)]
        # Blue: IndicationMaxPhase and MaxPhase = 4               [All drugs (both those being repurposed in trials and not)]
        drug_dict = {}
        for drug in drug_data:
            if drug[1] not in drug_dict.keys():
                drug_dict[drug[1]] = {'Outer1': 0, 'Outer2': 0, 'Outer3': 0, 'Outer4': 0, 'Inner': 0}
            #Approved Drug
            if drug[2] == 4:
                drug_dict[drug[1]]['Outer3'] += 1
            else:
                #Repurposed Drug
                if drug[0] in phase4.keys():
                    drug_dict[drug[1]]['Outer2'] += 1
                # New Agent
                else:
                    drug_dict[drug[1]]['Outer1'] += 1

        repurposed_tree, repurposed_tree_options, repurposed_circles, repurposed_receptors = LandingPage.generate_tree_plot(drug_dict)
        #Remove 0 circles
        for key, outer_dict in repurposed_circles.items():
            repurposed_circles[key] = {k: v for k, v in outer_dict.items() if v != 0}

        context['rep_tree'] = json.dumps(repurposed_tree)
        context['rep_tree_options'] = repurposed_tree_options
        context['rep_circles'] = json.dumps(repurposed_circles)
        context['rep_whole_dict'] = json.dumps(repurposed_receptors)

        return context

class DiseaseOverview(TemplateView):
    """
    Per class statistics of known ligands.
    """

    template_name = 'disease_overview.html'

    def get_context_data(self, **kwargs):

        context = super().get_context_data(**kwargs)

        def find_max_values(node, node_name):
            """
            Given a node (a dictionary) and its name (a string), return:
            (max_red, red_node, max_purple, purple_node, max_blue, blue_node)
            for the entire subtree starting at `node`.

            If a color key is missing at a node, treat its value as 0.
            """
            # Get current node's color values (default to 0 if missing)
            local_red = node.get('red', 0)
            local_purple = node.get('purple', 0)
            local_blue = node.get('blue', 0)

            # Initialize maxima with the current node's values
            max_red = local_red
            red_node = node_name
            max_purple = local_purple
            purple_node = node_name
            max_blue = local_blue
            blue_node = node_name

            # Iterate over children. Children are identified as dictionary values.
            for key, value in node.items():
                if isinstance(value, dict):
                    # Recurse into child
                    r_val, r_k, p_val, p_k, b_val, b_k = find_max_values(value, key)

                    # Compare and update max_red
                    if r_val > max_red:
                        max_red = r_val
                        red_node = r_k

                    # Compare and update max_purple
                    if p_val > max_purple:
                        max_purple = p_val
                        purple_node = p_k

                    # Compare and update max_blue
                    if b_val > max_blue:
                        max_blue = b_val
                        blue_node = b_k

            return max_red, red_node, max_purple, purple_node, max_blue, blue_node


        def find_node_by_name(nested_dict, target_name):
            """
            Recursively search nested_dict for a node keyed by target_name.
            Returns the dictionary for that node if found, else None.
            """
            for key, value in nested_dict.items():
                if key == target_name:
                    # Found the node
                    return value
                if isinstance(value, dict):
                    # Recurse into children
                    result = find_node_by_name(value, target_name)
                    if result is not None:
                        return result
            return None

        def increment_color_value(nested_dict, target_leaf, color):
            """
            Search through nested_dict for target_leaf. Once found, increment the
            value associated with 'color' by 1.
            Returns True if successful, False if not found.
            """
            for key, value in nested_dict.items():
                if key == target_leaf:
                    # Found the leaf node
                    if color in value:
                        # Handle case if value is None or not initialized
                        if value[color] is None:
                            value[color] = 0
                        value[color] += 1
                    else:
                        # If the color key doesn't exist, initialize it at 1
                        value[color] = 1
                    return True

                # If the value is a nested dict, recurse into it
                if isinstance(value, dict):
                    if increment_color_value(value, target_leaf, color):
                        return True

            # If not found in this branch, return False
            return False

        def gradient_color(value_fraction, start_color=(255,255,255), end_color=(0,0,0)):

            r = int(start_color[0] + (end_color[0] - start_color[0]) * value_fraction)
            g = int(start_color[1] + (end_color[1] - start_color[1]) * value_fraction)
            b = int(start_color[2] + (end_color[2] - start_color[2]) * value_fraction)

            return f"#{r:02x}{g:02x}{b:02x}"

        def find_global_maxima(node, global_max_values):
            """
            Recursively traverse 'node' and update global_max_values.
            Any node can have 'red', 'purple', 'blue' keys.
            Missing keys default to 0.
            """
            red_val = node.get('red', 0)
            purple_val = node.get('purple', 0)
            blue_val = node.get('blue', 0)

            if red_val > global_max_values['red']:
                global_max_values['red'] = red_val
            if purple_val > global_max_values['purple']:
                global_max_values['purple'] = purple_val
            if blue_val > global_max_values['blue']:
                global_max_values['blue'] = blue_val

            # Recurse into child dictionaries
            for k, v in node.items():
                if isinstance(v, dict):
                    find_global_maxima(v, global_max_values)


        def apply_colors(node, global_max_values, color_targets):
            """
            Recursively apply gradient colors to each node based on its values and global maxima.
            Replace 'red', 'purple', 'blue' values with their corresponding gradient colors.
            """
            # Extract node's own values (default 0 if missing)
            red_val = node.get('red', 0)
            purple_val = node.get('purple', 0)
            blue_val = node.get('blue', 0)

            red_frac = red_val / global_max_values['red'] if global_max_values['red'] > 0 else 0
            purple_frac = purple_val / global_max_values['purple'] if global_max_values['purple'] > 0 else 0
            blue_frac = blue_val / global_max_values['blue'] if global_max_values['blue'] > 0 else 0

            # Compute gradient colors
            # The function gradient_color(frac, start_color, end_color) must be defined beforehand.
            new_node = {}
            new_node['red'] = gradient_color(red_frac, (255,255,255), color_targets['red'])
            new_node['purple'] = gradient_color(purple_frac, (255,255,255), color_targets['purple'])
            new_node['blue'] = gradient_color(blue_frac, (255,255,255), color_targets['blue'])

            # Recurse into children
            for k, v in node.items():
                if isinstance(v, dict):
                    new_node[k] = apply_colors(v, global_max_values, color_targets)
                else:
                    # For non-dict, non-color keys, if you need to preserve them, reassign here
                    # Check that you're not overwriting the colors keys we just set
                    if k not in ['red', 'purple', 'blue']:
                        new_node[k] = v

            return new_node

        def convert_dict_to_colors(data_dict):
            # Define the target colors for each metric
            color_targets = {
                "red": (255,0,0),
                "purple": (128,0,128),
                "blue": (0,0,255)
            }

            # 1. Find global maxima for each color across the entire nested structure
            global_max_values = {"red": 0, "purple": 0, "blue": 0}
            find_global_maxima(data_dict, global_max_values)

            # 2. Apply gradient colors throughout the nested dict
            result = apply_colors(data_dict, global_max_values, color_targets)
            return result

        def color_style(color):
            # If the color is white, add a border style. Otherwise, no border.
            if color == '#ffffff':
                return f"background-color: {color}; border: 1px solid #ccc;"
            else:
                return f"background-color: {color};"

        def propagate_max_colors(node):
            """
            Given a nested dictionary node:
            - Keys that are colors (red, purple, blue) have integer values.
            - Other keys may lead to child dictionaries.
            This function updates each node's (red, purple, blue) to be the maximum
            found in that node or any of its descendants.

            Returns a tuple (max_red, max_purple, max_blue) for the subtree.
            """

            # Initialize current node's colors. If missing, default to 0.
            current_red = node.get('red', 0)
            current_purple = node.get('purple', 0)
            current_blue = node.get('blue', 0)

            # Check children
            for key, value in node.items():
                if isinstance(value, dict):
                    # Recursively propagate max colors down the subtree
                    r_val, p_val, b_val = propagate_max_colors(value)

                    # Update current node's max based on child
                    if r_val > current_red:
                        current_red = r_val
                    if p_val > current_purple:
                        current_purple = p_val
                    if b_val > current_blue:
                        current_blue = b_val

            # Update the node's own values after checking children
            node['red'] = current_red
            node['purple'] = current_purple
            node['blue'] = current_blue

            return current_red, current_purple, current_blue

        def render_nested_structure(data):
            html = []
            html.append('<ul class="no-bullet-indent">')
            for key, value in data.items():
                if key in ('red', 'purple', 'blue'):
                    continue

                # Safely extract colors
                red = value.get('red', '#ffffff') if isinstance(value, dict) else '#ffffff'
                purple = value.get('purple', '#ffffff') if isinstance(value, dict) else '#ffffff'
                blue = value.get('blue', '#ffffff') if isinstance(value, dict) else '#ffffff'

                # Check if node has children
                has_children = False
                if isinstance(value, dict):
                    for v in value.values():
                        if isinstance(v, dict):
                            has_children = True
                            break

                if isinstance(value, dict) and has_children:
                    # Node with children
                    html.append(f"""
                    <li>
                        <details>
                            <summary>
                                <span class="dot" style="{color_style(red)}"></span>
                                <span class="dot" style="{color_style(purple)}"></span>
                                <span class="dot" style="{color_style(blue)}"></span>
                                <b>{key}</b>
                            </summary>
                            {render_nested_structure(value)}
                        </details>
                    </li>
                    """)
                else:
                    # Leaf node
                    if isinstance(value, dict):
                        # Leaf node with colors
                        html.append(f"""
                        <li>
                            <span class="dot" style="{color_style(red)}"></span>
                            <span class="dot" style="{color_style(purple)}"></span>
                            <span class="dot" style="{color_style(blue)}"></span>
                            {key}
                        </li>
                        """)
                    else:
                        # Leaf node without colors (just a string or value)
                        # No colors defined, but if you want, you can handle differently
                        html.append(f"<li>{key}: {value}</li>")
            html.append('</ul>')
            return ''.join(html)

        lvl_0_1 = Indication.objects.filter(level__in=[0,1]).order_by('slug')
        all_indis = Indication.objects.all().order_by('slug')
        indication_data = {}
        listdata = {}
        listplot_data_variables = {}
        Label_Conversion = {}
        data_types_list = {'Col1': "Continuouos", 'Col2': "Continuouos", 'Col3': "Continuouos", 'Col4': "Continuouos",}
        #Gather data and setup master dict

        for indi in lvl_0_1:
           if indi.level == 0 and indi.title not in indication_data.keys():
               indication_data[indi.title] = {}
               listdata[indi.title] = []
           if indi.level == 1:
               lvl0 = indi.get_level_0().title
               if lvl0 not in indication_data.keys():
                   indication_data[lvl0] = {}
                   listdata[lvl0] = []
               if indi.title not in indication_data[lvl0].keys():
                   indication_data[lvl0][indi.title] = {'red': 0, 'purple': 0, 'blue': 0}
                   listdata[lvl0].append(indi.title)
                   listplot_data_variables[indi.title] = {'Value1': 'Circle', 'Value2': 0, 'Value3': 'Circle', 'Value4': 0, 'Value5': 'Circle','Value6': 0}
                   Label_Conversion[indi.title] = indi.title
        #with master dict set, now we add data from drugs
        drug_data = Drugs.objects.all().prefetch_related('ligand', 'indication', 'indication__parent',
                                                         'indication__parent__parent', 'indication__parent__parent__parent',
                                                         'indication__parent__parent__parent__parent',
                                                         'indication__parent__parent__parent__parent__parent')

        drug_dict = {}
        for item in drug_data:
            if item.ligand.name not in drug_dict.keys():
                drug_dict[item.ligand.name] = {'Max Phase': 0}
            if item.indication.get_level_1().title not in drug_dict[item.ligand.name].keys():
                drug_dict[item.ligand.name][item.indication.get_level_1().title] = [item.indication_max_phase, item.indication.get_level_0().title]
            if item.indication_max_phase > drug_dict[item.ligand.name]['Max Phase']:
                drug_dict[item.ligand.name]['Max Phase'] = item.indication_max_phase

        # Red: IndicationMaxPhase < 4 and MaxPhase (compound) < 4 [New agents in trial (lack approval)]
        # Purple: IndicationMaxPhase < 4 and MaxPhase = 4         [Drugs being repurposed in trials (have approval)]
        # Blue: IndicationMaxPhase and MaxPhase = 4               [All drugs (both those being repurposed in trials and not)]

        for drug in drug_dict.keys():
            max_phase = drug_dict[drug]['Max Phase']
            indications = list(drug_dict[drug].keys())[1:]
            for ind in indications:
                if max_phase < 4 and drug_dict[drug][ind][0] < 4:
                    listplot_data_variables[ind]['Value2'] += 1
                elif max_phase == 4 and drug_dict[drug][ind][0] < 4:
                    listplot_data_variables[ind]['Value4'] += 1
                elif max_phase == 4 and drug_dict[drug][ind][0] == 4:
                    listplot_data_variables[ind]['Value6'] += 1

        counter = 1
        labeled_listdata = {}
        for key in listdata.keys():
            new_key = f"{counter:02d} {key}"
            labeled_listdata[new_key] = listdata[key]
            counter += 1

        #Manually purging irrelevant classes
        to_purge = ['23 External causes of morbidity or mortality',
                    # 'Injury, poisoning or certain other consequences of external causes',
                    '26 Supplementary Chapter Traditional Medicine Conditions',
                    '27 Supplementary section for functioning assessment',
                    '24 Factors influencing health status or contact with health services',
                    '28 Extension Codes']

        for item in to_purge:
            del labeled_listdata[item]

        List_Data_result = {"category_array": [], "final_array": []}
        for key in labeled_listdata.keys():
            List_Data_result['category_array'].append('ReceptorFamily')
            List_Data_result['final_array'].append(key)
            for value in labeled_listdata[key]:
                List_Data_result['category_array'].append('Receptor')
                List_Data_result['final_array'].append(value)

        ####### TESTING ########

        records = []
        for indi in all_indis:
            records.append({'title':indi.title, 'level':indi.level})
        # Setup the master nested dict
        hierarchy = {}
        # Start with a "virtual root" at level -1
        stack = [(-1, hierarchy)]
        for record in records:
            level = record['level']
            title = record['title']
            # Pop until the top of the stack is strictly less than the current record's level
            while stack and stack[-1][0] >= level:
                stack.pop()
            # Now the stack top should be the parent level
            parent_dict = stack[-1][1]
            # Create a new entry for the current record
            parent_dict[title] = {}
            # Push current node
            stack.append((level, parent_dict[title]))

        counter = 1
        labeled_hierarchy = {}
        for key in hierarchy.keys():
            new_key = f"{counter:02d} {key}"
            labeled_hierarchy[new_key] = hierarchy[key]
            counter += 1

        for item in to_purge:
            del labeled_hierarchy[item]

        drug_dict_leaf = {}
        for item in drug_data:
            if item.ligand.name not in drug_dict_leaf.keys():
                drug_dict_leaf[item.ligand.name] = {'Max Phase': 0}
            if item.indication.title not in drug_dict_leaf[item.ligand.name].keys():
                drug_dict_leaf[item.ligand.name][item.indication.title] = [item.indication_max_phase, item.indication.title]
            if item.indication_max_phase > drug_dict_leaf[item.ligand.name]['Max Phase']:
                drug_dict_leaf[item.ligand.name]['Max Phase'] = item.indication_max_phase

        #now need to populated them
        for drug in drug_dict_leaf.keys():
            max_phase = drug_dict_leaf[drug]['Max Phase']
            indications = list(drug_dict_leaf[drug].keys())[1:]
            for ind in indications:
                if max_phase < 4 and drug_dict_leaf[drug][ind][0] < 4:
                    increment_color_value(labeled_hierarchy, ind, 'red')
                elif max_phase == 4 and drug_dict_leaf[drug][ind][0] < 4:
                    increment_color_value(labeled_hierarchy, ind, 'purple')
                elif max_phase == 4 and drug_dict_leaf[drug][ind][0] == 4:
                    increment_color_value(labeled_hierarchy, ind, 'blue')

        max_labels = {}
        for primary_key in labeled_hierarchy:
            r_val, r_k, p_val, p_k, b_val, b_k = find_max_values(labeled_hierarchy[primary_key], primary_key)
            max_labels[primary_key] = {
                'red': r_k if r_val != 0 else '#ffffff',
                'purple': p_k if p_val != 0 else '#ffffff',
                'blue': b_k if b_val != 0 else '#ffffff'
            }

        propagate_max_colors(labeled_hierarchy)
        colored = convert_dict_to_colors(labeled_hierarchy)

        # Now we update max_labels with actual colors from colored_result.
        for major_key, color_dict in max_labels.items():
            for color_key, node_name in color_dict.items():
                if node_name == '#ffffff':
                    continue

                # Find the node in colored_result that matches node_name
                node = find_node_by_name(colored, node_name)
                if node is not None and color_key in node:
                    # Replace the string in max_labels with the actual color tuple
                    max_labels[major_key][color_key] = node[color_key]

        del colored['red']
        del colored['purple']
        del colored['blue']
        # Combine data and labels into a single structure:
        combined = []
        for main_key, keys in colored.items():
            main_label = max_labels.get(main_key, {"red": "#ffffff", "purple": "#ffffff", "blue": "#ffffff"})
            rendered_html = mark_safe(render_nested_structure(keys))
            combined.append({
                "main_key": main_key,
                "main_label": main_label,
                "rendered_items": rendered_html
            })

        context['data'] = json.dumps(indication_data)
        context['listplot_data'] = json.dumps(labeled_listdata)
        context['listplot_data_variables'] = json.dumps(listplot_data_variables)
        context['listplot_datatypes'] = json.dumps(data_types_list)
        context['Label_Conversion'] = json.dumps(Label_Conversion)
        context['List_Data_result'] = json.dumps(List_Data_result)
        context['combined'] = combined
        context['data'] = colored
        context['labels'] = max_labels

        return context

class NewDrugsBrowser(TemplateView):
    # Template using this class #
    template_name = 'Drugs_Indications_Targets.html'
    # Get context for hmtl usage #
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Get data - server side - Queries #
        Drug_browser_data = Drugs.objects.all().prefetch_related('target','target__family__parent__parent','target__family__parent__parent__parent','indication','indication__code','ligand','ligand__ligand_type','moa')
        # Possible addons: 'target__family__parent','target__family__parent__parent','target__family__parent__parent__parent','indication','ligand','moa'
        # initialize context list for pushing data to html #
        context_data_browser = list()
        for entry in Drug_browser_data:
            ## For the drug browser table ##
            # Protein id, uniprot, and receptor name
            Protein_id = str(entry.target.id)
            Protein_uniprot = str(entry.target.accession)
            Protein_name = str(entry.target.entry_name).split("_")[0].upper()
            Protein_receptor_name = str(entry.target.name)
            Protein_family = str(entry.target.family.parent)
            Protein_class = str(entry.target.family.parent.parent.parent)
            Drug_id = str(entry.ligand.id)
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
                'DrugID': Drug_id,
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


#################################
####          Drugs          ####
#################################

class DrugsTable(TemplateView):
    # Template using this class #
    template_name = 'Drugs.html'
    # Get context for hmtl usage #
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Get data - server side - Queries #
        Drug_browser_data = Drugs.objects.all().prefetch_related('target','ligand')

        # initialize context list for pushing data to html #
        context_data_browser = list()

        # Drugs and drugs target (#number)
        Data_drugs_dict = {}
        for entry in Drug_browser_data:
            Drug_id = str(entry.ligand.id)
            Protein_id = str(entry.target.id)
            Drug_name = str(entry.ligand.name)
            if Drug_id not in Data_drugs_dict:
                Data_drugs_dict[Drug_id] = {}
                Data_drugs_dict[Drug_id]['Name'] = Drug_name
                Data_drugs_dict[Drug_id]['Targets'] = []
                Data_drugs_dict[Drug_id]['Targets'].append(Protein_id)
            else:
                if Protein_id not in Data_drugs_dict[Drug_id]['Targets']:
                    Data_drugs_dict[Drug_id]['Targets'].append(Protein_id)
                else:
                    pass

        for drug_id in Data_drugs_dict:
            ## For the drug browser table ##
            Drug_name = Data_drugs_dict[drug_id]['Name']
            target_number = len(Data_drugs_dict[drug_id]['Targets'])

            # create and append to context data
            jsondata_browser = {
                'Drug_id': drug_id,
                'Drugs': Drug_name,
                'Targets': target_number
            }
            context_data_browser.append(jsondata_browser)
        context['Data_Drugs'] = context_data_browser
        return context

#################################
####      Indications        ####
#################################

class Indications(TemplateView):
    # Template using this class #
    template_name = 'Indications.html'
    # Get context for hmtl usage #
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Get data - server side - Queries #
        Drug_browser_data = Drugs.objects.all().prefetch_related('target','indication','indication__code','ligand')
        # initialize context list for pushing data to html #
        context_data_browser = list()

        Indications_dict = {}
        for entry in Drug_browser_data:
            ## For the drug browser table ##
            Protein_id = str(entry.target.id)
            Drug_id = str(entry.ligand.id)
            Indication_id = str(entry.indication.id)
            Indication_name = str(entry.indication.name)
            Indication_code = str(entry.indication.code.index)

            if Indication_id not in Indications_dict:
                Indications_dict[Indication_id] = {}
                Indications_dict[Indication_id]['Name'] = Indication_name
                Indications_dict[Indication_id]['Code'] = Indication_code
                Indications_dict[Indication_id]['Drugs'] = []
                Indications_dict[Indication_id]['Targets'] = []
                # append drugs and targets to indications dict #
                Indications_dict[Indication_id]['Drugs'].append(Drug_id)
                Indications_dict[Indication_id]['Targets'].append(Protein_id)
            else:
                if Drug_id not in Indications_dict[Indication_id]['Drugs']:
                    Indications_dict[Indication_id]['Drugs'].append(Drug_id)
                else:
                    pass
                if Protein_id not in Indications_dict[Indication_id]['Targets']:
                    Indications_dict[Indication_id]['Targets'].append(Protein_id)
                else:
                    pass
        for indication in Indications_dict:

            Indication_name = str(Indications_dict[indication]['Name'])
            Indication_code = str(Indications_dict[indication]['Code'])
            Number_of_drugs = len(Indications_dict[indication]['Drugs'])
            Number_of_targets = len(Indications_dict[indication]['Targets'])

            # create and append to context data
            jsondata_browser = {
                'Indication_id': indication,
                'Indication_name': Indication_name,
                'Indication_code': Indication_code,
                'Drugs_number': Number_of_drugs,
                'Target_number': Number_of_targets
            }
            context_data_browser.append(jsondata_browser)
        context['Data_Indications'] = context_data_browser
        return context


#################################
####        Targets          ####
#################################

class Targets(TemplateView):
    # Template using this class #
    template_name = 'Targets.html'
    # Get context for hmtl usage #
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Get data - server side - Queries #
        Drug_browser_data = Drugs.objects.all().prefetch_related('target','target__family__parent__parent','target__family__parent__parent__parent')

        # initialize context list for pushing data to html #
        context_data_browser = list()

        for entry in Drug_browser_data:
            ## For the drug browser table ##
            # Protein id, uniprot, and receptor name
            Protein_id = str(entry.target.id)
            Protein_uniprot = str(entry.target.accession)
            Protein_name = str(entry.target.entry_name).split("_")[0].upper()
            Protein_receptor_name = str(entry.target.name)
            Protein_family = str(entry.target.family.parent)
            Protein_class = str(entry.target.family.parent.parent.parent)

            # create and append to context data
            jsondata_browser = {
                'Index_number': Protein_id,
                'Protein_uniprot': Protein_uniprot,
                'Protein_name': Protein_name,
                'Protein_receptor': Protein_receptor_name,
                'Protein_class': Protein_class,
                'Protein_family': Protein_family,
            }
            context_data_browser.append(jsondata_browser)
        context['Data_Targets'] = context_data_browser
        return context

class TargetSelectionTool(TemplateView):
    # Template using this class #
    template_name = 'TargetSelectionTool.html'
    # Get context for hmtl usage #
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # Get data - server side - Queries #
        TissueExp = TissueExpression.objects.all().prefetch_related('protein','tissue')
        Target_drug_data = Drugs.objects.all().prefetch_related('target','ligand','indication','indication__code')
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
        Drugs_target_ids = []
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
                # if Protein_id in structure_dict:
                #     target_indication_dict[Target_indication_pair]['Structures'] = {}
                #     target_indication_dict[Target_indication_pair]['Structures']['Total'] = structure_dict[Protein_id]['Total']
                #     target_indication_dict[Target_indication_pair]['Structures']['Active'] = structure_dict[Protein_id]['Active']
                #     target_indication_dict[Target_indication_pair]['Structures']['Inactive'] = structure_dict[Protein_id]['Inactive']
                #     target_indication_dict[Target_indication_pair]['Structures']['Intermediate'] = structure_dict[Protein_id]['Intermediate']
                # else:
                #     print("Something fishy")
                #     target_indication_dict[Target_indication_pair]['Structures'] = {}
                #     target_indication_dict[Target_indication_pair]['Structures']['Total'] = 0
                #     target_indication_dict[Target_indication_pair]['Structures']['Active'] = 0
                #     target_indication_dict[Target_indication_pair]['Structures']['Inactive'] = 0
                #     target_indication_dict[Target_indication_pair]['Structures']['Intermediate'] = 0
            else:
                # Drugs #
                target_indication_dict[Target_indication_pair]['Drugs'].append(Drug_name)
                target_indication_dict[Target_indication_pair]['Drug__status'][Drug_status] += 1

            if Protein_id not in target_selection_dict:
                target_selection_dict[Protein_id] = [Protein_uniprot,Protein_name]
        for key in target_indication_dict:
            key_id = str(target_indication_dict[key]['information'][0])
            Drugs_target_ids.append(key_id)
            jsondata_TargetSelectionTool = {
                    'Index_number': key_id,
                    'Target_name': target_selection_dict[key_id][1],
                    'Target_uniprot': target_selection_dict[key_id][0],
                    'Indication_name': target_indication_dict[key]['information'][2],
                    'Indication_id': target_indication_dict[key]['information'][1],
                    'Novelty_score': str(target_indication_dict[key]['Novelty_score']),
                    'IDG': "Coming soon",
                    'Drugs_approved_names': target_indication_dict[key]['Drugs'],
                    'Drugs_total': int(len(target_indication_dict[key]['Drugs'])),
                    'Drugs_approved': int(target_indication_dict[key]['Drug__status']['Approved']),
                    'Drugs_in_trial': int(target_indication_dict[key]['Drug__status']['Active']),
                    'Drugs_discontinued': int(target_indication_dict[key]['Drug__status']['Discontinued'])
                    # 'Structures_total': int(target_indication_dict[key]['Structures']['Total']),
                    # 'Structures_active' : int(target_indication_dict[key]['Structures']['Active']),
                    # 'Structures_inactive' : int(target_indication_dict[key]['Structures']['Inactive']),
                    # 'Structures_intermediate' : int(target_indication_dict[key]['Structures']['Intermediate'])
            }
            context_target_selection.append(jsondata_TargetSelectionTool)
        context['Drug_Targets'] = context_target_selection
        # Go through server side data and modulate into a dict #
        for entry in TissueExp:
            # string values for Tissue expression table #
            protein_id = entry.protein.entry_name
            value = entry.value
            ### Get rid of the nan values ### <--- should be removed through the build
            if value != value:
                value = None
            Tissue_id = entry.tissue.name
            # Index key #
            index_key = str(entry.protein.id)
            if index_key not in Drugs_target_ids:
                pass
            else:
                if protein_id not in index_dict:
                    index_dict[str(protein_id)] = index_key
                # Expression value linked to protein / target #
                if protein_id not in Tissue_expression_dict:
                    Tissue_expression_dict[str(protein_id)] = {}
                    Tissue_expression_dict[str(protein_id)][str(Tissue_id)] = float(value) if value else '-'
                else:
                    Tissue_expression_dict[str(protein_id)][str(Tissue_id)] = float(value) if value else '-'
        # Run through dict and assign the correct values into the context data #
        for key in Tissue_expression_dict:
            jsondata_tissue = {
                    'Index_number': index_dict[key],
                    'ProteinID': key,
                    'ProteinName': str(key).split("_")[0].upper(),
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
        context['Tissue_header_list'] = ['Adipose tissue',
                                       'Adrenal gland',
                                       'Amygdala',
                                       'Appendix',
                                       'Basal ganglia',
                                       'Bone marrow',
                                       'Breast',
                                       'Cerebellum',
                                       'Cerebral cortex',
                                       'Cervix',
                                       'Choroid plexus',
                                       'Colon',
                                       'Duodenum',
                                       'Endometrium',
                                       'Epididymis',
                                       'Esophagus',
                                       'Fallopian tube',
                                       'Gallbladder',
                                       'Heart muscle',
                                       'Hippocampal formation',
                                       'Hypothalamus',
                                       'Kidney',
                                       'Liver',
                                       'Lung',
                                       'Lymph node',
                                       'Midbrain',
                                       'Ovary',
                                       'Pancreas',
                                       'Parathyroid gland',
                                       'Pituitary gland',
                                       'Placenta',
                                       'Prostate',
                                       'Rectum',
                                       'Retina',
                                       'Salivary gland',
                                       'Seminal vesicle',
                                       'Skeletal muscle',
                                       'Skin',
                                       'Small intestine',
                                       'Smooth muscle',
                                       'Spinal cord',
                                       'Spleen',
                                       'Stomach',
                                       'Testis',
                                       'Thymus',
                                       'Thyroid gland',
                                       'Tongue',
                                       'Tonsil',
                                       'Urinary bladder',
                                       'Vagina']
        context['Tissue_header_list'] = [x.replace(' ', '<br>') for x in context['Tissue_header_list']]
        # Lastly return context for html usage #
        return context

# @cache_page(60 * 60 * 24 * 28)
def indication_detail(request, code):

    code = code.upper()
    context = dict()
    #code = '4A8Z'
    indication_data = Drugs.objects.filter(indication__code=code).prefetch_related('ligand',
                                                                                   'target',
                                                                                   'indication')

    indication_name = Indication.objects.filter(code=code).values_list('title', flat=True).distinct()[0]

    sankey = {"nodes": [],
              "links": []}
    caches = {'indication':[],
              'level_0': [],
              'ligands': [],
              'targets': [],
              'entries': []}

    node_counter = 0
    for record in indication_data:
        #assess the values for indication/ligand/protein
        indication_code = record.indication.title
        indication_0 = record.indication.get_level_0().title
        ligand_name = record.ligand.name.capitalize()
        uri = record.indication.uri.index
        ligand_id = record.ligand.id
        protein_name = record.target.name
        target_name = record.target.entry_name
        #check for each value if it exists and retrieve the source node value
        if indication_code not in caches['indication']:
            sankey['nodes'].append({"node": node_counter, "name": indication_code, "url":'https://icd.who.int/browse/2024-01/mms/en#'+uri})
            node_counter += 1
            caches['indication'].append(indication_code)
        indi_node = next((item['node'] for item in sankey['nodes'] if item['name'] == indication_code), None)

        if indication_0 not in caches['level_0']:
            sankey['nodes'].append({"node": node_counter, "name": indication_0, "url":'https://icd.who.int/browse/2024-01/mms/en#'+uri})
            node_counter += 1
            caches['level_0'].append(indication_0)
        level_0_node = next((item['node'] for item in sankey['nodes'] if item['name'] == indication_0), None)

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

        #append connection between indication and level 0
        sankey['links'].append({"source":indi_node, "target":level_0_node, "value":1, "ligtrace": ligand_name, "prottrace": None})
        #append connection between level 0 and ligand
        sankey['links'].append({"source":level_0_node, "target":lig_node, "value":1, "ligtrace": ligand_name, "prottrace": None})
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
    total_points = len(caches['targets']) + len(caches['targets']) + 1
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
