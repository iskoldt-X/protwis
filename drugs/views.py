from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect, JsonResponse
from django.conf import settings
from django.db.models import Count, Max, Q, F, Value, CharField, Case, When, IntegerField
from django.db.models import Count, Max
from django.core.cache import cache
from django.db import connection, reset_queries
from django.views.decorators.cache import cache_page
from django.views.generic import TemplateView
from django.utils.safestring import mark_safe
from common.views import AbsReferenceSelectionTable, getReferenceTable, getLigandTable, getLigandCountTable, AbsTargetSelection
from structure.models import Structure
from drugs.models import Drugs, Indication, ATCCodes, IndicationAssociation
from protein.views import get_sankey_data
from protein.models import Protein, ProteinFamily, Tissues, TissueExpression
from mapper.views import LandingPage
from ligand.models import AssayExperiment, LigandID

import re
import json
import numpy as np
from collections import OrderedDict, defaultdict
from copy import deepcopy
import pandas as pd
import os
import logging

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
            for key in list(phases_dict.keys()):
                # Calculate the count of items in the current phase
                count = len(phases_dict[key])
                
                # Update the key with the count in parentheses
                new_key = f"{key} {count}"
                
                # Update the dictionary: add a new key and remove the old one
                phases_dict[new_key] = '\n'.join(phases_dict[key])
                del phases_dict[key]  # Remove the old key
        else:
            # Call to get receptors in each maximum phase
            receptor_phases = Drugs.objects.all().values_list('indication_max_phase','target_id__entry_name').distinct()
            for item in receptor_phases:
                if item[0] not in phases_dict.keys():
                    phases_dict[item[0]] = []
                phases_dict[item[0]].append(item[1])
            phases_dict = {key_mapping[k]: phases_dict[k] for k in key_mapping if k in phases_dict}
            for key in list(phases_dict.keys()):
                # Calculate the count of items in the current phase
                count = len(phases_dict[key])
                
                # Update the key with the count in parentheses
                new_key = f"{key} ({count})"
                
                # Update the dictionary: add a new key and remove the old one
                phases_dict[new_key] = '\n'.join(phases_dict[key])
                del phases_dict[key]  # Remove the old key

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

            # Fetch distinct ligands and their names
            table_data = Drugs.objects.select_related('ligand').values(
                'ligand',  # Agent/Drug id
                'ligand__name',  # Agent/Drug name
            ).distinct()

            # Convert `table_data` to a DataFrame
            table_df = pd.DataFrame(table_data)

            # Ensure we have a unique list of ligand IDs
            unique_ligand_ids = table_df['ligand'].unique().tolist()

            # Fetch related resources: `web_resource__name` and `index` for the unique ligand IDs
            resource_data = LigandID.objects.filter(
                ligand__in=unique_ligand_ids
            ).values(
                'ligand',  # Ligand ID
                'web_resource__name',  # Web resource name (e.g., "PubChem", "ChEMBL")
                'index',  # Identifier ID
            )

            # Convert `resource_data` queryset to a list of dictionaries
            resource_data_list = list(resource_data)

            # Transform `resource_data_list` into a structure for merging
            resource_dict = defaultdict(dict)
            for row in resource_data_list:
                ligand_id = row['ligand']
                resource_name = row['web_resource__name']
                resource_index = row['index']
                resource_dict[ligand_id][resource_name] = resource_index

            # Add resource data to the table DataFrame by mapping from `resource_dict`
            table_df['Resources'] = table_df['ligand'].map(resource_dict)

            # Flatten the `Resources` dictionary into separate columns for each resource
            resource_df = table_df['Resources'].apply(pd.Series)

            # Merge `table_df` and `resource_df`
            final_df = pd.concat([table_df, resource_df], axis=1).drop(columns=['Resources'])

            # Rename columns for better readability
            final_df.rename(columns={
                'ligand': 'LigandID',
                'ligand__name': 'Ligand Name',
                'ChEMBL_compound_ids': 'ChEMBL',  # Rename for consistency
                'Guide To Pharmacology': 'GTP',  # Rename for consistency
                'Drug Central': 'DrugCentral',   # Rename for consistency
            }, inplace=True)

            # Convert the final DataFrame to JSON
            json_records = final_df.to_json(orient='records')

            # Pass the JSON data to the template context
            context['Drugs_identifier_table'] = json_records

        elif page == 'Targets':
            description = 'Search by target name'
            # Fetch distinct targets and create a dictionary of {target.name: target.id}
            search_data = Drugs.objects.all().prefetch_related('target').distinct('target')
            search_dict = {drug.target.name: drug.target.id for drug in search_data}

            # # Create sankey_dict_serialized using a dictionary comprehension
            # sankey_dict_serialized = {
            #     drug.target.id: get_sankey_data(self.get_entry_name_by_target_id(drug.target.id))
            #     for drug in search_data
            # }

            # Pass the serialized sankey_dict to the context
            # context['sankey_dict'] = json.dumps(sankey_dict_serialized)

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
            table_data = Drugs.objects.select_related(
                'target',
                'ligand__ligand_type',
                'indication',
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
            search_data = Drugs.objects.all().prefetch_related('indication').distinct('indication')
            search_dict = {drug.indication.title: drug.indication.id for drug in search_data}

            # ###########################
            # Single Data Query
            # ###########################
            # Fetch all data in a single query
            table_data = Drugs.objects.select_related(
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

        #     #### GPCRome Indication Stuff START ####
        #     data_dir = os.sep.join([settings.DATA_DIR, 'drug_data'])
        #     filepath = os.sep.join([data_dir, 'short_titles_ICD.csv'])
        #     titles = pd.read_csv(filepath, sep=';', low_memory=False)
        #     title_conversion = {key: [value] for key, value in zip(titles['title'], titles['title_short'])}

        #     indication_levels_01 = Indication.objects.filter(level__in=[0,1])
        #     indication_tree = {}
        #     conversion = {}
        #     wheel_data = {}
        #     wheel_slugs = {}
        #     crunch = {}

        #     for item in indication_levels_01:
        #         if item.title == 'Symptoms, signs or clinical findings, not elsewhere classified':
        #             item.title = 'Symptoms, signs or clinical findings'
        #         elif item.title == 'Certain conditions originating in the perinatal period':
        #             item.title = 'Certain conditions originating in perinatal period'
        #         elif item.title == 'Injury, poisoning or certain other consequences of external causes':
        #             item.title = 'Injury, poisoning or other external causes'
        #         elif item.title == 'Pregnancy, childbirth or the puerperium':
        #             item.title = 'Pregnancy, childbirth or puerperium'
        #         elif item.title == 'Diseases of the blood or blood-forming organs':
        #             item.title = 'Diseases of the blood or related organs'

                # if (item.level == 0) and (item.title.split(' ')[0] not in ['Supplementary', 'Extension', 'External', 'Factors']):
                #     indication_tree[item.slug] = []
                #     conversion[item.slug] = item.title
                # if (item.level == 1) and (item.parent.title.split(' ')[0] not in ['Supplementary', 'Extension', 'External', 'Factors']):
                #     root = item.slug[:4]
                #     if root not in indication_tree.keys():
                #         indication_tree[root] = []
                #     indication_tree[root].append(item.title)
                #     conversion[item.slug] = item.title
                #     wheel_data[item.title] = {'Value1': 0}
                #     wheel_slugs[item.slug] = {'Value1': 0}
                #     crunch[item.title] = {1: 0, 2: 0, 3: 0, 4: 0, 'unique': []}

        #     indication_tree2 = LandingPage.convert_keys(indication_tree, conversion)

        #     #Now get the drug data
        #     indication_drug_data = Drugs.objects.all().prefetch_related('indication')

        #     for item in indication_drug_data:
        #         try:
        #             title = item.indication.get_level_1().title
        #             slug = item.indication.get_level_1().slug
        #             phase = item.indication_max_phase
        #             wheel_data[title]['Value1'] +=1
        #             wheel_slugs[slug]['Value1'] +=1
        #             crunch[title][phase] += 1
        #             crunch[title]['unique'].append(item.ligand_id)
        #         except:
        #             continue

        #     for key in title_conversion.keys():
        #         title_conversion[key].append(wheel_data[key]['Value1'])
        #         title_conversion[key].append(crunch[key][1])
        #         title_conversion[key].append(crunch[key][2])
        #         title_conversion[key].append(crunch[key][3])
        #         title_conversion[key].append(crunch[key][4])
        #         uniq = len(set(crunch[key]['unique']))
        #         title_conversion[key].append(uniq)

        #     indication_full = {"NameList": indication_tree2, "DataPoints": wheel_data}
        #     context['GPCRome_data'] = json.dumps(indication_full["NameList"])
        #     context['GPCRome_data_variables'] = json.dumps(indication_full['DataPoints'])
        #     context['Title_conversion'] = json.dumps(title_conversion)
        #     #### GPCRome Indication Stuff END   ####

        # Convert to JSON string and pass to context
        if page == 'Drugs':
            pass
        else:
            context['search_data'] = json.dumps(search_dict)
        context['page'] = page
        context['title'] = title
        context['description'] = description
        # context['table_data'] = table_data

        return context

    def get_entry_name_by_target_id(self, target_id):
        """
        Retrieves the protein.entry_name based on the given drug.target.id.

        Args:
            target_id (int): The ID of the target (Protein).

        Returns:
            str or None: The entry_name of the associated Protein, or None if not found.
        """
        try:
            protein = Protein.objects.get(id=target_id)
            return protein.entry_name
        except Protein.DoesNotExist:
            return None

# Set up logging
logger = logging.getLogger(__name__)

def fetch_sankey_data_view(request):
    logger.info("fetch_sankey_data_view called")  # Log the function call
    entry_name = request.GET.get('entry_name', None)  # Extract `entry_name` from the request
    
    if not entry_name:
        logger.error("No entry_name provided in the request")
        return JsonResponse({'error': 'Missing entry_name'}, status=400)

    logger.info(f"Received entry_name: {entry_name}")  # Log the entry name

    try:
        sankey_data = get_sankey_data(entry_name)  # Call your helper function
        if sankey_data:
            logger.info(f"Sankey data successfully retrieved for entry_name: {entry_name}")
            return JsonResponse({'sankey_data': sankey_data}, status=200)
        else:
            logger.warning(f"No data found for entry_name: {entry_name}")
            return JsonResponse({'error': 'No data found for this entry'}, status=404)
    except Exception as e:
        logger.exception(f"An error occurred while fetching sankey data for entry_name: {entry_name}")
        return JsonResponse({'error': 'An internal server error occurred'}, status=500)

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
        drug_data = Drugs.objects.all().values_list('target__entry_name', 'ligand__name','indication_max_phase')
        
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
                # if target == 'drd2_human' and key == 'Outer4':
                #     print(len(unique_entries),sorted(unique_entries))
                drug_dict[target][key] = len(unique_entries)  # Replace the list with the count

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

class TargetSelectionTool(TemplateView):
    # Template using this class #
    template_name = 'TargetSelectionTool.html'
    # Get context for hmtl usage #
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)

        # Fetch all proteins
        all_proteins = Protein.objects.filter(
            species_id=1,
            parent_id__isnull=True,
            accession__isnull=False,
            family_id__slug__startswith='0'
        ).exclude(
            family_id__slug__startswith='007'
        ).exclude(
            family_id__slug__startswith='008'
        ).values(
            'id', 'entry_name', 'name', 'family__parent__name', 
            'family__parent__parent__name', 'family__parent__parent__parent__name'
        )

        # Convert to DataFrame
        proteins_df = pd.DataFrame(list(all_proteins))
        proteins_df.rename(columns={
            'id': 'Target ID',
            'entry_name': 'Gene name',
            'name': 'Protein name',
            'family__parent__name': 'Receptor family',
            'family__parent__parent__name': 'Ligand type',
            'family__parent__parent__parent__name': 'Class'
        }, inplace=True)
        
        # Fetch all data in a single query
        table_data = Drugs.objects.select_related(
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
                'indication_max_phase', # Max pahse
                'publication_count', # publication count
                'target_level' # IDG target level
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
            'ligand': "Ligand ID", # Ligand ID
            'novelty_score': 'Novelty (Pharos)', #Only novelty we have
            'drug_status': 'Status', #Approval
            'moa__name': 'Mode of action', #Modality
            'indication_max_phase': 'Phase', # phase
            'publication_count': 'Literature', # Pub count
            'target_level': 'IDG'

        }, inplace=True)

        # Merge the proteins DataFrame with the drugs DataFrame on 'Target ID'
        df = proteins_df.merge(df, on='Target ID', how='left', suffixes=('', '_dup'))

        # Drop duplicate columns introduced by the merge
        df = df.loc[:, ~df.columns.str.endswith('_dup')]

        # Fill missing values in the merged DataFrame with an empty string
        df.fillna("", inplace=True)

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
        df_first = df.merge(agg_data_targets, on=group_cols, how='left')

        # Keep only the specified columns in df_first
        keep_col_names = [
            'Target ID', 'Gene name', 'Protein name', 'Receptor family', 'Ligand type',
            'Class', 'Literature', 'Novelty (Pharos)', 'IDG',
            'Total', 'Active', 'Inactive', 'All_Max_Phase', 'All_Drugs', 'All_Agents',
            'Stimulatory_max_phase', 'Stimulatory_Drugs', 'Stimulatory_Agents',
            'Inhibitory_max_phase', 'Inhibitory_Drugs', 'Inhibitory_Agents',
        ]

        # Keep only the specified columns in df_first
        df_first = df_first[keep_col_names]

        # Drop duplicates based on the reduced set of columns
        df_first.drop_duplicates(inplace=True)

        # Fetch data from AssayExperiment with select_related for performance
        assay_data = (
            AssayExperiment.objects.select_related("protein", "ligand")
            .values("protein", "ligand", "value_type")
        )

        # Convert to DataFrame
        assay_df = pd.DataFrame(list(assay_data))

        # Step 1: Count occurrences of "pEC50" and "pIC50" for each target-ligand pair
        pair_counts = assay_df.groupby(["protein", "ligand", "value_type"]).size().unstack(fill_value=0)

        # Step 2: Identify and exclude pairs that have both "pEC50" and "pIC50"
        pairs_with_both = pair_counts[(pair_counts.get("pEC50", 0) > 0) & (pair_counts.get("pIC50", 0) > 0)].reset_index()
        remaining_pairs = assay_df[~assay_df.set_index(["protein", "ligand"]).index.isin(
            pairs_with_both.set_index(["protein", "ligand"]).index
        )]

        # Step 3: Count remaining unique target-ligand pairs and their "pEC50" and "pIC50"
        remaining_counts = (
            remaining_pairs.groupby(["protein", "ligand", "value_type"])
            .size()
            .unstack(fill_value=0)
            .reset_index()
        )

        # Step 4: Aggregate at the target level
        agg_results = remaining_counts.groupby("protein").agg(
            Total_Ligands=("ligand", "nunique"),
            pEC50_Count=("pEC50", "sum"),
            pIC50_Count=("pIC50", "sum"),
        ).reset_index()

        # Rename "protein" to "Target ID" to match df_first
        agg_results.rename(columns={"protein": "Target ID"}, inplace=True)

        # Merge the new aggregated data with df_first
        df_first = pd.merge(df_first, agg_results, on="Target ID", how="left")

        # Fill NaN values in the newly added columns with 0
        df_first[["Total_Ligands", "pEC50_Count", "pIC50_Count"]] = df_first[
            ["Total_Ligands", "pEC50_Count", "pIC50_Count"]
        ].fillna("")

        # Replace 0s with empty strings in relevant count columns
        df_first[["Total_Ligands", "pEC50_Count", "pIC50_Count", 
                "Total", "Active", "Inactive", "All_Drugs", "All_Agents", 
                "Stimulatory_max_phase", "Stimulatory_Drugs", "Stimulatory_Agents", 
                "Inhibitory_max_phase", "Inhibitory_Drugs", "Inhibitory_Agents"]] = \
            df_first[["Total_Ligands", "pEC50_Count", "pIC50_Count", 
                    "Total", "Active", "Inactive", "All_Drugs", "All_Agents", 
                    "Stimulatory_max_phase", "Stimulatory_Drugs", "Stimulatory_Agents", 
                    "Inhibitory_max_phase", "Inhibitory_Drugs", "Inhibitory_Agents"]].replace(0, "")
        
        # Step 1: Query the TissueExpression model
        tissue_datatable = TissueExpression.objects.select_related('tissue').values(
            'protein',        # Target ID
            'tissue__name',  # Tissue name
            'value'          # Expression value
        )

        # Step 2: Convert the QuerySet to a DataFrame
        df_tissue_table = pd.DataFrame(list(tissue_datatable))

        # Step 3: Rename columns for clarity
        df_tissue_table.rename(columns={
            'protein': 'Target ID',
            'tissue__name': 'Tissue Name',
            'value': 'Expression Value'
        }, inplace=True)

        # Step 4: Clean tissue names (e.g., replace underscores, capitalize)
        def clean_tissue_name(name):
            # Replace underscores with spaces, remove trailing "_X", and capitalize words
            name = name.replace("_", " ").rstrip(" 1234567890").title()
            return name

        df_tissue_table['Tissue Name'] = df_tissue_table['Tissue Name'].apply(clean_tissue_name)

        # Step 5: Pivot the DataFrame to make Tissue Names columns
        df_tissue_pivot = df_tissue_table.pivot_table(
            index='Target ID',         # Use Target ID as the index
            columns='Tissue Name',     # Use Tissue Name as columns
            values='Expression Value', # Use Expression Value as values
            aggfunc='first'            # If there are duplicates, take the first value
        )

        # Step 6: Reset index to turn Target ID back into a column
        df_tissue_pivot.reset_index(inplace=True)

        # Step 7: Fill NaN values with 0 (optional, based on your needs)
        df_tissue_pivot.fillna("", inplace=True)

        # Step 8: Generate a list of tissue column names (excluding Target ID)
        tissue_column_list = sorted([col for col in df_tissue_pivot.columns if col != 'Target ID'])

        # Merge df_first with df_tissue_pivot on 'Target ID'
        merged_df = pd.merge(df_first, df_tissue_pivot, on="Target ID", how="left")

        # Fill NaN values with empty strings for all columns
        merged_df.fillna("", inplace=True)

        # Add cancer to table 1

         # Add cancer data
        cancer_data = Protein.objects.select_related('Protein').values(
            'cancer__protein',
            'cancer__cancer__name',
            'cancer__expression__max_expression',
        )
        
        cancer_df = pd.DataFrame(list(cancer_data))
        cancer_df.rename(columns={
            'cancer__protein': 'Target ID',
            'cancer__cancer__name': 'Cancer',
            'cancer__expression__max_expression': 'Expression'
        }, inplace=True)

        # Filter only shared target IDs
        cancer_df = cancer_df[cancer_df['Target ID'].isin(target_ids)]

        # Ensure Target ID is an integer in cancer_df
        cancer_df['Target ID'] = cancer_df['Target ID'].astype(int)

       # Step 1: Calculate `NumberOfCancers` (sum of Medium and High for each Target ID)
        cancer_df['IsRelevantExpression'] = cancer_df['Expression'].isin(['Medium', 'High']).astype(int)
        number_of_cancers = (
            cancer_df.groupby('Target ID')['IsRelevantExpression']
            .sum()
            .reset_index(name='NumberOfCancers')
        )

        # Step 2: Pivot the cancer data
        cancer_pivot = cancer_df.pivot_table(
            index='Target ID',  # Rows: Target ID
            columns='Cancer',   # Columns: Cancer types
            values='Expression',  # Values: Expression levels
            aggfunc='first'      # Use the first value (only one row per Target ID + Cancer)
        )

        # Capitalize Cancer column names
        cancer_pivot.columns = [col.capitalize() for col in cancer_pivot.columns]

        # Reset index to make Target ID a column
        cancer_pivot.reset_index(inplace=True)

        # Step 3: Merge `NumberOfCancers` into the pivot table
        final_cancer_table = pd.merge(
            cancer_pivot,
            number_of_cancers,
            on='Target ID',
            how='left'
        )

        # Step 4: Merge cancer data with the main table
        Table1 = pd.merge(merged_df, final_cancer_table, on='Target ID', how='left')

        # Fill NaN values with empty strings
        Table1.fillna("", inplace=True)

        # Exclude 'Target ID' and sort the remaining columns
        cancer_column_list = sorted([col for col in cancer_pivot.columns if col != 'Target ID'])

        ########################################
        # Disease indications and associations #
        ########################################
        
        # Fetch the second data table with indication and master level
        table_data_2 = Drugs.objects.select_related(
            'target__family__parent__parent__parent',
            'moa',
            'indication'
        ).values(
            'target',  # Target ID
            'target__entry_name',  # Gene name
            'target__name',  # Protein name
            'target__family__parent__name',  # Receptor family
            'target__family__parent__parent__name',  # Ligand type
            'target__family__parent__parent__parent__name',  # Class
            'ligand',  # Ligand ID
            'indication__title',  # Indication title
            'indication__slug',  # Indication slug
            'disease_association__association_score',  # Association score
            'drug_status'  # Status
        )

        # Convert to DataFrame
        df_table_2 = pd.DataFrame(list(table_data_2))
        df_table_2.rename(columns={
            'target': 'Target ID',
            'target__entry_name': 'Gene name',
            'target__name': 'Protein name',
            'target__family__parent__name': 'Receptor family',
            'target__family__parent__parent__name': 'Ligand type',
            'target__family__parent__parent__parent__name': 'Class',
            'ligand': "Ligand ID",
            'indication__title': 'Indication',
            'indication__slug': 'Indication Slug',
            'disease_association__association_score': 'Association Score',
            'drug_status': 'Status'
        }, inplace=True)

        # Classify based on status
        df_table_2['Classification'] = df_table_2['Status'].apply(
            lambda x: 'Drug' if x == 'Approved' else 'Agent'
        )

        # Fetch top-level indications (master indications)
        master_indications = Indication.objects.filter(level=0).values('slug', 'title')

        # Create a mapping dictionary for master indications
        master_mapping = {entry['slug']: entry['title'] for entry in master_indications}

        # Map master titles to indications based on the top-level slug
        df_table_2['Master Indication'] = df_table_2['Indication Slug'].apply(
            lambda slug: master_mapping.get(slug.split("_")[0], "")
        )

        # Fetch ATC data
        atc_data = ATCCodes.objects.values(
            'ligand', 'code__index', 'name_0', 'name_1'
        )

        # Convert to DataFrame
        atc_df = pd.DataFrame(list(atc_data))
        atc_df.rename(columns={
            'ligand': 'Ligand ID',
            'code__index': 'ATC Code',
            'name_0': 'ATC Parent Name',
            'name_1': 'ATC Name'
        }, inplace=True)

        # Merge ATC data with the second table
        df_table_2 = df_table_2.merge(atc_df, on='Ligand ID', how='left')

        # Fill missing values with empty strings
        df_table_2.fillna("", inplace=True)

        # Aggregate the table efficiently
        def aggregate_table_fast(df):
            # Select unique rows for non-aggregated columns
            unique_df = df.drop_duplicates(subset=[
                'Target ID', 'Gene name', 'Protein name', 'Receptor family',
                'Ligand type', 'Class', 'Master Indication', 'Indication',
                'ATC Code', 'ATC Parent Name', 'ATC Name', 'Association Score'
            ])

            # Compute counts for 'Drug' and 'Agent' using crosstab
            classification_counts = pd.crosstab(
                index=[df['Target ID'], df['Indication']],  # Grouping keys
                columns=df['Classification']
            ).reset_index()

            # Rename columns for clarity
            classification_counts.rename(columns={
                'Drug': 'Drug Count',
                'Agent': 'Agent Count'
            }, inplace=True)

            # Merge classification counts back into the unique rows
            result = pd.merge(
                unique_df,
                classification_counts,
                on=['Target ID', 'Indication'],
                how='left'
            )

            # Fill NaN counts with 0
            result[['Drug Count', 'Agent Count']] = result[['Drug Count', 'Agent Count']].fillna(0).astype(int)

            return result

        # Aggregate the table to collapse unique rows and compute counts quickly
        df_table_2 = aggregate_table_fast(df_table_2)

        def remove_duplicates_and_fill_gaps(df):
            # Replace empty strings with NaN for better handling of missing values
            df.replace("", pd.NA, inplace=True)

            # Sort by completeness (count of non-missing values) and drop exact duplicates
            df = df.sort_values(by=df.columns.to_list(), key=lambda col: col.notna(), ascending=False)
            df = df.drop_duplicates(subset=df.columns.difference(['ATC Code', 'ATC Name', 'ATC Parent Name']), keep='first')

            # Fill gaps by replacing NaN back with empty strings for display purposes
            df.fillna("", inplace=True)

            return df
        # Remove duplicates and prioritize completeness
        df_table_2 = remove_duplicates_and_fill_gaps(df_table_2)

        # Retain only the specified columns
        keep_cols = [
            'Target ID', 'Gene name', 'Protein name', 'Receptor family', 'Ligand type',
            'Class', 'Master Indication', 'Indication', 'ATC Code',
            'ATC Parent Name', 'ATC Name', 'Association Score', 'Drug Count', 'Agent Count'
        ]
        df_table_2 = df_table_2[keep_cols]

        # Replace 0s with empty strings for display purposes
        df_table_2.replace(0, "", inplace=True)

        ####################
        # Untapped targets #
        ####################

        # Extract unique Target IDs from df_table_2
        existing_target_ids = df_table_2['Target ID'].unique()

        # Filter table_data_3 to exclude existing Target IDs
        filtered_table_data_3 = IndicationAssociation.objects.select_related(
            'target__family__parent__parent__parent',
            'indication'
        ).exclude(
            target__in=existing_target_ids  # Exclude targets that already exist in df_table_2
        ).values(
            'target',  # Target ID
            'target__entry_name',  # Gene name
            'target__name',  # Protein name
            'target__family__parent__name',  # Receptor family
            'target__family__parent__parent__name',  # Ligand type
            'target__family__parent__parent__parent__name',  # Class
            'indication__title',  # Indication title
            'indication__slug', # Indication slug
            'association_score',  # Association score
        )

        # Convert the filtered QuerySet to a DataFrame
        df_table_3 = pd.DataFrame(list(filtered_table_data_3))

        df_table_3.rename(columns={
            'target': 'Target ID',
            'target__entry_name': 'Gene name',
            'target__name': 'Protein name',
            'target__family__parent__name': 'Receptor family',
            'target__family__parent__parent__name': 'Ligand type',
            'target__family__parent__parent__parent__name': 'Class',
            'indication__title': 'Indication',
            'indication__slug': 'Indication Slug',
            'association_score': 'Association Score',
        }, inplace=True)

        # Fetch master indications and create mapping
        master_indications = Indication.objects.filter(level=0).values('slug', 'title')
        master_mapping = {entry['slug']: entry['title'] for entry in master_indications}

        # Add Master Indication to df_table_3
        df_table_3['Master Indication'] = df_table_3['Indication Slug'].apply(
            lambda slug: master_mapping.get(slug.split("_")[0], "")
        )

        # Convert the final DataFrame to JSON
        json_records_targets = Table1.to_json(orient='records')
        json_records_table_2 = df_table_2.to_json(orient='records')
        json_records_table_3 = df_table_3.to_json(orient='records')
        context['Full_data'] = json_records_targets
        context['Disease_table'] = json_records_table_2
        context['Untapped_table'] = json_records_table_3
        context['Tissue_cols'] = tissue_column_list
        context['cancer_cols'] = cancer_column_list

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
        


        # for key in sorted(IDC_hierarchy_drugs):
        #     print(key, IDC_hierarchy_drugs[key]['title'],len(IDC_hierarchy_drugs[key]['children']))
        
        # start=22
        # list_of_drugs = [start]
        # for key in sorted(IDC_hierarchy):
        #     last = list_of_drugs[-1]
        #     current = len(IDC_hierarchy[key]['children'])
        #     next = last+current+1
        #     list_of_drugs.append(next)
        #     print(key,IDC_hierarchy[key]['title'],current,next)
        # print(list_of_drugs)


        # start=706
        # list_of_drugs = [start]
        # for key in sorted(IDC_hierarchy_drugs):
        #     last = list_of_drugs[-1]
        #     current = len(IDC_hierarchy_drugs[key]['children'])
        #     next = last+current+1
        #     list_of_drugs.append(next)
        #     print(key,IDC_hierarchy_drugs[key]['title'],current,next)
        # print(list_of_drugs)

        # start=1363
        # list_of_drugs = [start]
        # for key in sorted(IDC_hierarchy_agents):
        #     last = list_of_drugs[-1]
        #     current = len(IDC_hierarchy_agents[key]['children'])
        #     next = last+current+1
        #     list_of_drugs.append(next)
        #     print(key,IDC_hierarchy_agents[key]['title'],current,next)
        # print(list_of_drugs)
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
    #code = '4A8Z'
    indication_data = Drugs.objects.filter(indication__code=code).prefetch_related('ligand',
                                                                                        'target',
                                                                                        'indication',
                                                                                        'indication__uri')

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
        indication_code = record.indication.title.capitalize()
        indication_uri = record.indication.uri.index
        ligand_name = record.ligand.name.capitalize()
        indication_0 = record.indication.get_level_0().title
        uri = record.indication.uri.index
        ligand_id = record.ligand.id
        protein_name = record.target.name
        target_name = record.target.entry_name
        #check for each value if it exists and retrieve the source node value
        if indication_code not in caches['indication']:
            sankey['nodes'].append({"node": node_counter, "name": indication_code, "url":'https://icd.who.int/browse/2024-01/mms/en#'+indication_uri})
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
    context['indication_uri'] = indication_uri
    context['indication'] = indication_name.capitalize()
    context['sankey'] = json.dumps(sankey)
    context['points'] = total_points
    context['targets'] = list(caches['entries'])
    context['ligands'] = list(caches['ligands'])
    return render(request, 'indication_detail.html', context)
