from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from django.db import IntegrityError
from django.utils.text import slugify

from common.models import WebResource, WebLink, Publication
from protein.models import Protein, TissueExpression, CancerType, CancerExpression, Tissues, ExpressionValue
from drugs.models import Drugs2024, Indication, IndicationAssociation
from ligand.models import Ligand, LigandID, LigandType, LigandRole
from common.tools import test_model_updates

import pandas as pd
import os
import django.apps
import logging

class Command(BaseCommand):
    help = 'Build Drug and NHS Data'

    publication_cache = {}

    def add_arguments(self, parser):
        parser.add_argument('--filename', action='store', dest='filename',
                            help='Filename to import. Can be used multiple times')

    logger = logging.getLogger(__name__)

    # source file directory
    data_dir = os.sep.join([settings.DATA_DIR, 'drug_data'])
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    test_model_updates(all_models, tracker, initialize=True)

    def add_arguments(self, parser):
        parser.add_argument("--test_run",
                            action="store_true",
                            help="Skip this during a test run",
                            default=False)
        parser.add_argument('-u', '--purge',
                            action='store_true',
                            dest='purge',
                            default=False,
                            help='Purge existing ligand records')

    def purge_data(self):
        try:
            Drugs2024.objects.all().delete()
            Indication.objects.all().delete()
            TissueExpression.objects.all().delete()
            Tissues.objects.all().delete()
            CancerType.objects.all().delete()
            CancerExpression.objects.all().delete()
            ExpressionValue.objects.all().delete()
            IndicationAssociation.objects.all().delete()
            # NHSPrescribings.objects.all().delete()
        except Exception as msg:
            print(msg)
            self.logger.warning('Existing data cannot be deleted')
            self.logger.warning('Drugs module not found: nothing to delete.')

    def handle(self, *args, **options):
        if options["test_run"]:
            print("Skipping in test run")
            return

        if options['purge']:
            print("Started purging bioactivity and ligand data")
            self.purge_data()
            self.tracker = {}
            test_model_updates(self.all_models, self.tracker, initialize=True)
            print("Ended purging data")

        print("\n\nWelcome to the Drugs2024 build process. Build steps will be printed.")
        print("##### STEP 0 START #####")
        print("\n\nStarted parsing data and setting up different dataframes")
        indication_df, tissue_df, cancer_df, drug_df, association_data = Command.setup_data()
        print("##### STEP 1 START #####")
        print("\n\nStarted parsing Indication data and building Indication Model")
        Command.generate_indications(indication_df)    #DONE
        print("\n\nIndication Model built. Performing checks")
        test_model_updates(self.all_models, self.tracker, check=True)
        print("##### STEP 2 START #####")
        print("\n\nStarted parsing Tissue Expression data and building TissueExpression Model")
        Command.generate_tissue_expression(tissue_df)  #DONE
        print("\n\TissueExpression Model built. Performing checks")
        test_model_updates(self.all_models, self.tracker, check=True)
        print("##### STEP 3 START #####")
        print("\n\nStarted parsing Cancer Prognostics data and building CancerPrognostics Model")
        Command.generate_cancer_prog(cancer_df)        #DONE
        print("\n\CancerPrognostics Model built. Performing checks")
        test_model_updates(self.all_models, self.tracker, check=True)
        print("##### STEP 4 START #####")
        print("\n\nStarted parsing Disease Association data and building IndicationAssociation Model")
        Command.generate_disease_associations(association_data)        #DONE
        print("\n\IndicationAssociation Model built. Performing checks")
        test_model_updates(self.all_models, self.tracker, check=True)
        print("##### STEP 5 START #####")
        print("\n\nStarted parsing Drug data and building Drug2024 Model")
        Command.create_drug_data(drug_df)
        print("\n\Drug Model built. Performing checks")
        test_model_updates(self.all_models, self.tracker, check=True)

    def setup_data():
        #Loading the two csv files
        all_data = Command.read_csv_data('03_FinalData.csv') #old one:03_FINAL_DATA_UPDATED.csv new one:03_FinalData.csv
        target_data = Command.read_csv_data('08_TargetPrioritazion_AllData.csv')
        atc_codes = Command.read_csv_data('03_ligand_ATCCodes.csv')
        opentarget_scores = Command.read_csv_data('08_TargetPrioritazion_Data_DiseaseAssociations.csv')
        #getting the cancer data for each protein
        cancer_data = target_data[['entry_name','Cancer','MaxExpression']]
        #Clean the cancer data from NaN data columns
        cancer_data = cancer_data.dropna(subset=['Cancer','MaxExpression'], how='all').drop_duplicates()
        #selecting tissues data columns plus the protein entry name
        tissues_cols = ['entry_name'] + [col for col in target_data.columns if col.startswith('Tissue')]
        tissues_data = target_data[tissues_cols]
        #Clean the tissues data from NaN data columns
        tissues_data = tissues_data.dropna(subset=[col for col in tissues_data.columns if col.startswith('Tissue')], how='all').drop_duplicates()
        #getting the indications data from the all_data file
        indication_data = all_data[['ICD_Title', 'ICD_Code']].drop_duplicates()
        #remove cancer and tissue columns from data
        columns_to_keep = ['entry_name', 'ICD_Code', 'genetic_association', 'affected_pathway', 'somatic_mutation', 'animal_model', 'novelty_score']
        #define a filtered version of the target_data dataframe
        shaved_data = target_data[columns_to_keep]
        shaved_data = shaved_data.dropna(subset=['genetic_association', 'affected_pathway', 'somatic_mutation', 'animal_model', 'novelty_score'], how='all').drop_duplicates()
        #merge dataframes to have everything connected for Drugs model
        drug_data = pd.merge(all_data, shaved_data, on=['entry_name', 'ICD_Code'], how='inner')
        drug_data = pd.merge(drug_data, atc_codes, on=['ligand_name'], how='inner')
        #Drop the duplicates
        drug_data = drug_data.drop_duplicates()
        opentarget_scores = opentarget_scores.drop_duplicates()
        opentarget_scores = opentarget_scores.where(pd.notnull(opentarget_scores), None)
        #Somehow, opentargets has new indications that we need to add to the database
        opentarget_ICD = opentarget_scores[['ICD_Title', 'ICD_Code']].drop_duplicates()
        indication_data = pd.concat([indication_data, opentarget_ICD]).drop_duplicates()
        return indication_data, tissues_data, cancer_data, drug_data, opentarget_scores

    @staticmethod
    def transform_column_name(col_name):
        # Split on " - " and take the second part, if it exists; otherwise, return the original
        tissue_name = col_name.split(" - ")[1] if " - " in col_name else col_name
        # Remove the "[nTPM]" or any other bracketed part
        tissue_name = tissue_name.split(' [')[0]
        # Replace spaces with underscores and convert to lowercase
        tissue_name = tissue_name.replace(" ", "_").lower()
        return tissue_name

    def generate_tissue_expression(tissues_data):
        #process the column headers
        tissues_data.columns = [Command.transform_column_name(col) for col in tissues_data.columns]
        tissues = list(tissues_data.columns)[1:]
        for i, row in tissues_data.iterrows():
            protein = Command.fetch_protein(row['entry_name'])
            for tissue in tissues:
                slug_tissue = slugify(tissue)
                t, _ = Tissues.objects.get_or_create(slug=slug_tissue, name=tissue)
                record, _ = TissueExpression.objects.get_or_create(value=row[tissue], protein=protein, tissue=t)

    def generate_disease_associations(association_data):
        #parse che association data to fill in the specific model
        #need to fetch protein and also indication
        for i, row in association_data.iterrows():
            protein = Command.fetch_protein(row['entry_name'])
            indication = Command.fetch_indication(row['ICD_Title'])
            association, _ = IndicationAssociation.objects.get_or_create(target = protein,
                                                                        indication = indication,
                                                                        association_score = row['Score'],
                                                                        ot_genetics_portal = row['ot_genetics_portal'],
                                                                        gene_burden = row['gene_burden'],
                                                                        clingen = row['clingen'],
                                                                        gene2phenotype = row['gene2phenotype'],
                                                                        orphanet = row['orphanet'],
                                                                        genomics_england = row['genomics_england'],
                                                                        uniprot_literature = row['uniprot_literature'],
                                                                        uniprot_variants = row['uniprot_variants'],
                                                                        sysbio = row['sysbio'],
                                                                        cancer_gene_census  = row['cancer_gene_census'],
                                                                        cancer_biomarkers = row['cancer_biomarkers'],
                                                                        intogen = row['intogen'],
                                                                        eva = row['eva'],
                                                                        eva_somatic = row['eva_somatic'],
                                                                        chembl = row['chembl'],
                                                                        slapenrich = row['slapenrich'],
                                                                        crispr = row['crispr'],
                                                                        crispr_screen = row['crispr_screen'],
                                                                        reactome = row['reactome'],
                                                                        europepmc = row['europepmc'],
                                                                        expression_atlas = row['expression_atlas'],
                                                                        impc = row['impc'])

    def generate_indications(indication_data):
        #Create the reference for the Ontology web resource
        EMBL = {
            'name': 'EMBL Ontology',
            'url': 'https://www.ebi.ac.uk/ols4/ontologies/efo/classes?short_form=$index'
        }
        #Create the reference for the ICD Ontology web resource
        ICD = {
            'name': 'ICD-11 Ontology',
            'url': 'https://icd.who.int/browse/2024-01/mms/en#$index'
        }
        #Create the reference for the ATC Ontology web resource
        ATC = {
            'name': 'ATC Ontology',
            'url': 'https://atcddd.fhi.no/atc_ddd_index/?code=$index'
        }
        wr_embl, created = WebResource.objects.get_or_create(slug='indication_embl', defaults=EMBL)
        wr_icd, created = WebResource.objects.get_or_create(slug='indication_icd', defaults=ICD)
        wr_atc, created = WebResource.objects.get_or_create(slug='indication_atc', defaults=ATC)
        #Define the Web Resource
        web_resource = WebResource.objects.get(slug='indication_icd')
        for i, row in indication_data.iterrows():
            indication = Indication()
            indication.name = row['ICD_Title']
            indication.code, created = WebLink.objects.get_or_create(index=row['ICD_Code'], web_resource=web_resource)
            indication.save()

    def generate_cancer_prog(cancer_data):
        #parse che cancer data and collate type and expression.
        #Then fill the ManyToMany
        for i, row in cancer_data.iterrows():
            slug_cancer = slugify(row['Cancer'])
            protein, _ = Protein.objects.get_or_create(entry_name=row['entry_name'])
            ct, _ = CancerType.objects.get_or_create(slug=slug_cancer, name=row['Cancer'])
            exp, _ = ExpressionValue.objects.get_or_create(max_expression = row['MaxExpression'])
            ce = CancerExpression()
            ce.cancer = ct
            ce.expression = exp
            ce.save()
            protein.cancer.add(ce)


    def create_drug_data(drug_data):
        #Parse the drug dataframe
        for i, row in drug_data.iterrows():
            #fetch the ligand or generate a new ligand record if there is no match
            ligand = Command.fetch_ligand(row)
            #Then add the different references (PubChem, DrugBank and ChEMBL)
            #TODO: add also UNII and CAS as values in the ManyToMany
            Command.add_drug_references(ligand, row)
            #Fetch the reference protein
            protein = Command.fetch_protein(row['entry_name'])
            #fetch the indication
            indication = Command.fetch_indication(row['ICD_Title'])
            #Fetch the ligand action (role)
            moa = Command.fetch_role(row['Action'])
            #Fetch the inidcation association
            association = Command.fetch_association(row['entry_name'], row['ICD_Title'])
            #to be human readable instead of numerical values (ask David)
            drug, _ = Drugs2024.objects.get_or_create(charge=row['Charge'],
                                                      complexity=row['Complexity'],
                                                      tpsa=row['TPSA'],
                                                      drug_status=row['Drug_Status'],
                                                      approval_year=row['Approval_Year'] if pd.notna(row['Approval_Year']) else None,
                                                      indication_max_phase=row['IndicationMaxPhase'],
                                                      affected_pathway=row['affected_pathway'] if pd.notna(row['affected_pathway']) else None,
                                                      somatic_mutation=row['somatic_mutation'] if pd.notna(row['somatic_mutation']) else None,
                                                      similarity_to_model=row['animal_model'] if pd.notna(row['animal_model']) else None,
                                                      novelty_score=row['novelty_score'],
                                                      genetic_association=row['genetic_association'] if pd.notna(row['genetic_association']) else None,
                                                      indication_status=row['IndicationStatus'],
                                                      atc_code=row['ATC_Code'],
                                                      moa=moa,
                                                      indication=indication,
                                                      ligand=ligand,
                                                      disease_association=association,
                                                      target=protein)

            #Commented for the sake of testing
            #since it's calculated to have more than 300k unique pubs to be added

            # if pd.notna(row['PubMedID']):
            #     ref = row['PubMedID'].split(';')
            #     try:
            #         for pmid in ref:
            #             if pmid != '':
            #                 publication = fetch_publication(pmid)
            #                 drug.reference.add(publication)
            #     except Exception as e:
            #         print(f'The Drugs {pmid} publication was not added to the data base'))

    @staticmethod
    def fetch_ligand(row):
        """
        fetch ligands with Ligand model
        requires: ligand id.
        """
        #will perform several checks
        mapper = {'chembl_ligand': str(row['ChEMBLID']).split(';'),
                  'pubchem': str(row['PubChemCID']).split(';'),
                  'drugbank': str(row['DrugBankID']).split(';')}
        check = None
        #Check for match of Inchikey.
        #If inchikey field has multiple inchi, split them and generate a list
        #if inchi is 'nan', apply a workaround
        if ';' in str(row['InChiKey']):
            inchi_list = str(row['InChiKey']).split(';')
        elif pd.notna(row['InChiKey']):
            inchi_list = ['NOT AVAILABLE']
        else:
            inchi_list = [str(row['InChiKey'])]
        for inchi in inchi_list:
            try:
                check = Ligand.objects.get(inchikey=inchi)
                return check
            except Ligand.DoesNotExist:
                for key, values in mapper.items():
                    for code in values:
                        if code != 'nan':
                            try:
                                check = LigandID.objects.filter(index=code, web_resource__slug=key)
                                if len(check) == 1:
                                    return check[0].ligand
                                else:
                                    check = None
                            except LigandID.DoesNotExist:
                                continue
            if check == None:
                type = Command.fetch_type(row['Drug_Type'])
                #TODO: adjust the length of float numbers
                check, _ = Ligand.objects.get_or_create(name=row['ligand_name'],
                                                        ambiguous_alias=False,
                                                        hacc=row['HBondAceptorCount'] if pd.notna(row['HBondAceptorCount']) else None,
                                                        hdon=row['HBondDonorCount'] if pd.notna(row['HBondDonorCount']) else None,
                                                        inchikey=inchi if inchi !='NOT AVAILABLE' else None,
                                                        ligand_type=type,
                                                        logp=row['XLogP'] if pd.notna(row['XLogP']) else None,
                                                        mw=row['MolecularWeight'] if pd.notna(row['MolecularWeight']) else None,
                                                        rotatable_bonds=row['RotableBondCount'] if pd.notna(row['RotableBondCount']) else None,
                                                        smiles=row['SMILES'])
                return check

    @staticmethod
    def fetch_association(target, indication):
        """
        fetch indication association with indication title and target entry name
        requires: indication title and target entry name
        """
        try:
            association = IndicationAssociation.objects.get(target_id__entry_name=target, indication_id__name=indication)
            return association
        except:
            print('No association found for this pais or target-indication')
            return None

    @staticmethod
    def fetch_type(record):
        """
        fetch ligand type with type through a conversion dict
        requires: type
        """
        conversion = {'Small molecule': 'small-molecule',
                      'Protein': 'protein',
                      'Antibody': 'protein',
                      'Oligonucleotide': 'peptide'}
        type_slug = conversion[record]
        lt, _ = LigandType.objects.get_or_create(slug=type_slug)
        return lt

    @staticmethod
    def add_drug_references(ligand, row):
        mapper = {'chembl_ligand': str(row['ChEMBLID']).split(';'),
                  'pubchem': str(row['PubChemCID']).split(';'),
                  'drugbank': str(row['DrugBankID']).split(';')}
        for key, values in mapper.items():
            for code in values:
                try:
                    check = LigandID.objects.get(index=code, ligand_id=ligand.id, web_resource__slug=key)
                except LigandID.DoesNotExist:
                    wr = WebResource.objects.get(slug=key)
                    LigandID(index=code, web_resource=wr, ligand_id=ligand.id).save()

    @staticmethod
    def fetch_role(action):
        """
        fetch ligand role based on action, from a conversion dict
        requires: action
        """
        conversion = {'Antagonist': 'Antagonist',
                      'Agonist': 'Agonist',
                      'Inverse agonist': 'Inverse agonist',
                      'Partial agonist': 'Partial agonist',
                      'partial agonist': 'Partial agonist',
                      'unknown': 'unknown',
                      'Modulator': 'NAM',
                      'NAM': 'NAM',
                      'Binding agent': 'Binding - unknown pharmacological activity',
                      'Inhibitor': 'Antagonist',
                      'PAM': 'PAM',
                      'partial antagonist': 'Antagonist',
                      'cross-linking agent': 'Agonist'
                      }
        lr = None
        if action in conversion.keys():
            query = conversion[action]
            role_slug = slugify(query)
            lr, _ = LigandRole.objects.get_or_create(slug=role_slug, defaults={'name': query})
        return lr

    @staticmethod
    def fetch_indication(name):
        """
        fetch indication with indication name
        requires: indication name
        """
        try:
            indication = Indication.objects.filter(name=name)[0]
            return indication
        except:
            print('No indication found for this entry name')
            return None

    @staticmethod
    def fetch_protein(target):
        """
        fetch receptor with Protein model
        requires: protein entry_name
        """
        try:
            protein = Protein.objects.get(entry_name=target)
            return protein
        except:
            print('No protein found for this entry name')
            return None

    @staticmethod
    def read_csv_data(filename):
        filepath = os.sep.join([Command.data_dir, filename])
        data = pd.read_csv(filepath, low_memory=False)
        return data

    @staticmethod
    def fetch_publication(publication_doi):
        """
        fetch publication with Publication model
        requires: publication doi or pmid
        """
        if pd.isna(publication_doi) is True:
            return None

        if ("ISBN" in publication_doi) or (publication_doi == '0'):
            return None

        try:
            float(publication_doi)
            publication_doi = str(int(publication_doi))
        except ValueError:
            pass

        if publication_doi.isdigit():  # assume pubmed
            pub_type = 'pubmed'
        else:  # assume doi
            pub_type = 'doi'

        if publication_doi not in Command.publication_cache:
            try:
                wl = WebLink.objects.get(
                    index=publication_doi, web_resource__slug=pub_type)
            except WebLink.DoesNotExist:
                try:
                    wl = WebLink.objects.create(
                        index=publication_doi, web_resource=WebResource.objects.get(slug=pub_type))
                except IntegrityError:
                    wl = WebLink.objects.get(
                        index=publication_doi, web_resource__slug=pub_type)

            try:
                pub = Publication.objects.get(web_link=wl)
            except Publication.DoesNotExist:
                pub = Publication()
                try:
                    pub.web_link = wl
                    pub.save()
                except IntegrityError:
                    pub = Publication.objects.get(web_link=wl)
                if pub_type == 'doi':
                    pub.update_from_doi(doi=publication_doi)
                elif pub_type == 'pubmed':
                    pub.update_from_pubmed_data(index=publication_doi)
                try:
                    pub.save()
                except Exception as e:
                    # if something off with publication, skip.
                    print("Build drugs Publication fetching error | module: fetch_publication. Row # is : " +
                          str(publication_doi) + ' ' + pub_type)
                    print(f'{type(e).__name__} {e}')

            Command.publication_cache[publication_doi] = pub
        else:
            pub = Command.publication_cache[publication_doi]

        return pub
