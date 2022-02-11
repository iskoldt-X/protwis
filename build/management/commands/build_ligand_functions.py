#from django.db import connection
#from django.db import IntegrityError
from django.db.models import Q
from ligand.models import Ligand, LigandID, LigandType, LigandRole
from common.models import WebResource
from common.tools import get_or_create_url_cache, fetch_from_web_api

import time
import os
import re

import datamol as dm

external_sources = ["pubchem", "gtoplig", "chembl_ligand", "drugbank", "drug_central"]
def get_or_create_ligand(name, ids = {}, lig_type = "small-molecule", unichem = False, extended_matching = True):
    """ This function tries to obtain a small molecule Ligand object.

        If the ligand already exists it will return the corresponding object. If
        not, it will processes the provided data and create the ligand object.
        When a ligand does not exist, or cannot be created, with the given input
        it will return None

        Parameters
        ----------
        name : str
            Name of the ligand
        ids : dict, optional
            Dictionary where the key is the id_type and the value is the id
            <id_type> is e.g. pubchem, gtp, chembl, etc. (slug of the web_resource)
            <id> is the Identifier itself
        ids : str, optional
            String indicating the ligand type (according to the LigandType slugs)
            By default this is set to "small-molecule".
        unichem : bool, optional
            Boolean indicating whether or not to match the compound to other
            sources via UniChem using the GtP ID. By default this is disabled.
        extended_matching : bool, optional
            Boolean indicating whether or not to match the compound
            sources via UniChem using the GtP ID. By default this is enabled.

        Returns
        -------
        obj
            Ligand object or None

        Raises
        ------
        ...Error
            If no ...
    """

    ligand = None

    # Check and filter IDs
    for type_id in list(ids.keys()):
        if ids[type_id] == None or ids[type_id] == "None" or ids[type_id] == "":
            del ids[type_id]
        elif ids[type_id].isnumeric():
            ids[type_id] = int(ids[type_id])
        elif is_float(ids[type_id]):
            ids[type_id] = int(float(ids[type_id]))

    # No IDs are provided, so the only way forward is a name match
    # Very short name = not unique - minimum is length of 3
    # TODO - add filter for numerical names and typical cpd X names
    if len(ids) == 0 and len(name) >= 3:
        results = Ligand.objects.filter(name__iexact=name)
        if results.count() == 1:
            ligand = results.first()
            # DEBUGGING
            print("Semi-ambiguous matching using name")
        # DEBUGGING
        elif results.count() > 1:
            print("Ambiguous name (", name,") as it has ", results.count(), " corresponding entries")
        else:
            print("No ligand found with name",name)
    else:
        # Try to find ligand via PubChem ID, Guide to Pharmacology ID, ChEMBL ID
        for type in external_sources:
            if ligand == None and type in ids:
                ligand = get_ligand_by_id(type, ids[type])

        # Peptide or protein entry Sequence - filter gtoplig as those entries are standalone and should not be merged
        if ligand == None and "sequence" in ids and len(ids["sequence"]) > 3 and "gtoplig" not in ids:
            result = Ligand.objects.filter(sequence = ids["sequence"])
            if result.count() > 0:
                ligand = result.first()
                # DEBUGGING
                if result.count() > 1:
                    print("Ambiguous sequence (", ids["sequence"],") as it has ", results.count(), " corresponding entries")

        # UniProt ID if there's no sequence or other IDs
        # TODO figure out best way of matching when sequence and UniProt ID are mixed as there can be multiple variants
        elif ligand == None and "uniprot" in ids and len(set.intersection(set(ids.keys()), set(external_sources)))==0:
            result = Ligand.objects.filter(uniprot = ids["uniprot"], sequence = None)
            if result.count() > 0:
                ligand = result.first()
                # DEBUGGING
                if result.count() > 1:
                    print("Ambiguous match - uniprot (", ids["uniprot"],") as it has ", results.count(), " corresponding entries")

        # InChiKey
        if ligand == None and "inchikey" in ids:
            ligand = get_ligand_by_inchikey(ids["inchikey"])

        # SMILES and calculated InChIKeys
        # clean_inchi = None
        if ligand == None and "smiles" in ids:
            result = Ligand.objects.filter(smiles = ids["smiles"])
            if result.count() > 0:
                ligand = result.first()
                # DEBUGGING
                if result.count() > 1:
                    print("Ambiguous SMILES (", ids["smiles"],") as it has ", results.count(), " corresponding entries")
            elif not "inchikey" in ids:
                # calculate inchikey from given SMILES and repeat inchikey check
                input_mol = dm.to_mol(ids["smiles"], sanitize=False)
                ligand = get_ligand_by_inchikey(dm.to_inchikey(input_mol))
            # elif ligand == None:
            #     # calculate cleaned InchiKey and repeate check
            #     clean_inchi = get_cleaned_inchikey(ids["smiles"])
            #     if clean_inchi != None:
            #         ligand = get_ligand_by_inchikey(clean_inchi)

        # Tried all default options => no ligand found => try UniChem IDs before creating new ligand
        if extended_matching and ligand == None:
            external_ids = list(set.intersection(set(ids.keys()), set(external_sources)))
            current_ids = ids.values()
            for type in external_ids:
                if ligand == None and type in ids:
                    unichem_ids = match_id_via_unichem(type, ids[type])
                    for row in unichem_ids:
                        if row["id"] not in current_ids:
                            ligand = get_ligand_by_id(row["type"], row["id"])
                            current_ids.append(row["id"])
                            if ligand != None:
                                # print("UNICHEM MATCH!!", name)
                                # print("MATCHING", type, ids[type], "via UniChem")
                                # print(ids)
                                # print(row)
                                # print("============")
                                break

        # If the ligand has not been found => create ligand using the provided IDs
        # TODO merge this into one function
        if ligand == None:
            # Creation of peptides and proteins
            if lig_type != "small-molecule" and "sequence" in ids:
                ligand = create_ligand_from_id(name, "sequence", ids["sequence"], lig_type)

            # Creation of small molecules and others without sequence/UniProt
            if ligand == None:
                for type in ["smiles", "pubchem", "gtoplig", "chembl_ligand"]:
                    if type in ids:
                        ligand = create_ligand_from_id(name, type, ids[type], lig_type)
                        if ligand != None:
                            break

            # Create empty ligand
            if ligand == None:
                print("Creating an empty ligand", name, ids)
                ligand = create_ligand_from_id(name, "", "", lig_type)
                ligand.ambiguous = True


        # Add missing IDs via (web)links to the ligand object
        if ligand != None:
            # Add provided entries to ligand when missing
            if "pdb" in ids and ligand.pdbe == None:
                ligand.pdbe = ids["pdb"]
            if "uniprot" in ids and ligand.uniprot == None:
                ligand.uniprot = ids["uniprot"]
            if "sequence" in ids and ligand.sequence == None and len(ids["sequence"]) < 1000:
                ligand.sequence = ids["sequence"]
            if "inchikey" in ids:
                ligand.inchikey = ids["inchikey"]

            # Create list of existing weblinks
            current_ids = []
            for wl in ligand.ids.all():
                current_ids.append(str(wl.index))

            ligand.save()
            for type_id in ids:
                # If the ID does not yet exist => add
                if ids[type_id] not in current_ids and type_id in external_sources:
                    wr = WebResource.objects.get(slug=type_id)
                    wl, created = LigandID.objects.get_or_create(ligand=ligand, index=ids[type_id], web_resource=wr)
                    ligand.ids.add(wl)
                    current_ids.append(str(ids[type_id]))

            # Match the GtP Identifier to other sources via UniChem
            if unichem and "gtoplig" in ids:
                unichem_ids = match_id_via_unichem("gtoplig", ids["gtoplig"])
                for row in unichem_ids:
                    if row["id"] not in current_ids:
                        if row["type"] == "pdb":
                            if ligand.pdbe == None:
                                ligand.pdbe = row["id"]
                                ligand.save()
                        else:
                            wr = WebResource.objects.get(slug=row["type"])
                            wl = LigandID.objects.get_or_create(ligand=ligand, index=row["id"], web_resource=wr)[0]
                            ligand.ids.add(wl)
    return ligand

unichem_src_types = {"1": "chembl_ligand", "2": "drugbank", "3": "pdb", "4": "gtoplig", "22": "pubchem", "34": "drug_central"}
unichem_src_types_inv = {v: k for k, v in unichem_src_types.items()}
def match_id_via_unichem(type, id):
    results = []
    cache_dir = ['unichem', 'id_match']
    if type in unichem_src_types_inv:
        type_id = unichem_src_types_inv[type]
        unichem_url = "https://www.ebi.ac.uk/unichem/rest/src_compound_id/$index/" + type_id
        cache_dir[1] = "id_match_" + type
        unichem = fetch_from_web_api(unichem_url, id, cache_dir)
        if unichem:
            for entry in unichem:
                if entry["src_id"] in unichem_src_types.keys():
                    results.append({"type" : unichem_src_types[entry["src_id"]], "id": entry["src_compound_id"]})
    return results


def get_ligand_by_id(type, id):
    result = Ligand.objects.filter(ids__index=id, ids__web_resource__slug=type)
    if result.count() > 0:
        if result.count() > 1:
            print("Multiple entries for the same ID - This should never happen - error", type, id)
        return result.first()
    else:
        return None

def create_ligand_from_id(name, type, id, lig_type):
    # create new ligand
    ligand = Ligand()
    ligand.name = name
    ligand.canonical = True # TODO: discuss if we should keep this
    ligand.ambiguous = False
    ligand.ligand_type = LigandType.objects.get_or_create(slug=lig_type, defaults={'name': lig_type})[0]

    # Obtain external data (if needed)
    try:

        # TODO - continue here and add peptide/protein support for web services

        if type == "sequence":
            ligand.sequence = sequence
        elif type == "uniprot":
            ligand.uniprot = id
        elif  type == "smiles":
            ligand.smiles = id
        elif type == "pubchem":
            # check cache
            cache_dir = ['pubchem', 'cid', 'property']
            url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/$index/property/IsomericSMILES,InChIKey/json'
            pubchem = fetch_from_web_api(url, id, cache_dir)

            if pubchem != False and pubchem['PropertyTable']['Properties'][0]:
                properties = pubchem['PropertyTable']['Properties'][0]
                if len(properties['IsomericSMILES'])>0:
                    ligand.smiles = properties['IsomericSMILES']
                if len(properties['InChIKey'])>0:
                    ligand.inchikey = properties['InChIKey']
        elif type == "gtoplig":
            # check cache
            cache_dir = ['guidetopharmacology', 'ligand_structure']
            url = 'http://www.guidetopharmacology.org/services/ligands/$index/structure'
            gtop = fetch_from_web_api(url, id, cache_dir)

            if gtop:
                if len(gtop["smiles"])>0:
                    ligand.smiles = gtop["smiles"]
                if len(gtop["inchiKey"])>0:
                    ligand.inchikey = gtop["inchiKey"]
                if len(gtop["oneLetterSeq"])>0 and len(gtop["oneLetterSeq"])<1000:
                    ligand.sequence = gtop["oneLetterSeq"]
                    url = 'https://www.guidetopharmacology.org/services/ligands/$index/databaseLinks?database=UniProtKB'
                    gtop_uniprot = fetch_from_web_api(url, id, cache_dir)
                    if gtop_uniprot and len(gtop["accession"])>0:
                        ligand.uniprot = gtop["accession"]
        elif type == "chembl_ligand":
            # check cache
            cache_dir = ['chembl', 'ligand_entry']
            url = 'https://www.ebi.ac.uk/chembl/api/data/molecule/$index/?format=json'
            chembl = fetch_from_web_api(url, id, cache_dir)

            # NOTE: canonical SMILES from ChEMBL are isomeric smiles
            if chembl and "molecule_structures" in chembl:
                if len(chembl["molecule_structures"]["canonical_smiles"]) > 0:
                    ligand.smiles = chembl["molecule_structures"]["canonical_smiles"]
                if len(chembl["molecule_structures"]["standard_inchi_key"]) > 0:
                    ligand.inchikey = chembl["molecule_structures"]["standard_inchi_key"]
    except:
        skip=True

    # perform RDkit calculations
    if lig_type == "small-molecule" and ligand.smiles != None and ligand.smiles != "":
        input_mol = dm.to_mol(ligand.smiles, sanitize=True)
        if input_mol:
            # Check if InChIKey has been set
            if ligand.inchikey == None:
                ligand.inchikey = dm.to_inchikey(input_mol)

            # Check if cleaned InChIKey has been set
            # if ligand.clean_inchikey == None:
                # ligand.clean_inchikey = get_cleaned_inchikey(ligand.smiles)

            # Calculate RDkit properties
            ligand.mw = dm.descriptors.mw(input_mol)
            ligand.rotatable_bonds = dm.descriptors.n_rotatable_bonds(input_mol)
            ligand.hacc = dm.descriptors.n_hba(input_mol)
            ligand.hdon = dm.descriptors.n_hbd(input_mol)
            ligand.logp = dm.descriptors.clogp(input_mol)
        else:
            print("Issues with molecule", name)
            # TODO - try again if we have other IDs

        # Verify InChiKey
        if ligand.inchikey == None and ligand.smiles != None:
            ligand.inchikey = dm.to_inchikey(input_mol)

        # Before storing - check one more time if we already have this ligand
        if ligand.inchikey != None:
            checkligand = get_ligand_by_inchikey(ligand.inchikey)
            if checkligand != None:
                return checkligand

        # if ligand.clean_inchikey != None and ligand.clean_inchikey != ligand.inchikey:
        #     checkligand = get_ligand_by_inchikey(ligand.clean_inchikey)
        #     if checkligand != None:
        #         return checkligand

        return ligand
    elif lig_type != "small-molecule" or type == "":
        # Before storing - check one more time if we already have this ligand based on sequence
        # TODO implement this double check
        # https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/71728433/XML
        # USE LXML for parsing => https://stackoverflow.com/questions/21628290/parsing-xml-with-xpath-in-python-3


        # https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL59?format=xml&_=1643906691331
        # https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL266481?format=xml
        #=> PEPTIDE parse HELM notation => no uniprot

        return ligand
    else:
        # Ligand could not be created
        return None

def get_ligand_by_inchikey(inchikey):
    result = Ligand.objects.filter(Q(inchikey = inchikey) | Q(clean_inchikey = inchikey))
    if result.count() > 0:
        if result.count() > 1:
            print("Multiple entries for the same InChIkey - This should never happen - error", inchikey)
        return result.first()
    else:
        return None

def get_cleaned_inchikey(smiles):
    return None
    # inchikey = None
    # try:
    #     input_mol = dm.to_mol(smiles, sanitize=True)
    #     if input_mol:
    #         # Sligthly odd set of steps to get optimal cleaned and dominant tautomer
    #         input_mol = dm.to_neutral(input_mol)
    #         input_mol = dm.fix_mol(input_mol)
    #         cleaned_smiles = dm.standardize_smiles(dm.to_smiles(input_mol), tautomer=True)
    #         cleaned_mol = dm.to_mol(cleaned_smiles, sanitize=False)
    #         inchikey = dm.to_inchikey(cleaned_mol)
    # except:
    #     skip = True
    # return inchikey

def is_float(element):
    try:
        float(element)
        return True
    except ValueError:
        return False
