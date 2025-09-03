# interaction/schrodinger_processor.py

from django.shortcuts import render, redirect
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings
from django.db.models import Count, Sum, Avg, Q
from django.utils.text import slugify
from django.conf import settings
import os
import yaml
from django.utils.text import slugify


# Import necessary Django models
from interaction.models import (
    ResidueFragmentInteraction,
    StructureLigandInteraction,
    ResidueFragmentInteractionType,
)

from structure.models import (
    Structure,
    PdbData,
    Rotamer,
    Fragment,
    StructureModel,
    StructureComplexModel,
    StructureExtraProteins,
    StructureVectors,
    StructureModelRMSD,
    StructureModelpLDDT,
    StructureAFScores,
)

# Note: Structure, Ligand, Residue, Rotamer are imported from their respective apps

from ligand.models import Ligand
from residue.models import Residue, Rotamer

import logging

logger = logging.getLogger(
    __name__
)  # Or 'gpcrdb.interaction.schrodinger' for more specific logging


def process_schrodinger_sm_interactions(
    current_structure_obj: Structure,
    current_ligand_db_obj: Ligand,
    ligand_pdb_het_code: str,
    pdb_code_str: str,
    schrodinger_interactions_dir_override: str = None,
) -> bool:
    """
    Processes Schrödinger pre-computed Small Molecule interaction data from a YAML file.
    Returns True if YAML was found and processed (regardless of new interactions created),
    Returns False if YAML was not found or a critical parsing error occurred,
    allowing the caller to decide on fallback.
    """
    logger.info(
        f"Processing Schrödinger SM interactions for {pdb_code_str} with ligand {ligand_pdb_het_code}"
    )

    schrodinger_interaction_dir = schrodinger_interactions_dir_override or getattr(
        settings, "SCHRODINGER_INTERACTIONS_DIR", None
    )

    if not schrodinger_interaction_dir:
        logger.error(
            "SCHRODINGER_INTERACTIONS_DIR is not set. Cannot process Schrödinger SM interactions."
        )
        return False  # Critical: cannot proceed

    # Ensure ligand_pdb_het_code and pdb_code_str are uppercase for filename consistency
    yaml_file_name = (
        f"{pdb_code_str.upper()}_{ligand_pdb_het_code.upper()}_interactions.yaml"
    )
    yaml_file_path = os.path.join(schrodinger_interaction_dir, yaml_file_name)

    if not os.path.exists(yaml_file_path):
        logger.warning(f"Schrödinger SM YAML file not found: {yaml_file_path}")
        return False  # File not found, indicates fallback is needed

    try:
        with open(yaml_file_path, "r") as f_yaml:
            schrodinger_data = yaml.safe_load(f_yaml)
    except yaml.YAMLError as e:
        logger.error(f"Error parsing YAML file {yaml_file_path}: {e}")
        return False  # Parsing error, indicates fallback is needed

    if (
        not schrodinger_data
        or "result" not in schrodinger_data
        or "interactions" not in schrodinger_data["result"]
    ):
        logger.warning(
            f"Invalid or empty YAML data in {yaml_file_path}. Expected 'result.interactions'."
        )
        return True  # File was found and "processed", but no data to add. Don't fallback if file implies "no interactions".

    # Get the main StructureLigandInteraction object.
    try:
        sli = StructureLigandInteraction.objects.get(
            structure=current_structure_obj,
            ligand=current_ligand_db_obj,
            pdb_reference=ligand_pdb_het_code.upper(),  # Ensure HET code matches casing
        )
    except StructureLigandInteraction.DoesNotExist:
        logger.error(
            f"StructureLigandInteraction not found for PDB {pdb_code_str}, "
            f"Ligand {current_ligand_db_obj.name} (HET: {ligand_pdb_het_code}). "
            "This SLI should have been created before calling this processor. Skipping."
        )
        return False  # SLI missing is a problem for associating interactions.
    except StructureLigandInteraction.MultipleObjectsReturned:
        logger.error(
            f"Multiple SLI objects for PDB {pdb_code_str}, Ligand {current_ligand_db_obj.name}. "
            "Ambiguous. Skipping."
        )
        return False

    interactions_created_count = 0
    for interaction_entry in schrodinger_data["result"]["interactions"]:
        receptor_res_info = interaction_entry["receptor_residue"]
        pdb_residue_number = int(receptor_res_info["pdb_number"])
        chain_id_from_yaml = receptor_res_info["chain_id"]  # Keep for logging/debugging
        residue_1_letter = receptor_res_info["name_1_letter"].upper()
        insertion_code = receptor_res_info.get("insertion_code")  # Handle if present

        # 1. Map receptor residue to Rotamer
        try:
            # The query for Residue should be specific enough.
            # PDB number is often not unique across chains in complex PDBs,
            # but within a single GPCR chain (preferred_chain), it usually is.
            # create_rotamers usually works on a specific protein_conformation.
            query_params_residue = {
                "protein_conformation": current_structure_obj.protein_conformation,
                "sequence_number": pdb_residue_number,
            }
            # If insertion_code is provided and not empty/null, add it to the query
            # The Residue model doesn't have an insertion_code field. This mapping
            # assumes that sequence_number is unique enough in the context of the
            # structure's protein_conformation. If PDBs have insertion codes that differentiate
            # residues with the same pdb_number, this needs more careful handling, potentially
            # at the PDB parsing/Residue creation stage in create_rotamers.
            # For now, we'll assume create_rotamers handles this disambiguation if necessary.

            receptor_residue_obj = Residue.objects.get(**query_params_residue)

            if receptor_residue_obj.amino_acid != residue_1_letter:
                logger.warning(
                    f"Residue AA mismatch for PDB {pdb_code_str}, StructID {current_structure_obj.id}, "
                    f"SeqNum {pdb_residue_number}, ChainFromYAML {chain_id_from_yaml} (Insertion: {insertion_code if insertion_code else 'N/A'}): "
                    f"DB has {receptor_residue_obj.amino_acid}, YAML has {residue_1_letter}. Skipping."
                )
                continue

            # Get Rotamer
            rotamer_obj = Rotamer.objects.get(
                structure=current_structure_obj, residue=receptor_residue_obj
            )

        except Residue.DoesNotExist:
            logger.warning(
                f"Residue not found in DB for {pdb_code_str} (StructID {current_structure_obj.id}), "
                f"SeqNum {pdb_residue_number}, ChainFromYAML {chain_id_from_yaml} (Insertion: {insertion_code if insertion_code else 'N/A'}). Skipping interaction."
            )
            continue
        except Rotamer.DoesNotExist:
            logger.warning(
                f"Rotamer not found for {pdb_code_str} (StructID {current_structure_obj.id}), "
                f"Residue {receptor_residue_obj} (SeqNum {pdb_residue_number}). Skipping interaction."
            )
            continue
        except (
            Residue.MultipleObjectsReturned,
            Rotamer.MultipleObjectsReturned,
        ) as e_multi:
            logger.error(
                f"Ambiguous mapping for residue/rotamer in {pdb_code_str}, SeqNum {pdb_residue_number}: {e_multi}. "
                "Skipping interaction."
            )
            continue

        # 2. Map interaction type
        feature = interaction_entry["feature"]
        feature_family = interaction_entry["feature_family"]

        interaction_type_slug = slugify(f"{feature_family}-{feature}")
        interaction_type_name = f"{feature_family} {feature}"
        interaction_category = feature_family  # e.g., "Acceptor", "Donor", "Aromatic"

        rfit_defaults = {
            "name": interaction_type_name,
            "type": interaction_category,
            # 'direction' can be inferred for some (e.g. Donor implies protein->ligand or vice-versa)
            # but let's keep it simple for now. It's nullable.
        }
        # If 'annotated_interaction' is provided and is a specific GPCRdb term, use it.
        # This field is 'null' in the YAML example.
        gpcrdb_specific_interaction = interaction_entry.get("annotated_interaction")
        if gpcrdb_specific_interaction:
            # Potentially override slug, name, type based on this.
            # Example: if it's 'polar_donor_protein', it gives more info.
            # This requires a mapping from Schrodinger's feature/family to GPCRdb's detailed types.
            # For now, we use feature/family.
            pass

        interaction_type_obj, created = (
            ResidueFragmentInteractionType.objects.get_or_create(
                slug=interaction_type_slug, defaults=rfit_defaults
            )
        )
        if created:
            logger.info(
                f"Created ResidueFragmentInteractionType: slug='{interaction_type_slug}'"
            )

        # 3. Get/Create Fragment object
        if not sli.pdb_file:
            logger.warning(
                f"SLI for PDB {pdb_code_str}, Ligand {ligand_pdb_het_code} "
                f"has no pdb_file. Cannot create Fragment. Skipping."
            )
            continue

        fragment_obj, created = Fragment.objects.get_or_create(
            ligand=current_ligand_db_obj,
            structure=current_structure_obj,
            residue=rotamer_obj.residue,  # Fragment is specific to the interacting receptor residue
            defaults={
                "pdbdata": sli.pdb_file
            },  # Use PDB data of the whole ligand in its bound pose
        )
        if created:
            logger.info(
                f"Created Fragment for ligand {current_ligand_db_obj.name} with residue {rotamer_obj.residue} in {pdb_code_str}"
            )

        # 4. Create ResidueFragmentInteraction
        rfi, created = ResidueFragmentInteraction.objects.get_or_create(
            structure_ligand_pair=sli,
            rotamer=rotamer_obj,
            fragment=fragment_obj,
            interaction_type=interaction_type_obj,
        )
        if created:
            interactions_created_count += 1

    logger.info(
        f"Successfully created {interactions_created_count} new SM interactions from {yaml_file_path} for {pdb_code_str}-{ligand_pdb_het_code}."
    )
    return True  # YAML found and processed.


# Placeholder for future peptide interaction processor
def process_schrodinger_peptide_interactions():
    # Similar logic, but YAML format and mapping might differ
    pass


# Placeholder for future protein-protein interaction processor
def process_schrodinger_protein_interactions():
    # Similar logic, but YAML format and mapping might differ
    pass
