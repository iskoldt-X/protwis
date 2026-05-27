# interaction/schrodinger_processor.py

from django.shortcuts import render, redirect
from django.http import HttpResponse, HttpResponseRedirect
from django.conf import settings
from django.db.models import Count, Sum, Avg, Q
from django.db import transaction
from django.utils.text import slugify
from django.conf import settings
import glob
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
from residue.models import Residue  # NB: Rotamer is in structure.models (imported above), not residue.models

import logging

logger = logging.getLogger(
    __name__
)  # Or 'gpcrdb.interaction.schrodinger' for more specific logging


def locate_interaction_yamls(base_dir, pdb_code, het_code):
    """Return the Engine 1 interaction YAML(s) for one (PDB, HET ligand).

    Engine 1 writes one YAML per ligand *instance* in a nested layout::

        {base}/{PDB}_{HET}/{PDB}_{HET}/{HET}_{chain}_{resnum}/{HET}_{chain}_{resnum}.yaml

    A single HET can occur multiple times in a PDB (different chain/resnum), so
    this returns a sorted list (possibly several instances). No DB access.
    """
    pdb = pdb_code.upper()
    het = het_code.upper()
    # Engine 1 doubles the {PDB}_{HET} directory, then one dir per ligand instance.
    instances_parent = os.path.join(base_dir, f"{pdb}_{het}", f"{pdb}_{het}")
    found = []
    for inst_dir in sorted(glob.glob(os.path.join(instances_parent, f"{het}_*"))):
        if not os.path.isdir(inst_dir):
            continue
        # The YAML inside an instance dir shares the dir's name: {HET}_{chain}_{resnum}.yaml
        candidate = os.path.join(inst_dir, os.path.basename(inst_dir) + ".yaml")
        if os.path.exists(candidate):
            found.append(candidate)
    return found


# The 20 standard amino acids (1-letter). Engine 1 writes "X" for anything
# non-standard (water, ions, glycans, covalent modifications fused into the
# receptor — Q-S1 §3.c / trap H5); those have no protwis Residue row.
STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def is_standard_residue(name_1_letter):
    """True iff ``name_1_letter`` is one of the 20 standard amino acids.

    Engine 1 marks non-standard receptor residues as ``"X"`` — they must be
    skipped (no matching protwis Residue) rather than queried (Step A5).
    """
    return name_1_letter.upper() in STANDARD_AMINO_ACIDS


_INTERACTION_TYPE_MAP = None  # lazily loaded {(family, direction): slug}


def _load_interaction_type_map():
    """Load and cache the Plan C (family, direction) -> slug routing table.

    The map lives next to this module (``interaction/interaction_type_map.yaml``);
    only its ``rules`` are loaded — the ``deferred`` section documents things the
    consumer (this module) handles itself (polar_backbone, dedup, dropped classes).
    """
    global _INTERACTION_TYPE_MAP
    if _INTERACTION_TYPE_MAP is None:
        map_path = os.path.join(os.path.dirname(__file__), "interaction_type_map.yaml")
        with open(map_path) as fh:
            doc = yaml.safe_load(fh)
        _INTERACTION_TYPE_MAP = {
            (rule["family"], rule["direction"]): rule["slug"]
            for rule in (doc.get("rules") or [])
        }
    return _INTERACTION_TYPE_MAP


def resolve_slug(feature_family, direction):
    """Map an Engine 1 ``(feature_family, direction)`` pair to a canonical protwis slug.

    Routing is by ``(family, direction)`` only — the chemical ``feature`` name is
    detail that does not change the slug (Plan C contract). Unknown combinations
    ``raise ValueError`` (fail-loud) so that non-default contacts
    (halogen/metal/water) surface immediately rather than silently dropping data.
    """
    lookup = _load_interaction_type_map()
    try:
        return lookup[(feature_family, direction)]
    except KeyError:
        raise ValueError(
            f"No interaction_type_map rule for (family={feature_family!r}, "
            f"direction={direction!r}). Known combos: {sorted(lookup)}"
        )


def passes_chain_filter(chain_id, preferred_chain):
    """True iff ``chain_id`` belongs to the structure's preferred receptor chain.

    protwis stores no per-residue chain; ``build_structures.py:1410`` drops any
    residue whose chain != ``Structure.preferred_chain``. Multi-chain receptors
    are reduced to their first chain (build_structures.py:248-249), which we
    replicate here. An empty preference keeps everything (Step A7).
    """
    if not preferred_chain:
        return True
    preferred_first = preferred_chain.split(",")[0]
    return chain_id == preferred_first


def get_receptor_pdb_block(interaction_entry):
    """Return the receptor-side PDB ATOM text for one interaction entry.

    Engine 1 attaches the interacting receptor residue's atoms as a PDB block.
    NB (trap H3 / §0.5 C4): this block uses a non-standard column layout (resname
    in col 17, alt-loc column swallowed) — store it as opaque text and NEVER feed
    it to BioPython (Step A8).
    """
    return interaction_entry.get("receptor_pdb_block") or ""


def parse_receptor_residue(receptor_res_info):
    """Pure (no-DB) parse of a YAML ``receptor_residue`` block into the fields
    the processor needs to locate the matching protwis Residue/Rotamer.

    Engine 1 YAML schema for ``receptor_residue``::

        name_1_letter: N
        pdb_residue_number: 406
        chain_id: A
        insertion_code: ''

    Unit-testable in isolation (Step A1).
    """
    return {
        "sequence_number": int(receptor_res_info["pdb_residue_number"]),
        "amino_acid": receptor_res_info["name_1_letter"].upper(),
        "chain_id": receptor_res_info["chain_id"],
        "insertion_code": receptor_res_info.get("insertion_code", "") or "",
    }


@transaction.atomic
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

    # Engine 1 nested layout — there may be several ligand instances per (PDB, HET).
    yaml_file_paths = locate_interaction_yamls(
        schrodinger_interaction_dir, pdb_code_str, ligand_pdb_het_code
    )
    if not yaml_file_paths:
        logger.warning(
            f"No Schrödinger SM YAML found for {pdb_code_str}_{ligand_pdb_het_code} "
            f"under {schrodinger_interaction_dir}"
        )
        return False  # File not found, indicates fallback is needed

    all_interactions = []
    for yaml_file_path in yaml_file_paths:
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
            continue  # This instance has nothing; keep scanning the others.

        all_interactions.extend(schrodinger_data["result"]["interactions"])

    if not all_interactions:
        # YAML(s) found but no usable interactions — return without touching the DB
        # so we never wipe existing rows on an empty/invalid file.
        logger.warning(
            f"No usable interactions in located YAML(s) for "
            f"{pdb_code_str}_{ligand_pdb_het_code}; leaving existing rows untouched."
        )
        return True

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

    # A9 idempotency (Q-C7): replace this ligand-pair's interactions. Because the
    # whole function is @transaction.atomic, a mid-write error rolls back BOTH this
    # delete and any partial writes — the original rows survive a failed run.
    ResidueFragmentInteraction.objects.filter(structure_ligand_pair=sli).delete()

    interactions_created_count = 0
    for interaction_entry in all_interactions:
        receptor_res_info = interaction_entry["receptor_residue"]
        parsed_residue = parse_receptor_residue(receptor_res_info)
        pdb_residue_number = parsed_residue["sequence_number"]
        chain_id_from_yaml = parsed_residue["chain_id"]  # Keep for logging/debugging
        residue_1_letter = parsed_residue["amino_acid"]
        insertion_code = parsed_residue["insertion_code"]

        # A5: skip non-standard residues (Engine 1 marks them "X") before any DB
        # query — they have no protwis Residue row.
        if not is_standard_residue(residue_1_letter):
            logger.warning(
                f"Skipping non-standard receptor residue '{residue_1_letter}' in "
                f"{pdb_code_str} at seq {pdb_residue_number} chain {chain_id_from_yaml}."
            )
            continue

        # A7: keep only the preferred receptor chain (protwis has no per-residue
        # chain; mirrors build_structures.py:1410).
        if not passes_chain_filter(
            chain_id_from_yaml, current_structure_obj.preferred_chain
        ):
            logger.info(
                f"Skipping non-preferred chain '{chain_id_from_yaml}' in {pdb_code_str} "
                f"at seq {pdb_residue_number} (preferred={current_structure_obj.preferred_chain})."
            )
            continue

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

        # 2. Map interaction type via the Plan C map: (family, direction) -> canonical slug.
        # The `feature` name is chemical detail and intentionally NOT used for routing
        # (so a comma-joined feature string can never corrupt the slug — trap H2 / A4).
        feature_family = interaction_entry["feature_family"]
        direction = interaction_entry.get("direction")
        interaction_type_slug = resolve_slug(feature_family, direction)  # fail-loud on unknown

        # A6: use .get() (NOT get_or_create) — the 18 canonical slugs are pre-seeded
        # via the interaction_types.json fixture (Plan B); never mint new slugs here.
        try:
            interaction_type_obj = ResidueFragmentInteractionType.objects.get(
                slug=interaction_type_slug
            )
        except ResidueFragmentInteractionType.DoesNotExist:
            logger.error(
                f"Canonical slug '{interaction_type_slug}' not in DB — "
                "run `loaddata interaction_types.json` (Plan B fixture) before importing."
            )
            raise

        # 3. Get/Create Fragment, backed by the YAML's receptor_pdb_block (A8 / §0.5 B4).
        # The old code keyed Fragment.pdbdata off the always-None SLI pdb file and
        # skipped creation when it was missing, so Fragments were never created. We
        # store the block as opaque text in a PdbData (trap H3: non-standard columns —
        # never parse with BioPython). One Fragment per (ligand, structure, residue),
        # reused across multiple interactions on the same residue.
        fragment_obj = Fragment.objects.filter(
            ligand=current_ligand_db_obj,
            structure=current_structure_obj,
            residue=rotamer_obj.residue,
        ).first()
        if fragment_obj is None:
            pdb_data_obj = PdbData.objects.create(
                pdb=get_receptor_pdb_block(interaction_entry)
            )
            fragment_obj = Fragment.objects.create(
                ligand=current_ligand_db_obj,
                structure=current_structure_obj,
                residue=rotamer_obj.residue,
                pdbdata=pdb_data_obj,
            )
            logger.info(
                f"Created Fragment for ligand {current_ligand_db_obj.name} "
                f"with residue {rotamer_obj.residue} in {pdb_code_str}"
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
