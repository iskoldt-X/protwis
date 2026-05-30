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

    Engine 1 writes one YAML per ligand *instance*, but the depth varies between
    exports. Two layouts are seen in the wild and both are supported here:

    * **nested** (6LN2)::

        {base}/{PDB}_{HET}/{PDB}_{HET}/{HET}_{chain}_{resnum}/{HET}_{chain}_{resnum}.yaml

    * **shallow** (2Y02), one level fewer — instance dirs sit directly under the
      top ``{PDB}_{HET}`` dir::

        {base}/{PDB}_{HET}/{HET}_{chain}_{resnum}/{HET}_{chain}_{resnum}.yaml

    A single HET can occur multiple times in a PDB (different chain/resnum), so
    this returns a sorted list (possibly several instances). The non-preferred
    chains are dropped later by ``passes_chain_filter``. No DB access.
    """
    pdb = pdb_code.upper()
    het = het_code.upper()
    top = os.path.join(base_dir, f"{pdb}_{het}")
    # Look for instance dirs ({HET}_*) under BOTH the doubled (nested) parent and
    # the top (shallow) parent. The doubled dir is named {PDB}_{HET} and the
    # instance dirs are named {HET}_* — they never collide, so each layout only
    # matches its own instances.
    instance_parents = [os.path.join(top, f"{pdb}_{het}"), top]
    found = []
    seen = set()
    for parent in instance_parents:
        for inst_dir in sorted(glob.glob(os.path.join(parent, f"{het}_*"))):
            if not os.path.isdir(inst_dir):
                continue
            # The YAML inside an instance dir shares the dir's name: {HET}_{chain}_{resnum}.yaml
            candidate = os.path.join(inst_dir, os.path.basename(inst_dir) + ".yaml")
            if os.path.exists(candidate) and candidate not in seen:
                seen.add(candidate)
                found.append(candidate)
    return sorted(found)


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


# A6b (ADR-005) — main-chain N/O H-bonds reclassify to ``polar_backbone``.
#
# views.py:1008-1019 puts backbone first in the slug elif-chain (it wins over
# donor/acceptor). resolve_slug above can't know about it: backbone vs sidechain
# is a property of which receptor atom the H-bond lands on, and the (family,
# direction) pair alone doesn't carry that. Pre-ADR-005 the wrapper YAML didn't
# expose the receptor atom name either — the only signal was the verbose
# receptor_pdb_block, which would have forced fragile re-parsing here.
#
# ADR-005 closed that gap: the wrapper now emits ``receptor_atom_name`` per
# interaction (a `.pdbname.strip()` off the actual Schrödinger atom that
# resolve_slug already pointed at). With that field present, A6b reduces to a
# tiny post-routing override applied at each resolve_slug call site.
_BACKBONE_PROMOTABLE_SLUGS = frozenset({
    "polar_donor_protein",     # protein donates H from main-chain N
    "polar_acceptor_protein",  # protein accepts H at main-chain carbonyl O
})
_BACKBONE_ATOM_NAMES = frozenset({"N", "O"})


def apply_backbone_override(slug, receptor_atom_name):
    """Promote a polar H-bond slug to ``polar_backbone`` when the protein
    partner is a main-chain N or O. No-op for any other slug or atom.

    ``receptor_atom_name`` may be ``None`` or empty (old YAML fixtures without
    the ADR-005 field) — those silently pass through unchanged, so adding A6b
    cannot break replays of older inputs.
    """
    if slug not in _BACKBONE_PROMOTABLE_SLUGS:
        return slug
    name = (receptor_atom_name or "").strip()
    if name in _BACKBONE_ATOM_NAMES:
        return "polar_backbone"
    return slug


# ---------------------------------------------------------------------------
# A6c — de-duplication only (ADR-008 / ADR-009; reverts the earlier 2c-1
# "charge suppresses h-bond" branch).
#
# protwis views.py:998/1002 forces ``hydrogenmatch = False`` on charged residues
# so that a charge-assisted h-bond is recorded as ONE charge row, not two. The
# first A6c implementation mirrored that at the slug level. ADR-008 (2026-05-28)
# reverses the call: charge-assisted h-bonds carry two independent physical
# forces (Coulomb + directional N-H...O dipole) and Schrödinger's two detectors
# fire independently — keeping both rows is a faithful record, not double
# counting. ADR-009 generalises this: the DB layer preserves information; UI
# compactness is solved by front-end views, not by dropping rows here.
#
# What remains in this pass is the only loss-free reduction: collapsing
# identical ``(sequence_number, slug)`` rows (e.g. 6LN2 N407's two Acceptor
# entries that both map to ``polar_donor_protein``, or 2Y02 N2's twin Donor
# entries to the same residue). Aromatic slugs are unaffected.
# ---------------------------------------------------------------------------


def apply_priority_dedup(records):
    """De-duplicate resolved (residue, slug) interactions; no suppression.

    ``records`` is a list of dicts, each with at least ``sequence_number`` (int),
    ``amino_acid`` (1-letter str) and ``slug`` (canonical slug str). Returns the
    surviving records as a list (input order preserved), deduped on
    ``(sequence_number, slug)``.

    Semantics (ADR-008 / ADR-009): identical ``(sequence_number, slug)`` rows
    collapse to one. Nothing else is dropped — distinct slugs on the same
    residue (e.g. a charge slug + an h-bond slug on a charge-assisted contact)
    are all preserved as faithful records of independent Schrödinger detectors.
    """
    survivors = []
    seen = set()
    for r in records:
        key = (r["sequence_number"], r["slug"])
        if key in seen:
            continue
        seen.add(key)
        survivors.append(r)
    return survivors


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

    # A6c (ADR-008 / ADR-009): compute the surviving (seq, slug) set up front
    # (pure / no DB) to collapse identical-(seq, slug) duplicates. No
    # suppression — distinct slugs on the same residue (e.g. charge + h-bond on
    # a charge-assisted contact) all survive as independent Schrödinger
    # detector readings. resolve_slug stays fail-loud here too (unknown combos
    # raise before any write).
    prelim_records = []
    for interaction_entry in all_interactions:
        parsed = parse_receptor_residue(interaction_entry["receptor_residue"])
        if not is_standard_residue(parsed["amino_acid"]):
            continue
        if not passes_chain_filter(
            parsed["chain_id"], current_structure_obj.preferred_chain
        ):
            continue
        prelim_slug = resolve_slug(
            interaction_entry["feature_family"],
            interaction_entry.get("direction"),
        )
        prelim_slug = apply_backbone_override(
            prelim_slug, interaction_entry.get("receptor_atom_name")
        )
        prelim_records.append(
            {
                "sequence_number": parsed["sequence_number"],
                "amino_acid": parsed["amino_acid"],
                "slug": prelim_slug,
            }
        )
    survivor_keys = {
        (r["sequence_number"], r["slug"]) for r in apply_priority_dedup(prelim_records)
    }

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
        # A6b (ADR-005): promote main-chain N/O H-bonds to polar_backbone before
        # the dedup keys are checked — survivor_keys above was computed with the
        # same override, so the two passes stay in lockstep.
        feature_family = interaction_entry["feature_family"]
        direction = interaction_entry.get("direction")
        interaction_type_slug = resolve_slug(feature_family, direction)  # fail-loud on unknown
        interaction_type_slug = apply_backbone_override(
            interaction_type_slug, interaction_entry.get("receptor_atom_name")
        )

        # A6c (ADR-008 / ADR-009): drop entries the dedup pass collapsed —
        # only identical (seq, slug) duplicates are removed; distinct slugs on
        # the same residue all survive. survivor_keys was computed up front.
        if (pdb_residue_number, interaction_type_slug) not in survivor_keys:
            logger.info(
                "A6c dedup: duplicate (seq, slug) (%s on %s residue %s).",
                interaction_type_slug,
                pdb_code_str,
                pdb_residue_number,
            )
            continue

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


# ===========================================================================
# Engine 2 (peptide / protein-protein) consumer — D3.
#
# This is the THIRD Strategy implementation. Engine 1 (small molecule) above
# writes RFI rows (residue ↔ ligand-fragment). Engine 2's chemistry is a
# residue-residue interface (receptor residue ↔ peptide/protein-partner
# residue), so it consumes a *different* schema (ADR-014, engine2/1.0) and
# writes to *different* models.
#
# TARGET MODELS (investigation finding — D3): protwis already ships
# residue-residue interface tables in the contactnetwork app:
#
#   contactnetwork.InteractingPeptideResiduePair
#       receptor_residue (FK residue.Residue) ↔ a peptide residue stored as
#       three loose fields (peptide_amino_acid / _three_letter /
#       _sequence_number) + FK ligand.LigandPeptideStructure (the peptide
#       chain). The peptide residue is NOT a protwis Residue row (peptides are
#       not in the Residue table), which matches Engine 2's data: only the
#       receptor side resolves to a Residue/generic number.
#
#   contactnetwork.InteractionPeptide
#       per-pair interaction detail rows: peptide_atom / receptor_atom /
#       interaction_type / specific_type / interaction_level.
#
# These are populated today by the *legacy* BioPython path
# (contactnetwork/cube.py compute_interactions → InteractingPair
# .save_peptide_interactions, run by build_crystal_interactions). The API
# (StructurePeptideLigandInteractionSerializer) + front-end already read them.
# So Engine 2 reuses the same DB contract; no new model/migration for the
# *core* interface interactions.
#
# ADR-014 / ADR-009 — interface metrics (RESOLVED by ADR-015, 2026-05-30):
#   Engine 2 produces buried-SASA, surface complementarity (Sc) and the
#   backbone/sidechain hbond 4-class — protein-protein interface metrics the
#   legacy RDKit path never had. ADR-009 says the DB layer must preserve
#   information. This consumer therefore:
#     * writes the core interactions to the existing models (fail-loud),
#     * preserves the hbond 4-class + face/edge subtype in ``specific_type``
#       (a free-text field, lossless for that signal),
#     * (ADR-015) writes per-residue-pair buried-SASA + Sc into three NULLABLE
#       fields added to InteractingPeptideResiduePair (buried_sasa_receptor /
#       buried_sasa_peptide / surface_complementarity). These are additive: the
#       legacy BioPython path leaves them NULL (it never computed them). The
#       per-pair grain matches the engine2/1.0 residue_pair_summaries schema.
# ===========================================================================


# Engine 2 interaction-type vocabulary, mapped to the legacy InteractionPeptide
# vocabulary the API + front-end already understand (contactnetwork/interaction.py
# CI subclasses). The legacy interaction_type values are: van-der-waals /
# hydrophobic / ionic / polar / aromatic. specific_type carries the finer detail.
#
# Engine 2 type        -> (legacy interaction_type, specific_type seed)
_ENGINE2_TYPE_MAP = {
    "hydrogen_bond": ("polar", "h-bond"),
    "salt_bridge": ("ionic", "salt-bridge"),
    "pi_pi_stacking": ("aromatic", "pi-stacking"),
    "pi_cation": ("aromatic", "pi-cation"),
    "hydrophobic_contact": ("hydrophobic", ""),
    "steric_clash": ("steric-clash", "clash"),
}

# Engine 2 instance YAMLs live at:
#   {base}/{PDB}/{PDB}_{recvChain}_{ligChain}/{PDB}_{recvChain}_{ligChain}.yaml
# (the adapter-emitted, ADR-014 schema file — NOT the *_worker.yaml sidecar,
#  which is the raw worker schema). Mirrors the D4 pioneer layout in /tmp.


def locate_engine2_yamls(base_dir, pdb_code):
    """Return the Engine 2 (engine2/1.0) instance YAML(s) for one PDB.

    Layout (D4 pioneer, schema_adapter output)::

        {base}/{PDB}/{PDB}_{recv}_{lig}/{PDB}_{recv}_{lig}.yaml

    A single PDB can have several interface instances (multiple peptide chains,
    or peptide + protein-protein), so this returns a sorted list. The raw
    ``*_worker.yaml`` sidecar (worker schema, not ADR-014) is excluded. No DB
    access.
    """
    pdb = pdb_code.upper()
    top = os.path.join(base_dir, pdb)
    found = []
    for inst_dir in sorted(glob.glob(os.path.join(top, f"{pdb}_*"))):
        if not os.path.isdir(inst_dir):
            continue
        candidate = os.path.join(inst_dir, os.path.basename(inst_dir) + ".yaml")
        if os.path.exists(candidate):
            found.append(candidate)
    return sorted(found)


def map_engine2_interaction(interaction_entry):
    """Map one Engine 2 ``interface_interactions`` entry to the legacy
    ``(interaction_type, specific_type, interaction_level)`` triple plus the
    receptor/peptide atom names.

    Returns a dict (pure, no DB)::

        {interaction_type, specific_type, interaction_level,
         receptor_atom, peptide_atom}

    fail-loud (ADR-009): an unknown Engine 2 ``type`` raises ValueError so a
    new interaction class surfaces immediately instead of silently dropping.

    ``specific_type`` is enriched losslessly with the Engine 2-only signal that
    the legacy models have no column for — the hbond 4-class (hb_bb/bs/sb/ss)
    for H-bonds, and the pi subtype (Face-to-Face / Edge-to-Face) for stacking.
    interaction_level is always 0 (Engine 2 emits a single strict definition,
    unlike the legacy strict/loose split).
    """
    itype = interaction_entry.get("type")
    if itype not in _ENGINE2_TYPE_MAP:
        raise ValueError(
            f"Unknown Engine 2 interaction type {itype!r}. "
            f"Known: {sorted(_ENGINE2_TYPE_MAP)}"
        )
    legacy_type, specific = _ENGINE2_TYPE_MAP[itype]

    receptor_atom = ""
    peptide_atom = ""

    if itype == "hydrogen_bond":
        # Preserve the backbone/sidechain 4-class (ADR-014 interface signal).
        hbond_class = interaction_entry.get("hbond_class")
        if hbond_class:
            specific = f"{specific}:{hbond_class}"
        donor = interaction_entry.get("donor") or {}
        acceptor = interaction_entry.get("acceptor") or {}
        receptor_atom, peptide_atom = _engine2_atoms_by_side(donor, acceptor)
    elif itype == "salt_bridge":
        anion = interaction_entry.get("anion") or {}
        cation = interaction_entry.get("cation") or {}
        receptor_atom, peptide_atom = _engine2_atoms_by_side(anion, cation)
    elif itype == "pi_pi_stacking":
        subtype = interaction_entry.get("subtype")
        if subtype:
            specific = f"{specific}:{subtype}"
        # ring-ring: residue-level, no single atom
        r1 = interaction_entry.get("residue1") or {}
        r2 = interaction_entry.get("residue2") or {}
        receptor_atom, peptide_atom = _engine2_atoms_by_side(r1, r2)
    elif itype == "pi_cation":
        cat = interaction_entry.get("cation_residue") or {}
        pi = interaction_entry.get("pi_residue") or {}
        receptor_atom, peptide_atom = _engine2_atoms_by_side(cat, pi)
    else:  # hydrophobic_contact / steric_clash — atom1/atom2
        a1 = interaction_entry.get("atom1") or {}
        a2 = interaction_entry.get("atom2") or {}
        receptor_atom, peptide_atom = _engine2_atoms_by_side(a1, a2)

    return {
        "interaction_type": legacy_type,
        "specific_type": specific,
        "interaction_level": 0,
        "receptor_atom": receptor_atom,
        "peptide_atom": peptide_atom,
    }


def _engine2_atoms_by_side(part_a, part_b):
    """Given two interaction partners (each a dict with a ``side`` field of
    ``selection1`` = receptor / ``selection2`` = peptide), return
    ``(receptor_atom_name, peptide_atom_name)`` regardless of which partner is
    which. ``atom_name`` may be absent (ring/residue-level) → empty string.

    Engine 2 always tags every partner with ``side`` (D2 adapter contract). If
    a partner is missing its side tag, we fall back to the textual order
    (a=receptor, b=peptide) but the per-pair side resolution in the processor
    re-derives the residues from ``side`` anyway, so this is only for atom names.
    """
    def name(d):
        return (d.get("atom_name") or "").strip()

    if part_a.get("side") == "selection2" or part_b.get("side") == "selection1":
        # a is peptide, b is receptor
        return name(part_b), name(part_a)
    return name(part_a), name(part_b)


def _engine2_partner_side(part):
    """``selection1`` -> 'receptor', ``selection2`` -> 'peptide', else None."""
    side = part.get("side")
    if side == "selection1":
        return "receptor"
    if side == "selection2":
        return "peptide"
    return None


def _engine2_pair_partners(interaction_entry):
    """Return ``(receptor_partner_dict, peptide_partner_dict)`` for one Engine 2
    entry, using each partner's ``side`` tag (selection1=receptor,
    selection2=peptide). Returns ``(None, None)`` if the two partners are not on
    opposite sides (intra-chain pair — out of interface scope, skipped by the
    caller).
    """
    itype = interaction_entry.get("type")
    pair_keys = {
        "hydrogen_bond": ("donor", "acceptor"),
        "salt_bridge": ("anion", "cation"),
        "pi_pi_stacking": ("residue1", "residue2"),
        "pi_cation": ("cation_residue", "pi_residue"),
        "hydrophobic_contact": ("atom1", "atom2"),
        "steric_clash": ("atom1", "atom2"),
    }[itype]
    p1 = interaction_entry.get(pair_keys[0]) or {}
    p2 = interaction_entry.get(pair_keys[1]) or {}

    s1, s2 = _engine2_partner_side(p1), _engine2_partner_side(p2)
    if s1 == "receptor" and s2 == "peptide":
        return p1, p2
    if s1 == "peptide" and s2 == "receptor":
        return p2, p1
    return None, None


def build_residue_pair_metrics(residue_pair_summaries):
    """Index ``residue_pair_summaries`` by ``(recv_seq, pep_seq)`` -> interface
    metrics dict (ADR-015).

    The engine2/1.0 schema reports buried-SASA + surface complementarity (Sc)
    once per residue pair (under ``residue_pair_summaries[].properties``), NOT
    per individual interaction — so these are stored on
    ``InteractingPeptideResiduePair`` (the pair-level model), not on the
    per-interaction ``InteractionPeptide`` rows.

    Each summary entry tags its two residues with ``side`` (selection1=receptor,
    selection2=peptide), and ``properties.set_1_buried_sasa`` / ``set_2_buried_sasa``
    correspond to selection1 / selection2 respectively (worker contract). We
    resolve receptor vs peptide by ``side`` (not positional order) so the mapping
    holds even if residue1/residue2 ordering ever flips.

    Returns ``{(recv_seq, pep_seq): {"buried_sasa_receptor", "buried_sasa_peptide",
    "surface_complementarity"}}``. Pure (no DB). Missing metric keys map to None
    (nullable fields — additive, ADR-010/015). Entries whose two residues are not
    on opposite sides are skipped (out of interface scope).
    """
    metrics = {}
    for summary in residue_pair_summaries or []:
        r1 = summary.get("residue1") or {}
        r2 = summary.get("residue2") or {}
        s1, s2 = _engine2_partner_side(r1), _engine2_partner_side(r2)
        # set_1 == selection1, set_2 == selection2 (worker contract).
        props = summary.get("properties") or {}
        set1 = props.get("set_1_buried_sasa")
        set2 = props.get("set_2_buried_sasa")
        if s1 == "receptor" and s2 == "peptide":
            recv, pep = r1, r2
            buried_receptor, buried_peptide = set1, set2
        elif s1 == "peptide" and s2 == "receptor":
            recv, pep = r2, r1
            # residue1 is the peptide here, so set_1 is the peptide-side SASA.
            buried_receptor, buried_peptide = set2, set1
        else:
            # intra-chain / untagged — not an interface pair.
            continue
        key = (int(recv["resid"]), int(pep["resid"]))
        metrics[key] = {
            "buried_sasa_receptor": buried_receptor,
            "buried_sasa_peptide": buried_peptide,
            "surface_complementarity": props.get("surface_complementarity"),
        }
    return metrics


@transaction.atomic
def process_schrodinger_peptide_interactions(
    current_structure_obj: Structure,
    ligand_chain: str,
    pdb_code_str: str,
    schrodinger_interactions_dir_override: str = None,
) -> bool:
    """Process Engine 2 (engine2/1.0) peptide interface interactions into the
    contactnetwork InteractingPeptideResiduePair / InteractionPeptide models.

    Mirrors the Engine 1 small-molecule processor's good habits:
      * fail-loud on unknown interaction types / missing receptor Residue;
      * per-PDB (per ligand_chain) boundary;
      * @transaction.atomic delete-then-insert idempotency — a mid-write error
        rolls back this chain's rows, leaving the prior state intact.

    Returns True if a matching YAML was found and processed (even if it yielded
    no interactions), False if no YAML was found or a parse error occurred, so
    the caller can decide on fallback (same contract as the SM processor).

    ``ligand_chain`` selects which LigandPeptideStructure / interface instance
    to write (a PDB can have several peptide chains). When several Engine 2
    YAML instances exist for this PDB, only the one whose ligand_chain matches
    is consumed.
    """
    logger.info(
        f"Processing Schrödinger Engine 2 peptide interactions for {pdb_code_str} "
        f"chain {ligand_chain}"
    )

    base_dir = schrodinger_interactions_dir_override or getattr(
        settings, "SCHRODINGER_INTERACTIONS_DIR", None
    )
    if not base_dir:
        logger.error(
            "SCHRODINGER_INTERACTIONS_DIR is not set. Cannot process Engine 2 "
            "peptide interactions."
        )
        return False

    yaml_paths = locate_engine2_yamls(base_dir, pdb_code_str)
    if not yaml_paths:
        logger.warning(
            f"No Engine 2 YAML found for {pdb_code_str} under {base_dir}"
        )
        return False

    # Locate the LigandPeptideStructure for this (structure, chain). build_structures
    # creates it at the ligand-type decision point for type∈{peptide,protein}.
    from ligand.models import LigandPeptideStructure
    from contactnetwork.models import (
        InteractingPeptideResiduePair,
        InteractionPeptide,
    )

    lps_qs = LigandPeptideStructure.objects.filter(
        structure=current_structure_obj, chain=ligand_chain
    )
    lps = lps_qs.first()
    if lps is None:
        logger.error(
            f"No LigandPeptideStructure for {pdb_code_str} chain {ligand_chain}; "
            "build_structures must create it before importing Engine 2 interactions."
        )
        return False

    # Pick the YAML instance whose metadata.ligand_chain matches this chain.
    chosen = None
    for path in yaml_paths:
        try:
            with open(path) as fh:
                doc = yaml.safe_load(fh)
        except yaml.YAMLError as e:
            logger.error(f"Error parsing Engine 2 YAML {path}: {e}")
            return False
        if not doc or doc.get("schema_version") != "engine2/1.0":
            logger.warning(
                f"Skipping non-engine2/1.0 YAML {path} "
                f"(schema_version={doc.get('schema_version') if doc else None})."
            )
            continue
        meta = doc.get("metadata") or {}
        if (meta.get("ligand_chain") or "") == ligand_chain:
            chosen = doc
            break

    if chosen is None:
        logger.warning(
            f"No engine2/1.0 instance for {pdb_code_str} chain {ligand_chain} "
            f"among {len(yaml_paths)} YAML(s)."
        )
        return False

    interactions = chosen.get("interface_interactions") or []

    # ADR-015: per-residue-pair interface metrics (buried-SASA / Sc). The
    # engine2/1.0 schema carries these once per residue pair under
    # residue_pair_summaries, not per interaction — so they attach to the
    # pair-level model (InteractingPeptideResiduePair), populated when the pair
    # is first created below.
    pair_metrics = build_residue_pair_metrics(chosen.get("residue_pair_summaries"))

    # Idempotency: drop any prior peptide interface rows for THIS peptide
    # instance before re-inserting. CASCADE on the pair removes its
    # InteractionPeptide children.
    InteractingPeptideResiduePair.objects.filter(peptide=lps).delete()

    if not interactions:
        logger.warning(
            f"Engine 2 YAML for {pdb_code_str} chain {ligand_chain} has no "
            "interface_interactions; leaving zero rows (delete already ran)."
        )
        return True

    # Group InteractionPeptide rows per (receptor_seq, peptide_seq) pair.
    pairs_cache = {}  # (recv_seq, pep_seq) -> InteractingPeptideResiduePair
    created_count = 0

    for entry in interactions:
        recv_partner, pep_partner = _engine2_pair_partners(entry)
        if recv_partner is None or pep_partner is None:
            # Not an interface (receptor↔peptide) pair — skip (intra-chain or
            # untagged). The interface worker should not emit these, but be safe.
            logger.info(
                f"Skipping non-interface Engine 2 entry (type={entry.get('type')}) "
                f"in {pdb_code_str} chain {ligand_chain}."
            )
            continue

        recv_seq = int(recv_partner["resid"])
        pep_seq = int(pep_partner["resid"])
        recv_resname = (recv_partner.get("resname") or "").upper()
        pep_resname = (pep_partner.get("resname") or "").upper()

        mapped = map_engine2_interaction(entry)  # fail-loud on unknown type

        # Resolve the receptor Residue (same lookup as the SM processor).
        try:
            receptor_residue = Residue.objects.get(
                protein_conformation=current_structure_obj.protein_conformation,
                sequence_number=recv_seq,
            )
        except Residue.DoesNotExist:
            logger.warning(
                f"Receptor Residue not found for {pdb_code_str} seq {recv_seq}; "
                "skipping Engine 2 interaction."
            )
            continue
        except Residue.MultipleObjectsReturned:
            logger.error(
                f"Ambiguous receptor Residue for {pdb_code_str} seq {recv_seq}; "
                "skipping Engine 2 interaction."
            )
            continue

        # AA-guard (Gate-2, ADR-018 / fix #1): the receptor join above keys
        # purely on author seq-number, so a -N register offset between the PDB
        # author numbering and GPCRdb's Residue numbering would silently map
        # this interaction onto the WRONG receptor residue. The Engine 2 YAML
        # already carries the receptor residue's three-letter name, so we
        # cross-check it against the DB residue's amino acid (the same guard the
        # SM processor enforces at its receptor lookup, :448). _three_to_one
        # normalises PrepWizard protonation aliases (HIE/HID/HIP/ASH/…) so they
        # are not flagged as false mismatches. On disagreement we warn + skip
        # rather than write a mis-registered row (北极星: do not silently
        # corrupt data — surface it).
        recv_one_letter = _three_to_one(recv_resname)
        if recv_resname and recv_one_letter != receptor_residue.amino_acid:
            logger.warning(
                f"Receptor AA mismatch for {pdb_code_str} seq {recv_seq}: "
                f"DB={receptor_residue.amino_acid} YAML={recv_resname} "
                f"({recv_one_letter}); likely author-seqnum register offset. "
                "Skipping Engine 2 interaction (no mis-registered row written)."
            )
            continue

        pep_one_letter = _three_to_one(pep_resname)

        key = (recv_seq, pep_seq)
        pair = pairs_cache.get(key)
        if pair is None:
            # ADR-015: attach the pair-level interface metrics if the
            # residue_pair_summaries listed this pair. Pairs absent from the
            # summaries (or rows from the legacy BioPython path) keep NULL —
            # the fields are nullable / additive.
            metrics = pair_metrics.get(key, {})
            pair = InteractingPeptideResiduePair.objects.create(
                peptide_amino_acid_three_letter=pep_resname[:3],
                peptide_amino_acid=pep_one_letter,
                peptide_sequence_number=pep_seq,
                peptide=lps,
                receptor_residue=receptor_residue,
                buried_sasa_receptor=metrics.get("buried_sasa_receptor"),
                buried_sasa_peptide=metrics.get("buried_sasa_peptide"),
                surface_complementarity=metrics.get("surface_complementarity"),
            )
            pairs_cache[key] = pair

        InteractionPeptide.objects.create(
            interacting_peptide_pair=pair,
            peptide_atom=mapped["peptide_atom"],
            receptor_atom=mapped["receptor_atom"],
            interaction_type=mapped["interaction_type"],
            specific_type=mapped["specific_type"],
            interaction_level=mapped["interaction_level"],
        )
        created_count += 1

    logger.info(
        f"Created {created_count} Engine 2 peptide interaction rows across "
        f"{len(pairs_cache)} residue pairs for {pdb_code_str} chain {ligand_chain}."
    )
    return True


# A protein-protein interface (G protein / arrestin / receptor dimer) is the
# SAME residue-residue interface shape as a peptide. ADR-014 keeps it on the
# Engine 2 schema. The contactnetwork peptide models are keyed off
# LigandPeptideStructure, which build_structures ALSO creates for
# ligand['type'] == 'protein' (build_structures.py:1791) — so the peptide
# processor handles both. This thin alias documents the Strategy slot and keeps
# the two-tier dispatch explicit; if protein-protein later needs a distinct
# target model, it diverges here.
def process_schrodinger_protein_interactions(
    current_structure_obj: Structure,
    ligand_chain: str,
    pdb_code_str: str,
    schrodinger_interactions_dir_override: str = None,
) -> bool:
    """Process Engine 2 protein-protein interface interactions.

    Currently identical to the peptide path: both are receptor-residue ↔
    partner-residue interfaces stored on the contactnetwork peptide models
    (build_structures creates a LigandPeptideStructure for type == 'protein'
    too). Kept as a separate Strategy entry point so a future divergence (e.g.
    a dedicated protein-protein target model) lands here without touching the
    peptide path.
    """
    return process_schrodinger_peptide_interactions(
        current_structure_obj=current_structure_obj,
        ligand_chain=ligand_chain,
        pdb_code_str=pdb_code_str,
        schrodinger_interactions_dir_override=schrodinger_interactions_dir_override,
    )


# Three-letter -> one-letter for peptide residues. Engine 2 peptides include
# non-standard / unnatural residues (DPN, NLE, HIE, ...) that have no protwis
# Residue row (peptides are not in the Residue table) — so we cannot rely on the
# DB. Fall back to "X" for anything unrecognised (the legacy contactnetwork
# path consults unnatural_amino_acids.yaml; mirroring that fully is out of D3
# scope — "X" is a faithful, non-crashing placeholder, flagged in the Inbox).
_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Common Schrödinger/PrepWizard protonation-state aliases:
    "HIE": "H", "HID": "H", "HIP": "H", "ASH": "D", "GLH": "E",
    "LYN": "K", "CYX": "C", "ARN": "R",
}


def _three_to_one(resname):
    """Map a (possibly non-standard) three-letter residue name to one letter,
    falling back to 'X' for unknown / unnatural residues."""
    return _THREE_TO_ONE.get((resname or "").upper(), "X")
