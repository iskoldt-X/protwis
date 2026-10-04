"""Stored interactions of a structure GPCRdb has, for the online calculation page.

interaction.views.calculate computes interactions for a PDB file a user
uploads or a PDB code a user types. For a typed code that GPCRdb already has,
the answer is the interactions stored in the database -- imported from the
Schrodinger deliveries -- rather than a fresh run of the legacy calculation on
a file downloaded from RCSB, which can disagree with the structure page. An
uploaded file is still calculated: it may not be the deposited structure.

The results come in the shape interaction.views.calculate_interactions returns,
so the page needs no other change.
"""

import collections

from interaction.models import ResidueFragmentInteraction
from structure.models import Structure

PEPTIDE_REFERENCE = "PEP"

THREE_LETTER = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
    "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
    "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}


def ligand_key(pdb_reference, ligand_name):
    """The results key of an anchor: its HET code, or the ligand's name for a "pep" chain."""
    reference = (pdb_reference or "").strip().upper()
    return ligand_name if reference == PEPTIDE_REFERENCE else reference


def build_results(rows, chain):
    """Results from (ligand key, one-letter amino acid, residue number, slug, name, type,
    direction) rows, in the shape calculate_interactions returns.

    Each interaction is [residue, fragment file, slug, name, type, direction]
    with residue as three-letter name, number and chain (ASP113A, what
    interaction.views.regexaa reads); the fragment file is left empty, the page
    does not read it. A residue that is not a standard amino acid is left out.
    The ligand with most rows comes first, ties in key order: calculate takes the
    first ligand as the main one. The score is the number of rows.
    """
    per = collections.OrderedDict()
    for key, amino_acid, number, slug, name, type_, direction in rows:
        three = THREE_LETTER.get((amino_acid or "").upper())
        if three is None:
            continue
        per.setdefault(key, []).append(
            ["{}{}{}".format(three, number, chain), "", slug, name, type_ or "", direction or ""])
    ordered = sorted(per.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    return collections.OrderedDict(
        (key, {"score": len(interactions), "interactions": interactions})
        for key, interactions in ordered)


def stored_results(pdbname):
    """(results, stored structure text) for an experimental structure GPCRdb has, else None."""
    structure = (Structure.objects
                 .filter(pdb_code__index__iexact=(pdbname or "").strip(),
                         structure_type__origin="experiment")
                 .select_related("pdb_data")
                 .first())
    if structure is None or structure.pdb_data is None:
        return None
    chain = (structure.preferred_chain or "").split(",")[0].strip()
    rows = (ResidueFragmentInteraction.objects
            .filter(structure_ligand_pair__structure=structure)
            .order_by("structure_ligand_pair_id", "rotamer__residue__sequence_number",
                      "interaction_type__slug", "id")
            .values_list("structure_ligand_pair__pdb_reference",
                         "structure_ligand_pair__ligand__name",
                         "rotamer__residue__amino_acid",
                         "rotamer__residue__sequence_number",
                         "interaction_type__slug", "interaction_type__name",
                         "interaction_type__type", "interaction_type__direction"))
    return (build_results(((ligand_key(ref, name), aa, number, slug, tname, ttype, direction)
                           for ref, name, aa, number, slug, tname, ttype, direction in rows), chain),
            structure.pdb_data.pdb)
