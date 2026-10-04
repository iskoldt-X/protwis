"""The ligand 3D file of an anchor (StructureLigandInteraction.pdb_file).

The interaction page loads it into each ligand's viewer (through
interaction.views.download) and structure/pdb/<pdb>/ligand/<lig> serves it.
It used to be a by-product of the legacy calculation. The imports write it
now: the ligand and the receptor residues the import wrote interaction rows
for, as lines of GPCRdb's own stored structure text, so the viewer shows the
residues the table lists, with the names and coordinates the rest of the site
uses.

The ligand is found in that text by coordinates. Product coordinates come
through structure preparation and are written again with three decimals; on
the corpus they differ from the stored text by at most 0.007 A (measured
2026-10-04 on 533 ligand atoms of 60 instances), far below any bond length.
"""

from interaction import schrodinger_chain_map as cm
from structure.models import PdbData

# A product ligand atom names the residue of the stored text that has an atom
# within this distance of it.
LIGAND_MATCH_TOL = 0.05  # A


def _is_hydrogen(line):
    element = line[76:78].strip() if len(line) >= 78 else ""
    return element in ("H", "D") or (not element and line[12:16].strip().startswith("H"))


def text_residues(gpcrdb_text):
    """[(residue id, lines, heavy atom xyz)] of the stored text, in text order.

    First model only, waters left out; a residue id is
    (chain, residue number, insertion code, residue name), read from the
    columns parse_gpcrdb_pdb reads.
    """
    order, lines, atoms = [], {}, {}
    for line in (gpcrdb_text or "").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            continue
        resname = line[17:20].strip()
        if resname in cm.WATER:
            continue
        rid = (line[21], line[22:26].strip(), line[26].strip() if len(line) > 26 else "", resname)
        if rid not in lines:
            order.append(rid)
            lines[rid], atoms[rid] = [], []
        lines[rid].append(line.rstrip("\n"))
        if not _is_hydrogen(line):
            atoms[rid].append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return [(rid, lines[rid], atoms[rid]) for rid in order]


def ligand_line_xyz(lines):
    """Heavy-atom coordinates of standard-column ligand lines."""
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54 and not _is_hydrogen(line):
            out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return out


def _cell(xyz):
    return tuple(int(c // LIGAND_MATCH_TOL) for c in xyz)


def complex_text(gpcrdb_text, receptor_chain, receptor_seqs, ligand_xyz=(), ligand_chain=""):
    """The anchor's 3D file, or "" when no residue of the text is the ligand.

    The ligand is every residue of the text with a heavy atom within
    LIGAND_MATCH_TOL of a ``ligand_xyz`` atom, or, for an anchor that is a
    chain, every residue on ``ligand_chain``. The receptor part is every
    residue on ``receptor_chain`` whose number is in ``receptor_seqs`` (no
    insertion code). Lines keep the text's order and columns.
    """
    grid = {}
    for xyz in ligand_xyz:
        grid.setdefault(_cell(xyz), []).append(xyz)
    tol2 = LIGAND_MATCH_TOL ** 2

    def matches(xyz):
        cx, cy, cz = _cell(xyz)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for other in grid.get((cx + dx, cy + dy, cz + dz), ()):
                        if sum((a - b) ** 2 for a, b in zip(xyz, other)) <= tol2:
                            return True
        return False

    residues = text_residues(gpcrdb_text)
    ligand = {rid for rid, _lines, atoms in residues
              if (ligand_chain and rid[0] == ligand_chain) or any(matches(a) for a in atoms)}
    if not ligand:
        return ""
    seqs = {int(s) for s in receptor_seqs}
    out = []
    for rid, lines, _atoms in residues:
        if rid in ligand or (rid[0] == receptor_chain and not rid[2]
                             and rid[1].lstrip("-").isdigit() and int(rid[1]) in seqs):
            out.extend(lines)
    return "\n".join(out + ["END"]) + "\n"


def write_complex_file(sli, text):
    """Point the anchor at ``text`` as its 3D file; "" leaves it with none.

    Returns (status, replaced PdbData id or None). status is "unchanged"
    (the anchor's file already holds exactly this text, as on a second run),
    "written", "cleared" or "none". The replaced PdbData is not deleted here:
    StructureLigandInteraction.pdb_file cascades, so the caller deletes it only
    once nothing references it (schrodinger_import.delete_unreferenced_pdbdata).
    """
    old = sli.pdb_file_id
    if not text:
        if old is None:
            return "none", None
        sli.pdb_file = None
        sli.save(update_fields=["pdb_file"])
        return "cleared", old
    if old is not None and sli.pdb_file.pdb == text:
        return "unchanged", None
    sli.pdb_file = PdbData.objects.create(pdb=text)
    sli.save(update_fields=["pdb_file"])
    return "written", old
