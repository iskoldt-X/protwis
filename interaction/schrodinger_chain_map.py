"""Reconcile chain names between Schrodinger products and GPCRdb.

The products name chains as the RCSB mmCIF does (author chain, up to four
characters). GPCRdb names them as its stored PDB-format structure text does
(one character, sometimes renamed, split or hand-edited by curators). This
module builds two maps offline, once per (GPCRdb dump, product run):

* anchor map: one row per GPCRdb ligand-anchor copy (pdb, HET, chain_res
  token) naming the product instance that is the same ligand copy;
* receptor map: one row per structure naming the product (author) chain whose
  atoms are GPCRdb's preferred chain.

Exact coordinate identity decides. The annotation's label_asym_id is a second,
independent witness: it confirms coordinate answers, flags annotation errors
when it disagrees, and is the only thing that can confirm a name-based
fallback where GPCRdb stores an older model whose coordinates no longer match.

No database access and no Django imports: the command in
interaction/management/commands/build_schrodinger_chain_map.py feeds this
module with text it reads from the database and the file system.
"""

import collections
import re

# ---------------------------------------------------------------------------
# mmCIF _atom_site (author side)
# ---------------------------------------------------------------------------

_CIF_TOKEN = re.compile(r"""'(?:[^']|'(?=\S))*'(?=\s|$)|"(?:[^"]|"(?=\S))*"(?=\s|$)|\S+""")

WATER = frozenset({"HOH", "DOD", "WAT"})


class ParseError(ValueError):
    """An input file does not have the shape this module relies on."""


def _cif_tokens(line):
    out = []
    for m in _CIF_TOKEN.finditer(line):
        tok = m.group(0)
        if len(tok) > 1 and tok[0] == tok[-1] and tok[0] in "'\"":
            tok = tok[1:-1]
        out.append(tok)
    return out


def coord_key(x, y, z):
    """Coordinates as a hashable key at the 3-decimal precision both formats use."""
    return "%.3f %.3f %.3f" % (float(x), float(y), float(z))


def parse_mmcif_atoms(text):
    """Return first-model, non-hydrogen, non-water atoms of an mmCIF _atom_site loop.

    Columns are read by name; every data row must have exactly as many tokens
    as the loop declares, otherwise ParseError. Each atom is a dict with
    label_asym, auth_asym, comp, auth_seq, icode, atom, group, key.
    """
    lines = text.splitlines()
    i = 0
    cols = []
    while i < len(lines):
        if lines[i].startswith("loop_") and i + 1 < len(lines) and lines[i + 1].startswith("_atom_site."):
            i += 1
            while i < len(lines) and lines[i].startswith("_atom_site."):
                cols.append(lines[i].split(".", 1)[1].strip())
                i += 1
            break
        i += 1
    rows = []
    while i < len(lines) and lines[i].startswith(("ATOM", "HETATM")):
        tokens = _cif_tokens(lines[i])
        if len(tokens) != len(cols):
            raise ParseError("atom_site row has %d fields, loop declares %d" % (len(tokens), len(cols)))
        rows.append(dict(zip(cols, tokens)))
        i += 1
    if not rows:
        raise ParseError("no _atom_site rows found")
    first_model = rows[0].get("pdbx_PDB_model_num")
    atoms = []
    for r in rows:
        if r.get("pdbx_PDB_model_num") != first_model:
            continue
        comp = r.get("auth_comp_id", r["label_comp_id"])
        if comp in WATER or r.get("type_symbol") in ("H", "D"):
            continue
        icode = r.get("pdbx_PDB_ins_code", "?")
        atoms.append({
            "label_asym": r["label_asym_id"],
            "auth_asym": r["auth_asym_id"],
            "comp": comp,
            "auth_seq": r["auth_seq_id"],
            "icode": "" if icode in ("?", ".") else icode,
            "atom": r.get("auth_atom_id", r["label_atom_id"]),
            "group": r["group_PDB"],
            "key": coord_key(r["Cartn_x"], r["Cartn_y"], r["Cartn_z"]),
        })
    return atoms


# ---------------------------------------------------------------------------
# GPCRdb stored PDB-format text (GPCRdb side)
# ---------------------------------------------------------------------------

def parse_gpcrdb_pdb(text):
    """Return first-model, non-hydrogen, non-water atoms of GPCRdb's stored text.

    Fields are read exactly as build_structures reads them: chain = column 22
    (line[21]), residue number = columns 23-26, residue name = columns 18-20
    (so five-character CCD codes appear truncated to three characters).
    """
    atoms = []
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            continue
        resname = line[17:20].strip()
        element = line[76:78].strip() if len(line) >= 78 else ""
        atom = line[12:16].strip()
        if resname in WATER or element in ("H", "D") or (not element and atom.startswith("H")):
            continue
        atoms.append({
            "chain": line[21],
            "resnum": line[22:26].strip(),
            "icode": line[26].strip() if len(line) > 26 else "",
            "resname": resname,
            "atom": atom,
            "group": line[:6].strip(),
            "key": coord_key(line[30:38], line[38:46], line[46:54]),
        })
    return atoms


# ---------------------------------------------------------------------------
# Annotation (upstream ligands.tsv) -> label_asym_id per copy
# ---------------------------------------------------------------------------

def annotation_labels(rows):
    """Map (PDB, HET, token) -> label_asym_id from ligands.tsv rows.

    Residue_seq_id and label_asym_id are comma lists aligned copy for copy.
    Rows whose two lists differ in length contribute nothing.
    """
    out = {}
    for r in rows:
        tokens = [t.strip() for t in (r.get("Residue_seq_id") or "").split(",") if t.strip()]
        labels = [t.strip() for t in (r.get("label_asym_id") or "").split(",") if t.strip()]
        if not tokens or len(tokens) != len(labels):
            continue
        key = ((r.get("PDB") or "").strip().upper(), (r.get("Name") or "").strip().upper())
        for tok, lab in zip(tokens, labels):
            out[key + (tok.replace(" ", ""),)] = lab
    return out


# ---------------------------------------------------------------------------
# Anchor map
# ---------------------------------------------------------------------------

TOKEN_RE = re.compile(r"^(?P<chain>[A-Za-z0-9]):(?P<resnum>-?\d+)(?P<icode>[A-Za-z]?)$")

ANCHOR_COLUMNS = ("pdb", "het", "token", "instance", "status", "source",
                  "n_exact", "n_gpcrdb", "label", "label_instance", "note")


def instance_name(het, auth_asym, auth_seq, icode):
    return "{}_{}_{}{}".format(het, auth_asym, auth_seq, icode)


def split_tokens(chain_res):
    """chain_res -> list of 'C:NUM' tokens, or [] when empty / chain only."""
    items = [t.strip() for t in (chain_res or "").split(",") if t.strip()]
    if items and all(TOKEN_RE.match(t) for t in items):
        return items
    return []


def resolve_anchor(pdb, het, token, cif_atoms, gpcrdb_atoms, product_instances, label):
    """Decide which product instance is the ligand copy GPCRdb calls `token`.

    Returns a dict with the ANCHOR_COLUMNS fields. Rules:

    1. exact coordinates: the GPCRdb residue (chain, number, name[:3]) is
       compared atom by atom with every product instance of the HET; exactly
       one instance with the most shared atoms decides (source coord_exact).
    2. the label route names an instance via label_asym_id; it confirms rule 1
       (source coord+label), or disagrees (status errata, rule 1 still wins),
       or is absent (source coord_only).
    3. no shared atoms at all (GPCRdb stores an older model): the name-based
       candidate HET_<chain>_<number> is accepted only when the label route
       names the same instance (source fallback+label); otherwise unresolved.
    4. the product has no instance of this HET: no_product.
    """
    het = het.upper()
    copies = sorted(i for i in product_instances if i.split("_", 1)[0].upper() == het)
    row = {"pdb": pdb, "het": het, "token": token, "instance": "", "status": "",
           "source": "", "n_exact": 0, "n_gpcrdb": 0, "label": label or "",
           "label_instance": "", "note": ""}
    m = TOKEN_RE.match(token)
    chain, resnum, icode = m.group("chain"), m.group("resnum"), m.group("icode")

    by_instance = collections.defaultdict(set)
    label_candidates = set()
    for a in cif_atoms:
        if a["comp"].upper() != het:
            continue
        name = instance_name(het, a["auth_asym"], a["auth_seq"], a["icode"])
        by_instance[name].add(a["key"])
        if label and a["label_asym"] == label:
            label_candidates.add(name)
    if len(label_candidates) == 1:
        row["label_instance"] = next(iter(label_candidates))
    elif len(label_candidates) > 1:
        row["note"] = "label names several residues"

    gkeys = {a["key"] for a in gpcrdb_atoms
             if a["chain"] == chain and a["resnum"] == resnum and a["icode"] == icode
             and a["resname"].upper() == het[:3]}
    row["n_gpcrdb"] = len(gkeys)

    if not copies:
        row["status"], row["note"] = "no_product", (row["note"] or "product has no instance of this HET")
        return row

    scores = sorted(((len(gkeys & by_instance.get(c, set())), c) for c in copies), reverse=True)
    best, best_name = scores[0]
    if best > 0:
        if len(scores) > 1 and scores[1][0] == best:
            row["status"], row["note"] = "unresolved", "two instances share the same coordinates"
            return row
        row["instance"], row["n_exact"] = best_name, best
        if not row["label_instance"]:
            row["status"], row["source"] = "ok", "coord_only"
        elif row["label_instance"] == best_name:
            row["status"], row["source"] = "ok", "coord+label"
        else:
            row["status"], row["source"] = "errata", "coord_exact"
            row["note"] = "annotation label_asym_id points at {}".format(row["label_instance"])
        return row

    fallback = instance_name(het, chain, resnum, icode)
    if fallback in copies and row["label_instance"] == fallback:
        row["instance"], row["status"], row["source"] = fallback, "ok", "fallback+label"
        row["note"] = "no exact coordinates (older model in GPCRdb?)"
        return row
    row["status"] = "unresolved"
    row["note"] = "no exact coordinates; name candidate {} {}; label candidate {}".format(
        fallback, "exists" if fallback in copies else "absent", row["label_instance"] or "-")
    return row


# ---------------------------------------------------------------------------
# Receptor map
# ---------------------------------------------------------------------------

RECEPTOR_COLUMNS = ("pdb", "preferred_chain", "auth_chain", "status", "method",
                    "n_ca_gpcrdb", "n_ca_matched", "renumbered", "note")


def resolve_receptor(pdb, preferred_chain, cif_atoms, gpcrdb_atoms):
    """Name the product (author) chain whose CA atoms are GPCRdb's preferred chain.

    exact: the author chain sharing the most CA coordinates with the GPCRdb
    preferred chain; every matched CA must keep its residue number, otherwise
    status renumbered (the importer must refuse the structure).
    identity_drift: no CA matches at all (older model in GPCRdb); accepted only
    if an author chain of the same name carries every GPCRdb (number, name)
    pair of the preferred chain.
    """
    pref = (preferred_chain or "").split(",")[0].strip()
    row = {"pdb": pdb, "preferred_chain": pref, "auth_chain": "", "status": "",
           "method": "", "n_ca_gpcrdb": 0, "n_ca_matched": 0, "renumbered": 0, "note": ""}
    g_ca = {a["key"]: (a["resnum"], a["icode"], a["resname"]) for a in gpcrdb_atoms
            if a["group"] == "ATOM" and a["atom"] == "CA" and a["chain"] == pref}
    row["n_ca_gpcrdb"] = len(g_ca)
    if not g_ca:
        row["status"], row["note"] = "unresolved", "GPCRdb text has no CA on the preferred chain"
        return row
    per_chain = collections.Counter()
    renumbered = collections.Counter()
    for a in cif_atoms:
        if a["group"] == "ATOM" and a["atom"] == "CA" and a["key"] in g_ca:
            per_chain[a["auth_asym"]] += 1
            if (a["auth_seq"], a["icode"]) != g_ca[a["key"]][:2]:
                renumbered[a["auth_asym"]] += 1
    if per_chain:
        (chain, n), = per_chain.most_common(1)
        row.update(auth_chain=chain, n_ca_matched=n, method="exact",
                   renumbered=renumbered[chain])
        row["status"] = "renumbered" if renumbered[chain] else "ok"
        if len(per_chain) > 1:
            row["note"] = "CA also matched on " + ",".join(
                "%s:%d" % kv for kv in sorted(per_chain.items()) if kv[0] != chain)
        return row
    wanted = {(v[0], v[1], v[2].upper()) for v in g_ca.values()}
    have = {(a["auth_seq"], a["icode"], a["comp"].upper()) for a in cif_atoms
            if a["group"] == "ATOM" and a["atom"] == "CA" and a["auth_asym"] == pref}
    if have and wanted <= have:
        row.update(auth_chain=pref, method="identity_drift", status="ok",
                   note="no CA coordinate matches; every GPCRdb (number, name) found on the same-named chain")
    else:
        row["status"] = "unresolved"
        row["note"] = "no CA coordinate matches and the same-named chain does not carry the GPCRdb residues"
    return row
