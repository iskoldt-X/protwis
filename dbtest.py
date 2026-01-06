from protein.models import ProteinFamily, Protein
from ligand.models import LigandType, Ligand
from django.db.models import Prefetch
import random

print("\n=== Root ProteinFamily (parent is NULL) — top 5 ===")
roots = list(ProteinFamily.objects.filter(parent__isnull=True).order_by("id")[:5])
for pf in roots:
    print(f"- id={pf.id}  name={pf.name!r}  slug={pf.slug!r}")

print("\n=== Pick one root family to inspect children/grandchildren ===")
# Prefer something that looks like a GPCR root (often slug '000'), otherwise first root
root = ProteinFamily.objects.filter(slug="000").first() or (roots[0] if roots else None)
if not root:
    print("No root ProteinFamily rows found.")
else:
    print(f"Using root: id={root.id}  name={root.name!r}  slug={root.slug!r}")

    children = list(ProteinFamily.objects.filter(parent=root).order_by("slug", "id"))
    print(f"\nChildren of {root.slug!r} ({len(children)}):")
    for c in children:
        print(f"  - id={c.id}  name={c.name!r}  slug={c.slug!r}")

    # Show grandchildren for first few children (so output stays readable)
    for c in children[:5]:
        gkids = list(ProteinFamily.objects.filter(parent=c).order_by("slug", "id"))
        print(f"\nGrandchildren of child {c.slug!r} ({len(gkids)}):")
        for g in gkids[:30]:
            print(f"    - id={g.id}  name={g.name!r}  slug={g.slug!r}")
        if len(gkids) > 30:
            print("    ... (truncated)")

print("\n=== LigandType slugs — all ===")
lts = list(LigandType.objects.all().order_by("slug", "id"))
if not lts:
    print("No LigandType rows found.")
else:
    for lt in lts:
        print(f"- id={lt.id}  name={lt.name!r}  slug={lt.slug!r}")

print("\n=== One random Protein and its family lineage ===")
p = Protein.objects.select_related("family", "family__parent", "family__parent__parent", "family__parent__parent__parent").order_by("?").first()
if not p:
    print("No Protein rows found.")
else:
    fam = p.family
    print(f"Protein: id={p.id} entry_name={p.entry_name!r} accession={p.accession!r} name={p.name!r}")
    if fam:
        p1 = fam.parent
        p2 = p1.parent if p1 else None
        p3 = p2.parent if p2 else None
        print(f"  Family:        id={fam.id} name={fam.name!r} slug={fam.slug!r}")
        print(f"  Parent:        {None if not p1 else (p1.name, p1.slug)}")
        print(f"  Grandparent:   {None if not p2 else (p2.name, p2.slug)}")
        print(f"  Great-grandp.: {None if not p3 else (p3.name, p3.slug)}")
    else:
        print("  Protein has no family? (unexpected)")

    # Check whether this Protein has any direct relation to LigandType:
    # It does NOT in the schema; LigandType is linked via Ligand, and via Endogenous_GTP/ProteinCouplings, etc.
    # Still, we can show a couple of LigandTypes observed among endogenous ligands for this receptor, if any.
    try:
        from ligand.models import Endogenous_GTP
        lt_names = (
            Endogenous_GTP.objects
            .filter(receptor=p)
            .values_list("ligand__ligand_type__name", flat=True)
            .distinct()
        )
        lt_names = [x for x in lt_names if x]
        print(f"  Endogenous_GTP ligand_type names (distinct): {lt_names[:20]}")
        if len(lt_names) > 20:
            print("  ... (truncated)")
    except Exception as e:
        print(f"  (Could not query Endogenous_GTP ligand types: {e})")

print("\nDone.")