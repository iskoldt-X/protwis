"""Add the three interaction types Engine 1 writes that older databases lack.

aro_ion_protein (pi-cation), halogen_protein and metal_coordination_protein
are the targets of interaction_type_map.yaml rules that no earlier build
created. Existing rows are never modified: a slug that is already present is
left as it is, except that an existing halogen_protein row is renamed to
"halogen contact".

halogen_protein is named "halogen contact", not "halogen bond": the producer
criterion is a halogen within 3.5 A of a polar atom, and its angle limits
remove no rows, so the rows are contacts rather than proven halogen bonds.

The reverse operation is a deliberate no-op. Deleting interaction types would
cascade to every ResidueFragmentInteraction that uses them.
"""

from django.db import migrations


ENGINE1_TYPES = (
    # slug, name, type, direction
    ("aro_ion_protein", "aromatic (pi-cation)", "aromatic", "protein"),
    ("halogen_protein", "halogen contact", "polar", ""),
    ("metal_coordination_protein", "metal coordination", "polar", ""),
)


def seed(apps, schema_editor):
    InteractionType = apps.get_model("interaction", "ResidueFragmentInteractionType")
    for slug, name, type_, direction in ENGINE1_TYPES:
        row, created = InteractionType.objects.get_or_create(
            slug=slug, defaults={"name": name, "type": type_, "direction": direction})
        # A database seeded from the earlier fixture calls this slug
        # "halogen bond"; the name is the only field corrected here.
        if not created and slug == "halogen_protein" and row.name != name:
            row.name = name
            row.save(update_fields=["name"])


class Migration(migrations.Migration):

    dependencies = [
        ("interaction", "0007_structureligandinteraction_site_chain_res"),
    ]

    operations = [
        migrations.RunPython(seed, migrations.RunPython.noop),
    ]
