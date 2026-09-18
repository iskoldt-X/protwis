"""Add the three interaction types Engine 1 writes that older databases lack.

aro_ion_protein (pi-cation), halogen_protein and metal_coordination_protein
are the targets of interaction_type_map.yaml rules that no earlier build
created. Existing rows are never modified: a slug that is already present is
left exactly as it is.

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
        InteractionType.objects.get_or_create(
            slug=slug, defaults={"name": name, "type": type_, "direction": direction})


class Migration(migrations.Migration):

    dependencies = [
        ("interaction", "0007_structureligandinteraction_site_chain_res"),
    ]

    operations = [
        migrations.RunPython(seed, migrations.RunPython.noop),
    ]
