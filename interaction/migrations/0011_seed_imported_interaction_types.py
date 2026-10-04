"""Seed every interaction type the Schrodinger imports write.

build_structures no longer computes ligand interactions, so nothing in a build
creates the types the legacy calculation used to create on the fly. These ten
complete schrodinger_import.required_slugs() -- the targets of
interaction_type_map.yaml outside the excluded families, plus polar_backbone,
the backbone override -- with 0009 and 0010 seeding the other four. Each keeps
the slug, name, type and direction the legacy calculation gave it, so a
database built either way names them alike.

Existing rows are never modified. The reverse operation is a deliberate no-op:
deleting interaction types would cascade to every ResidueFragmentInteraction
that uses them.
"""

from django.db import migrations


IMPORTED_TYPES = (
    # slug, name, type, direction
    ("acc", "accessible", "hidden", ""),
    ("aro_ef_protein", "aromatic (edge-to-face)", "aromatic", "protein"),
    ("aro_ff", "aromatic (face-to-face)", "aromatic", "none"),
    ("hyd", "hydrophobic", "hydrophobic", ""),
    ("polar_acceptor_protein", "polar (hydrogen bond)", "polar", "protein"),
    ("polar_backbone", "polar (hydrogen bond with backbone)", "polar", "protein"),
    ("polar_donor_protein", "polar (hydrogen bond)", "polar", "protein"),
    ("polar_double_neg_protein", "polar (charge-charge)", "polar", ""),
    ("polar_double_pos_protein", "polar (charge-charge)", "polar", ""),
    ("Van der Waals", "Van der Waals", "waals", ""),
)


def seed(apps, schema_editor):
    InteractionType = apps.get_model("interaction", "ResidueFragmentInteractionType")
    for slug, name, type_, direction in IMPORTED_TYPES:
        InteractionType.objects.get_or_create(
            slug=slug, defaults={"name": name, "type": type_, "direction": direction})


class Migration(migrations.Migration):

    dependencies = [
        ("interaction", "0010_seed_covalent_interaction_type"),
    ]

    operations = [
        migrations.RunPython(seed, migrations.RunPython.noop),
    ]
