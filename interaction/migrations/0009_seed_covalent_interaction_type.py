"""Add the interaction type Engine 1's Covalent rows route to.

interaction_type_map.yaml routes (Covalent, '') to the slug ``covalent``, which
no earlier build created. A Covalent row is a bond of order >= 1 from a ligand
atom to a receptor atom, as the Suite drew it from the input file's
_struct_conn rows or as one of its bond builders added it during preparation;
it says a bond exists, not its order.

The type is ``covalent`` and not ``hidden``: pages leave out hidden types, and
so does the scorecard (scorecard/score.py counts only types that are not
hidden). The rows are imported so that both show them (Binghan 2026-09-28).

An existing row with this slug is left as it is. The reverse operation is a
deliberate no-op: deleting the type would cascade to every
ResidueFragmentInteraction that uses it.
"""

from django.db import migrations


COVALENT_TYPE = ("covalent", "covalent bond", "covalent", "")  # slug, name, type, direction


def seed(apps, schema_editor):
    InteractionType = apps.get_model("interaction", "ResidueFragmentInteractionType")
    slug, name, type_, direction = COVALENT_TYPE
    InteractionType.objects.get_or_create(
        slug=slug, defaults={"name": name, "type": type_, "direction": direction})


class Migration(migrations.Migration):

    dependencies = [
        ("interaction", "0008_seed_engine1_interaction_types"),
    ]

    operations = [
        migrations.RunPython(seed, migrations.RunPython.noop),
    ]
