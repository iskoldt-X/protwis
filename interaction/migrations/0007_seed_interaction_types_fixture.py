"""Seed ResidueFragmentInteractionType with the 21-slug fixture.

Phase 1c axis 1D: the import_schrodinger_interactions command relies on
the 21-slug catalog (ADR-011/012/013) being present. Build-time loaddata
was the legacy path; a data migration makes the seed automatic on
manage.py migrate so a fresh deployment is not silently underseeded.

The migration is idempotent: it upserts by slug rather than inserting
blindly, so re-running on an already-seeded DB is a no-op. Removing a
row from the fixture does NOT delete it from the DB (delete is an
explicit operational decision, not a side effect of a fixture edit).
"""

import json
import os

from django.db import migrations


FIXTURE_PATH = os.path.join(
    os.path.dirname(__file__), '..', 'fixtures', 'interaction_types.json')


def seed(apps, schema_editor):
    Model = apps.get_model('interaction', 'ResidueFragmentInteractionType')
    with open(FIXTURE_PATH, 'r') as f:
        rows = json.load(f)
    for row in rows:
        if row.get('model') != 'interaction.residuefragmentinteractiontype':
            continue
        fields = row['fields']
        Model.objects.update_or_create(
            slug=fields['slug'],
            defaults={
                'name': fields['name'],
                'type': fields['type'],
                'direction': fields.get('direction', ''),
            },
        )


def unseed(apps, schema_editor):
    Model = apps.get_model('interaction', 'ResidueFragmentInteractionType')
    with open(FIXTURE_PATH, 'r') as f:
        rows = json.load(f)
    slugs = [r['fields']['slug'] for r in rows
             if r.get('model') == 'interaction.residuefragmentinteractiontype']
    Model.objects.filter(slug__in=slugs).delete()


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0006_auto_20241031_1314'),
    ]

    operations = [
        migrations.RunPython(seed, unseed),
    ]
