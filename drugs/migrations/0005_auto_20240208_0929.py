# Generated by Django 3.2.18 on 2024-02-08 08:29

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0006_releasestatistics_database'),
        ('drugs', '0004_auto_20240208_0927'),
    ]

    operations = [
        migrations.RenameField(
            model_name='drugs2024',
            old_name='publication',
            new_name='reference',
        ),
        migrations.AddField(
            model_name='drugs',
            name='publication',
            field=models.ManyToManyField(to='common.Publication'),
        ),
    ]
