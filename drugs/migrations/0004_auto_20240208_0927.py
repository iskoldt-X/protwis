# Generated by Django 3.2.18 on 2024-02-08 08:27

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0025_auto_20230831_1052'),
        ('common', '0006_releasestatistics_database'),
        ('drugs', '0003_drugs_publication'),
    ]

    operations = [
        migrations.CreateModel(
            name='Indication',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=55)),
                ('code', models.CharField(max_length=10)),
            ],
            options={
                'db_table': 'indication',
            },
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='publication',
        ),
        migrations.AlterModelTable(
            name='drugs',
            table='drugs_old',
        ),
        migrations.CreateModel(
            name='Drugs2024',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('charge', models.CharField(max_length=2)),
                ('complexity', models.FloatField(max_length=4)),
                ('tpsa', models.CharField(max_length=10)),
                ('drug_status', models.CharField(max_length=15)),
                ('approval_year', models.IntegerField(max_length=4)),
                ('indication_max_phase', models.IntegerField(max_length=1)),
                ('moa', models.CharField(max_length=30, null=True)),
                ('affected_pathway', models.CharField(max_length=30)),
                ('somatic_mutation', models.CharField(max_length=30)),
                ('similarity_to_model', models.FloatField(max_length=4)),
                ('novelty_score', models.FloatField(max_length=4)),
                ('indication', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='drugs.indication')),
                ('ligand', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='ligand.ligand')),
                ('publication', models.ManyToManyField(to='common.Publication')),
                ('target', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='protein.protein')),
            ],
            options={
                'db_table': 'drugs',
            },
        ),
    ]
