# Generated by Django 3.2.25 on 2024-12-04 11:28

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('ligand', '0025_auto_20230831_1052'),
        ('common', '0006_releasestatistics_database'),
        ('protein', '0021_auto_20240223_1135'),
        ('drugs', '0003_drugs_publication'),
    ]

    operations = [
        migrations.CreateModel(
            name='Indication',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('title', models.CharField(max_length=255)),
                ('code', models.CharField(max_length=8, null=True)),
                ('slug', models.CharField(max_length=255, null=True, unique=True)),
                ('level', models.IntegerField(null=True)),
                ('parent', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='drugs.indication')),
                ('uri', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='common.weblink')),
            ],
            options={
                'db_table': 'indication',
            },
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='indication',
        ),
        migrations.AddField(
            model_name='drugs',
            name='indication',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='drugs.indication'),
        ),
        migrations.RenameField(
            model_name='drugs',
            old_name='publication',
            new_name='reference',
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='approval',
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='clinicalstatus',
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='drugtype',
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='externallink',
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='name',
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='novelty',
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='phase',
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='phasedate',
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='references',
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='status',
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='synonym',
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='targetlevel',
        ),
        migrations.AddField(
            model_name='drugs',
            name='affected_pathway',
            field=models.CharField(max_length=30, null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='approval_year',
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='charge',
            field=models.CharField(max_length=5, null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='complexity',
            field=models.FloatField(max_length=4, null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='drug_status',
            field=models.CharField(max_length=15, null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='genetic_association',
            field=models.CharField(max_length=30, null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='indication_max_phase',
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='indication_status',
            field=models.CharField(max_length=15, null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='ligand',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='ligand.ligand'),
        ),
        migrations.AddField(
            model_name='drugs',
            name='novelty_score',
            field=models.FloatField(max_length=4, null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='publication_count',
            field=models.IntegerField(null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='similarity_to_model',
            field=models.FloatField(max_length=4, null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='somatic_mutation',
            field=models.CharField(max_length=30, null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='target_level',
            field=models.CharField(max_length=30, null=True),
        ),
        migrations.AddField(
            model_name='drugs',
            name='tpsa',
            field=models.CharField(max_length=10, null=True),
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='moa',
        ),
        migrations.AddField(
            model_name='drugs',
            name='moa',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='ligand.ligandrole'),
        ),
        migrations.RemoveField(
            model_name='drugs',
            name='target',
        ),
        migrations.AddField(
            model_name='drugs',
            name='target',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='protein.protein'),
        ),
        migrations.CreateModel(
            name='IndicationAssociation',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('association_score', models.FloatField(max_length=7, null=True)),
                ('ot_genetics_portal', models.FloatField(max_length=7, null=True)),
                ('gene_burden', models.FloatField(max_length=7, null=True)),
                ('clingen', models.FloatField(max_length=7, null=True)),
                ('gene2phenotype', models.FloatField(max_length=7, null=True)),
                ('orphanet', models.FloatField(max_length=7, null=True)),
                ('genomics_england', models.FloatField(max_length=7, null=True)),
                ('uniprot_literature', models.FloatField(max_length=7, null=True)),
                ('uniprot_variants', models.FloatField(max_length=7, null=True)),
                ('sysbio', models.FloatField(max_length=7, null=True)),
                ('cancer_gene_census', models.FloatField(max_length=7, null=True)),
                ('cancer_biomarkers', models.FloatField(max_length=7, null=True)),
                ('intogen', models.FloatField(max_length=7, null=True)),
                ('eva', models.FloatField(max_length=7, null=True)),
                ('eva_somatic', models.FloatField(max_length=7, null=True)),
                ('chembl', models.FloatField(max_length=7, null=True)),
                ('slapenrich', models.FloatField(max_length=7, null=True)),
                ('crispr', models.FloatField(max_length=7, null=True)),
                ('crispr_screen', models.FloatField(max_length=7, null=True)),
                ('reactome', models.FloatField(max_length=7, null=True)),
                ('europepmc', models.FloatField(max_length=7, null=True)),
                ('expression_atlas', models.FloatField(max_length=7, null=True)),
                ('impc', models.FloatField(max_length=7, null=True)),
                ('indication', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='drugs.indication')),
                ('target', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='protein.protein')),
            ],
            options={
                'db_table': 'indication_association',
            },
        ),
        migrations.CreateModel(
            name='ATCCodes',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('code', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='common.weblink')),
                ('ligand', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='ligand.ligand')),
            ],
            options={
                'db_table': 'atc_codes',
            },
        ),
        migrations.AddField(
            model_name='drugs',
            name='disease_association',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='drugs.indicationassociation'),
        ),
    ]