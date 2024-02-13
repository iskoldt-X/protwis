# Generated by Django 3.2.18 on 2024-02-12 14:44

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0018_auto_20231026_2107'),
    ]

    operations = [
        migrations.CreateModel(
            name='CancerType',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('slug', models.SlugField(max_length=100, unique=True)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'cancer',
            },
        ),
        migrations.CreateModel(
            name='ExpressionValue',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('max_expression', models.CharField(max_length=30)),
            ],
            options={
                'db_table': 'expression',
            },
        ),
        migrations.CreateModel(
            name='Tissues',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('slug', models.SlugField(max_length=100, unique=True)),
                ('name', models.CharField(max_length=100)),
            ],
            options={
                'db_table': 'tissues',
            },
        ),
        migrations.CreateModel(
            name='TissueExpression',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('value', models.FloatField(blank=True, null=True)),
                ('protein', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='protein.protein')),
                ('tissue', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='protein.tissues')),
            ],
            options={
                'db_table': 'tissueexpression',
            },
        ),
        migrations.CreateModel(
            name='CancerExpression',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('cancer', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='protein.cancertype')),
                ('expression', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='protein.expressionvalue')),
            ],
            options={
                'db_table': 'cancer_expression',
            },
        ),
        migrations.AddField(
            model_name='protein',
            name='cancer',
            field=models.ManyToManyField(to='protein.CancerExpression'),
        ),
    ]
