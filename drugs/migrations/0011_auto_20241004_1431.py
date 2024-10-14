# Generated by Django 3.2.25 on 2024-10-04 12:31

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('drugs', '0010_auto_20240220_2019'),
    ]

    operations = [
        migrations.AddField(
            model_name='drugs2024',
            name='indication_status',
            field=models.CharField(max_length=15, null=True),
        ),
        migrations.AlterField(
            model_name='drugs2024',
            name='approval_year',
            field=models.IntegerField(null=True),
        ),
        migrations.AlterField(
            model_name='drugs2024',
            name='indication_max_phase',
            field=models.IntegerField(null=True),
        ),
    ]
