from django.db import models


class Drugs(models.Model):
    target = models.ManyToManyField('protein.Protein')
    name = models.CharField(max_length=100)
    synonym = models.CharField(max_length=300, null=True) # addName as list, make as ManyField
    drugtype = models.CharField(max_length=40, null=True)
    status = models.CharField(max_length=15)
    indication = models.CharField(max_length=150, null=True)
    approval = models.CharField(max_length=5)
    clinicalstatus = models.CharField(max_length=30, null=True, default=False)
    phasedate = models.CharField(max_length=15, null=True, default=True)
    phase = models.CharField(max_length=1, null=True, default=True)
    moa = models.CharField(max_length=30, null=True, default=True)
    targetlevel = models.CharField(max_length=10, null=True, default='NA')
    externallink = models.CharField(max_length=150, null=True) # Link Framework
    novelty = models.CharField(max_length=15) #Boolean
    references = models.CharField(max_length=180) #Boolean
    publication = models.ManyToManyField('common.Publication')


    def __str__(self):
        return self.name

    class Meta():
        db_table = 'drugs'


class Drugs2024(models.Model):
    target = models.ForeignKey('protein.Protein', on_delete=models.CASCADE)
    indication = models.ForeignKey('Indication', on_delete=models.CASCADE)
    # indication = models.ManyToManyField('Indication')
    ligand = models.ForeignKey('ligand.Ligand', on_delete=models.CASCADE)
    charge = models.CharField(max_length=5, null=True)
    complexity = models.FloatField(max_length=4, null=True)
    tpsa = models.CharField(max_length=10, null=True)
    drug_status = models.CharField(max_length=15, null=True)
    approval_year = models.IntegerField(null=True)
    indication_max_phase = models.IntegerField(null=True)
    indication_status = models.CharField(max_length=15, null=True)
    moa = models.ForeignKey('ligand.LigandRole', on_delete=models.CASCADE, null=True)
    genetic_association = models.CharField(max_length=30, null=True)
    affected_pathway = models.CharField(max_length=30, null=True)
    somatic_mutation = models.CharField(max_length=30, null=True)
    similarity_to_model = models.FloatField(max_length=4, null=True)
    novelty_score = models.FloatField(max_length=4, null=True)
    reference = models.ManyToManyField('common.Publication')

    # TO BE ADDED
    # IndicationStatus
    # ICD name
    # ICD code (website https://icd.who.int/browse/2024-01/mms/en#821852937)
    # ATC code (website https://atcddd.fhi.no/atc_ddd_index/?code=N06B)

    def __str__(self):
        return self.ligand.name

    class Meta():
        db_table = 'drugs_new'

#Do we need to fix this model structure?
class Indication(models.Model):
    name =  models.CharField(max_length=70)
    web_link =  models.ForeignKey('common.WebLink', on_delete=models.CASCADE, null=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'indication'

# class CancerExpression(models.Model):
#     cancer = models.ForeignKey('CancerType', on_delete=models.CASCADE)
#     expression = models.ForeignKey('ExpressionValue', on_delete=models.CASCADE)
#
#     def __str__(self):
#         return self.cancer.slug
#
#     class Meta():
#         db_table = 'cancer_expression'
