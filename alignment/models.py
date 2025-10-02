from django.db import models, connection
from enum import IntEnum

# Create your models here.

class AlignmentConsensus(models.Model):
    slug = models.SlugField(max_length=100, unique=True)
    alignment = models.BinaryField()
    gn_consensus = models.BinaryField(blank=True) # Store conservation calculation for each GN


class CustomClassSimilarityManager(models.Manager):
    def truncate_table(self):
        cursor = connection.cursor()
        table_name = self.model._meta.db_table
        sql = 'TRUNCATE TABLE "{0}" CASCADE'.format(table_name)
        cursor.execute(sql)

class ClassSimilarity(models.Model):
    protein_family1 = models.ForeignKey('protein.ProteinFamily',null=False, on_delete=models.CASCADE, related_name='class_similarity_protein_family1')
    protein_family2 = models.ForeignKey('protein.ProteinFamily',null=False, on_delete=models.CASCADE, related_name='class_similarity_protein_family2')
    similar_protein1 = models.ForeignKey('protein.Protein',null=True, on_delete=models.CASCADE, related_name='class_similarity_similar_protein1')
    similar_protein2 = models.ForeignKey('protein.Protein',null=True, on_delete=models.CASCADE, related_name='class_similarity_similar_protein2')
    ident_protein1 = models.ForeignKey('protein.Protein',null=True, on_delete=models.CASCADE, related_name='class_similarity_ident_protein1')
    ident_protein2 = models.ForeignKey('protein.Protein',null=True, on_delete=models.CASCADE, related_name='class_similarity_ident_protein2')
    similarity = models.IntegerField(null=False)
    identity = models.IntegerField(null=False)

    objects = models.Manager()  # The default manager.
    custom_objects = CustomClassSimilarityManager()  # The custom manager.
    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['protein_family1', 'protein_family2'], name='unique_class_similarity_protein_family'),
            models.UniqueConstraint(fields=['similar_protein1', 'similar_protein2'], name='unique_class_similarity_similar_proteins'),
            models.UniqueConstraint(fields=['ident_protein1', 'ident_protein2'], name='unique_class_similarity_ident_proteins'),
        ]
    def __str__(self):
        return str(self.protein_family1)+" - "+str(self.protein_family2)+": "+str(self.similarity)

class CustomClassRepresentativeSpeciesManager(models.Manager):
    def truncate_table(self):
        cursor = connection.cursor()
        table_name = self.model._meta.db_table
        sql = 'TRUNCATE TABLE "{0}"'.format(table_name)
        cursor.execute(sql)

class ClassRepresentativeSpecies(models.Model):
    protein_family = models.OneToOneField('protein.ProteinFamily',null=False, on_delete=models.CASCADE, related_name='class_representative_protein_family')
    species = models.ForeignKey('protein.Species', on_delete=models.CASCADE, related_name='class_representative_species')

    objects = models.Manager()  # The default manager.
    custom_objects = CustomClassRepresentativeSpeciesManager()  # The custom manager.

class ClassSimilarityType(IntEnum):
    IDENTITY = 0
    SIMILARITY = 1
    @classmethod
    def choices(cls):
        return [(key.value, key.name) for key in cls]

class ClassSimilarityTie(models.Model):
    class_similarity = models.ForeignKey('alignment.ClassSimilarity',null=False, on_delete=models.CASCADE, related_name='class_similarity_tie_class_similarity')
    protein1 = models.ForeignKey('protein.Protein',null=True, on_delete=models.CASCADE, related_name='class_similarity_tie_protein1')
    protein2 = models.ForeignKey('protein.Protein',null=True, on_delete=models.CASCADE, related_name='class_similarity_tie_protein2')
    type = models.IntegerField(choices=ClassSimilarityType.choices(), default=int(ClassSimilarityType.SIMILARITY))
    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['protein1', 'protein2','type'], name='unique_class_similarity_tie_proteins_type'),   
        ]
    def __str__(self):
        key2value = {}
        for key,value in ClassSimilarityType.choices():
            key2value[key] = value
        return str(self.protein1)+" vs "+str(self.protein2)+": "+key2value[self.type]

class CustomReceptorSimilarityManager(models.Manager):
    def truncate_table(self):
        with connection.cursor() as cursor:
            cursor.execute(f'TRUNCATE TABLE "{self.model._meta.db_table}" CASCADE')

class ReceptorSimilarity(models.Model):
    # existing
    protein_ref    = models.ForeignKey(
        'protein.Protein',
        on_delete=models.CASCADE,
        related_name='receptor_similarity_as_ref',
        db_column='ref',
        db_index=True,
    )
    protein_target = models.ForeignKey(
        'protein.Protein',
        on_delete=models.CASCADE,
        related_name='receptor_similarity_as_target',
        db_column='target',
        db_index=True,
    )
    identity   = models.PositiveSmallIntegerField()
    similarity = models.PositiveSmallIntegerField()

    # NEW: top-level class FKs (nullable for backfill)
    ref_class = models.ForeignKey(
        'protein.ProteinFamily',
        null=True, blank=True,
        on_delete=models.CASCADE,
        related_name='sim_as_ref_class',
        db_column='ref_class',
        db_index=True,
    )
    target_class = models.ForeignKey(
        'protein.ProteinFamily',
        null=True, blank=True,
        on_delete=models.CASCADE,
        related_name='sim_as_target_class',
        db_column='target_class',
        db_index=True,
    )

    objects = models.Manager()
    custom_objects = CustomReceptorSimilarityManager()

    class Meta:
        db_table = 'alignment_receptorsimilarity'
        constraints = [
            models.UniqueConstraint(fields=['protein_ref', 'protein_target'],
                                    name='uniq_receptor_similarity_pair'),
            models.CheckConstraint(check=~models.Q(protein_ref=models.F('protein_target')),
                                   name='check_ref_ne_target'),
        ]
        indexes = [
            # Forward order
            models.Index(fields=['ref_class', 'target_class', 'identity'],   name='rs_cls_id_idx'),
            models.Index(fields=['ref_class', 'target_class', 'similarity'], name='rs_cls_sim_idx'),
            # Reverse order (helps the (b,a) half of your OR)
            models.Index(fields=['target_class', 'ref_class', 'identity'],   name='rs_cls_id_rev_idx'),
            models.Index(fields=['target_class', 'ref_class', 'similarity'], name='rs_cls_sim_rev_idx'),
        ]

    def __str__(self):
        return f'{self.protein_ref} vs {self.protein_target}: sim={self.similarity} id={self.identity}'


