from django.db import models, connection

# Create your models here.

class AlignmentConsensus(models.Model):
    slug = models.SlugField(max_length=100, unique=True)
    alignment = models.BinaryField()
    gn_consensus = models.BinaryField(blank=True) # Store conservation calculation for each GN


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


