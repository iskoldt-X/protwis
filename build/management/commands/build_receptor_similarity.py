from build.management.commands.base_build import Command as BaseBuild

from django.db.models import F, Q
from django.db.models.functions import Substr

from protein.models import Protein, ProteinSegment, ProteinFamily, Species, CLASSLESS_PARENT_GPCR_SLUGS
from alignment.models import ReceptorSimilarity

from common.alignment import Alignment

from collections import OrderedDict
import logging, os, sys, time, re
from datetime import datetime

import warnings, numpy as np
warnings.filterwarnings("ignore", category=RuntimeWarning, module="numpy")

logger = logging.getLogger('receptor_similarity')
hdlr = logging.FileHandler('./logs/receptor_similarity_to_db.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)

class_prefix_re = re.compile(r'^(Class)\s+', flags=re.I)

DEFAULT_BATCH_SIZE = 2000


class Command(BaseBuild):
    help = 'Build receptor similarity and identity and store them in alignment_receptorsimilarity (no CSV).'

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--verbose', default=False, action='store_true',
                            help='Print progress in stdout')
        parser.add_argument('--limit', type=int, default=False, action='store',
                            help='Use only this many GPCRs per class (testing).')
        # default None so we can apply DEFAULT_BATCH_SIZE in code
        parser.add_argument('--batch-size', type=int, default=None, action='store',
                            help='bulk_create batch size.')

    # ---------- helpers ----------
    def get_parent_gpcr_families(self, exclude_classless_artificial_class=True, include_classless_natural_classes=True):
        parent_family = ProteinFamily.objects.get(slug='000')
        parent_gpcr_families = ProteinFamily.objects.filter(
            parent_id=parent_family.pk, slug__startswith='0').exclude(pk=parent_family.pk)
        if exclude_classless_artificial_class:
            for slug in CLASSLESS_PARENT_GPCR_SLUGS:
                parent_gpcr_families = parent_gpcr_families.exclude(slug__startswith=slug)
        classless_protein_families = []
        if include_classless_natural_classes:
            classless_protein_families = self.get_classless_bottom_protein_families()
        parent_gpcr_families = list(parent_gpcr_families) + classless_protein_families
        return sorted(parent_gpcr_families, key=lambda f: (int(f.slug.split('_')[0])))

    def get_human_species(self):
        return Species.objects.get(common_name__iexact='Human')

    def get_yeast_species(self):
        return Protein.objects.filter(entry_name__iendswith='_yeast')[0].species

    def __slug_tree_branch(self, slug_parts, slug_tree_dict):
        if len(slug_parts) == 1:
            slug_tree_dict[slug_parts[0]] = None
            return slug_tree_dict
        slug_subtree_dict = slug_tree_dict.get(slug_parts[0]) or {}
        slug_tree_dict[slug_parts[0]] = self.__slug_tree_branch(slug_parts[1:], slug_subtree_dict)
        return slug_tree_dict

    def __sort_slug_tree_branch(self, slug_tree_dict):
        slug_tree_ordered_dict = OrderedDict()
        for slug in sorted(sorted(slug_tree_dict.keys(), key=lambda x: int(x[1:])), key=lambda x: x[0]):
            subtree = slug_tree_dict[slug]
            slug_tree_ordered_dict[slug] = self.__sort_slug_tree_branch(subtree) if subtree is not None else None
        return slug_tree_ordered_dict

    def __parse_slug_tree_(self, slug_tree_dict, slug_list_list, slug_list):
        for slug, subtree in slug_tree_dict.items():
            slug_list.append(slug)
            if subtree is not None:
                self.__parse_slug_tree_(subtree, slug_list_list, slug_list)
            else:
                slug_list_list.append(slug_list.copy())
            slug_list.pop()

    def get_classless_bottom_protein_families(self):
        classless_parent_gpcrs_slugs_list = sorted(sorted(CLASSLESS_PARENT_GPCR_SLUGS, key=lambda x: int(x[1:])), key=lambda x: x[0])
        parent_gpcr_families = ProteinFamily.objects.filter(slug__startswith=classless_parent_gpcrs_slugs_list[0])
        for slug in classless_parent_gpcrs_slugs_list[1:]:
            parent_gpcr_families = parent_gpcr_families.filter(slug__startswith=slug)

        slug_2_family = {f.slug: f for f in parent_gpcr_families}
        tree = {}
        for slug in slug_2_family.keys():
            self.__slug_tree_branch(slug.split('_'), tree)
        tree = self.__sort_slug_tree_branch(tree)
        slug_list_list, tmp = [], []
        self.__parse_slug_tree_(tree, slug_list_list, tmp)
        return [slug_2_family['_'.join(parts)] for parts in slug_list_list]

    def filter_out_non_species_parent_gpcr_families(self, parent_gpcr_families, species):
        new_slugs = set()
        species_list = species if isinstance(species, (list, tuple)) else [species]
        slugs = [f.slug for f in parent_gpcr_families]
        for slug in slugs:
            q = Protein.objects.annotate(family_slug=F('family__slug')).filter(
                family_slug__startswith=slug, species__in=species_list)
            if q.exists():
                new_slugs.add(slug)
        return [f for f in parent_gpcr_families if f.slug in new_slugs]

    def filter_out_non_human_parent_gpcr_families(self, parent_gpcr_families):
        return self.filter_out_non_species_parent_gpcr_families(parent_gpcr_families, self.get_human_species())

    def filter_out_non_yeast_parent_gpcr_families(self, parent_gpcr_families):
        return self.filter_out_non_species_parent_gpcr_families(parent_gpcr_families, self.get_yeast_species())
    # -----------------------------------------------------------

    def handle(self, *args, **options):
        verbose    = options['verbose']
        batch_size = options.get('batch_size') or DEFAULT_BATCH_SIZE
        initial_step1 = 380
        initial_step2 = 380

        start_time = time.time()
        if verbose: print('Truncating alignment_receptorsimilarity...')
        ReceptorSimilarity.custom_objects.truncate_table()

        parent_families = self.get_parent_gpcr_families(
            exclude_classless_artificial_class=True,
            include_classless_natural_classes=True
        )
        human_parent_gpcr_families = self.filter_out_non_human_parent_gpcr_families(parent_families)
        human_species = self.get_human_species()

        yeast_parent_gpcr_families = self.filter_out_non_yeast_parent_gpcr_families(parent_families)
        yeast_non_human = [f for f in yeast_parent_gpcr_families if f not in set(human_parent_gpcr_families)]
        yeast_species = self.get_yeast_species()

        gpcr_segments = ProteinSegment.objects.filter(
            Q(proteinfamily='GPCR') & (Q(slug__regex='TM[1-7]') | Q(slug='H8'))
        )

        # ===== cache top-level class families by slug prefix ('001', '002', ...) =====
        root = ProteinFamily.objects.get(slug='000')
        top_class_fams = ProteinFamily.objects.filter(parent_id=root.id, slug__regex=r'^\d{3}$')
        TOP_CLASS_BY_PREFIX = {pf.slug: pf for pf in top_class_fams}  # e.g. {'001': <PF Class A>, ...}

        def top_class_family_for(protein):
            """
            Return the top-level class ProteinFamily for a Protein,
            using the first chunk of its family.slug.
            """
            fam = protein.family
            if not fam or not fam.slug:
                return None
            prefix = fam.slug.split('_', 1)[0]  # '001_...' -> '001'
            return TOP_CLASS_BY_PREFIX.get(prefix)

        # ===== IMPORTANT: when building your per-class protein lists, add select_related('family') =====
        human_map, human_counts = {}, {}
        for fam in human_parent_gpcr_families:
            qs = (Protein.objects
                    .annotate(family_slug=F('family__slug'))
                    .filter(species=human_species, family_slug__startswith=fam.slug)
                    .exclude(accession=None)
                    .order_by('family_slug', 'entry_name')
                    .select_related('family'))  # <--- so we can use family.slug without extra queries
            human_map[fam] = list(qs)
            human_counts[fam] = len(qs)

        yeast_map, yeast_counts = {}, {}
        for fam in yeast_non_human:
            qs = (Protein.objects
                    .annotate(family_slug=F('family__slug'))
                    .filter(species=yeast_species, family_slug__startswith=fam.slug)
                    .exclude(accession=None)
                    .order_by('family_slug', 'entry_name')
                    .select_related('family'))  # <--- same here
            yeast_map[fam] = list(qs)
            yeast_counts[fam] = len(qs)

        selected_families = human_parent_gpcr_families + yeast_non_human
        selected_map, selected_counts = {}, {}
        for fam in selected_families:
            if fam in human_map and fam in yeast_map:
                selected_map[fam] = human_map[fam] + yeast_map[fam]
                selected_counts[fam] = human_counts[fam] + yeast_counts[fam]
            elif fam in human_map:
                selected_map[fam] = human_map[fam]
                selected_counts[fam] = human_counts[fam]
            else:
                selected_map[fam] = yeast_map[fam]
                selected_counts[fam] = yeast_counts[fam]

        # batching + duplicate prevention across all loops
        to_create = []
        seen_pairs = set()          # (min_id, max_id) to avoid duplicate rows
        unique_class_pairs = set()  # track cross-class pairs done, e.g., "A@B1"

        def flush_batch():
            nonlocal to_create
            if to_create:
                ReceptorSimilarity.objects.bulk_create(
                    to_create, batch_size=batch_size, ignore_conflicts=True
                )
                to_create.clear()

        step1 = int(initial_step1)
        step2 = int(initial_step2)
        step_halved = False

        for fam1 in selected_families:
            fam1_total = selected_counts[fam1]
            while_loop_continue = False
            while True:
                if options['limit'] and fam1_total > options['limit']:
                    fam1_total = options['limit']

                fam1_name = class_prefix_re.sub(r'', fam1.name.replace('<i>', '').replace('</i>', ''))

                for clim in range(0, fam1_total, step1):
                    block1 = selected_map[fam1][clim:clim+step1]
                    if options['limit']:
                        block1 = block1[:options['limit']]
                    # defensive clean
                    block1 = [p for p in block1 if p is not None]
                    if not block1:
                        continue

                    for fam2 in selected_families:
                        fam2_total = selected_counts[fam2]
                        if options['limit'] and fam2_total > options['limit']:
                            fam2_total = options['limit']

                        fam2_name = class_prefix_re.sub(r'', fam2.name.replace('<i>', '').replace('</i>', ''))
                        key = '@'.join(sorted([fam1_name, fam2_name]))
                        # only skip reverse if this is cross-class
                        if fam1 != fam2 and key in unique_class_pairs:
                            continue

                        for clim2 in range(0, fam2_total, step2):
                            # within-class: skip lower triangle block pairs
                            if fam1 == fam2 and clim2 < clim:
                                continue

                            if options['verbose']:
                                print(fam1, f"from:{clim+1} to:{clim+min(step1, fam1_total-clim)} (of:{fam1_total})",
                                      'vs', fam2, f"from:{clim2+1} to:{clim2+min(step2, fam2_total-clim2)} (of:{fam2_total})")

                            block2 = selected_map[fam2][clim2:clim2+step2]
                            if options['limit']:
                                block2 = block2[:options['limit']]
                            block2 = [p for p in block2 if p is not None]
                            if not block2:
                                continue

                            # avoid duplicate identical block work
                            if fam1 == fam2 and clim == clim2 and len(block1) == len(block2):
                                proteins = block1
                            else:
                                proteins = block1 + block2

                            cs_alignment = Alignment()
                            cs_alignment.load_proteins(proteins)
                            cs_alignment.load_segments(gpcr_segments)
                            r = cs_alignment.build_alignment()
                            if r == "Too large":
                                print('Alignment too large. Retrying...', file=sys.stderr)
                                while_loop_continue = True
                                break

                            cs_alignment.remove_non_generic_numbers_from_alignment()
                            cs_alignment.calculate_similarity_matrix()

                            pos_of = {p.protein.entry_name: i for i, p in enumerate(cs_alignment.proteins)}

                            for p1 in block1:
                                for p2 in block2:
                                    if p1.entry_name == p2.entry_name:
                                        continue
                                    a, b = p1.id, p2.id
                                    pair_key = (a, b) if a < b else (b, a)
                                    if pair_key in seen_pairs:
                                        continue

                                    # similarity (p1 vs p2) lives in [p2][pos_of[p1]]
                                    pos_s = pos_of[p1.entry_name]
                                    sim_val = int(cs_alignment.similarity_matrix[p2.entry_name]['values'][pos_s][0])
                                    # identity (p1 vs p2) lives in [p1][pos_of[p2]]
                                    pos_i = pos_of[p2.entry_name]
                                    id_val = int(cs_alignment.similarity_matrix[p1.entry_name]['values'][pos_i][0])

                                    ref, tgt = (p1, p2) if a < b else (p2, p1)

                                    # compute top-level classes from family.slug prefix map
                                    ref_cls = top_class_family_for(ref)
                                    tgt_cls = top_class_family_for(tgt)

                                    to_create.append(ReceptorSimilarity(
                                        protein_ref=ref,
                                        protein_target=tgt,
                                        identity=id_val,
                                        similarity=sim_val,
                                        ref_class=ref_cls,
                                        target_class=tgt_cls,
                                    ))
                                    seen_pairs.add(pair_key)

                                    if len(to_create) >= batch_size:
                                        flush_batch()

                            flush_batch()

                        # finished all blocks for this cross-class pair → mark it done
                        if not while_loop_continue and fam1 != fam2:
                            unique_class_pairs.add(key)

                        if while_loop_continue:
                            break

                    # memory trim only if not retrying
                    if not while_loop_continue:
                        try:
                            selected_map[fam1][clim:clim+step1] = [None] * (min(step1, fam1_total - clim))
                        except Exception:
                            pass

                    if while_loop_continue:
                        break

                if while_loop_continue:
                    if options['verbose']: print("Halving step1 and step2...")
                    step1 //= 2
                    step2 //= 2
                    step_halved = True
                    while_loop_continue = False
                    if options['verbose']: print("Retrying last class similarity computation with the new steps...")
                    continue
                elif step_halved:
                    if options['verbose']: print("Restoring initial step1 and step2...")
                    step1 = initial_step1
                    step2 = initial_step2
                    step_halved = False

                # free whole fam1 list—safe now because reverse pairs are skipped
                del selected_map[fam1]
                break

        flush_batch()

        elapsed = time.time() - start_time
        if verbose:
            print(f'Inserted receptor similarities in {elapsed:.2f}s')
        logger.info(f'Inserted receptor similarities in {elapsed:.2f}s')
