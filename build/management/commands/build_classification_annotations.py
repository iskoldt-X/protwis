from collections import defaultdict

from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from django.utils.text import slugify

import xlrd

from protein.models import (
    Gene,
    Protein,
    ProteinFamilyClassification,
    ProteinFamilyClassificationChemotype,
    ProteinFamilyClassificationModality,
    ProteinFamilyClassificationSense,
)


class Command(BaseCommand):
    help = "Import GPCR family classification annotations from Excel."

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._level4_cache = {}

    def add_arguments(self, parser):
        parser.add_argument(
            "-f",
            "--filename",
            action="store",
            dest="filename",
            help="Path to Classification.xlsx",
        )
        parser.add_argument(
            "-s",
            "--sheet",
            action="store",
            dest="sheet",
            default="Classification",
            help="Worksheet name to import",
        )

    def handle(self, *args, **options):
        filename = options.get("filename") or self.default_filename()
        sheet_name = options.get("sheet") or "Classification"

        data = self.parse_excel(filename, remove_header_linebreak=True)
        if sheet_name not in data:
            raise CommandError("Worksheet not found: {}".format(sheet_name))

        rows = data[sheet_name]
        self.stdout.write("Using file: {} (sheet: {})".format(filename, sheet_name))

        family_data, stats, seen_family_ids = self.aggregate_rows(rows)
        self.apply_updates(family_data)
        self.prune_relationships(seen_family_ids, stats)

        self.stdout.write("Rows processed: {}".format(stats["rows_total"]))
        self.stdout.write("Gene matches: {}".format(stats["gene_matches"]))
        self.stdout.write("Entry short matches: {}".format(stats["entry_matches"]))
        self.stdout.write("Ambiguous skipped: {}".format(stats["ambiguous_skipped"]))
        self.stdout.write("Unmatched rows: {}".format(stats["unmatched"]))
        self.stdout.write("Families updated: {}".format(len(family_data)))
        self.stdout.write("Warnings: {}".format(stats["warnings"]))

    def default_filename(self):
        return settings.DATA_DIR + "/protein_data/Classification.xlsx"

    def aggregate_rows(self, rows):
        family_data = {}
        stats = defaultdict(int)
        seen_family_ids = set()

        for row in rows.values():
            stats["rows_total"] += 1
            gene_name = self.value_from_row(
                row, ["GPCRs (Gene name)", "GPCRs (Gene name)"]
            )
            uniprot_short = self.value_from_row(
                row, ["GPCRs (UniProt)", "GPCRs (UniProt)"]
            )

            protein, method, ambiguous, candidates = self.resolve_protein(
                gene_name, uniprot_short
            )
            if ambiguous:
                stats["warnings"] += 1
                stats["ambiguous_skipped"] += 1
                self.stdout.write(
                    "WARNING: ambiguous gene mapping for {} (candidates: {}).".format(
                        gene_name or "n/a", ", ".join(candidates) or "n/a"
                    )
                )
                continue
            if not protein:
                stats["unmatched"] += 1
                continue
            if method == "gene":
                stats["gene_matches"] += 1
            else:
                stats["entry_matches"] += 1

            family = self.resolve_level4_family(protein.family)
            if not family:
                stats["unmatched"] += 1
                continue
            seen_family_ids.add(family.id)

            if family.id not in family_data:
                family_data[family.id] = {
                    "family": family,
                    "senses": set(),
                    "chemotypes": {},
                    "modalities": {},
                }

            sense_value = self.value_from_row(row, ["Sense"])
            if sense_value:
                family_data[family.id]["senses"].add(sense_value)

            chemotype_value = self.value_from_row(row, ["Chemotype"])
            chemotype_order = self.value_from_row(row, ["Chemotype order"])
            self.apply_ordered_value(
                family_data[family.id]["chemotypes"],
                chemotype_order,
                chemotype_value,
                family,
                "chemotype",
                stats,
            )

            modality_value = self.value_from_row(row, ["Modality"])
            modality_order = self.value_from_row(row, ["Modality order"])
            self.apply_ordered_value(
                family_data[family.id]["modalities"],
                modality_order,
                modality_value,
                family,
                "modality",
                stats,
            )

        return family_data, stats, seen_family_ids

    def apply_updates(self, family_data):
        with transaction.atomic():
            for data in family_data.values():
                family = data["family"]

                # Sense logic: Non-sensory overrides all; Unknown is ignored if other values exist
                senses_to_apply = data["senses"]
                if "Non-sensory" in senses_to_apply:
                    senses_to_apply = {"Non-sensory"}
                elif "Unknown" in senses_to_apply and len(senses_to_apply) > 1:
                    senses_to_apply.remove("Unknown")

                sense_obj = None
                if senses_to_apply:
                    sense_name = sorted(senses_to_apply)[0]
                    sense_obj = self.get_or_create_vocab(
                        ProteinFamilyClassificationSense, sense_name
                    )

                orders = set(data["chemotypes"].keys()) | set(data["modalities"].keys())

                # Replace rows for this family to keep grouped-by-order layout.
                ProteinFamilyClassification.objects.filter(
                    protein_family=family
                ).delete()

                if orders:
                    for order_num in sorted(orders):
                        chemotype_obj = None
                        modality_obj = None
                        if order_num in data["chemotypes"]:
                            chemotype_obj = self.get_or_create_vocab(
                                ProteinFamilyClassificationChemotype,
                                data["chemotypes"][order_num],
                            )
                        if order_num in data["modalities"]:
                            modality_obj = self.get_or_create_vocab(
                                ProteinFamilyClassificationModality,
                                data["modalities"][order_num],
                            )
                        ProteinFamilyClassification.objects.get_or_create(
                            protein_family=family,
                            sense=sense_obj,
                            chemotype=chemotype_obj,
                            chemotype_order=order_num if chemotype_obj else None,
                            modality=modality_obj,
                            modality_order=order_num if modality_obj else None,
                        )
                elif sense_obj:
                    ProteinFamilyClassification.objects.get_or_create(
                        protein_family=family,
                        sense=sense_obj,
                    )

    def prune_relationships(self, seen_family_ids, stats):
        stale_rows = (
            ProteinFamilyClassification.objects.exclude(
                protein_family_id__in=seen_family_ids
            )
            if seen_family_ids
            else ProteinFamilyClassification.objects.all()
        )

        removed_count = stale_rows.count()

        if removed_count:
            stale_rows.delete()

        stats["pruned_total"] = removed_count

        self.stdout.write("Pruned stale annotations: total={}".format(removed_count))

    def apply_ordered_value(
        self, storage, order_value, name_value, family, label, stats
    ):
        if not name_value:
            return
        try:
            order = int(order_value)
        except (TypeError, ValueError):
            stats["warnings"] += 1
            self.stdout.write(
                "WARNING: invalid {} order for {}.".format(label, family.slug)
            )
            return
        if order not in (1, 2):
            stats["warnings"] += 1
            self.stdout.write(
                "WARNING: {} order out of range for {}.".format(label, family.slug)
            )
            return

        if order in storage and storage[order] != name_value:
            stats["warnings"] += 1
            self.stdout.write(
                "WARNING: conflicting {} order {} for {} (kept {}).".format(
                    label, order, family.slug, storage[order]
                )
            )
            return
        storage[order] = name_value

    def resolve_protein(self, gene_name, uniprot_short):
        if gene_name:
            genes = Gene.objects.filter(
                name__iexact=gene_name, species__common_name="Human"
            )
            proteins = (
                Protein.objects.filter(genes__in=genes, species__common_name="Human")
                .distinct()
                .order_by("id")
            )
            protein_count = proteins.count()
            if protein_count == 1:
                return proteins.first(), "gene", False, []
            if protein_count > 1:
                if uniprot_short:
                    short = uniprot_short.strip().lower()
                    short_matches = proteins.filter(
                        entry_name__istartswith="{}_".format(short)
                    )
                    if short_matches.count() == 1:
                        return short_matches.first(), "gene+entry", False, []
                primary_genes = genes.filter(position=0)
                primary_proteins = (
                    Protein.objects.filter(
                        genes__in=primary_genes, species__common_name="Human"
                    )
                    .distinct()
                    .order_by("id")
                )
                if primary_proteins.count() == 1:
                    return primary_proteins.first(), "gene+primary", False, []
                return (
                    None,
                    "gene",
                    True,
                    list(proteins.values_list("entry_name", flat=True)),
                )
            # No human match, fall back to any species
            genes_any = Gene.objects.filter(name__iexact=gene_name)
            proteins_any = (
                Protein.objects.filter(genes__in=genes_any).distinct().order_by("id")
            )
            any_count = proteins_any.count()
            if any_count == 1:
                return proteins_any.first(), "gene_any", False, []
            if any_count > 1:
                if uniprot_short:
                    short = uniprot_short.strip().lower()
                    short_matches = proteins_any.filter(
                        entry_name__istartswith="{}_".format(short)
                    )
                    if short_matches.count() == 1:
                        return short_matches.first(), "gene_any+entry", False, []
                primary_genes_any = genes_any.filter(position=0)
                primary_proteins_any = (
                    Protein.objects.filter(genes__in=primary_genes_any)
                    .distinct()
                    .order_by("id")
                )
                if primary_proteins_any.count() == 1:
                    return primary_proteins_any.first(), "gene_any+primary", False, []
                return (
                    None,
                    "gene_any",
                    True,
                    list(proteins_any.values_list("entry_name", flat=True)),
                )

        if uniprot_short:
            short = uniprot_short.strip().lower()
            proteins = Protein.objects.filter(
                entry_name__istartswith="{}_".format(short),
                species__common_name="Human",
            ).order_by("id")
            if proteins.exists():
                return proteins.first(), "entry", False, []

        return None, None, False, []

    def resolve_level4_family(self, family):
        if not family:
            return None
        if family.id in self._level4_cache:
            return self._level4_cache[family.id]

        lineage = []
        current = family
        while current:
            lineage.append(current)
            current = current.parent

        chain = list(reversed(lineage))
        if len(chain) < 5:
            resolved = None
        else:
            resolved = chain[4]

        self._level4_cache[family.id] = resolved
        return resolved

    def get_or_create_vocab(self, model, name):
        slug = slugify(name)
        obj, created = model.objects.get_or_create(slug=slug, defaults={"name": name})
        if not created and obj.name != name:
            obj.name = name
            obj.save(update_fields=["name"])
        return obj

    def value_from_row(self, row, keys):
        for key in keys:
            if key in row and row[key] != "":
                return row[key]
        return None

    def parse_excel(self, path, remove_header_linebreak=False):
        workbook = xlrd.open_workbook(path)
        worksheets = workbook.sheet_names()
        data = {}
        for worksheet_name in worksheets:
            if worksheet_name in data:
                continue

            data[worksheet_name] = {}
            worksheet = workbook.sheet_by_name(worksheet_name)
            num_rows = worksheet.nrows - 1
            num_cells = worksheet.ncols

            headers = []
            for col in range(num_cells):
                header = worksheet.cell_value(0, col)
                if header == "":
                    header = "i_{}".format(col)
                if header in headers:
                    header += "_{}".format(col)
                if remove_header_linebreak:
                    header = header.replace("\n", " ")
                header = header.strip()
                headers.append(header)

            for row in range(1, num_rows + 1):
                key = worksheet.cell_value(row, 0)
                if key == "":
                    continue
                row_dict = {}
                for col in range(num_cells):
                    row_dict[headers[col]] = worksheet.cell_value(row, col)
                data[worksheet_name][key] = row_dict
        return data
