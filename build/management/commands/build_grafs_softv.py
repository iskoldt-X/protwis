from django.core.management.base import BaseCommand
from django.conf import settings
from django.db import transaction
from django.utils.text import slugify

import os
import pandas as pd

from common.models import Modality
from protein.models import Chemotype, Sense, Protein


# =============================================================================
# GRAFS-SOFTV "Allowed Vocabulary" (Source of Truth)
# =============================================================================
# This is a RECEPTOR-centric classification (Protein sidecar), NOT ligand.LigandType.
# The build script is the enforcer: it will NOT create unknown Chemotypes from files.

ALLOWED_MODALITIES = [
    # 5 types (broad)
    {"name": "Small molecule", "slug": "small-molecule", "description": ""},
    {"name": "Peptide", "slug": "peptide", "description": ""},
    {"name": "Protein", "slug": "protein", "description": ""},
    {"name": "Ion", "slug": "ion", "description": ""},
    {"name": "Orphan", "slug": "orphan", "description": ""},
]

ALLOWED_SENSES = [
    # 6 types
    {"name": "Non-sensory", "slug": "non-sensory", "description": ""},
    {"name": "Odorant", "slug": "odorant", "description": ""},
    {"name": "Taste", "slug": "taste", "description": ""},
    {"name": "Vision", "slug": "vision", "description": ""},
    {"name": "Other sensory", "slug": "other-sensory", "description": ""},
    {"name": "Unknown", "slug": "unknown", "description": ""},
]

# 15 types (specific). Keys: chemotype slug; values include modality slug.
ALLOWED_CHEMOTYPES = [
    {"name": "Aminergic", "slug": "aminergic", "modality": "small-molecule", "description": ""},
    {"name": "Lipid", "slug": "lipid", "modality": "small-molecule", "description": ""},
    {"name": "Melatonin", "slug": "melatonin", "modality": "small-molecule", "description": ""},
    {"name": "Nucleotide", "slug": "nucleotide", "modality": "small-molecule", "description": ""},
    {"name": "Steroid", "slug": "steroid", "modality": "small-molecule", "description": ""},
    {"name": "Carboxylic acid", "slug": "carboxylic-acid", "modality": "small-molecule", "description": ""},
    {"name": "Amino acid", "slug": "amino-acid", "modality": "small-molecule", "description": ""},
    {"name": "Ion", "slug": "ion", "modality": "ion", "description": ""},
    {"name": "Peptide", "slug": "peptide", "modality": "peptide", "description": ""},
    {"name": "Protein", "slug": "protein", "modality": "protein", "description": ""},
    {"name": "Adhesion", "slug": "adhesion", "modality": "protein", "description": ""},
    {"name": "Sensory", "slug": "sensory", "modality": "small-molecule", "description": ""},
    {"name": "Orphan", "slug": "orphan", "modality": "orphan", "description": ""},
    {"name": "Other", "slug": "other", "modality": "small-molecule", "description": ""},
    {"name": "Fungal pheromone", "slug": "fungal-pheromone", "modality": "peptide", "description": ""},
]


def _norm(v) -> str:
    """Normalize values from spreadsheets for safe, strict matching."""
    if v is None:
        return ""
    s = str(v).strip()
    if s.lower() in {"nan", "none", ""}:
        return ""
    return s


def _build_lookup_by_name_and_slug(items):
    """
    Returns dict mapping normalized name/slug -> canonical slug.
    Accepts both exact names and slugs (case-insensitive) from ingestion files.
    """
    out = {}
    for it in items:
        nm = it["name"].strip()
        slug = it["slug"].strip()
        out[nm.casefold()] = slug
        out[slug.casefold()] = slug
        # also accept slugified name for convenience (still maps to canonical slug)
        out[slugify(nm).casefold()] = slug
    return out


class Command(BaseCommand):
    help = "Build receptor-centric GRAFS-SOFTV sidecar classification (Modality/Chemotype/Sense)"

    def add_arguments(self, parser):
        parser.add_argument(
            "--filename",
            action="store",
            dest="filename",
            default=None,
            help="Optional explicit source file path (xlsx/csv). Defaults to DATA_DIR/protein_data/GRAFS_SOFTV_Source.xlsx",
        )

    def handle(self, *args, **options):
        data_dir = os.path.join(settings.DATA_DIR, "protein_data")
        default_xlsx = os.path.join(data_dir, "GRAFS_SOFTV_Source.xlsx")
        default_csv = os.path.join(data_dir, "GRAFS_SOFTV_Source.csv")

        filename = options.get("filename") or (default_xlsx if os.path.exists(default_xlsx) else default_csv)
        file_exists = filename and os.path.exists(filename)

        # -------- Step 1: Ensure framework exists (Always run) --------
        self.stdout.write("##### STEP 1: Initializing allowed GRAFS-SOFTV vocabulary (strict schema) #####")

        modality_by_slug = {}
        sense_by_slug = {}
        chemotype_by_slug = {}

        with transaction.atomic():
            # Modalities
            for m in ALLOWED_MODALITIES:
                obj, _ = Modality.objects.update_or_create(
                    slug=m["slug"],
                    defaults={"name": m["name"], "description": m.get("description", "")},
                )
                modality_by_slug[obj.slug] = obj

            # Senses
            for s in ALLOWED_SENSES:
                obj, _ = Sense.objects.update_or_create(
                    slug=s["slug"],
                    defaults={"name": s["name"], "description": s.get("description", "")},
                )
                sense_by_slug[obj.slug] = obj

            # Chemotypes (and their modality mapping)
            for c in ALLOWED_CHEMOTYPES:
                modality_slug = c["modality"]
                modality = modality_by_slug.get(modality_slug)
                if modality is None:
                    raise RuntimeError(f"Invalid ALLOWED_CHEMOTYPES mapping: modality '{modality_slug}' not defined")

                obj, _ = Chemotype.objects.update_or_create(
                    slug=c["slug"],
                    defaults={
                        "name": c["name"],
                        "description": c.get("description", ""),
                        "modality": modality,
                    },
                )
                chemotype_by_slug[obj.slug] = obj

        self.stdout.write(
            f"Initialized/validated: {len(modality_by_slug)} modalities, {len(sense_by_slug)} senses, {len(chemotype_by_slug)} chemotypes."
        )

        # -------- Step 2/3: File ingestion + strict validation --------
        if not file_exists:
            self.stdout.write(
                "WARNING: No source file found. Framework initialized with default vocabulary. No proteins linked."
            )
            self.stdout.write(f"Expected: {default_xlsx} (or {default_csv})")
            return

        self.stdout.write("##### STEP 2: Ingesting source file (strict validation, no unknown vocabulary) #####")
        self.stdout.write(f"Reading: {filename}")

        try:
            if filename.lower().endswith(".xlsx"):
                df = pd.read_excel(filename)
            else:
                df = pd.read_csv(filename)
        except Exception as e:
            self.stderr.write(f"ERROR: Failed to read source file: {e}")
            self.stderr.write(
                "Framework was initialized, but no proteins were linked due to read failure."
            )
            return

        # Column resolution (robust)
        def pick_col(candidates):
            cols = {str(c).strip().casefold(): c for c in df.columns}
            for cand in candidates:
                key = cand.casefold()
                if key in cols:
                    return cols[key]
            return None

        col_entry = pick_col(["Protein Entry", "Entry", "entry_name", "Entry name", "Uniprot", "UniProt"])
        col_chemotype = pick_col(["Chemotype", "Receptor Chemotype"])
        col_sense = pick_col(["Sense", "Receptor Sense"])

        if col_entry is None:
            self.stderr.write("ERROR: Could not find a protein identifier column (e.g., 'Protein Entry' or 'entry_name').")
            self.stderr.write("Framework was initialized, but no proteins were linked.")
            return

        chemotype_lookup = _build_lookup_by_name_and_slug(ALLOWED_CHEMOTYPES)
        sense_lookup = _build_lookup_by_name_and_slug(ALLOWED_SENSES)

        linked = 0
        skipped = 0

        for idx, row in df.iterrows():
            entry = _norm(row.get(col_entry))
            if not entry:
                skipped += 1
                continue

            chemotype_raw = _norm(row.get(col_chemotype)) if col_chemotype else ""
            sense_raw = _norm(row.get(col_sense)) if col_sense else ""

            chemotype_slug = chemotype_lookup.get(chemotype_raw.casefold()) if chemotype_raw else None
            sense_slug = sense_lookup.get(sense_raw.casefold()) if sense_raw else None

            if chemotype_raw and chemotype_slug is None:
                self.stdout.write(
                    f"WARNING: row {idx}: unknown Chemotype '{chemotype_raw}' for protein '{entry}' (skipping; not creating)."
                )
                skipped += 1
                continue

            if sense_raw and sense_slug is None:
                self.stdout.write(
                    f"WARNING: row {idx}: unknown Sense '{sense_raw}' for protein '{entry}' (skipping; not creating)."
                )
                skipped += 1
                continue

            try:
                protein = Protein.objects.get(entry_name__iexact=entry)
            except Protein.DoesNotExist:
                self.stdout.write(f"WARNING: row {idx}: protein '{entry}' not found (skipping).")
                skipped += 1
                continue

            updated_fields = []
            if chemotype_slug:
                protein.chemotype = chemotype_by_slug[chemotype_slug]
                updated_fields.append("chemotype")
            if sense_slug:
                protein.sense = sense_by_slug[sense_slug]
                updated_fields.append("sense")

            if updated_fields:
                protein.save(update_fields=updated_fields)
                linked += 1
            else:
                # No classification info in row
                skipped += 1

        self.stdout.write(f"Done. Proteins linked: {linked}. Rows skipped: {skipped}.")

        # -------- Step 3 (Future): placeholder for richer mapping logic --------
        # TODO (Future):
        # - When Gaspar's curated file is available, extend ingestion to:
        #   - support additional columns (e.g., confidence, notes, overrides)
        #   - link receptor subfamilies / families at scale
        #   - track provenance (source file version, timestamp)


