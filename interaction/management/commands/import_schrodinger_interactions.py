"""Import Schrodinger Engine 1 interaction YAMLs, one transaction per structure.

Usage::

    python manage.py import_schrodinger_interactions \\
        --data-dir /app/data/schrodinger --pdb 2RH1 --pdb 6CM4 \\
        --anomaly-csv /runs/anomalies.csv

Each structure is imported in its own transaction. A structure that fails
(unreadable YAML, a row the type map cannot route, or any unexpected error)
is rolled back and reported; the others are unaffected. The command exits
non-zero when any structure failed, after all structures have been attempted.

An anchor with no matching product instance keeps its existing rows and is
reported as a WARNING (anchor_not_replaced).

The anomaly CSV is written outside the transactions and flushed per row, so
it survives any rollback. Every row that was read but not written is
accounted for in it.
"""

import csv
import datetime
import json
import os

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from interaction import schrodinger_import as si
from interaction.models import ResidueFragmentInteractionType
from structure.models import Structure


class _Rollback(Exception):
    """Raised inside a dry-run transaction to undo it."""


class AnomalyLog(object):
    """CSV sink for per-anchor accounting; independent of the database."""

    COLUMNS = ["timestamp", "pdb", "sli_id", "het", "level", "category", "count", "detail"]

    def __init__(self, path):
        self._fh = open(path, "w", newline="")
        self._writer = csv.writer(self._fh)
        self._writer.writerow(self.COLUMNS)
        self._fh.flush()
        self.levels = {"INFO": 0, "WARNING": 0, "ERROR": 0}

    def log(self, pdb, level, category, sli_id="", het="", count=1, detail=""):
        self.levels[level] += 1
        self._writer.writerow([
            datetime.datetime.utcnow().isoformat(timespec="seconds"),
            pdb, sli_id, het, level, category, count, detail,
        ])
        self._fh.flush()

    def close(self):
        self._fh.close()


class Command(BaseCommand):
    help = "Import Schrodinger Engine 1 interactions (one transaction per structure)."

    def add_arguments(self, parser):
        parser.add_argument("--data-dir", required=True,
                            help="Root of the product tree: {data_dir}/{PDB}/{instance}/.")
        parser.add_argument("--pdb", action="append", default=[],
                            help="PDB code to import; repeatable.")
        parser.add_argument("--pdb-list", default=None,
                            help="File with one PDB code per line (# comments allowed).")
        parser.add_argument("--anomaly-csv", required=True,
                            help="Where to write the per-anchor accounting CSV.")
        parser.add_argument("--report-json", default=None,
                            help="Optional path for a machine-readable per-anchor report.")
        parser.add_argument("--dry-run", action="store_true",
                            help="Run every structure and roll each one back.")

    def _pdb_codes(self, options):
        codes = [c.strip().upper() for c in options["pdb"] if c.strip()]
        if options["pdb_list"]:
            with open(options["pdb_list"]) as fh:
                for line in fh:
                    line = line.split("#", 1)[0].strip()
                    if line:
                        codes.append(line.upper())
        if not codes:
            raise CommandError("no PDB codes given (use --pdb or --pdb-list)")
        return list(dict.fromkeys(codes))

    def _check_slugs(self):
        present = set(ResidueFragmentInteractionType.objects.values_list("slug", flat=True))
        missing = sorted(si.required_slugs() - present)
        if missing:
            raise CommandError(
                "interaction types missing from the database: {} "
                "(run migrate; interaction 0008 seeds them)".format(", ".join(missing)))

    def handle(self, *args, **options):
        codes = self._pdb_codes(options)
        if not os.path.isdir(options["data_dir"]):
            raise CommandError("--data-dir {!r} is not a directory".format(options["data_dir"]))
        self._check_slugs()
        log = AnomalyLog(options["anomaly_csv"])
        report = []
        failed = []
        totals = {"anchors": 0, "untouched": 0, "deleted": 0, "written": 0}
        try:
            for pdb in codes:
                entry = {"pdb": pdb}
                report.append(entry)
                structure = (Structure.objects
                             .filter(pdb_code__index__iexact=pdb)
                             .select_related("structure_type", "pdb_code", "protein_conformation")
                             .first())
                if structure is None:
                    log.log(pdb, "WARNING", "structure_not_in_db")
                    entry["status"] = "structure_not_in_db"
                    continue
                if structure.structure_type.origin != si.STRUCTURE_ORIGIN:
                    log.log(pdb, "INFO", "not_experimental",
                            detail=structure.structure_type.slug)
                    entry["status"] = "not_experimental"
                    continue
                if not os.path.isdir(os.path.join(options["data_dir"], pdb)):
                    log.log(pdb, "WARNING", "no_product_dir",
                            detail="every in-scope anchor of this structure is left untouched")
                try:
                    with transaction.atomic():
                        outcomes, out_of_scope = si.import_structure(
                            structure, options["data_dir"])
                        if options["dry_run"]:
                            raise _Rollback()
                except _Rollback:
                    pass
                except Exception as exc:
                    # The structure's transaction has been rolled back; record
                    # the failure and go on with the next structure.
                    message = "{}: {}".format(type(exc).__name__, exc)[:300]
                    log.log(pdb, "ERROR", type(exc).__name__, detail=message)
                    entry["status"] = "failed"
                    entry["error"] = message
                    failed.append(pdb)
                    continue
                entry["status"] = "rolled_back" if options["dry_run"] else "imported"
                entry["out_of_scope_anchors"] = out_of_scope
                entry["anchors"] = []
                for o in outcomes:
                    self._log_outcome(log, pdb, o)
                    entry["anchors"].append({
                        "sli_id": o.sli_id, "het": o.het, "mode": o.mode,
                        "instances": o.instances, "missing": o.missing,
                        "deleted": o.deleted, "written": o.written,
                        "fragments_created": o.fragments_created, "counts": dict(o.counts),
                        "other_chain_by_chain": o.other_chain_by_chain,
                        "dropped": dict(o.dropped),
                    })
                    totals["anchors"] += 1
                    totals["untouched"] += o.mode in si.UNTOUCHED_MODES
                    totals["deleted"] += o.deleted
                    totals["written"] += o.written
        finally:
            log.close()
            if options["report_json"]:
                with open(options["report_json"], "w") as fh:
                    json.dump({"dry_run": options["dry_run"], "totals": totals,
                               "failed": failed, "structures": report}, fh, indent=1)

        self.stdout.write(
            "{} structures, {} in-scope anchors ({} left untouched), {} RFI rows deleted, "
            "{} written{}; anomalies INFO={} WARNING={} ERROR={}".format(
                len(codes), totals["anchors"], totals["untouched"], totals["deleted"],
                totals["written"],
                " (dry run: rolled back)" if options["dry_run"] else "",
                log.levels["INFO"], log.levels["WARNING"], log.levels["ERROR"]))
        if failed:
            raise CommandError("{} structure(s) failed and were rolled back: {}".format(
                len(failed), ", ".join(failed)))

    @staticmethod
    def _log_outcome(log, pdb, o):
        c = o.counts
        if o.mode in si.UNTOUCHED_MODES:
            log.log(pdb, "WARNING", "anchor_not_replaced", o.sli_id, o.het,
                    detail="{}; existing rows kept; missing={}".format(
                        o.mode, ",".join(o.missing)))
            return
        if o.mode == "exact_partial":
            log.log(pdb, "WARNING", "anchor_instances_partial", o.sli_id, o.het,
                    len(o.missing), detail="missing={}".format(",".join(o.missing)))
        if c["rows_in"] == 0:
            log.log(pdb, "INFO", "product_has_zero_rows", o.sli_id, o.het,
                    detail=",".join(o.instances))
        for category, level in (("excluded_family", "INFO"),
                                ("nonstandard_residue", "INFO"),
                                ("duplicate", "INFO"),
                                ("other_chain", "WARNING")):
            if c[category]:
                detail = ""
                if category == "other_chain":
                    detail = ",".join("{}:{}".format(k, v)
                                      for k, v in sorted(o.other_chain_by_chain.items()))
                log.log(pdb, level, category, o.sli_id, o.het, c[category], detail)
        for category, n in sorted(o.dropped.items()):
            log.log(pdb, "WARNING", category, o.sli_id, o.het, n)
