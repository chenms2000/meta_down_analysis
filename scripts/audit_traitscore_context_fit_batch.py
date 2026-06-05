"""Batch audit TraitScore CSV context-fit behavior.

This script runs the read-only evidence service over TraitScore group-difference
CSV files, using structured cancer/tissue/cell/comparison/goal context inferred
from the filename. Outputs are derived validation reports only; they do not
write graph facts or change release artifacts.
"""

from __future__ import annotations

import argparse
import csv
import html
import importlib.util
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SERVICE_PATH = ROOT / "scripts" / "metabo_service.py"
SPEC = importlib.util.spec_from_file_location("metabo_service_batch", SERVICE_PATH)
metabo_service = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = metabo_service
SPEC.loader.exec_module(metabo_service)


CONTEXTS: dict[str, dict[str, str]] = {
    "CESC": {
        "cancer_type": "cervical cancer",
        "tissue": "cervix",
        "cell_type": "epithelial",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "CESC",
    },
    "CRC": {
        "cancer_type": "colorectal cancer",
        "tissue": "colon and rectum",
        "cell_type": "epithelial",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "CRC COAD READ",
    },
    "cSCC": {
        "cancer_type": "cutaneous squamous cell carcinoma",
        "tissue": "skin",
        "cell_type": "epithelial keratinocyte",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "cSCC",
    },
    "ESCA": {
        "cancer_type": "esophageal cancer",
        "tissue": "esophagus",
        "cell_type": "epithelial",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "ESCA",
    },
    "HCC": {
        "cancer_type": "hepatocellular carcinoma",
        "tissue": "liver",
        "cell_type": "epithelial hepatocyte",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "HCC",
    },
    "HNSC": {
        "cancer_type": "head and neck squamous cell carcinoma",
        "tissue": "head and neck",
        "cell_type": "epithelial",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "HNSC",
    },
    "ICC": {
        "cancer_type": "intrahepatic cholangiocarcinoma",
        "tissue": "liver bile duct",
        "cell_type": "epithelial cholangiocyte",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "ICC cholangiocarcinoma",
    },
    "KICH": {
        "cancer_type": "chromophobe renal cell carcinoma",
        "tissue": "kidney",
        "cell_type": "renal epithelial",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "KICH",
    },
    "KIRC": {
        "cancer_type": "clear cell renal cell carcinoma",
        "tissue": "kidney",
        "cell_type": "renal epithelial",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "KIRC ccRCC",
    },
    "LUAD": {
        "cancer_type": "lung adenocarcinoma",
        "tissue": "lung",
        "cell_type": "epithelial",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "LUAD",
    },
    "OV": {
        "cancer_type": "ovarian cancer",
        "tissue": "ovary",
        "cell_type": "epithelial",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "OV ovarian",
    },
    "STAD": {
        "cancer_type": "gastric cancer",
        "tissue": "stomach",
        "cell_type": "epithelial",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "STAD",
    },
    "THCA": {
        "cancer_type": "thyroid cancer",
        "tissue": "thyroid",
        "cell_type": "epithelial follicular",
        "comparison": "tumor vs adjacent",
        "analysis_goal": "tumor metabolism",
        "context_terms": "THCA",
    },
}


def cancer_from_name(path: Path) -> str:
    parts = path.name.split("_")
    return parts[3] if len(parts) > 3 else ""


def brief_rows(rows: list[dict[str, Any]] | None, limit: int = 8) -> list[dict[str, Any]]:
    result = []
    for row in (rows or [])[:limit]:
        result.append(
            {
                "rank": row.get("rank"),
                "name": row.get("display_name") or row.get("drug_name") or row.get("pathway_uid") or row.get("target_uid") or row.get("disease_uid"),
                "score": row.get("score"),
                "fit": row.get("context_fit_score"),
                "tier": row.get("context_fit_tier"),
                "appendix": bool(row.get("context_fit_appendix") or row.get("appendix")),
                "reasons": row.get("context_fit_reasons") or [],
                "penalties": row.get("context_fit_penalties") or [],
            }
        )
    return result


def audit_flags(service: Any, pack: dict[str, Any], pred: dict[str, Any], context: dict[str, str]) -> tuple[list[str], dict[str, Any]]:
    summary = pack.get("input_summary") or {}
    input_count = int(summary.get("input_count") or 0)
    matched = int(summary.get("matched_count") or 0)
    ambiguous = int(summary.get("ambiguous_count") or 0)
    unresolved = ambiguous + int(summary.get("unmatched_count") or 0) + int(summary.get("invalid_count") or 0)
    matched_ratio = matched / input_count if input_count else 0.0
    unresolved_ratio = unresolved / input_count if input_count else 0.0
    pathways = pack.get("pathway_rankings") or []
    diseases = pack.get("disease_rankings") or []
    targets = pack.get("target_rankings") or []
    lit = pack.get("literature_evidence_pack") or {}
    context_groups = set(service.context_cancer_groups(context.values()))
    top_pathway = pathways[0] if pathways else {}
    top_disease = diseases[0] if diseases else {}
    top_target = targets[0] if targets else {}
    top10 = diseases[:10]
    top3_group_match = any(set(row.get("context_fit_display_groups") or []).intersection(context_groups) for row in diseases[:3])
    top10_all_low = all(row.get("context_fit_tier") in {"low_fit", "appendix_fit"} for row in top10) if top10 else False
    wrong_specific = [
        row.get("display_name")
        for row in top10
        if set(row.get("context_fit_display_groups") or [])
        and context_groups
        and not set(row.get("context_fit_display_groups") or []).intersection(context_groups)
        and row.get("context_fit_tier") in {"high_fit", "medium_fit"}
    ]
    flags: list[str] = []
    if matched_ratio < 0.25:
        flags.append("low_matched_ratio")
    if unresolved_ratio > 0.75:
        flags.append("high_ambiguous_unmatched_ratio")
    if not pathways or float(top_pathway.get("context_fit_score") or 0.0) < 0.25:
        flags.append("low_top_pathway_context_fit")
    if not top3_group_match and not top10_all_low:
        flags.append("no_context_matching_disease_in_top3")
    elif not top3_group_match:
        flags.append("current_context_disease_not_recovered_low_fit_background_only")
    if wrong_specific:
        flags.append("wrong_cancer_specific_disease_high_or_medium")
    if float(lit.get("max_context_relevance") or 0.0) < 0.45:
        flags.append("low_literature_context_relevance")
    if not targets:
        flags.append("no_targets")
    if not pred.get("drug_rankings"):
        flags.append("no_drug_rankings")
    metrics = {
        "input_count": input_count,
        "matched": matched,
        "ambiguous": ambiguous,
        "unresolved": unresolved,
        "matched_ratio": round(matched_ratio, 6),
        "unresolved_ratio": round(unresolved_ratio, 6),
        "context_groups": sorted(context_groups),
        "top_pathway": top_pathway.get("display_name", ""),
        "top_pathway_tier": top_pathway.get("context_fit_tier", ""),
        "top_pathway_fit": top_pathway.get("context_fit_score", 0),
        "top_disease": top_disease.get("display_name", ""),
        "top_disease_tier": top_disease.get("context_fit_tier", ""),
        "top_disease_fit": top_disease.get("context_fit_score", 0),
        "top_target": top_target.get("display_name", ""),
        "top_target_fit": top_target.get("context_fit_score", 0),
        "max_context_relevance": lit.get("max_context_relevance", 0),
        "max_biomedbert_relevance": lit.get("max_biomedbert_relevance", 0),
        "top3_group_match": top3_group_match,
        "wrong_specific_diseases": wrong_specific[:5],
    }
    return flags, metrics


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_html(path: Path, rows: list[dict[str, Any]], detail: dict[str, Any], started: str) -> None:
    cols = [
        "status",
        "flags",
        "set",
        "cancer_code",
        "input_count",
        "selected_count",
        "matched",
        "ambiguous",
        "unresolved",
        "matched_ratio",
        "top_pathway",
        "top_pathway_tier",
        "top_disease",
        "top_disease_tier",
        "max_context_relevance",
        "top_explanation_paths",
        "elapsed_sec",
        "file",
    ]
    flagged = [row for row in rows if row.get("status") != "ok"]
    lines = [
        '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><title>TraitScore Batch Context-Fit Audit</title>',
        "<style>body{font-family:system-ui,-apple-system,Segoe UI,sans-serif;max-width:1280px;margin:28px auto;padding:0 18px;line-height:1.5;color:#17202a}table{border-collapse:collapse;width:100%;font-size:12px}th,td{border:1px solid #dbe3ea;padding:6px 7px;vertical-align:top}th{background:#f4f7fa;text-align:left}code{background:#f4f7fa;padding:1px 4px;border-radius:4px}.ok{color:#157347}.warn{color:#a15c00}.err{color:#b42318}pre{white-space:pre-wrap;background:#f6f8fa;border:1px solid #e5e7eb;border-radius:6px;padding:8px}</style></head><body>",
        "<h1>TraitScore_GroupDiff_allCelltypes batch context-fit audit</h1>",
        f"<p>Run: <code>{html.escape(started)}</code>. Files: <code>{len(rows)}</code>. Needs review/errors: <code>{len(flagged)}</code>.</p>",
        '<p><a href="traitscore_batch_context_fit_audit.csv">CSV summary</a> <a href="traitscore_batch_context_fit_audit.json">JSON detail</a></p>',
        "<h2>Summary</h2><table><thead><tr>",
    ]
    lines.extend(f"<th>{html.escape(col)}</th>" for col in cols)
    lines.append("</tr></thead><tbody>")
    for row in rows:
        cls = "ok" if row.get("status") == "ok" else ("err" if row.get("status") == "error" else "warn")
        lines.append("<tr>")
        for col in cols:
            value = row.get(col, "")
            lines.append(f'<td class="{cls if col == "status" else ""}">{html.escape(str(value))}</td>')
        lines.append("</tr>")
    lines.append("</tbody></table><h2>Flag Details</h2>")
    for row in flagged:
        rel = str(row.get("file") or "")
        item = detail.get(rel, {})
        lines.append(f"<h3>{html.escape(rel)}</h3>")
        lines.append(f'<p>Status: <code>{html.escape(str(row.get("status")))}</code>; flags: <code>{html.escape(str(row.get("flags")))}</code></p>')
        lines.append("<pre>" + html.escape(json.dumps({k: item.get(k) for k in ["context", "metrics", "pathways", "diseases", "targets", "drugs"]}, ensure_ascii=False, indent=2)) + "</pre>")
    lines.append("</body></html>")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Batch audit TraitScore CSV context-fit behavior.")
    parser.add_argument("--workspace", type=Path, default=ROOT)
    parser.add_argument("--input-root", type=Path, default=ROOT / "TraitScore_GroupDiff_allCelltypes")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "validation_reports" / "new_evidence_traitscore" / "batch_traitscore_context_fit_audit")
    parser.add_argument("--normalized-root", default="normalized_store_scispacy_abstract_full_20260604T172741")
    parser.add_argument("--literature-root", default="literature_evidence_biomedbert_full_20260604T172741")
    parser.add_argument("--release-id", default="mvp_20260513T002254")
    parser.add_argument("--max-records", type=int, default=300)
    parser.add_argument("--top-per-group", type=int, default=120)
    args = parser.parse_args()

    workspace = args.workspace.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    service = metabo_service.MetaboService(
        workspace=workspace,
        normalized_root=(workspace / args.normalized_root).resolve(),
        literature_root=(workspace / args.literature_root).resolve(),
        release_id=args.release_id,
    )
    files = sorted(path for path in args.input_root.rglob("*.csv") if path.name.startswith("trait_score_diff_"))
    rows: list[dict[str, Any]] = []
    detail: dict[str, Any] = {}
    started = datetime.now(timezone.utc).isoformat()
    for index, path in enumerate(files, start=1):
        rel = path.resolve().relative_to(workspace).as_posix() if path.resolve().is_relative_to(workspace) else str(path.resolve())
        cancer = cancer_from_name(path)
        context = CONTEXTS.get(cancer)
        row: dict[str, Any] = {"file": rel, "set": path.parent.name, "cancer_code": cancer, "status": "started"}
        t0 = time.time()
        print(f"[{index}/{len(files)}] {rel}", flush=True)
        try:
            if not context:
                raise RuntimeError(f"Missing context mapping for {cancer}")
            records = metabo_service.load_records_file(path)
            analyzed = service.analyze_differential_table(
                records,
                max_paths=25,
                max_hops=4,
                context=context,
                context_mode="soft",
                q_threshold=0.05,
                min_abs_effect=0.0,
                top_per_group=args.top_per_group,
                max_records=args.max_records,
            )
            pack = analyzed.get("analysis_pack") or {}
            pred = analyzed.get("predictions") or {}
            flags, metrics = audit_flags(service, pack, pred, context)
            review_flags = [flag for flag in flags if flag != "current_context_disease_not_recovered_low_fit_background_only"]
            status = "ok" if not flags else ("ok_with_boundary" if not review_flags else "needs_review")
            row.update(metrics)
            row.update(
                {
                    "status": status,
                    "flags": ";".join(flags),
                    "elapsed_sec": round(time.time() - t0, 3),
                    "selected_count": (analyzed.get("differential_table_selection") or {}).get("selected_count", len(analyzed.get("selected_records") or [])),
                    "support_count": (pack.get("literature_evidence_pack") or {}).get("support_count", 0),
                    "top_explanation_paths": len(pack.get("top_explanation_paths") or []),
                    "fallback_chains": len(pack.get("fallback_explanation_chains") or []),
                }
            )
            detail[rel] = {
                "context": context,
                "status": row["status"],
                "flags": flags,
                "metrics": metrics,
                "quality_warnings": pack.get("quality_warnings") or [],
                "pathways": brief_rows(pack.get("pathway_rankings"), 12),
                "diseases": brief_rows(pack.get("disease_rankings"), 12),
                "targets": brief_rows(pack.get("target_rankings"), 10),
                "drugs": brief_rows(pred.get("drug_rankings"), 10),
            }
        except Exception as exc:  # pragma: no cover - batch diagnostics
            row.update({"status": "error", "flags": "exception", "error": str(exc), "elapsed_sec": round(time.time() - t0, 3)})
            detail[rel] = {"context": context, "status": "error", "error": str(exc)}
        rows.append(row)
        write_csv(output_dir / "traitscore_batch_context_fit_audit.csv", rows)
        (output_dir / "traitscore_batch_context_fit_audit.json").write_text(
            json.dumps({"started": started, "updated": datetime.now(timezone.utc).isoformat(), "rows": rows, "detail": detail}, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
    write_html(output_dir / "traitscore_batch_context_fit_audit.html", rows, detail, started)
    print(output_dir / "traitscore_batch_context_fit_audit.html")


if __name__ == "__main__":
    main()
