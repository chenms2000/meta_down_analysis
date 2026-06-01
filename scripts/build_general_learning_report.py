"""Build a general-purpose learning report from a learning run.

This is the default report for generalized discovery and prioritization. It
does not require a cancer, tissue, cell-line, or metabolite-table context.
Context-specific reports can be generated separately as case studies.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq


DEFAULT_OUTPUT_ROOT = "learning_runs"

DISEASE_MODEL_PATHWAY_TERMS = (
    "cancer",
    "carcinoma",
    "melanoma",
    "glioma",
    "glioblastoma",
    "leukemia",
    "leukaemia",
    "lymphoma",
    "myeloma",
    "sarcoma",
    "adenoma",
    "tumor",
    "tumour",
    "metastasis",
    "metastatic",
    "disease",
    "disorder",
    "syndrome",
    "deficiency",
    "defective",
    "cell line",
    "cancer cells",
    "carcinoma cells",
    "tumor cells",
    "tumour cells",
)

PRIMARY_METABOLIC_ALLOW_TERMS = (
    "glycolysis",
    "gluconeogenesis",
    "tca",
    "citric acid",
    "tricarboxylic",
    "oxidative phosphorylation",
    "fatty acid",
    "beta oxidation",
    "carnitine",
    "amino acid",
    "glutamine",
    "glutamate",
    "arginine",
    "nitric oxide",
    "purine",
    "pyrimidine",
    "nucleotide",
    "one carbon",
    "methionine",
    "folate",
    "tryptophan",
    "kynurenine",
    "glutathione",
    "redox",
    "ferroptosis",
    "lipid",
    "phospholipid",
    "sphingolipid",
    "glycosylation",
    "glycan",
    "uremic",
    "organic anion",
    "osmolyte",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_parquet(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pq.read_table(path).to_pandas()


def evidence_tier(row: dict[str, Any]) -> str:
    label = str(row.get("weak_label", ""))
    score = float(row.get("priority_score", 0.0) or 0.0)
    if label == "strong_positive" and score >= 0.8:
        return "strong"
    if label in {"strong_positive", "medium_positive"} and score >= 0.45:
        return "medium"
    if label == "weak_positive" and score >= 0.25:
        return "weak"
    return "hypothesis_or_negative"


def top_frame(frame: pd.DataFrame, n: int) -> pd.DataFrame:
    if frame.empty:
        return frame
    return frame.head(n).copy()


def pathway_appendix_classification(display_name: Any) -> tuple[bool, str]:
    text = str(display_name or "").casefold()
    if not text:
        return False, ""
    has_disease_model = any(term in text for term in DISEASE_MODEL_PATHWAY_TERMS)
    if not has_disease_model:
        return False, ""
    has_primary_metabolic = any(term in text for term in PRIMARY_METABOLIC_ALLOW_TERMS)
    if has_primary_metabolic and not any(term in text for term in (" in cancer", " cancer cell", "tumor cell", "cell line")):
        return False, ""
    return True, "disease_or_model_specific_pathway_appendix_in_generalized_report"


def apply_pathway_appendix_policy(pathways: pd.DataFrame) -> pd.DataFrame:
    if pathways.empty:
        return pathways
    frame = pathways.copy()
    if "appendix" not in frame.columns or "downgrade_reason" not in frame.columns:
        classifications = [pathway_appendix_classification(value) for value in frame.get("display_name", pd.Series("", index=frame.index))]
        frame["appendix"] = [flag for flag, _reason in classifications]
        frame["downgrade_reason"] = [reason for _flag, reason in classifications]
    frame["appendix"] = frame["appendix"].fillna(False).astype(bool)
    if "priority_score_raw" not in frame.columns and "priority_score" in frame.columns:
        frame["priority_score_raw"] = frame["priority_score"]
    for column in ("priority_score", "mean_priority_score"):
        if column in frame.columns:
            values = pd.to_numeric(frame[column], errors="coerce").fillna(0.0)
            frame[column] = values.mask(frame["appendix"], values.clip(upper=0.25))
    if "display_tier" not in frame.columns:
        frame["display_tier"] = frame["appendix"].map(lambda flag: "appendix" if flag else "primary")
    return frame.sort_values(["priority_score", "positive_signal_count"], ascending=False)


def filter_excluded(frame: pd.DataFrame, exclude_terms: list[str]) -> pd.DataFrame:
    if frame.empty or not exclude_terms:
        return frame
    terms = [term.casefold() for term in exclude_terms if term]
    if not terms:
        return frame
    string_columns = [column for column in frame.columns if frame[column].dtype == object]
    if not string_columns:
        return frame
    mask = pd.Series(False, index=frame.index)
    for column in string_columns:
        text = frame[column].fillna("").astype(str).str.casefold()
        for term in terms:
            mask = mask | text.str.contains(term, regex=False)
    return frame.loc[~mask].copy()


def md_table(frame: pd.DataFrame, columns: list[str], max_rows: int = 12) -> str:
    if frame.empty:
        return "_No rows._\n"
    subset = frame[[column for column in columns if column in frame.columns]].head(max_rows).copy()
    cols = subset.columns.tolist()
    headers = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    rows = []
    for item in subset.to_dict(orient="records"):
        values = [str(item.get(column, "")).replace("\n", " ")[:180] for column in cols]
        rows.append("| " + " | ".join(values) + " |")
    return "\n".join([headers, sep, *rows]) + "\n"


def relation_summary(relations: pd.DataFrame) -> pd.DataFrame:
    if relations.empty:
        return pd.DataFrame()
    grouped = (
        relations.groupby(["source_kind", "weak_label"], dropna=False)
        .agg(
            rows=("candidate_uid", "count"),
            mean_priority=("priority_score", "mean"),
            max_priority=("priority_score", "max"),
        )
        .reset_index()
        .sort_values(["source_kind", "rows"], ascending=[True, False])
    )
    grouped["mean_priority"] = grouped["mean_priority"].round(6)
    grouped["max_priority"] = grouped["max_priority"].round(6)
    return grouped


def cluster_summary(clusters: pd.DataFrame) -> pd.DataFrame:
    if clusters.empty:
        return pd.DataFrame()
    frame = clusters.copy()
    keep = ["cluster_id", "entity_count", "node_type_counts_json", "examples_json", "mean_centroid_similarity"]
    frame = frame[[column for column in keep if column in frame.columns]].copy()
    if "mean_centroid_similarity" in frame.columns:
        frame["mean_centroid_similarity"] = pd.to_numeric(frame["mean_centroid_similarity"], errors="coerce").round(6)
    return frame.sort_values("entity_count", ascending=False)


def write_report(
    path: Path,
    run_id: str,
    pathway: pd.DataFrame,
    pathway_appendix: pd.DataFrame,
    target: pd.DataFrame,
    drug: pd.DataFrame,
    relation: pd.DataFrame,
    summary: pd.DataFrame,
    clusters: pd.DataFrame,
) -> None:
    lines = [
        f"# General Learning Report: {run_id}",
        "",
        f"- Created: `{utc_now()}`",
        "- Scope: generalized discovery and research-prioritization",
        "- Context: none; no cancer, tissue, cell line, or metabolite table filter is applied",
        "",
        "This report is the generalized view of the learning run. It should be used to inspect broad patterns, candidate relations, and global priority surfaces. Context-specific reports are optional downstream slices for explicitly selected questions.",
        "",
        "## Top Global Pathways",
        md_table(pathway, ["pathway_uid", "display_name", "priority_score", "mean_priority_score", "candidate_count", "positive_signal_count", "top_counterpart_name", "top_label_reason", "display_tier"]),
        "## Appendix: Disease/Model-Specific Pathways",
        "These rows keep their raw evidence but are downgraded out of the generalized primary pathway table because they encode a disease, cell model, or phenotype-specific context.",
        "",
        md_table(pathway_appendix, ["pathway_uid", "display_name", "priority_score_raw", "priority_score", "candidate_count", "positive_signal_count", "top_counterpart_name", "downgrade_reason"], max_rows=12),
        "## Top Global Targets",
        md_table(target, ["target_uid", "display_name", "priority_score", "mean_priority_score", "candidate_count", "positive_signal_count", "top_counterpart_name", "top_label_reason"]),
        "## Top Global Drugs",
        md_table(drug, ["drug_uid", "display_name", "priority_score", "mean_priority_score", "candidate_count", "positive_signal_count", "top_counterpart_name", "top_label_reason"]),
        "## Top Global Relations",
        md_table(relation, ["candidate_uid", "subject_name", "predicate", "object_name", "source_kind", "weak_label", "priority_score", "evidence_tier", "label_reason"]),
        "## Weak Label Distribution",
        md_table(summary, ["source_kind", "weak_label", "rows", "mean_priority", "max_priority"], max_rows=20),
        "## Largest Discovery Clusters",
        md_table(clusters, ["cluster_id", "entity_count", "node_type_counts_json", "mean_centroid_similarity"], max_rows=15),
        "## Interpretation",
        "- Priority scores are global research-ranking signals, not truth probabilities.",
        "- Strong and medium tiers mean the candidate has stronger weak-supervision evidence, not manual validation.",
        "- Disease/model-specific pathways are retained in appendix outputs and should not be treated as primary generalized metabolic conclusions.",
        "- This global report is the right default for generalized learning. Use context-specific reports only when asking a specific biological question.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a generalized learning report for a learning run.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--output-prefix", default="global")
    parser.add_argument(
        "--exclude-terms",
        default="",
        help="Comma-separated terms to omit from generated report rows. This only filters report outputs, not source data.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    run_dir = workspace / args.output_root / args.run_id
    ranking_dir = run_dir / "rankings"
    report_dir = run_dir / "reports"

    pathways = read_parquet(ranking_dir / "pathway_priorities.parquet")
    targets = read_parquet(ranking_dir / "target_priorities.parquet")
    drugs = read_parquet(ranking_dir / "drug_priorities.parquet")
    relations = read_parquet(ranking_dir / "relation_priorities.parquet")
    clusters = cluster_summary(read_parquet(ranking_dir / "entity_clusters.parquet"))
    exclude_terms = [term.strip() for term in str(args.exclude_terms or "").split(",") if term.strip()]

    if not relations.empty:
        relations = relations.copy()
        relations["evidence_tier"] = [evidence_tier(row) for row in relations.to_dict(orient="records")]

    pathways = apply_pathway_appendix_policy(pathways)
    pathways = filter_excluded(pathways, exclude_terms)
    targets = filter_excluded(targets, exclude_terms)
    drugs = filter_excluded(drugs, exclude_terms)
    relations = filter_excluded(relations, exclude_terms)
    clusters = filter_excluded(clusters, exclude_terms)

    primary_pathways = pathways[~pathways.get("appendix", pd.Series(False, index=pathways.index)).fillna(False).astype(bool)].copy() if not pathways.empty else pathways
    appendix_pathways = pathways[pathways.get("appendix", pd.Series(False, index=pathways.index)).fillna(False).astype(bool)].copy() if not pathways.empty else pathways
    top_pathways = top_frame(primary_pathways, args.top_n)
    top_appendix_pathways = top_frame(appendix_pathways, args.top_n)
    top_targets = top_frame(targets, args.top_n)
    top_drugs = top_frame(drugs, args.top_n)
    top_relations = top_frame(relations, args.top_n)
    weak_summary = relation_summary(relations)

    report_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_prefix
    top_pathways.to_csv(report_dir / f"{prefix}_top_pathways.csv", index=False, encoding="utf-8-sig")
    top_appendix_pathways.to_csv(report_dir / f"{prefix}_appendix_pathways.csv", index=False, encoding="utf-8-sig")
    top_targets.to_csv(report_dir / f"{prefix}_top_targets.csv", index=False, encoding="utf-8-sig")
    top_drugs.to_csv(report_dir / f"{prefix}_top_drugs.csv", index=False, encoding="utf-8-sig")
    top_relations.to_csv(report_dir / f"{prefix}_top_relations.csv", index=False, encoding="utf-8-sig")
    weak_summary.to_csv(report_dir / f"{prefix}_weak_label_summary.csv", index=False, encoding="utf-8-sig")
    clusters.head(args.top_n).to_csv(report_dir / f"{prefix}_cluster_summary.csv", index=False, encoding="utf-8-sig")
    write_report(
        report_dir / f"{prefix}_learning_report.md",
        args.run_id,
        top_pathways,
        top_appendix_pathways,
        top_targets,
        top_drugs,
        top_relations,
        weak_summary,
        clusters,
    )

    manifest = {
        "created_at_utc": utc_now(),
        "run_id": args.run_id,
        "scope": "generalized_no_context_filter",
        "exclude_term_count": len(exclude_terms),
        "outputs": {
            "learning_report": str(report_dir / f"{prefix}_learning_report.md"),
            "top_pathways": str(report_dir / f"{prefix}_top_pathways.csv"),
            "appendix_pathways": str(report_dir / f"{prefix}_appendix_pathways.csv"),
            "top_targets": str(report_dir / f"{prefix}_top_targets.csv"),
            "top_drugs": str(report_dir / f"{prefix}_top_drugs.csv"),
            "top_relations": str(report_dir / f"{prefix}_top_relations.csv"),
            "weak_label_summary": str(report_dir / f"{prefix}_weak_label_summary.csv"),
            "cluster_summary": str(report_dir / f"{prefix}_cluster_summary.csv"),
        },
        "notes": [
            "No cancer, tissue, cell-line, or metabolite-table context was applied.",
            "Use context-specific reports only as optional slices, not as the generalized model output.",
        ],
    }
    (report_dir / f"{prefix}_learning_report_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest["outputs"], ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
