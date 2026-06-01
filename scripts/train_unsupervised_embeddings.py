"""Train lightweight unsupervised entity embeddings from learning views."""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize


DEFAULT_OUTPUT_ROOT = "learning_runs"


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def write_parquet(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), path)


def signal_score(frame: pd.DataFrame) -> pd.Series:
    return (
        np.log1p(pd.to_numeric(frame.get("graph_degree", 0), errors="coerce").fillna(0.0))
        + 2.0 * np.log1p(pd.to_numeric(frame.get("literature_mention_count", 0), errors="coerce").fillna(0.0))
        + 3.0 * np.log1p(pd.to_numeric(frame.get("literature_support_count", 0), errors="coerce").fillna(0.0))
        + 2.0 * np.log1p(pd.to_numeric(frame.get("prediction_overlay_count", 0), errors="coerce").fillna(0.0))
    )


def select_entities(entities: pd.DataFrame, max_entities: int) -> pd.DataFrame:
    frame = entities.copy()
    frame["signal_score"] = signal_score(frame)
    frame["document"] = frame["document"].fillna("").astype(str)
    frame = frame[frame["document"].str.len() > 0].copy()
    if max_entities > 0 and len(frame) > max_entities:
        frame = frame.sort_values(["signal_score", "graph_degree"], ascending=False).head(max_entities).copy()
    return frame.reset_index(drop=True)


def make_embeddings(frame: pd.DataFrame, max_features: int, components: int) -> tuple[np.ndarray, dict[str, Any]]:
    docs = frame["document"].fillna("").astype(str).tolist()
    vectorizer = TfidfVectorizer(
        max_features=max_features,
        min_df=2,
        max_df=0.98,
        ngram_range=(1, 2),
        token_pattern=r"(?u)\b[a-zA-Z][\w\-]{1,}\b",
        sublinear_tf=True,
    )
    matrix = vectorizer.fit_transform(docs)
    usable_components = max(2, min(components, matrix.shape[0] - 1, matrix.shape[1] - 1))
    svd = TruncatedSVD(n_components=usable_components, random_state=13)
    dense = svd.fit_transform(matrix)
    dense = normalize(dense, norm="l2")
    metadata = {
        "tfidf_features": int(matrix.shape[1]),
        "tfidf_rows": int(matrix.shape[0]),
        "svd_components": int(usable_components),
        "svd_explained_variance_ratio_sum": float(np.sum(svd.explained_variance_ratio_)),
    }
    return dense.astype(np.float32), metadata


def build_embedding_frame(frame: pd.DataFrame, embeddings: np.ndarray) -> pd.DataFrame:
    base = frame[
        [
            "entity_uid",
            "node_type",
            "display_name",
            "canonical_name",
            "primary_external_id",
            "signal_score",
            "graph_degree",
            "literature_mention_count",
            "literature_support_count",
            "prediction_overlay_count",
        ]
    ].reset_index(drop=True)
    dims = pd.DataFrame(
        embeddings,
        columns=[f"dim_{idx:03d}" for idx in range(embeddings.shape[1])],
    )
    return pd.concat([base, dims], axis=1)


def build_neighbors(frame: pd.DataFrame, embeddings: np.ndarray, neighbors: int, query_limit: int) -> pd.DataFrame:
    n_neighbors = min(neighbors + 1, len(frame))
    model = NearestNeighbors(n_neighbors=n_neighbors, metric="cosine", algorithm="brute")
    model.fit(embeddings)
    query_count = len(frame) if query_limit <= 0 else min(query_limit, len(frame))
    distances, indices = model.kneighbors(embeddings[:query_count], return_distance=True)
    rows: list[dict[str, Any]] = []
    for query_idx in range(query_count):
        query = frame.iloc[query_idx]
        rank = 0
        for distance, neighbor_idx in zip(distances[query_idx], indices[query_idx]):
            if int(neighbor_idx) == query_idx:
                continue
            rank += 1
            neighbor = frame.iloc[int(neighbor_idx)]
            rows.append(
                {
                    "query_uid": query.entity_uid,
                    "query_type": query.node_type,
                    "query_name": query.display_name,
                    "neighbor_uid": neighbor.entity_uid,
                    "neighbor_type": neighbor.node_type,
                    "neighbor_name": neighbor.display_name,
                    "rank": rank,
                    "cosine_similarity": float(1.0 - distance),
                }
            )
            if rank >= neighbors:
                break
    return pd.DataFrame(rows)


def build_clusters(frame: pd.DataFrame, embeddings: np.ndarray, cluster_count: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    usable_clusters = max(2, min(cluster_count, len(frame)))
    model = MiniBatchKMeans(n_clusters=usable_clusters, random_state=13, batch_size=4096, n_init="auto")
    labels = model.fit_predict(embeddings)
    assignments = frame[["entity_uid", "node_type", "display_name", "signal_score"]].copy()
    assignments["cluster_id"] = labels.astype(int)

    rows: list[dict[str, Any]] = []
    for cluster_id, group in assignments.groupby("cluster_id"):
        type_counts = group["node_type"].value_counts().head(8).to_dict()
        examples = (
            group.sort_values("signal_score", ascending=False)
            .head(12)[["entity_uid", "node_type", "display_name"]]
            .to_dict(orient="records")
        )
        centroid = model.cluster_centers_[int(cluster_id)].reshape(1, -1)
        member_embeddings = embeddings[group.index.to_numpy()]
        mean_similarity = float(np.mean(cosine_similarity(member_embeddings, centroid))) if len(group) else 0.0
        rows.append(
            {
                "cluster_id": int(cluster_id),
                "entity_count": int(len(group)),
                "node_type_counts_json": json.dumps(type_counts, ensure_ascii=False, sort_keys=True),
                "examples_json": json.dumps(examples, ensure_ascii=False),
                "mean_centroid_similarity": mean_similarity,
            }
        )
    return assignments, pd.DataFrame(rows).sort_values("entity_count", ascending=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train lightweight unsupervised embeddings from learning views.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-root", default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--max-entities", type=int, default=120_000)
    parser.add_argument("--max-features", type=int, default=60_000)
    parser.add_argument("--components", type=int, default=64)
    parser.add_argument("--neighbors", type=int, default=10)
    parser.add_argument("--neighbor-query-limit", type=int, default=25_000)
    parser.add_argument("--clusters", type=int, default=80)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(args.workspace).resolve()
    run_dir = workspace / args.output_root / args.run_id
    views_dir = run_dir / "views"
    embedding_dir = run_dir / "embeddings"
    ranking_dir = run_dir / "rankings"

    entities = pq.read_table(views_dir / "entities.parquet").to_pandas()
    selected = select_entities(entities, args.max_entities)
    embeddings, embedding_metadata = make_embeddings(selected, args.max_features, args.components)
    embedding_frame = build_embedding_frame(selected, embeddings)
    assignments, cluster_summary = build_clusters(selected, embeddings, args.clusters)
    embedding_frame = embedding_frame.merge(assignments[["entity_uid", "cluster_id"]], on="entity_uid", how="left")
    neighbors = build_neighbors(selected, embeddings, args.neighbors, args.neighbor_query_limit)

    write_parquet(embedding_dir / "entity_embeddings.parquet", embedding_frame)
    write_parquet(ranking_dir / "entity_neighbors.parquet", neighbors)
    write_parquet(ranking_dir / "entity_clusters.parquet", cluster_summary)
    write_parquet(ranking_dir / "entity_cluster_assignments.parquet", assignments)

    manifest = {
        "created_at_utc": utc_now(),
        "run_id": args.run_id,
        "input_entities": str(views_dir / "entities.parquet"),
        "selected_entities": int(len(selected)),
        "embedding_metadata": embedding_metadata,
        "outputs": {
            "entity_embeddings": {"path": str(embedding_dir / "entity_embeddings.parquet"), "rows": int(len(embedding_frame))},
            "entity_neighbors": {"path": str(ranking_dir / "entity_neighbors.parquet"), "rows": int(len(neighbors))},
            "entity_clusters": {"path": str(ranking_dir / "entity_clusters.parquet"), "rows": int(len(cluster_summary))},
        },
        "notes": [
            "Embeddings are unsupervised discovery features, not factual assertions.",
            "Selection favors entities with graph, literature, or prediction-overlay signal.",
        ],
    }
    embedding_dir.mkdir(parents=True, exist_ok=True)
    (embedding_dir / "embedding_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest["outputs"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
