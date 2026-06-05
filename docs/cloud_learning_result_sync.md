# 云端学习结果下载与复现说明

本文记录云服务器训练结果如何同步回本地，以及哪些内容应该进入 Git。

当前示例对应上一轮云端学习 run：

```text
run_id: learn_scibio_20260605T165940Z
cloud_workspace: /disk3/ms/meta_down_analysis_learning_scibio
cloud_run_dir: /disk3/ms/meta_down_analysis_learning_scibio/learning_runs/learn_scibio_20260605T165940Z
```

这些输出属于研究优先级和解释辅助层。它们可以用于排序、聚类、邻居推荐和验证优先级提示，但不能直接升级为事实，也不能作为临床或治疗建议。

## 一、本地下载

在本地 PowerShell 中进入仓库目录：

```powershell
cd D:\meta_down_analysis
```

如果云服务器 SSH 用户、IP 和 run id 与默认值一致，直接运行：

```powershell
.\scripts\sync_cloud_learning_results.ps1
```

默认参数等价于：

```powershell
.\scripts\sync_cloud_learning_results.ps1 `
  -RemoteUser ms `
  -RemoteHost 10.64.2.166 `
  -RemoteWorkspace /disk3/ms/meta_down_analysis_learning_scibio `
  -RunId learn_scibio_20260605T165940Z `
  -LocalWorkspace D:\meta_down_analysis
```

下载后，本地结果会落在：

```text
D:\meta_down_analysis\learning_runs\learn_scibio_20260605T165940Z
```

同步 manifest 会写入：

```text
D:\meta_down_analysis\manifests\cloud_learning_sync\learn_scibio_20260605T165940Z.sync.json
```

`learning_runs/` 和 `manifests/` 已被 `.gitignore` 排除，不会被误提交到 GitHub。

如果本地已经有同名 run，需要覆盖时加：

```powershell
.\scripts\sync_cloud_learning_results.ps1 -Overwrite
```

如果希望保留下载的压缩包以便另存或校验：

```powershell
.\scripts\sync_cloud_learning_results.ps1 -KeepArchive
```

## 二、应检查的关键输出

脚本会检查以下核心文件是否存在：

```text
rankings/drug_priorities.parquet
rankings/pathway_priorities.parquet
rankings/relation_priorities.parquet
rankings/target_priorities.parquet
rankings/entity_embeddings.parquet
rankings/entity_neighbors.parquet
rankings/entity_clusters.parquet
rankings/entity_cluster_assignments.parquet
weak_labels/weak_labels.parquet
models/priority_ranker.joblib
reports/priority_ranker_report.json
reports/global_learning_report.md
reports/global_top_pathways.csv
reports/global_top_targets.csv
reports/global_top_drugs.csv
reports/global_top_relations.csv
```

如果脚本退出码为 `2`，表示下载完成但缺少某些预期文件，应先查看云端 run 是否完整。

## 三、云端复现命令

本次 run 使用隔离工作区，不覆盖原始发布目录：

```bash
cd /disk3/ms/meta_down_analysis

PY=/disk3/ms/miniconda3/envs/metabo_3.11/bin/python3
RUN=learn_scibio_$(date -u +%Y%m%dT%H%M%SZ)
W=/disk3/ms/meta_down_analysis_learning_scibio

mkdir -p "$W"
ln -sfn /disk3/ms/meta_down_analysis/graph_projection "$W/graph_projection"
ln -sfn /disk3/ms/meta_down_analysis/normalized_store_scispacy_abstract_full_20260604T172741 "$W/normalized_store"
ln -sfn /disk3/ms/meta_down_analysis/literature_evidence_biomedbert_full_20260604T172741 "$W/literature_evidence"
ln -sfn /disk3/ms/meta_down_analysis/manual_sources "$W/manual_sources"

$PY scripts/build_learning_views.py --workspace "$W" --release-id mvp_20260513T002254 --run-id "$RUN"
$PY scripts/train_unsupervised_embeddings.py --workspace "$W" --run-id "$RUN" --max-entities 60000 --max-features 50000 --components 48 --neighbors 10 --neighbor-query-limit 10000 --clusters 60
$PY scripts/build_weak_labels.py --workspace "$W" --run-id "$RUN" --negatives-per-positive 1.0
$PY scripts/train_priority_ranker.py --workspace "$W" --run-id "$RUN"
$PY scripts/build_general_learning_report.py --workspace "$W" --run-id "$RUN" --output-prefix global
$PY scripts/run_temporal_holdout_validation.py --workspace "$W" --run-id "$RUN"
$PY scripts/run_source_heldout_validation.py --workspace "$W" --run-id "$RUN" --heldout-source label_reason=curated_and_literature
```

后台运行版本：

```bash
nohup bash -lc '
set -euo pipefail
cd /disk3/ms/meta_down_analysis
PY=/disk3/ms/miniconda3/envs/metabo_3.11/bin/python3
RUN=learn_scibio_$(date -u +%Y%m%dT%H%M%SZ)
W=/disk3/ms/meta_down_analysis_learning_scibio

mkdir -p "$W"
ln -sfn /disk3/ms/meta_down_analysis/graph_projection "$W/graph_projection"
ln -sfn /disk3/ms/meta_down_analysis/normalized_store_scispacy_abstract_full_20260604T172741 "$W/normalized_store"
ln -sfn /disk3/ms/meta_down_analysis/literature_evidence_biomedbert_full_20260604T172741 "$W/literature_evidence"
ln -sfn /disk3/ms/meta_down_analysis/manual_sources "$W/manual_sources"

echo "RUN=$RUN"
$PY scripts/build_learning_views.py --workspace "$W" --release-id mvp_20260513T002254 --run-id "$RUN"
$PY scripts/train_unsupervised_embeddings.py --workspace "$W" --run-id "$RUN" --max-entities 60000 --max-features 50000 --components 48 --neighbors 10 --neighbor-query-limit 10000 --clusters 60
$PY scripts/build_weak_labels.py --workspace "$W" --run-id "$RUN" --negatives-per-positive 1.0
$PY scripts/train_priority_ranker.py --workspace "$W" --run-id "$RUN"
$PY scripts/build_general_learning_report.py --workspace "$W" --run-id "$RUN" --output-prefix global
$PY scripts/run_temporal_holdout_validation.py --workspace "$W" --run-id "$RUN" || true
$PY scripts/run_source_heldout_validation.py --workspace "$W" --run-id "$RUN" --heldout-source label_reason=curated_and_literature || true
echo "DONE $RUN"
' > logs/learning_scibio_latest.log 2>&1 &
```

查看进度：

```bash
tail -f /disk3/ms/meta_down_analysis/logs/learning_scibio_latest.log
```

## 四、Git 同步原则

提交到 GitHub 的内容应包括：

```text
scripts/sync_cloud_learning_results.ps1
docs/cloud_learning_result_sync.md
docs/learning_pipeline_mvp.md
README.md
```

不要提交：

```text
learning_runs/
normalized_store*/
literature_evidence*/
graph_projection/
manual_sources/
raw_lake/
manifests/
*.tgz
*.parquet
*.joblib
```

推荐提交流程：

```powershell
git status --short
git add scripts\sync_cloud_learning_results.ps1 docs\cloud_learning_result_sync.md docs\learning_pipeline_mvp.md README.md
git commit -m "Document cloud learning result sync workflow"
git push origin main
```

如果已有未提交的本地改动，先用 `git status --short` 分辨哪些是本次文档和同步脚本，避免把临时数据或旧实验输出一起提交。

## 五、解释边界

下载回来的 learning run 必须按以下边界使用：

- `rankings/*.parquet` 是研究优先级排序，不是生物学真值概率。
- `entity_embeddings`、`entity_neighbors` 和 `entity_clusters` 是探索层，只能用于相似性、聚类和候选推荐。
- weak labels 是规则和来源构造的训练信号，不等于人工金标准。
- temporal/source-held-out validation 是代理验证；没有人工 adjudication 时，应明确写出验证局限。
- 任何自然语言结论都要回链到输入代谢物、图谱边、文献证据、排序字段或显式模型输出，并标注置信度、降级原因和解释边界。
