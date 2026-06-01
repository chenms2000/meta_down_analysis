# 代谢组知识图谱与智能解释工作台

这个目录现在不只是下载器，也包含一个本地可运行的研究级代谢组知识图谱服务：它能把代谢物列表或差异代谢表解析成规范实体，基于冻结 release 做通路/靶点/疾病排序，回链文献证据，并通过 Web 界面给出可审计解释。

快速打开互动界面：

```powershell
python .\scripts\metabo_service.py --workspace . --release-id mvp_20260513T002254 serve --port 8765
```

浏览器访问：

```text
http://127.0.0.1:8765/
```

流程测试快捷入口：

```text
http://127.0.0.1:8765/flow-test
```

## TraitScore_GroupDiff_allCelltypes 能不能分析？

可以。`TraitScore_GroupDiff_allCelltypes/` 里的 `trait_score_diff_*.csv`
和 `.csv.gz` 属于已经算好的两组差异表，不是普通 LC-MS 原始丰度表。
推荐用通用差异表入口 `analyze-differential-table`：

```powershell
python .\scripts\metabo_service.py --workspace . --release-id mvp_20260513T002254 analyze-differential-table .\TraitScore_GroupDiff_allCelltypes\significant\trait_score_diff_LUAD_Epi_LUAD_Tumor_vs_Adjacent_significant_q0.05.csv --max-records 80 --top-per-group 80
```

总表也能作为输入，例如：

```powershell
python .\scripts\metabo_service.py --workspace . --release-id mvp_20260513T002254 analyze-differential-table .\TraitScore_GroupDiff_allCelltypes\trait_score_diff_ALLCELLTYPES_Tumor_vs_Adjacent_significant_q0.05.csv.gz --max-records 120 --top-per-group 10
```

这个入口只读取差异结果列，例如 `trait/GCST`、`group1/group2`、`mean_diff`、
`pseudo_log2FC_shifted`、`z_wilcoxon`、`q_wilcoxon_global` 等，不需要原始丰度矩阵。
它会把预计算差异转换成研究解释用的代谢 trait seed。结果里会明确区分：
严格匹配的代谢物、低权重扩展候选、歧义项和未匹配项。

注意：TraitScore 或其他差异统计量是“预计算差异信号”，不等同于直接代谢物丰度。
当没有标准 `log2FC` 时，系统会把 `mean_diff`、`cohen_d`、`z_wilcoxon` 等标成
非丰度效应量 surrogate。输出适合做研究候选通路/靶点解释，不能写成“某代谢物实测上调/下调”的结论。
全量总表较慢，GitHub 演示建议先跑 `significant/` 里的单个癌种 CSV 或用 `--max-records`
限制输入规模。旧入口 `analyze-trait-score` 仍可用，但新入口语义更贴近差异表。

详细流程、输入格式、界面说明、API 和重建命令见 [docs/user_guide.md](docs/user_guide.md)。

## 当前主要功能

- 原始数据下载与 manifest 审计：`download_all.ps1`、`scripts/download_databases.py`
- 规范化数据仓：`scripts/build_normalized_store.py`
- 图谱投影与解析索引：`scripts/build_graph_projection.py`
- 化合物匹配索引：`scripts/build_compound_match_index.py`
- 文献证据 overlay：`scripts/build_literature_evidence.py`
- 只读分析服务与 Web 工作台：`scripts/metabo_service.py`、`web/chat.html`
- 证据绑定的安全解释层：`scripts/llm_safe_adapter.py`
- 外部 LLM 连接自检：页面内“测试 LLM 连接”按钮和 `/llm/test` 接口可诊断鉴权、代理、DNS、模型名和超时问题
- 差异表 / TraitScore 差异表可通过 `/chat` 包装路径进入本地或外部 LLM 叙述层；LLM 只读取冻结分析包和证据，不参与解析或评分
- 大输入可用“外部 LLM（分层文字）”：按解析质量、ratio/class/identity、排名、证据和低置信附录分块叙述，逐块 guard 后再汇总
- 结果可信度与交互浏览：Web 工作台提供可信度阅读卡、可点击证据图谱、JSON/Markdown 导出
- 研究优先级 overlay：`manual_sources/prediction_overlays/<release_id>/`
- 发布验证与回归测试：`scripts/run_phase15_validation.py`、`tests/`

> 注意：本系统是研究工具，不是临床决策系统。歧义匹配不会被强行纳入评分，外部 LLM 也不能创建事实、修改图谱或改写评分。

## 原始数据下载层

这个目录现在包含一个可审计的原始数据下载层，用于把开放核心层数据源下载到版本化 `raw_lake`，并为每次运行生成 manifest。

## 文件

- `config/source_catalog.toml`：数据源清单、官方入口、许可层级、下载规则。
- `scripts/download_databases.py`：纯标准库 Python 下载器，支持目录索引解析、断点续传、sha256、`.md5` sidecar 校验、manifest。
- `download_all.ps1`：Windows 一键入口。
- `manifests/`：运行后生成的审计清单。
- `raw_lake/`：运行后生成的原始数据湖。

## 快速使用

先做 dry run，不下载文件：

```powershell
.\download_all.ps1
```

解析远程目录，查看这次会下载哪些动态文件：

```powershell
.\download_all.ps1 -ResolveIndexes
```

正式下载开放核心层：

```powershell
.\download_all.ps1 -Run -ResolveIndexes
```

只下载某几个源：

```powershell
.\download_all.ps1 -Run -ResolveIndexes -Sources chebi,reactome,wikipathways
```

包含默认关闭的重型开放源，例如 PubTator、Europe PMC、PMC OA 或完整 PubChem SDF：

```powershell
.\download_all.ps1 -Run -ResolveIndexes -IncludeDisabled -Sources pubtator3,europe_pmc_oa
```

许可增强层默认不会下载。确实完成许可复核后再显式开启：

```powershell
.\download_all.ps1 -Run -ResolveIndexes -Layers licensed_enhancement -IncludeDisabled -AcceptLicensed
```

## Python CLI

```powershell
python .\scripts\download_databases.py --list-sources
python .\scripts\download_databases.py --dry-run --resolve-indexes --sources chebi
python .\scripts\download_databases.py --yes --resolve-indexes --sources chebi --jobs 4
```

测试目录解析时可以限制每个源的文件数：

```powershell
python .\scripts\download_databases.py --dry-run --resolve-indexes --limit-per-source 3
```

## 许可分层

默认层是 `open_core`。其中 PubMed baseline、Open Targets evidence、PubTator、Europe PMC/PMC OA、完整 PubChem SDF 可能非常大，部分被标记为 `enabled_by_default = false`，需要 `-IncludeDisabled` 才会进入计划。

`licensed_enhancement` 和 `controlled_access` 不会因为 `-IncludeDisabled` 被误抓取，还必须额外加 `-AcceptLicensed`。HMDB、KEGG、DisGeNET 目前保留为 manual 规则，避免在许可未确认时自动化下载。

## 输出结构

正式运行后，每个源写入：

```text
raw_lake/<source_id>/<release_id>/...
```

manifest 写入：

```text
manifests/<release_id>.manifest.json
```

manifest 记录官方入口、许可策略、下载 URL、本地路径、字节数、sha256、md5 校验状态和失败原因。这个文件后续可以直接作为发布审计、解析任务输入或回滚依据。

## 调整版本

大部分源使用 `current/` 或目录索引动态解析。少数季度发布源需要人工更新目录，例如 Open Targets 当前配置为：

```toml
https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.03/output/
```

下一次 Open Targets 发布后，只需要在 `config/source_catalog.toml` 中把 `26.03` 改成新版本，再 dry run 验证。
