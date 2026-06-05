# meta_down_analysis

本仓库提供一个研究原型，用于对上游分析产生的 GCST/trait 差异表和普通代谢物表进行结构化解析与证据组织。该原型旨在将输入记录分层连接到本地注释、代谢物候选、知识图谱关系、文献证据、排序结果和置信度边界，形成可审计的研究解释。该工具仍需进一步方法学验证和工程完善，不应视为临床决策系统。

## 范围与定位

- **输入对象**：GCST/trait 差异表、普通代谢物差异表、带 HMDB/ChEBI/PubChem/InChIKey 的代谢物列表。
- **主要用途**：将上游表格转换为可追溯的研究解释线索，而不是重新执行上游统计分析。
- **输出形式**：实体匹配、歧义项、低置信候选、通路/靶点/疾病排序、证据引用、解释边界和导出 JSON/Markdown。
- **当前边界**：GCST 解析依赖本地 annotation 文件；图谱传播和模型排序仅作为研究优先级信号，不支持临床或治疗结论。

## 本地数据依赖

如果输入表只有 GCST accession 而没有代谢物名称，本地运行时需要额外准备 GCST annotation 文件。该文件作为数据依赖单独管理，不随代码仓库分发，应放在：

```text
raw_lake/European/European_trait_annotations.csv
```

或：

```text
raw_lake/European_point/European_trait_annotations.csv
```

字段模板见 [config/european_trait_annotations.template.csv](config/european_trait_annotations.template.csv)，详细说明见 [docs/gcst_annotation_dependency.md](docs/gcst_annotation_dependency.md)。

普通代谢物表不需要该 GCST annotation 文件；包含 `name`、`metabolite`、`HMDB`、`ChEBI`、`PubChem CID`、`InChIKey`、`log2FC`、`pvalue`、`padj`、`direction` 等常见列时，可使用常规代谢物解析路径。

文献证据 overlay 使用本地重建的文献语料。GitHub 仓库只记录检索式、时间范围、筛选条件和重建边界，不分发文献记录、摘要、全文、影响因子表或句级抽取结果。检索策略见 [config/literature_search_strategy.json](config/literature_search_strategy.json)，说明见 [docs/literature_corpus_dependency.md](docs/literature_corpus_dependency.md)。

## 方法学概览

1. **输入识别**：系统先判断输入是普通代谢物表、GCST/trait 表，还是已经计算好的两组差异结果表。差异表只读取标识符、效应量、显著性、方向和分组元数据，不读取原始丰度矩阵。
2. **GCST 注释映射**：对 GCST-only 输入，系统用本地 `European_trait_annotations.csv` 将 accession 映射到 reported trait、候选代谢物、ratio component、class/pool 线索或稳定化合物 ID。没有 annotation 时，GCST 只能作为 trait-level 记录保留。
3. **实体解析**：代谢物名称和稳定 ID 会进入冻结的 compound match index。严格匹配保留为较高可信度种子；歧义名称、ratio trait、class/pool 或 analog candidate 仅作为低权重扩展种子。
4. **差异信号处理**：`log2FC`、`mean_diff`、`cohen_d`、`z_wilcoxon`、FDR/q 值等被保留为带标签的效应量。若输入并非直接丰度测量，输出会标记为 surrogate differential signal，避免将其表述为实测代谢物丰度变化。
5. **知识图谱与证据连接**：解析后的种子连接到代谢物、通路、反应、基因/靶点、疾病和文献证据 overlay。排名结果保留 edge、source record、literature evidence 或模型输出引用。
6. **置信度与边界**：结论按严格匹配、扩展候选、图传播、文献支持和模型-only 信号分层。弱证据、歧义匹配、ratio component、未复核 GCST 和纯模型排序会被降级。
7. **叙述层**：本地模板或外部 LLM 只读取冻结 analysis pack。外部 LLM 不能创建实体、改写评分、修改图谱或补充无来源事实。

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

详细流程、输入格式、界面说明、API 和重建命令见 [docs/user_guide.md](docs/user_guide.md)。

## 当前主要功能

- 原始数据下载与 manifest 审计：`download_all.ps1`、`scripts/download_databases.py`
- 规范化数据仓：`scripts/build_normalized_store.py`
- 图谱投影与解析索引：`scripts/build_graph_projection.py`
- 化合物匹配索引：`scripts/build_compound_match_index.py`
- 文献证据 overlay：`scripts/build_literature_evidence.py`
- 文献语料重建策略：`config/literature_search_strategy.json`、`docs/literature_corpus_dependency.md`
- 只读分析服务与 Web 工作台：`scripts/metabo_service.py`、`web/chat.html`
- 证据绑定的安全解释层：`scripts/llm_safe_adapter.py`
- 外部 LLM 连接自检：页面内“测试 LLM 连接”按钮和 `/llm/test` 接口可诊断鉴权、代理、DNS、模型名和超时问题
- GCST 差异结果表可通过 `/chat` 包装路径进入本地或外部 LLM 叙述层；LLM 只读取冻结分析包和证据，不参与解析或评分
- 大输入可用“外部 LLM（分层文字）”：按解析质量、ratio/class/identity、排名、证据和低置信附录分块叙述，逐块 guard 后再汇总
- 结果可信度与交互浏览：Web 工作台提供可信度阅读卡、可点击证据图谱、JSON/Markdown 导出
- 研究优先级 overlay：`manual_sources/prediction_overlays/<release_id>/`
- 云端学习结果下载与复现：`scripts/sync_cloud_learning_results.ps1`、`docs/cloud_learning_result_sync.md`
- 发布验证与回归测试：`scripts/run_phase15_validation.py`、`tests/`

> 注意：本系统用于研究解释和验证优先级排序，不构成临床决策系统。歧义匹配不会被强行提升为高置信证据，外部 LLM 也不能创建事实、修改图谱或改写评分。

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

## 引用与开源许可

本仓库代码以 MIT License 开源，见 [LICENSE](LICENSE)。如果在论文、报告或复现实验中使用本仓库，请优先引用冻结的 GitHub Release / Zenodo DOI；在 DOI 生成前，可引用本仓库地址：

```text
https://github.com/chenms2000/meta_down_analysis
```

仓库根目录提供 [CITATION.cff](CITATION.cff)，GitHub 会据此显示引用入口。论文 DOI 或 Zenodo DOI 生成后，应同步更新 `CITATION.cff`、README 和论文方法/数据可用性声明。

许可边界：MIT License 适用于本仓库中的原创代码、脚本和项目文档。通过 `download_all.ps1`、`scripts/download_databases.py` 或配置文件引用、下载、重建的第三方数据库、文献、标识符、外部知识库和派生数据，应遵循其各自来源、许可条款和引用要求。本仓库不重新授权第三方数据源。
