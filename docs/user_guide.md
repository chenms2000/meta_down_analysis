# 代谢组知识图谱工作流与界面使用说明

当前可用发布：`mvp_20260513T002254`

这个项目是研究级代谢组知识图谱与解释服务。它把公开数据库、规范化实体、图谱投影、文献证据和可审计解释层组合起来，用于分析用户上传的代谢物列表或差异代谢表。它不是临床决策系统。

## 一、最快使用流程

如果当前目录已经有 `normalized_store/`、`graph_projection/`、`compound_match_index/` 和 `literature_evidence/`，可以直接启动本地界面：

```powershell
python .\scripts\metabo_service.py --workspace . --release-id mvp_20260513T002254 serve --port 8765
```

然后打开：

```text
http://127.0.0.1:8765/
```

如需直接运行内置样例，可打开流程测试快捷入口，它会自动载入样例并提交完整分析链路：

```text
http://127.0.0.1:8765/flow-test
```

界面使用顺序：

1. 上传 CSV/TSV/TXT/JSON，或直接粘贴代谢物表。
2. 选择一个问题，必要时填写癌种、组织、细胞类型等背景。
3. 选择解释后端。默认 `本地安全模板` 最稳；外部 LLM 可在页面里填写一次性 API 配置，且只负责叙述，不参与评分或改图。
4. 点击“开始分析”。
5. 先看“输入检查”和“实体匹配队列”，再看解释、文献证据、通路/靶点/疾病排序和研究优先级提示。

## 二、输入格式

推荐上传带表头的 CSV/TSV：

```csv
metabolite,log2FC,pvalue,padj,direction
glucose,1.2,0.01,0.03,up
lactate,1.1,0.03,0.04,up
glutamine,0.6,0.04,0.08,up
```

也可以使用 JSON：

```json
{
  "records": [
    {"HMDB": "HMDB0000122", "log2FC": 1.2, "padj": 0.01, "direction": "up"},
    {"ChEBI": "CHEBI:16651", "log2FC": 1.1, "padj": 0.02, "direction": "up"}
  ]
}
```

支持的常用列：

| 类别 | 列名示例 | 用途 |
| --- | --- | --- |
| 名称 | `name`, `metabolite`, `compound` | 按名称/同义词匹配代谢物 |
| 外部 ID | `HMDB`, `ChEBI`, `PubChem CID`, `KEGG`, `InChIKey` | 更稳定的实体规范化 |
| 质谱特征 | `formula`, `mz`, `adduct`, `ppm_tolerance` | 用于质量/分子式辅助匹配 |
| 差异信息 | `log2FC`, `pvalue`, `padj`, `direction` | 用于方向性和种子权重 |
| 背景 | `cancer_context`, `cancer_type`, `tissue`, `cell_type`, `cell_line` | 用于 overlay 排名加权 |

如果输入只有自然语言，服务会尽量解析类似 `glucose 上调 log2FC=1.2 p=0.01` 的片段。正式分析建议使用表格或 JSON。

## 三、界面结果怎么看

| 区块 | 含义 | 建议动作 |
| --- | --- | --- |
| 输入检查 | 输入数、匹配数、歧义数、未匹配数、文献证据数 | 如果歧义或未匹配很多，先修正输入 |
| 可信度阅读卡 | 将输入匹配、核心通路和扩展假设分层 | 先看总体可信度，再决定是否只作探索性解释 |
| 提示与阻断原因 | 质量警告、解释阻断、安全 guard 信息 | 高风险提示需要人工复核 |
| 实体匹配队列 | 已匹配、歧义、未匹配三类输入 | 歧义项不会进入评分 |
| 解释摘要 | 本地或外部解释器生成的证据绑定叙述 | 每段都应有来源引用 |
| 文献证据包 | 支持关系数、PMID、支持类别、最高文献概率 | 用于回查原始文献证据 |
| 交互证据图谱 | 可点击浏览输入代谢物、通路、靶点和疾病/表型关联 | 用于理解分析流，不代表相邻节点都是直接反应边 |
| 通路/靶点/疾病排序 | 图传播、通路富集和证据加权后的排序 | 看分数，也看“证据可追溯” |
| 解释路径 | 从输入代谢物到终端节点的图谱路径 | 用于解释为什么排到前面 |
| 研究优先级提示 | 药物/细胞背景 overlay 的研究提示 | 只能作为研究优先级，不是临床建议 |
| 原始 JSON | 完整接口返回 | 调试、复现、下游集成使用 |

## 四、主要功能

- 数据源下载：`download_all.ps1` 和 `scripts/download_databases.py` 根据 `config/source_catalog.toml` 抓取开放核心数据源，并生成 manifest。
- 规范化数据仓：`scripts/build_normalized_store.py` 生成代谢物、基因、通路、疾病、靶点、文章和句子等规范表。
- 图谱投影：`scripts/build_graph_projection.py` 生成节点、边、解析索引和稀疏图。
- 化合物匹配索引：`scripts/build_compound_match_index.py` 支持名称、外部 ID、InChIKey、分子式、m/z、RT/MS2 扩展匹配。
- 文献证据层：`scripts/build_literature_evidence.py` 生成句级提及、候选关系和边支持证据。
- 分析服务：`scripts/metabo_service.py` 提供只读 API、CLI 和 Web 工作台。
- 安全解释层：`scripts/llm_safe_adapter.py` 把结构化分析改写为有来源引用的解释，禁止创建事实、修改评分或写入图谱。
- 预测 overlay：`manual_sources/prediction_overlays/<release_id>/` 可放入药物靶点和细胞背景表，作为研究优先级提示。
- 回归验证：`scripts/run_phase15_validation.py`、`scripts/run_llm_safe_adapter_validation.py` 等用于发布前验收。

## 五、API 与 CLI

启动服务后可用端点：

| 端点 | 方法 | 用途 |
| --- | --- | --- |
| `/` 或 `/chat` | GET | Web 工作台 |
| `/releases` | GET | 当前 release、表行数、图谱规模 |
| `/resolve` | POST | 单实体解析 |
| `/precheck/metabolites` | POST | 只做代谢物输入预检查 |
| `/analyze/metabolites` | POST | 分析代谢物表，返回 analysis pack |
| `/analyze/differential-table` | POST | 分析已计算好的差异结果表，不读取原始丰度矩阵 |
| `/explain` | POST | 分析并生成安全解释 |
| `/chat` | POST | Web 工作台使用的综合入口 |
| `/llm/status` | GET | 查看服务级外部 LLM 配置状态 |
| `/llm/test` | POST | 用一次性或服务级配置发送极小测试请求，诊断鉴权、代理、DNS、超时等问题 |
| `/entity/<id>` | GET | 实体详情 |
| `/pubchem/<cid>` | GET | PubChem CID 缓存详情 |
| `/subgraph` | GET/POST | 返回局部子图 |
| `/evidence` | GET/POST | 查询边、实体、句子的证据 |

CLI 示例：

```powershell
python .\scripts\metabo_service.py --workspace . --release-id mvp_20260513T002254 releases
python .\scripts\metabo_service.py --workspace . --release-id mvp_20260513T002254 resolve glucose --entity-type metabolite
python .\scripts\metabo_service.py --workspace . --release-id mvp_20260513T002254 analyze-metabolites .\records.json --max-paths 25 --max-hops 4
python .\scripts\metabo_service.py --workspace . --release-id mvp_20260513T002254 analyze-differential-table .\TraitScore_GroupDiff_allCelltypes\significant\trait_score_diff_LUAD_Epi_LUAD_Tumor_vs_Adjacent_significant_q0.05.csv --max-records 80 --top-per-group 80
python .\scripts\metabo_service.py --workspace . --release-id mvp_20260513T002254 explain .\records.json --question "请解释主要通路和证据"
```

`TraitScore_GroupDiff_allCelltypes/` 里的 CSV 推荐使用 `analyze-differential-table`，
不要用普通 `analyze-metabolites`。该入口会把两组比较表中的
`trait/GCST`、`group1/group2`、`mean_diff`、`pseudo_log2FC_shifted`、
`z_wilcoxon` 和 q 值列转换成 TraitScore 分析输入，并在结果里标出严格
代谢物匹配、低权重扩展候选、歧义项和未匹配项。TraitScore 结果是研究候选
解释，不是直接 LC-MS 丰度结论。

更推荐的新入口是 `analyze-differential-table`。它把输入明确视为“已经计算好的差异结果表”，
只读取代谢物/trait 标识、效应量、P 值、FDR/q 值、方向和分组元数据，不需要原始丰度矩阵。
如果没有标准 `log2FC`，会把 `mean_diff`、`cohen_d`、`z_wilcoxon` 等保留为带标签的方向性效应量，
仅作为图谱种子权重和研究优先级信号，不改写成实测丰度结论。

Web 工作台选择 `差异结果表` 或 `TraitScore 差异表` 时，提交入口仍是 `/chat`：
服务会先调用差异表选择/解析逻辑生成冻结 `analysis_pack`，再把该结果交给 `/explain`
的本地模板或外部 LLM 叙述层。

## 六、重建或更新发布

常规顺序：

```powershell
.\download_all.ps1 -Run -ResolveIndexes
.\build_mvp.ps1 -ReleaseId mvp_20260513T002254
python .\scripts\build_graph_projection.py --release-id mvp_20260513T002254
python .\scripts\build_pubchem_cid_cache.py --release-id mvp_20260513T002254
python .\scripts\build_compound_match_index.py --workspace . --release-id mvp_20260513T002254
python .\scripts\build_literature_evidence.py --workspace . --release-id mvp_20260513T002254
python .\scripts\run_phase15_validation.py --workspace . --release-id mvp_20260513T002254
```

可选优化：

- 有 RT/MS2 手工特征时，放入 `manual_sources/compound_features/<release_id>/` 后重建 compound match index。
- 有药物靶点或细胞背景 overlay 时，放入 `manual_sources/prediction_overlays/<release_id>/`，或运行 `scripts/build_prediction_overlays.py`。
- 文献证据规则变更后，先运行 `scripts/run_evidence_precision_qa.py`，再决定是否进入发布。

## 七、外部 LLM

外部 LLM 默认关闭。启用时仍然只允许读 `/analyze/metabolites`、`/analyze/differential-table`、`/evidence`、`/subgraph`、`/releases` 的冻结结果，输出必须通过本地 guard。

### 方式 A：页面内一次性配置

在 Web 工作台里选择 `外部 LLM（文字）`、`外部 LLM（分层文字）` 或 `外部 LLM（JSON）`，展开“外部 LLM API 设置”，填写：

- `Endpoint`：OpenAI 兼容 chat completions 地址，例如 `https://api.openai.com/v1/chat/completions`
- `Model`：模型名
- `API key`：本次请求临时使用
- `请求超时（秒）`：网络慢或模型响应慢时可提高到 120-300
- `代理 URL`：本机代理地址，例如 `http://127.0.0.1:7890`

这种方式不需要重启服务。API key 只随本次请求发送到本机服务进程，不写入 manifest、返回 JSON 或审计哈希。建议先点“测试 LLM 连接”，确认连接、鉴权、模型名和代理都正常，再用 `外部 LLM（文字）` 正式分析；文字模式对模型返回格式要求更低。大差异表或 TraitScore 表推荐 `外部 LLM（分层文字）`：系统会把解析质量、ratio/class/identity、排名、证据和低置信附录分块送入模型，每块先过本地 guard，再用通过校验的分块生成总述。

可检查当前服务级 LLM 状态：

```text
http://127.0.0.1:8765/llm/status
```

也可以用接口自检一次性配置：

```powershell
Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8765/llm/test -ContentType 'application/json' -Body (@{
  llm_config = @{
    endpoint = "https://api.deepseek.com/chat/completions"
    model = "deepseek-v4-pro"
    api_key = "<secret>"
    proxy_url = "http://127.0.0.1:7890"
    timeout_seconds = 120
  }
} | ConvertTo-Json -Depth 4)
```

### 方式 B：环境变量启动

如果希望每次请求默认复用同一套外部 LLM 配置，可以在启动服务前设置环境变量：

```powershell
$env:LLM_SAFE_ADAPTER_ENABLE_EXTERNAL = "1"
$env:LLM_SAFE_ADAPTER_MODEL = "your-narrator-model"
$env:LLM_SAFE_ADAPTER_API_KEY = "<secret>"
python .\scripts\metabo_service.py --workspace . --release-id mvp_20260513T002254 --enable-external-llm --llm-model your-narrator-model serve --port 8765
```

如果没有配置外部模型或 API key，界面选择外部后端会返回 guard 阻断信息，这是预期行为。

## 八、常见问题

- 页面显示 `release: unavailable`：确认服务是通过 `metabo_service.py serve` 启动，而不是直接打开 HTML 文件。
- 结果里待复核很多：优先使用 HMDB、ChEBI、PubChem CID、InChIKey 等稳定 ID，减少只用通用名称。
- 通路或靶点为空：检查输入是否匹配成功，以及 `graph_projection/<release_id>/` 是否存在。
- 文献证据为空：检查 `literature_evidence/<release_id>/` 是否存在，或降低分析样本对证据覆盖的期待。
- 外部 LLM 被阻断：先点页面里的“测试 LLM 连接”或调用 `/llm/test`。`401/403` 通常是 API key 问题，`404` 多半是 endpoint 或 model 问题，`transport_error` 通常是 DNS、代理、TLS 或超时。
