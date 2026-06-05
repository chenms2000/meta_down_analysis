# 优化后可选代谢组分析插件论文说明草案

> 使用建议：本稿适合放在论文“方法学/可选分析模块”或补充方法中。文中的 release ID、GitHub URL、验证集规模和性能指标请在最终投稿前替换为正式版本。

## 论文正文推荐版本

### 可选分析插件：证据约束的代谢组解释、候选机制排序与可复核输出

为便于读者对代谢算法输出进行进一步生物学解释，本研究提供一个可选的代谢组分析插件。该插件作为主算法之后的独立后处理模块使用，不参与主算法训练、参数选择或核心性能评估；读者可根据研究目的选择是否运行该插件，未运行插件不影响本文主算法结果的复现。插件的目标是将代谢物列表或差异代谢表映射到冻结版本的代谢知识图谱中，并生成带证据来源、置信等级和自动化质量状态的通路、靶点、疾病及文献支持候选结果。

插件以版本化知识图谱为背景运行。知识图谱由公开数据库、规范化实体表、代谢物-反应-基因/蛋白-通路-疾病/靶点关系、可选谱学特征索引以及句级文献证据层组成。当前文献层采用 PubMed/PMC 文本资源的 title/abstract 句级抽取，并保留本地全文文本路径、字节数和 checksum 作为审计线索；若需要全文句级挖掘，应以 `--article-sentence-scope full` 单独重建。输入可以是代谢物名称、外部数据库 ID、InChIKey、分子式、m/z、保留时间、MS/MS 特征，或包含 log2 fold change、P 值、校正 P 值和变化方向的差异代谢表。插件首先执行输入质控和实体规范化，将每个输入条目分为已匹配、歧义、未匹配或无效四类。只有达到预设匹配分数和候选间隔阈值的代谢物进入后续图评分；歧义和未匹配条目以自动化审计队列形式输出，不会被强制纳入传播分析。

实体解析采用多证据加权策略，综合外部 ID、结构键、名称/同义词、分子式/质量、RT、MS/MS 和图谱支持信息。针对同分异构体、泛化名称、常见高频同义词和质量近似候选，插件采用降权或阻断规则，以降低错误匹配对后续解释的影响。由于当前阶段暂时无法进行人工 adjudication，实体解析暂不报告人工精标准确率，而报告固定 gold/ambiguity fixture、compound-name 风险审计、abstention 数量和错误进入评分面的阻断情况；人工 adjudication 将作为后续优化环节补充，正式发布时同时保留解析阈值、配置哈希和验证集版本。

对于大规模 compound-name 风险项，插件不要求逐条人工确认，而是采用确定性 triage 策略：高风险 name-only 表面自动 abstain，脂质简称、异构体或多实体映射名称要求稳定数据库 ID、InChIKey、精确结构、RT/MS2 或其他正交证据，分隔符敏感名称仅作为格式风险提示，低优先级同义词表面进入 monitor/downweight 队列。该策略不等同于人工判定正确性，只是在人工 adjudication 暂不可用时阻止高风险名称直接进入确定性解析和高置信评分面。

在分析阶段，插件将已匹配代谢物转换为图谱种子节点，并根据差异方向、显著性和重复输入情况生成种子权重。随后在类型约束的异构知识图谱上进行传播：

\[
\pi=(1-\alpha)y+\alpha W\pi
\]

其中 \(y\) 表示输入代谢物种子向量，\(W\) 表示经边类型门控、归一化和度惩罚后的转移矩阵，\(\pi\) 为传播后的节点得分。插件同时计算通路覆盖、传播得分、文献支持、方向一致性和路径解释成本。解释路径采用边置信度负对数和定义：

\[
C(path)=\sum_{e\in path}-\log(p_{final}(e))
\]

其中 \(p_{final}(e)\) 由结构化数据库证据、文献证据、拓扑支持和用户输入支持共同构成。为减少高度连接节点和重复路径带来的解释冗余，插件对传播图采用类型 beam 压缩，并对解释路径进行 path signature 聚类，仅将代表性路径返回给用户。

文献证据层作为 overlay 使用，而非规范图谱事实来源。该层从本地文献语料中抽取句级实体提及和候选关系，并根据句内距离、提及密度、实验语境、背景语境、否定/不确定表达、综述性描述和文章内重复证据进行加权。多个来自同一 PMID/PMCID 的句子证据按文章聚合，避免重复句子导致证据膨胀。文献证据被标记为 `confirm`、`support_direction`、`novel_candidate`、`conflict_candidate` 或 `unresolved`。其中，`novel_candidate` 和 `conflict_candidate` 只能作为自动化假设提示，不会自动写入规范图谱，也不会被表述为已验证事实。OpenAI-compatible 或 DeepSeek 等外部模型仅作为可选叙述层，不参与证据抽取、实体决策或图谱写入。

插件输出结构化的 `analysis_pack` 和 `interpretation_report`，包括输入检查摘要、已匹配/歧义/未匹配队列、机制事实候选、通路排序、靶点排序、疾病排序、解释路径、文献证据包、质量警告、证据引用、上下文降级信息以及可复现性哈希。核心候选结论优先来自 `database_accuracy_store.v2` 中的 `mechanism_ready_facts`；通路、靶点、疾病和药物 ranking 作为 supporting hypotheses 或附录式候选保留。所有候选排序均带有 `claim_refs` 和 `evidence_refs`，至少回链到图边、来源记录或句级证据。若候选结果缺少证据回链，或上下文与输入背景存在明显不匹配，则结果会被降级为探索性输出或进入附录式候选队列。

为提高候选排序的可解释性和可比较性，插件将输出划分为 high、medium、exploratory 和 low 四个置信等级。该等级不是生物学真值概率，而是研究优先级分层，依据输入支持数、实体匹配质量、机制事实类型、方向一致性、图距离、文献支持和校准状态确定。`role_unknown_reaction_fact` 只支持“参与可追溯反应”的表述；只有 `directional_reaction_fact` 等方向明确事实才可作为方向相关候选线索，且仍需正交验证。当前阶段因条件限制暂未进行人工 adjudication，排序和置信等级先通过时间切分验证、来源留一验证、证据可追溯性检查、LLM guard 回归和保守降级规则进行代理评估；论文中应明确说明人工证据判读是后续优化环节，并报告 held-out 指标、证据回链通过率和 unsupported-claim rate。

该插件仅用于研究解释和假设生成，不用于临床诊断、治疗决策或因果证明。可选语言模型解释层只允许读取冻结的 `analysis_pack`、`interpretation_report`、证据包、子图和 release 元数据，并将结构化结果改写为带来源引用的自然语言摘要。语言模型不能参与实体解析最终裁决、证据选择、打分、图谱写入或事实创建；若生成内容未能通过本地证据绑定和禁止字段检查，则解释文本会被阻断。

## 具体流程说明

### 1. 输入与预检查

插件支持三类输入：

1. 代谢物列表：用于快速评估通路和机制覆盖。
2. 差异代谢表：用于方向性解释、通路排序和候选机制优先级排序。
3. 带谱学特征的候选表：用于辅助区分同分异构体或低特异性名称。

推荐输入字段如下：

| 字段类别 | 示例列名 | 用途 |
|---|---|---|
| 名称 | `name`, `metabolite`, `compound` | 名称或同义词匹配 |
| 外部 ID | `HMDB`, `ChEBI`, `PubChem CID`, `KEGG`, `InChIKey` | 稳定实体解析 |
| 谱学特征 | `formula`, `mz`, `adduct`, `rt`, `ms2_peaks` | 质量、保留时间和碎片辅助匹配 |
| 差异信息 | `log2FC`, `pvalue`, `padj`, `direction` | 种子权重和方向一致性 |
| 背景信息 | `cancer_type`, `tissue`, `cell_type`, `cell_line`, `cell_state` | 上下文加权或不匹配降级 |

预检查阶段输出输入数量、有效行数、重复行、缺失字段、已匹配数、歧义数、未匹配数和无效行数。若输入主要由泛化名称、无 ID 条目或低置信质量特征组成，插件会给出质量警告，并将排序解释降级为探索性。

### 2. Release 与数据层

每次分析均绑定到一个冻结 release。推荐 release 结构如下：

| 层级 | 内容 | 作用 |
|---|---|---|
| `raw_lake` | 原始数据库文件、下载 URL、许可、checksum | 审计和重建 |
| `normalized_store` | 规范实体、交叉引用、来源记录、基础关系 | 事实来源层 |
| `graph_projection` | 节点、边、边类型、稀疏图和解析索引 | 图传播和解释路径 |
| `compound_match_index` | 名称、ID、InChIKey、分子式/质量、RT、MS/MS 索引 | 实体解析 |
| `literature_evidence` | 句级提及、候选关系、文献边支持 | 证据 overlay |
| `validation_reports` | 回归测试、验证集结果、校准报告 | 发布门槛 |

正式论文中应报告 release ID、构建日期、主要数据源版本、节点数、边数、文献语料规模、验证报告哈希和代码 commit hash。

### 3. 实体规范化优化

每条输入记录生成多个候选代谢物，并计算综合解析分数。候选证据包括：

1. 外部 ID 精确命中。
2. InChIKey 完整或骨架命中。
3. 名称、同义词和常见生化名称命中。
4. 分子式、精确质量、m/z 和 adduct 匹配。
5. RT 和 MS/MS 谱图辅助匹配。
6. 图谱中与代谢通路、反应或文献证据的支持程度。

默认只有 top1 分数达到阈值且 top1-top2 margin 足够时才判定为已匹配。优化后建议增加以下控制：

| 风险类型 | 处理策略 |
|---|---|
| 同分异构体 | 需要 RT/MS2 或稳定 ID 支持，否则保留为歧义 |
| 泛化名称，如 lipid、hexose、fatty acid | 不直接映射为单一实体，可映射到 biological pool 或自动化审计队列 |
| 高频别名，如 common ion、generic synonym | 降权或阻断 |
| 仅质量匹配且候选过多 | 不进入评分，仅输出候选范围 |
| 名称与结构/质量冲突 | 降级或阻断 |

建议验证指标：

| 指标 | 说明 |
|---|---|
| top-1 accuracy | 已匹配条目的首位候选正确率 |
| abstention precision | 被标为歧义/未匹配的条目中确实不应强制匹配的比例 |
| ambiguous recall | 高风险输入被成功拦截为歧义的比例 |
| false-match rate | 错误实体进入图评分的比例 |
| abstention burden | 每 100 个输入中被自动保留为歧义/未匹配、未进入评分的条目数 |

### 4. 图传播与候选排序

已匹配代谢物被转换为种子节点。若输入包含差异信息，种子权重综合 log2FC、显著性和方向；若输入为普通列表，则使用列表模式并降低方向性解释权重。

排序分为三类：

1. 通路排序：由输入覆盖、通路富集、传播分数、文献支持和方向一致性共同决定。
2. 靶点排序：由代谢物-基因/靶点关系、传播分数、文献支持和图距离共同决定。
3. 疾病/表型排序：由代谢物、通路、靶点到疾病节点的证据路径和传播分数共同决定。

所有候选均带有：

| 字段 | 含义 |
|---|---|
| `score_components` | 富集、传播、文献、方向和上下文分量 |
| `confidence_tier` | high、medium、exploratory 或 low |
| `calibrated_confidence` | 校准后的研究优先级分数 |
| `boundary` | 该结果可如何被解释的边界说明 |
| `claim_refs` | 图边、来源记录或文献支持引用 |
| `evidence_refs` | 具体证据对象、PMID/PMCID、句子或边支持 |
| `context_mismatch` | 是否与显式上下文不匹配 |
| `needs_validation` | 是否需要额外实验、外部证据或未来人工复核；当前不会被表述为已验证事实 |

### 5. 文献证据精度控制

文献证据层建议使用以下精度控制：

| 控制项 | 目的 |
|---|---|
| 句内距离 | 降低远距离共现噪声 |
| 提及密度 | 降低实体堆叠句造成的假阳性 |
| 否定/不确定识别 | 区分支持、反对和未确定关系 |
| 背景/综述降权 | 避免把背景知识当成直接实验证据 |
| 直接实验语境加权 | 提高 assay、knockdown、inhibition 等实验句权重 |
| PMID/PMCID 聚合 | 避免同一文章多句重复放大 |
| novel/conflict 隔离 | 新关系和冲突关系不进入规范事实层 |

无人工 adjudication 时，建议使用以下自动化代理指标：

| 指标 | 说明 |
|---|---|
| evidence traceability pass rate | 证据是否能回链到句子、PMID/PMCID、关系候选和来源记录 |
| overlay isolation pass rate | `novel_candidate` 和 `conflict_candidate` 是否保持 overlay-only、不进入规范评分 |
| precision-filter stability | 高风险 surface block/downweight 后，既有边支持是否没有异常塌陷 |
| duplicate inflation rate | 同一文章重复句子导致的证据膨胀比例 |
| unsupported-claim rate | 自然语言解释中缺少来源支持的声明比例 |

### 6. 上下文建模与降级规则

当用户提供癌种、组织、细胞类型、细胞系或细胞状态时，插件采用上下文加权策略：

1. `soft` 模式：匹配上下文的候选加权，不匹配候选降权但保留。
2. `hard` 模式：仅保留显式匹配上下文的 overlay 候选。

优化后建议显式区分：

| 证据层级 | 解释策略 |
|---|---|
| 泛生物学证据 | 可作为背景机制，不能直接声称特定癌种相关 |
| 泛癌证据 | 可作为肿瘤相关候选，但需标明非特定癌种 |
| 癌种特异证据 | 可提高上下文置信度 |
| 组织/细胞类型特异证据 | 可用于上下文排序，但需标明来源 |
| 细胞系证据 | 作为模型系统证据，不等同于患者组织 |

若结果与显式上下文不符，例如神经元、免疫细胞或造血特异路径出现在上皮肿瘤背景中，插件应保留该候选但标记 `context_mismatch = true`，并将其降级为探索性或附录候选。

### 7. 暂无人工 adjudication 时的代理验证模式

如果当前研究条件暂时不允许人工判读文献或逐条校验实体，插件应采用代理验证模式。该模式不试图给出人工精标准确率，而是通过更保守的自动化边界降低误用风险；人工 adjudication 应作为后续优化环节继续补充：

1. 歧义实体不强行解析，保留为 `ambiguous`、`unmatched` 或 biological pool。
2. `novel_candidate` 和 `conflict_candidate` 文献证据保持 overlay-only，不进入规范图谱事实和强评分面。
3. 高置信结果必须具备结构化数据库边、清洁输入支持、证据回链和校准字段；长距离图传播、上下文不匹配、疾病/药物候选默认降级。
4. 使用 temporal holdout、source-held-out、compound name audit、compound-name deterministic triage、evidence contract、LLM guard regression 作为代理验证。
5. 论文中明确声明：本插件当前阶段暂未进行人工 adjudication，后续将以人工抽样复核或专家判读进一步优化；当前结果用于研究优先级排序和假设生成，不作为已验证生物学事实。

### 8. 结果校准与验证

优化后版本应将验证作为插件说明的一部分。推荐至少包含：

| 验证类型 | 目的 |
|---|---|
| 实体解析 gold set | 评估输入代谢物匹配准确性和歧义拦截能力 |
| 文献证据代理验证 | 在人工 adjudication 暂不可用时，评估证据回链、来源记录、overlay 隔离和精度过滤稳定性 |
| 时间切分验证 | 检验高分候选是否更容易被后续证据支持 |
| 来源留一验证 | 检验模型能否恢复被暂时移除的数据源关系 |
| 回归测试 | 保证相同输入、release 和配置返回相同结果 |
| 安全解释测试 | 保证自然语言解释不创建事实、不修改分数、不脱离证据 |

论文中可使用如下表格占位：

| 模块 | 验证集/方法 | 指标 | 结果 |
|---|---|---|---|
| 实体解析 | [gold set 名称，n=] | top-1 accuracy / false-match rate | [填入] |
| 歧义拦截 | [ambiguous set，n=] | ambiguous recall / review burden | [填入] |
| 文献证据 | [evidence contract + precision filters] | traceability / overlay isolation / unsupported-claim rate | [填入] |
| 排序校准 | [temporal holdout] | high-vs-low enrichment trend | [填入] |
| 来源泛化 | [source-held-out] | held-out recovery / top-fraction hit rate | [填入] |
| 解释安全 | [adversarial fixtures] | unsupported-claim rate | [填入] |

### 9. 输出契约

插件核心输出为 `analysis_pack.v1`。建议在论文或补充材料中说明以下字段：

| 字段 | 内容 |
|---|---|
| `input_summary` | 输入规模、匹配状态、重复项、种子权重 |
| `matched` | 成功解析并进入评分的代谢物 |
| `ambiguous` | 自动判定为证据不足、不进入评分的候选；可供未来复核 |
| `unmatched` | 当前 release 无法解析的输入 |
| `biological_entity_pools` | 泛化或池化生物实体解释层 |
| `quality_warnings` | 输入、解析、传播、上下文和证据警告 |
| `pathway_rankings` | 通路候选及其置信等级和证据 |
| `target_rankings` | 靶点候选及其置信等级和证据 |
| `disease_rankings` | 疾病/表型候选及其置信等级和证据 |
| `prediction_model` | 校准状态、置信分层、上下文降级和附录候选 |
| `top_explanation_paths` | 代表性解释路径 |
| `literature_evidence_pack` | 文献支持摘要、PMID/PMCID、句级证据 |
| `evidence_refs` | 所有引用的图边、来源记录和文献证据 |
| `determinism` | 输入哈希、release ID、配置哈希和结果哈希 |
| `blocked_reasons` | 阻断评分或解释的原因 |

### 10. 限制说明

建议在论文中明确列出以下限制：

1. 插件依赖当前 release 中已有的数据源，缺失数据库或许可限制会影响覆盖率。
2. 名称解析不能替代实验级化合物鉴定；缺少 MS/MS、RT 或稳定 ID 时，同分异构体解释应谨慎。
3. 文献证据可能受到发表偏倚、综述重复描述和共现噪声影响。
4. 图传播结果表示研究优先级，不等于因果关系或临床可操作性。
5. 上下文 overlay 依赖可用的癌种、组织、细胞类型或细胞系证据；在人工 adjudication 尚未完成前，具体样本外推应降级为探索性。
6. 语言模型输出只是证据绑定的叙述层，不能创建事实或替代结构化结果。
7. 当前阶段暂未进行人工 adjudication，因此不能声称文献抽取精度等同于人工精标准确率；现阶段只能报告自动化代理验证和保守阻断结果，人工 adjudication 作为后续优化方向。

### 11. GitHub 发布说明

若将插件放到 GitHub，建议 README 中明确：

1. 本插件是 optional analysis plugin，不是主算法复现的必要步骤。
2. 插件用于研究解释、证据回链和候选机制优先级排序，不用于临床决策。
3. 数据 release 与代码 release 分离；受许可限制的数据不直接提交到仓库。
4. 每个 release 记录数据源、版本、许可、下载日期、checksum、配置哈希和验证报告。
5. 仓库提供最小示例输入、输出 JSON schema、验证命令和已知限制。
6. 论文引用固定 tag、commit hash 或 Zenodo DOI。

可在论文数据可用性部分写：

> The optional analysis plugin, input templates, output schema, validation fixtures, and reproducibility scripts will be released at [GitHub URL] under version [tag/commit/DOI]. The plugin is not required to reproduce the core algorithmic results of this study. It is provided as a research-oriented interpretation layer for readers who wish to inspect pathway, target, disease, and literature-evidence links derived from their own metabolite tables.
