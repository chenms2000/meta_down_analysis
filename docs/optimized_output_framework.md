# 代谢库输出端优化框架

## 一、目标定位

输出端从“知识库浏览器”升级为“研究解释报告”。

核心目标不是展示更多字段，而是把已有的代谢物匹配、图谱关系、文献证据、排序结果和置信度组织成用户能读懂的结论链。

输出端必须遵守三条边界：

- 只做研究解释和候选优先级提示，不做临床决策。
- 每个自然语言结论都应能回到输入代谢物、图谱边、文献证据或排序字段。
- 对歧义、弱证据、图传播和模型排序结果必须明确降级，不把假设写成事实。

## 二、第一屏信息结构

默认第一屏按以下顺序展示：

1. 核心结论
2. 结论链路
3. 关键节点解释
4. 证据与置信度原因
5. 不确定性与复核队列
6. 下一步验证建议

详细节点、原始 JSON、完整子图和证据表继续保留，但作为可展开内容，不作为第一阅读入口。

## 三、输出生成流程

### 1. 输入整理

输入来源：

- 用户上传或粘贴的代谢物表
- `/precheck/metabolites` 的匹配结果
- `/analyze/metabolites` 的 `analysis_pack`
- `/evidence` 的句级证据
- `/subgraph` 的局部图谱
- 可选 context，例如癌种、组织、细胞类型、细胞系

处理动作：

- 统计 matched、ambiguous、unmatched、invalid。
- 标记稳定 ID 输入、名称输入、质量/分子式输入、RT/MS2 输入。
- 把输入代谢物按方向、显著性、重复命中和生物主题分组。
- 把 ambiguous 和 unmatched 进入复核队列，不直接参与强结论。

### 2. 结论候选生成

从以下来源生成候选结论：

- pathway rankings
- target rankings
- disease rankings
- prediction model metabolic themes
- top explanation paths
- literature evidence pack
- biological entity pools

每个候选结论必须形成统一结构：

```json
{
  "claim_id": "claim_001",
  "claim_type": "pathway_theme",
  "headline": "输入代谢物显示糖酵解/乳酸代谢相关信号增强",
  "supporting_inputs": [],
  "supporting_nodes": [],
  "supporting_edges": [],
  "evidence_refs": [],
  "confidence_tier": "medium",
  "confidence_reasons": {
    "positive_factors": [],
    "downgrade_factors": []
  },
  "boundary": "支持代谢主题富集，不等同于证明通路活性增强",
  "next_validation": []
}
```

### 3. 结论链组织

每条结论链按固定模板展示：

```text
主要结论：
输入代谢物显示某一代谢主题或通路相关信号。

依据：
哪些输入代谢物、方向和显著性共同支持该结论。

关系链：
输入代谢物 -> 通路/反应 -> 靶点/疾病/候选机制。

可能机制：
用谨慎语言解释可能的生物学含义。

相关候选：
列出相关通路、基因/靶点、疾病背景或药物候选。

可信度：
高/中/探索性/低，并解释加分和降级原因。

解释边界：
说明当前证据能支持什么，不能支持什么。

建议验证：
给出 2-4 个可执行的实验或数据验证方向。
```

### 4. 结论排序

排序优先级建议：

1. 输入覆盖度高的代谢主题
2. 有稳定数据库 ID 支持的结论
3. 有 curated graph edge 的结论
4. 有句级文献证据支持的结论
5. 与用户给定癌种、组织、细胞类型一致的结论
6. 图距离较短、解释路径更直接的结论

降级规则：

- 只由名称匹配支持：降级。
- 只由图传播支持：降级。
- 只由模型排序支持：降级。
- 与 context 不匹配：降级。
- novel/conflict literature evidence：作为探索性提示，不进入强结论。
- ambiguous 输入只能支持主题覆盖，不能支持精确化学身份结论。

## 四、节点说明卡

每个关键节点生成说明卡，用于替代生硬的节点详情。

### 节点说明卡字段

```json
{
  "node_uid": "",
  "display_name": "",
  "node_type": "metabolite | pathway | reaction | gene | target | disease | drug",
  "plain_language_definition": "",
  "role_in_network": "",
  "why_in_this_analysis": "",
  "relation_to_input": "",
  "evidence_source_summary": "",
  "confidence_tier": "",
  "confidence_reason": "",
  "boundary": "",
  "open_details": {
    "xrefs": [],
    "edges": [],
    "evidence_refs": []
  }
}
```

### 节点类型解释规则

代谢物节点：

- 说明它是否来自用户输入。
- 说明匹配方式：稳定 ID、名称、InChIKey、formula/mass、RT/MS2。
- 说明方向和显著性。
- 标记是否存在歧义或重复命中。

通路节点：

- 说明通路的简明生物学含义。
- 说明命中的输入代谢物。
- 说明是直接通路命中、富集结果还是传播结果。
- 明确“通路富集”不等于“通路活性已被证明”。

靶点/基因节点：

- 说明它通过哪些代谢物、通路或文献证据出现。
- 标记直接证据、图传播证据、文献 overlay 或 prediction overlay。
- 若来自长距离传播，默认写成候选机制，不写成确定机制。

疾病节点：

- 说明它是疾病背景关联，不是诊断结论。
- 如果用户给了癌种 context，要标记是否匹配。
- 如果是泛癌或疾病库关联，应降级为背景提示。

药物节点：

- 只能作为研究候选或靶点相关药物提示。
- 必须继承上游靶点置信度，不能比靶点更高。
- 禁止表达为治疗建议。

## 五、置信度解释器

置信度不只显示分数，还要显示为什么。

### 推荐分层

高：

- 输入匹配可靠，多个稳定 ID 支持。
- 多个输入代谢物指向同一通路或机制。
- 有 curated graph edge。
- 有句级文献证据。
- 与用户 context 一致。

中：

- 输入覆盖较好，但部分是名称匹配。
- 有 curated graph 或文献之一。
- 图路径较短，但缺少当前样本的独立验证。

探索性：

- 主要来自图传播、弱监督排序或 overlay。
- 证据不直接，或 context 不完全匹配。
- 可作为候选，但不应写成结论。

低：

- 输入匹配不稳定。
- 证据来自歧义名称、长距离传播或泛化关系。
- 缺少可追溯证据。

### 加分因素

- 多个输入代谢物命中同一通路。
- 输入有 HMDB、ChEBI、PubChem CID、InChIKey 等稳定 ID。
- 方向一致，例如相关代谢物共同上调或下调。
- 有 curated graph edge。
- 有 sentence-level literature evidence。
- 与癌种、组织、细胞类型 context 匹配。
- 解释路径短，图距离清楚。

### 降级因素

- 只靠名称匹配。
- 输入存在异构体、脂质、复合名、pool、ratio 或多实体歧义。
- 只由图传播或模型排序支持。
- 文献证据是 novel_candidate 或 conflict_candidate。
- 缺少当前癌种或细胞背景证据。
- 解释路径太长或中间节点过于泛化。

## 六、页面组件框架

### 1. 核心结论区

展示 3-5 条最重要结论。

每条结论包含：

- 一句话结论
- 可信度标签
- 输入支持数
- 证据数
- 是否需要复核
- 展开按钮

### 2. 结论链详情区

点击核心结论后展示：

- 依据代谢物
- 关系链
- 可能机制
- 关键节点
- 文献证据
- 置信度原因
- 解释边界
- 下一步验证

### 3. 关键节点说明区

按“输入代谢物、通路、靶点、疾病、药物候选”分组展示节点卡。

每张卡默认显示：

- 节点名称和类型
- 它为什么出现在本次分析中
- 与输入的关系
- 可信度和边界

详细 ID、边列表、xref、原始证据放在展开层。

### 4. 证据区

证据按用途组织，而不是按原始表格堆叠：

- 支持核心结论的证据
- 支持关键节点的证据
- 支持候选机制的证据
- 冲突或不确定证据

每条证据显示：

- 来源类型：curated graph、literature sentence、overlay、prediction
- 支持对象
- 支持强度
- PMID/PMCID 或 source record
- 原始句子或摘要

### 5. 复核队列

单独展示：

- ambiguous 输入
- unmatched 输入
- invalid 输入
- 只靠名称匹配的高风险输入
- 需要稳定 ID 或 RT/MS2 辅助确认的输入

复核队列不应藏在底部，因为它直接决定结论可信度。

### 6. 下一步验证区

根据结论类型生成建议：

通路结论：

- 检查同通路更多代谢物。
- 结合通路关键酶基因或蛋白表达。
- 做靶向代谢验证。

靶点候选：

- 验证相关基因/蛋白表达。
- 检查 knockdown、inhibitor 或公开扰动数据。
- 对照癌种或细胞类型 context。

疾病或药物候选：

- 仅作为文献检索和候选排序。
- 需要独立队列、实验或公开数据验证。

## 七、后端组装模块

建议新增一个解释组装层，而不是直接改底层图谱。

模块名称建议：

```text
scripts/build_output_report.py
```

或在服务层增加：

```text
MetaboService.build_interpretation_report()
```

输入：

```text
analysis_pack
evidence
subgraph
release metadata
user question/context
```

输出：

```json
{
  "report_version": "interpretation_report.v1",
  "executive_summary": [],
  "conclusion_chains": [],
  "node_cards": [],
  "confidence_explanations": [],
  "review_queue": [],
  "evidence_sections": [],
  "next_validation": [],
  "appendix": {
    "raw_rankings": {},
    "raw_subgraph": {},
    "raw_evidence": {}
  },
  "determinism": {}
}
```

## 八、LLM 使用边界

LLM 只做语言组织，不做事实生成。

允许：

- 把结构化结论改写成自然语言。
- 把证据和置信度原因组织成清晰段落。
- 给出基于已有结论的验证建议。

禁止：

- 新增 graph 中不存在的关系。
- 替用户解析 ambiguous 输入为确定实体。
- 修改 ranking score 或 confidence tier。
- 引用不存在的文献或证据。
- 把探索性候选写成确定机制。

所有 LLM 输出必须通过现有 safe adapter 或等价 guard。

## 九、最小可行版本

第一阶段只做 4 个能力：

1. 核心结论卡
2. 结论链详情
3. 节点说明卡
4. 置信度原因解释

暂不重做全部图谱浏览器。

验收标准：

- 用户第一屏能看懂“主要发现是什么”。
- 每个结论都有输入依据、关键节点、证据和边界。
- 每个关键节点都能解释“为什么出现在本次分析中”。
- 每个置信度都有加分因素和降级因素。
- ambiguous/unmatched 输入不会被隐藏。

## 十、推荐落地顺序

1. 定义 `interpretation_report.v1` JSON contract。
2. 在后端从 `analysis_pack` 组装 conclusion chains。
3. 为 pathway、target、disease、drug、metabolite 生成 node cards。
4. 建立 confidence reason 规则表。
5. 在 Web 第一屏渲染核心结论和结论链。
6. 把原有节点详情、证据表和 JSON 移入可展开 appendix。
7. 加入验证 fixture，确保每条结论都有 traceable evidence refs。
8. 再接入 LLM 改写，但只允许基于 report JSON 叙述。

## 十一、数据库准确度重构优先级

本文件只定义输出端组织方式。若目标是提高系统本身准确度，应优先执行：

```text
docs/database_accuracy_reform.md
docs/database_schema_v2.md
```

输出端不得绕过 v2 事实层直接把名称匹配、ratio component、class/pool 或 pathway membership 包装成机制结论。`interpretation_report.v1` 后续应消费 `database_accuracy_store.v2` 生成的 `mechanism_ready_facts`、`evidence_assertions` 和 `trait_entities`。
