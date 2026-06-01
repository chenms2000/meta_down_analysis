# 自动高置信机制结论工具：根本整改流程

> 状态说明：本文保留为机制结论层的远期设计。当前主线整改已调整为数据库准确度优先，见 `docs/database_accuracy_reform.md` 和 `docs/database_schema_v2.md`。不要把本文理解为“提高门槛即可解决问题”；必须先完成实体、关系、证据和上下文事实层重构。

## 目标重新定义

当前系统不能继续定位为“把代谢物列表映射到图谱后生成解释报告”。如果目标是自动高置信机制结论工具，系统必须回答一个更严格的问题：

> 在给定癌种、组织、细胞类型、实验设计和代谢差异数据的条件下，哪些机制结论可以被稳定化学身份、方向一致的代谢模式、通路/反应结构、背景匹配证据和独立验证信号共同支持？

因此，核心产品不应是 pathway ranking、target ranking 或 disease ranking，而应是 `mechanism_claim`。

旧逻辑：

```text
输入代谢物 -> 实体匹配 -> 图谱传播/排名 -> 输出解释
```

新逻辑：

```text
实验问题 -> 输入资格审查 -> 化学身份锁定 -> 方向性代谢模式 -> 机制断言生成
-> 多源证据门控 -> 背景一致性验证 -> 高置信结论发布
```

没有通过门控的内容不进入“核心结论”，只能进入“未通过原因/待补证据”。

## 一、输入端必须重建

高置信机制结论不能接受“任意代谢物名称列表”作为充分输入。输入必须先分成四条轨道：

### 1. Direct metabolomics track

用于强机制结论的主轨。

必须字段：

- 稳定化学 ID：HMDB、ChEBI、PubChem CID、InChIKey 至少一种
- 代谢物名称
- 方向：up/down
- 效应量：log2FC、mean difference 或模型系数
- 显著性：pvalue/padj
- 实验背景：癌种、组织、细胞类型、分组、样本量、患者/样本 ID

可加分字段：

- RT
- MS2
- adduct
- platform
- batch
- patient_id
- replicate_id

### 2. Trait score track

例如 GCST trait score、代谢物 GWAS score、signature score。

这类输入不能直接生成“代谢物丰度变化机制”。它只能生成：

```text
trait score 层机制关联
```

必须保留：

- trait accession
- reported trait
- trait type：metabolite_level / metabolite_ratio / unknown_trait / class_trait
- score difference
- score effect label，例如 pseudo_log2FC_shifted、mean_diff、cohen_d

禁止：

- 把 ratio trait 拆成两个普通代谢物并作为同向丰度输入
- 把 pseudo_log2FC_shifted 显示为标准 log2FC
- 把 GCST trait 直接等同于 LC-MS 代谢物丰度

### 3. Ratio track

代谢物比值必须作为独立实体：

```text
ratio_trait(A / B)
```

高置信解释只能写成：

```text
A 相对 B 的比例发生变化
```

不能自动写成：

```text
A 上升
B 下降
```

除非原始输入同时提供 A 和 B 的独立丰度证据。

### 4. Review-only track

以下输入默认不参与高置信机制结论：

- ambiguous name
- unmatched
- score = 0 candidate
- lipid class name
- pool name
- X-unknown feature
- only name match without stable ID
- conflicting top candidates

这些输入只能进入复核队列。

## 二、建立机制断言对象

系统需要新增核心结构：

```json
{
  "mechanism_claim_id": "mc_001",
  "claim_level": "mechanism",
  "context": {
    "cancer_type": "cSCC",
    "tissue": "skin",
    "cell_type": "epithelial cell",
    "comparison": "Tumor vs Adjacent"
  },
  "claim_text": "cSCC tumor epithelial cells show increased glycolytic pyruvate/lactate handling.",
  "mechanism_type": "pathway_activity | enzyme_axis | transport_axis | redox_state | mitochondrial_state | nutrient_dependency",
  "input_pattern": [],
  "required_evidence": [],
  "supporting_evidence": [],
  "contradictory_evidence": [],
  "confidence_gate": {
    "status": "passed | failed",
    "tier": "high | not_high",
    "failed_reasons": []
  },
  "allowed_wording": "",
  "forbidden_wording": []
}
```

它不是排名结果，而是可审计的机制断言。

## 三、高置信机制结论的硬门槛

一条机制结论必须同时通过 7 个门槛。

### Gate 1：输入身份门槛

必须满足：

- 支持该结论的核心代谢物中，至少 70% 是 strict identity
- 至少 3 个独立 strict metabolites 支持同一机制轴
- ratio、pool、class、ambiguous 输入不能作为核心支持

失败则输出：

```text
未形成高置信机制结论：化学身份支持不足。
```

### Gate 2：方向一致性门槛

必须满足：

- 同一机制轴内代谢物方向符合预定义模式
- 不能只因为多个代谢物命中同一通路就判定机制增强

例：

糖酵解增强不能只看 glucose、pyruvate、lactate 出现；需要判断：

```text
lactate up
pyruvate up or downstream lactate/pyruvate ratio up
glycolytic intermediates direction consistent
TCA/mitochondrial branch not contradictory or explicitly解释
```

### Gate 3：统计稳健性门槛

必须满足至少一项：

- 患者层面 pseudobulk 显著
- mixed model 显著
- 样本/患者重复一致
- 外部队列复现

单细胞级 p 值极小不能单独支持高置信机制。

### Gate 4：通路结构门槛

机制必须落到 reaction / enzyme / transporter / pathway module，而不是只落到大通路名称。

可接受：

```text
glucose uptake -> glycolysis -> lactate production axis
glutamine uptake -> glutamate/alpha-ketoglutarate anaplerosis axis
cysteine/glutathione -> redox buffering axis
acylcarnitine accumulation -> mitochondrial FAO transport/pressure axis
```

不可接受：

```text
Metabolism
Amino acid metabolism
Cancer
melanoma
```

### Gate 5：背景一致性门槛

必须匹配：

- 癌种或泛鳞癌/上皮癌背景
- 组织或细胞类型
- 肿瘤 vs 邻近/正常的比较方向

不匹配结果只能进 appendix。

### Gate 6：独立证据门槛

高置信至少需要两类独立证据：

- curated reaction/pathway evidence
- 输入数据统计证据
- 文献中同癌种/相近背景直接支持
- 转录/蛋白/酶活/flux/DepMap 等正交证据

只有图传播或只有文献泛背景不够。

### Gate 7：反证检查门槛

系统必须主动检查：

- 同机制内是否存在方向相反的关键代谢物
- 是否主要由 ambiguous 输入驱动
- 是否由 ratio 拆分造成假方向
- 是否存在 context mismatch
- 是否有 conflict literature evidence

有强反证则不能发布高置信结论。

## 四、机制库要从“图谱边”升级为“机制模板库”

需要新增 `mechanism_templates`，每个模板定义可验证模式。

示例：

```json
{
  "template_id": "glycolysis_lactate_axis",
  "name": "Glycolysis-lactate axis",
  "mechanism_type": "pathway_activity",
  "required_context": ["tumor_or_proliferating_cell"],
  "core_metabolites": ["glucose", "glucose-6-phosphate", "pyruvate", "lactate"],
  "supporting_genes": ["SLC2A1", "HK2", "PFKP", "PKM", "LDHA"],
  "positive_pattern": {
    "lactate": "up",
    "pyruvate": "up_or_context_dependent",
    "glycolytic_intermediates": "up_or_coherent"
  },
  "contradiction_rules": [
    "only_ratio_without_independent_abundance",
    "only_one_metabolite_support",
    "no_context_match"
  ],
  "minimum_high_confidence_support": {
    "strict_metabolites": 3,
    "independent_evidence_types": 2,
    "patient_level_validation": true
  },
  "allowed_wording": "supports increased glycolysis/lactate handling",
  "forbidden_wording": ["proves Warburg effect", "causes cancer"]
}
```

第一批建议只做 6 个机制模板：

1. glycolysis/lactate axis
2. glutamine-glutamate anaplerosis axis
3. glutathione/cysteine redox buffering axis
4. mitochondrial TCA pressure
5. acylcarnitine/FAO transport pressure
6. nucleotide synthesis/proliferation axis

不要一开始覆盖所有通路。高置信工具要先窄后准。

## 五、分析流程重构

### Step 1：输入资格审查

输出：

```json
{
  "input_readiness": {
    "can_generate_high_confidence_mechanism": true,
    "blocking_reasons": [],
    "usable_track_counts": {
      "direct_metabolomics": 0,
      "trait_score": 0,
      "ratio": 0,
      "review_only": 0
    }
  }
}
```

如果输入主要是 GCST trait score，系统应直接声明：

```text
当前输入不足以生成代谢丰度层高置信机制结论；可生成 trait-score 机制关联，但需要原始代谢物丰度或患者层面验证补证。
```

这不是降低置信度，而是拒绝错误任务。

### Step 2：化学身份锁定

只允许 strict identity 进入高置信机制判断。

要求：

- ID 优先于名称
- InChIKey / PubChem CID / ChEBI / HMDB 一致
- lipid species 必须到具体分子种，class name 不进入精确机制
- unknown feature 不进入机制

### Step 3：方向性代谢模式构建

把输入转为机制模式表：

```json
{
  "metabolite_uid": "",
  "canonical_name": "",
  "direction": "up",
  "effect_value": 0.42,
  "effect_label": "log2FC",
  "padj": 0.01,
  "patient_consistency": 0.8,
  "mechanism_template_hits": []
}
```

### Step 4：机制模板匹配

每个模板计算：

- core metabolite coverage
- direction coherence
- statistical robustness
- context match
- contradiction count
- orthogonal evidence count

只输出通过模板门槛的机制断言。

### Step 5：证据包组装

每条机制断言必须带证据包：

```json
{
  "data_evidence": [],
  "identity_evidence": [],
  "reaction_evidence": [],
  "literature_evidence": [],
  "orthogonal_evidence": [],
  "contradiction_checks": []
}
```

### Step 6：高置信发布门

不是算一个分数，而是规则门控：

```text
high = Gate1 AND Gate2 AND Gate3 AND Gate4 AND Gate5 AND Gate6 AND NOT strong_contradiction
```

没过门的机制不进入核心结论。

### Step 7：输出端只渲染机制断言

第一屏改成：

1. 高置信机制结论
2. 每条结论的证据链
3. 未形成结论的原因
4. 需要补充的数据
5. 附录：探索性图谱排名

## 六、对当前项目的整改模块

### 新增模块 1：`mechanism_templates`

建议文件：

```text
configs/mechanism_templates/*.json
```

职责：

- 定义机制轴
- 定义核心代谢物
- 定义方向模式
- 定义反证规则
- 定义允许/禁止表述

### 新增模块 2：`MechanismClaimBuilder`

建议位置：

```text
scripts/metabo_service.py
```

或拆分：

```text
scripts/mechanism_claims.py
```

职责：

- 从 strict inputs 构建 mechanism evidence matrix
- 匹配 mechanism_templates
- 生成 `mechanism_claim`
- 输出 gate status

### 新增模块 3：`HighConfidenceGate`

职责：

- 执行 7 个硬门槛
- 记录 failed reasons
- 禁止未过门结果进入核心结论

### 新增模块 4：`OrthogonalEvidenceAdapter`

职责：

- 接入转录/蛋白/酶活/DepMap/flux/外部队列
- 没有正交证据时，高置信机制必须阻断或降级为“数据层强模式，机制待证”

### 新增模块 5：`MechanismReportRenderer`

职责：

- 只渲染通过 gate 的 mechanism_claim
- 对未通过的结果显示“为什么不能形成高置信结论”
- pathway ranking、disease ranking、drug overlay 全部进入附录

## 七、验收标准

自动高置信机制工具的验收不能看界面好不好看，而要看错误率。

最低验收：

- 20 个 gold positive mechanism cases
- 20 个 gold negative / trap cases
- ratio trap：Tryptophan/pyruvate ratio 不得推出 tryptophan up 或 pyruvate down
- class trap：Sphingomyelin 不得匹配到无关单体并进入机制
- context trap：神经元/星形胶质通路不得成为 cSCC epithelial 核心结论
- p-value trap：单细胞极小 p 值但患者不一致时不得高置信
- graph trap：只有图传播路径不得高置信

通过标准：

- high-confidence precision >= 0.9
- trap false-positive rate <= 0.05
- 每条 high claim 都有完整 evidence pack
- 每条 high claim 都有 allowed wording 和 forbidden wording 检查

## 八、对你这次 cSCC 结果的正确处理

这组数据目前最多能形成“待补证据的机制候选”，不能直接形成自动高置信机制结论，原因不是置信度调低，而是输入类型不满足高置信机制工具的前提：

- 输入是 trait score diff，不是直接 LC-MS 代谢物丰度
- 大量 GCST 是 ratio trait
- 部分实体匹配错误
- 缺少患者层面验证
- 缺少正交基因/蛋白/酶活/flux 证据

如果要让它变成高置信结论，需要补齐：

1. 每个 GCST 的 level/ratio 类型和稳定化学身份
2. 原始代谢物或 trait score 的患者层面 pseudobulk/mixed model
3. cSCC epithelial 背景下关键酶/转运体表达或蛋白证据
4. 每个机制模板的方向一致性检查
5. 外部队列或靶向 LC-MS/MS 复核

补齐后系统才能输出类似：

```text
高置信机制结论：
cSCC 肿瘤上皮细胞存在 glutathione/cysteine redox buffering 增强证据。

依据：
strict identity 的 cysteine、cystine、glutathione 相关代谢物方向一致；
患者层面 pseudobulk 复现；
Reactome/Rhea 反应链支持；
cSCC epithelial 中 SLC7A11/GCLC/GSS 表达上升作为正交证据；
未发现强方向反证。

边界：
支持氧化还原缓冲相关机制增强，不等同于证明该机制是肿瘤发生原因。
```

## 九、落地顺序

第一阶段：机制工具内核

1. 新增 `mechanism_claim.v1` contract
2. 新增 6 个机制模板
3. 构建 strict input evidence matrix
4. 实现 7 gate high-confidence 判定
5. 阻断 pathway ranking 直接进入核心结论

第二阶段：输入补证

1. 支持 patient_id / sample_id
2. 支持 pseudobulk summary
3. 支持 transcript/protein evidence upload
4. 支持 LC-MS/MS identity evidence

第三阶段：验证集

1. 建 gold positive cases
2. 建 trap negative cases
3. 建 cSCC epithelial 专项 fixture
4. 发布 precision / false-positive report

第四阶段：输出端

1. 第一屏只展示 high mechanism claims
2. 没有 high claim 时显示阻断原因和补证清单
3. 探索性图谱排名放入附录
4. LLM 只能改写 mechanism_claim，不允许生成新机制
