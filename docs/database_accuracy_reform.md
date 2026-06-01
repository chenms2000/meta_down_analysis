# 数据库准确度重构方案

## 目标

本项目要提高的不是置信度显示，而是系统事实层的准确度。核心问题应从：

```text
怎样让输出更谨慎？
```

改成：

```text
怎样让数据库中的实体、关系、证据和上下文本身更准确？
```

如果数据库把 trait、ratio、lipid class、unknown feature、精确代谢物、通路共现、文献泛关联都压成同一种节点和同一种边，再好的输出端也只能包装不确定性。根本整改必须先重构数据层。

## 当前结构的主要准确度问题

### 1. 化学实体粒度混杂

当前 `metabolites` 和匹配索引更像“可被查询的化学名集合”，但高准确度分析需要区分：

- 精确化学实体：有 InChIKey / PubChem CID / ChEBI / HMDB 的 compound
- 脂质分子种：如 PC(16:0/18:1)
- 脂质类别：如 sphingomyelin
- 代谢物池：如 glutathione redox pair
- 平台未知峰：如 X-21383
- GWAS trait：如 GCST90200911
- ratio trait：如 tryptophan/pyruvate ratio

这些不能都当成 metabolite seed。

### 2. 匹配结果缺少“身份决策层”

当前 resolver 给出 matched / ambiguous / unmatched，但没有完整保存：

- 为什么匹配
- 哪些候选被排除
- 候选之间的冲突类型
- 该匹配属于 exact identity、class-level、ratio component、analog 还是 review-only
- 是否经过人工或规则审定

所以错误匹配会直接流入图谱。

### 3. 通路关系过粗

`metabolite_pathway_edges` 和 pathway ranking 容易把“通路共现”误读成机制关系。

例如：

```text
metabolite -> pathway
```

不等于：

```text
该通路活性增强
该代谢物流量增强
该反应被驱动
```

数据库需要 reaction / enzyme / transporter / compartment / participant role 层，而不是只靠 pathway membership。

### 4. 文献证据缺少断言结构

文献句子应该支持一个具体断言，而不是简单支持一条泛边。

错误做法：

```text
metabolite associated_with disease
```

正确做法：

```text
在某癌种/组织/实验条件下，某代谢物上升/下降/调控某过程
```

需要把文献证据拆成：

- 主体
- 谓词
- 客体
- 方向
- 上下文
- 方法
- 支持/反证/不确定/背景
- 证据句

### 5. context 不是一等公民

癌种、组织、细胞类型、实验模型、样本来源目前更像筛选条件，而不是事实表中的上下文坐标。结果是 melanoma、astrocytic、taste perception 这类背景不匹配内容容易进入 cSCC epithelial 结果。

## 新数据库分层

建议把数据库从现在的三层：

```text
normalized_store
compound_match_index
graph_projection
```

升级为六层：

```text
source_registry
entity_store
identity_resolution_store
relation_store
evidence_store
analysis_view / graph_projection
```

`graph_projection` 不再是事实源，只是为了检索、传播和界面浏览而生成的投影。

## 一、source_registry：来源和版本层

每条事实必须能回到来源。

### 表：`source_records`

字段：

```text
source_record_uid
source_name
source_type              curated_db | literature | user_upload | computed_overlay
source_version
download_date
license_id
file_path
checksum
parser_name
parser_hash
raw_record_id
raw_payload_json
```

用途：

- 防止不同版本数据混用
- 支持追踪错误来源
- 支持按来源重建和回滚

## 二、entity_store：实体层

实体必须拆成不同类型，不能统一塞进 metabolite。

### 表：`chemical_entities`

只存精确或近精确化学实体。

字段：

```text
chemical_uid
entity_granularity       exact_compound | stereoisomer | lipid_species | salt_or_adduct | unknown
canonical_name
formula
monoisotopic_mass
inchi
inchikey
inchikey14
smiles
charge
source_priority
identity_status          canonical | merged | deprecated | review
source_record_uid
```

### 表：`chemical_xrefs`

字段：

```text
chemical_uid
xref_source              HMDB | ChEBI | PubChem | KEGG | LipidMaps | CAS
xref_id
xref_type                primary | secondary | synonym | obsolete
source_record_uid
```

### 表：`chemical_names`

名称需要带来源和风险等级。

```text
chemical_uid
name
name_type                canonical | synonym | abbreviation | class_name | ambiguous_alias
language
source_record_uid
name_risk                low | medium | high
```

例如 `sphingomyelin` 应是 `class_name/high risk`，不能支持精确身份。

### 表：`metabolite_classes`

用于 lipid class、pool、family。

```text
class_uid
class_name
class_type               lipid_class | metabolite_family | redox_pair | pathway_pool
parent_class_uid
definition
source_record_uid
```

### 表：`chemical_class_members`

```text
class_uid
chemical_uid
membership_type          exact_member | possible_member | inferred_member
source_record_uid
```

这能允许 `sphingomyelin` 支持“鞘脂类别主题”，但不能错误匹配到某个无关单体。

### 表：`trait_entities`

GCST、signature、score、ratio 都进 trait，不进 metabolite。

```text
trait_uid
accession_id
trait_name
trait_type               metabolite_level | metabolite_ratio | class_trait | xenobiotic_trait | unknown_trait
trait_source             GWAS | signature_score | user_score
reported_trait
summary_statistics_url
source_record_uid
identity_scope           exact_chemical | ratio | class_level | unknown
```

### 表：`trait_components`

ratio 或复合 trait 的组件。

```text
trait_uid
component_role           numerator | denominator | component | class_component
component_uid
component_entity_type    chemical | class | unknown
direction_semantics      same_direction | inverse_direction | undefined
curation_status          reviewed | rule_inferred | review_required
source_record_uid
```

关键规则：

- ratio trait 可以连接到组件，但不能自动生成组件丰度方向。
- 组件连接必须有 `direction_semantics`。

## 三、identity_resolution_store：身份解析层

resolver 不能只返回一个结果，必须保存完整决策。

### 表：`input_features`

用户输入先成为观测特征，而不是直接成为 metabolite。

```text
input_feature_uid
input_row_id
feature_label
feature_type             direct_metabolite | trait_score | ratio | lipid_class | unknown_peak
effect_value
effect_label             log2FC | mean_diff | pseudo_log2FC_shifted | cohen_d
direction
pvalue
padj
sample_context_uid
raw_input_json
```

### 表：`identity_candidates`

```text
input_feature_uid
candidate_uid
candidate_entity_type    chemical | class | trait
candidate_name
match_basis              exact_id | inchikey | exact_name | synonym | mass | rt_ms2 | class_name | trait_annotation
score
margin
rank
false_match_risk         low | medium | high
blocking_reason
source_record_uid
```

### 表：`identity_decisions`

```text
input_feature_uid
decision_status          accepted_exact | accepted_class | accepted_trait | ambiguous | unmatched | rejected
accepted_entity_uid
accepted_entity_type
decision_rule
decision_confidence
review_status            auto | needs_review | human_accepted | human_rejected
reviewer
reviewed_at
decision_notes
```

准确度提升点：

- Tryptophan -> pyruvate 这类错误会停在 candidate/decision 层。
- score=0 candidate 只能是 rejected。
- ratio/class/unknown 不会污染 exact metabolite 图谱。

## 四、relation_store：关系事实层

图谱边要拆成事实关系，不再只有一个泛化 edge。

### 表：`reaction_entities`

```text
reaction_uid
reaction_name
reaction_source          Reactome | Rhea | KEGG | WikiPathways
equation
compartment_uid
source_record_uid
```

### 表：`reaction_participants_v2`

```text
reaction_uid
chemical_uid
participant_role         substrate | product | cofactor | inhibitor | activator
stoichiometry
compartment_uid
directionality           left_to_right | right_to_left | reversible | unknown
source_record_uid
```

### 表：`reaction_catalysts`

```text
reaction_uid
gene_uid
protein_uid
catalyst_role            enzyme | transporter | complex_member
source_record_uid
```

### 表：`pathway_modules`

大通路要切成机制模块。

```text
module_uid
module_name
parent_pathway_uid
module_type              pathway_module | transport_axis | redox_axis | biosynthesis_axis
definition
source_record_uid
```

### 表：`module_members`

```text
module_uid
member_uid
member_type              reaction | chemical | gene | class
member_role              core | supporting | context | contradiction_marker
source_record_uid
```

准确度提升点：

- 系统可以判断“命中反应链”还是“只是通路共现”。
- 机制结论可以落在 module，而不是泛 pathway。

## 五、evidence_store：证据断言层

证据不是边的备注，而是一等表。

### 表：`evidence_assertions`

```text
assertion_uid
subject_uid
subject_type
predicate                increases | decreases | associated_with | participates_in | catalyzes | transports | biomarker_of
object_uid
object_type
polarity                 positive | negative | neutral
direction                up | down | mixed | not_applicable
evidence_class           direct_assay | omics_association | curated_pathway | literature_background | model_prediction
support_status           support | contradict | uncertain | background
context_uid
method                   LC-MS | RNA-seq | proteomics | flux | knockdown | curated_db | text_mined
source_record_uid
sentence_uid
confidence_components_json
```

### 表：`evidence_sentences`

```text
sentence_uid
pmid
pmcid
section
sentence_text
sentence_hash
source_record_uid
```

### 表：`evidence_contexts`

```text
context_uid
cancer_type_uid
tissue_uid
cell_type_uid
cell_state_uid
species
model_system             patient | cell_line | organoid | mouse | in_vitro
comparison
context_text
```

准确度提升点：

- 文献中的泛癌背景不会自动支持 cSCC epithelial。
- support 和 contradict 可以同时保存，并在分析时显式冲突。

## 六、analysis_view / graph_projection：分析投影层

`graph_projection` 应由上面事实表生成，而不是承担事实语义。

投影时要生成多张图：

### 1. identity graph

用于实体解析。

包含：

- chemical_entities
- chemical_xrefs
- chemical_names
- metabolite_classes

### 2. reaction graph

用于机制解释。

包含：

- chemical -> reaction -> chemical
- reaction -> enzyme/transporter
- reaction -> pathway_module

### 3. evidence graph

用于证据追踪。

包含：

- assertion -> sentence/source/context
- assertion -> supported relation

### 4. context graph

用于背景匹配。

包含：

- cancer type ontology
- tissue ontology
- cell type ontology
- model system

### 5. report graph

用于前端浏览。

这是压缩图，只能展示，不作为高准确度事实源。

## 七、从数据库层提升准确度的关键规则

### Rule 1：实体类型不跨层晋升

```text
trait 不自动变 metabolite
ratio 不自动变 metabolite
class 不自动变 exact compound
unknown feature 不自动变 named metabolite
```

### Rule 2：名称只能作为候选，不作为事实

只有以下情况可以 accepted_exact：

- stable ID 命中
- InChIKey 命中
- PubChem/ChEBI/HMDB 多源一致
- 名称命中且唯一、低风险、无冲突，并有额外证据支持

### Rule 3：pathway membership 不能代表机制

`chemical in pathway` 只能生成候选，不生成机制结论。

机制必须来自：

- reaction chain
- enzyme/transporter
- module pattern
- context-matched evidence assertion

### Rule 4：证据必须带上下文

没有 context 的文献证据只能作为 background，不支持强机制。

### Rule 5：用户输入是 observation，不是 canonical fact

用户输入表只产生 `input_features` 和 `identity_decisions`，不能写入 canonical entity。

### Rule 6：所有分析结果必须可回放

分析结果需要记录：

- release_id
- source_record_uid
- identity_decision_uid
- relation/assertion UID
- context_uid
- scoring version

## 八、重构后的分析流程

### Step 1：输入进入 `input_features`

保留原始字段和输入类型，不直接匹配成 metabolite。

### Step 2：生成 `identity_candidates`

从 ID、名称、结构、mass、RT/MS2、trait annotation 多通道产生候选。

### Step 3：写入 `identity_decisions`

用规则或人工审定决定 accepted_exact / accepted_class / accepted_trait / ambiguous。

### Step 4：构建 observation matrix

只把 accepted_exact chemical 进入精确代谢物矩阵。

### Step 5：映射到 reaction/module

优先使用 reaction_participants_v2 和 pathway_modules，不再直接用 pathway membership 生成机制。

### Step 6：聚合 context-matched evidence assertions

证据必须按 context 过滤和加权。

### Step 7：输出 mechanism-ready facts

输出给上层算法的是：

```json
{
  "exact_observations": [],
  "class_observations": [],
  "trait_observations": [],
  "reaction_hits": [],
  "module_hits": [],
  "context_matched_assertions": [],
  "contradictory_assertions": [],
  "blocked_inputs": []
}
```

而不是直接输出 pathway ranking。

## 九、当前项目的落地顺序

### 第一阶段：不破坏现有服务，新增事实层

新增文档和 schema：

```text
docs/database_accuracy_reform.md
docs/database_schema_v2.md
```

新增构建脚本：

```text
scripts/build_database_accuracy_store_v2.py
scripts/database_accuracy_v2.py
```

该入口一次性构建 `entity_store`、`relation_store`、`evidence_store` 和 `analysis_view`。后续如果需要更细的流水线，可再拆成：

```text
scripts/build_entity_store_v2.py
scripts/build_identity_resolution_store.py
scripts/build_relation_store_v2.py
scripts/build_evidence_store_v2.py
scripts/build_analysis_views_v2.py
```

现有 `normalized_store` 保留，但逐步迁移为 v2 来源。

### 第二阶段：实体解析先替换

优先改 resolver，因为它是最大假阳性来源。

输出从：

```text
matched / ambiguous / unmatched
```

升级为：

```text
input_feature -> identity_candidates -> identity_decision
```

并把 resolver 的每一步写入表。

### 第三阶段：reaction/module 替换 pathway-first 解释

新增 `pathway_modules`，先人工定义 20 个高价值模块：

- glycolysis/lactate
- glutamine/glutamate
- TCA
- glutathione/cysteine
- acylcarnitine/FAO
- nucleotide synthesis
- one-carbon/serine
- sphingolipid membrane
- bile acid
- urea/arginine

然后把通路解释从 pathway ranking 改为 module hit。

### 第四阶段：文献证据断言化

把 `literature_edge_support` 升级成 `evidence_assertions`。

尤其要区分：

- direct assay
- omics association
- background/review
- uncertainty
- contradiction
- context mismatch

### 第五阶段：前端只消费 analysis_view

前端不直接读旧 ranking，而是读：

```text
mechanism_ready_facts
blocked_inputs
evidence_assertions
context_match_summary
```

## 十、验收指标

数据库准确度验收不是看输出好不好听，而是看基础事实错误率。

### 实体解析指标

- exact identity precision >= 0.95
- high-risk name false match rate <= 0.02
- ratio/class/unknown false promotion rate = 0
- ambiguous recall >= 0.9

### 关系准确度指标

- reaction participant role accuracy >= 0.95
- pathway module membership precision >= 0.9
- context mismatch false inclusion <= 0.05

### 证据准确度指标

- direct evidence precision >= 0.9
- background evidence 不进入 direct support
- contradiction evidence 可被检出
- 每条证据可回到 sentence/source/context

### 系统级指标

- trap cases false positive <= 0.05
- 每条机制结果都能追踪到 identity_decision、reaction/module、evidence_assertion
- pathway membership-only 结果不能进入核心机制

## 十一、对当前 cSCC 数据的数据库层解释

这次 cSCC 结果暴露的不是输出问题，而是数据库结构问题：

- GCST trait 被当成 metabolite 使用
- ratio trait 没有独立实体
- class/unknown/name-risk 没有独立身份层
- pathway membership 被过早解释成机制
- 文献证据没有严格按 context 和 evidence class 分层

因此正确整改不是“调低置信度”，而是：

1. 给 GCST 建 `trait_entities`
2. 给 ratio 建 `trait_components`
3. 给每个输入建 `input_features`
4. 把匹配过程写入 `identity_candidates` 和 `identity_decisions`
5. 用 `reaction_participants_v2` 和 `pathway_modules` 替代 pathway-first 解释
6. 用 `evidence_assertions` 替代泛化 literature edge
7. 最后再让输出端渲染这些更准确的事实

这样系统准确度才会真实提高。
