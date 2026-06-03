# 文献语料重建依赖

文献语料是文献证据 overlay 的本地数据依赖，不随 GitHub 代码仓库分发。本仓库仅记录检索式、时间范围、筛选条件和重建边界，使用户可在自己的环境中重建兼容语料。

## 检索策略

- 检索数据库/语法：PubMed。
- 本地获取环境：项目使用的 PubMed/PMC 文献记录及 Amazon Web Services Open Data 环境。
- 发表时间范围：2015-01-01 至 2026-05-05。
- 初始检索规模：约 140,000 条记录。
- 检索后纳入条件：期刊影响因子 > 5。
- 本地纳入语料规模：约 50,000 条记录。

完整检索式记录在 [../config/literature_search_strategy.json](../config/literature_search_strategy.json)。

## 检索式

```text
(
  cancer[tiab] OR cancers[tiab] OR tumor[tiab] OR tumors[tiab] OR tumour[tiab] OR tumours[tiab]
  OR neoplasm*[tiab] OR carcinoma*[tiab] OR malignan*[tiab] OR oncolog*[tiab]
  OR "Neoplasms"[MeSH Terms]
)
AND
(
  metabolism[tiab] OR metabolic[tiab] OR metabolite*[tiab] OR metabolomic*[tiab] OR metabolome[tiab]
  OR "cancer metabolism"[tiab] OR "tumor metabolism"[tiab] OR "tumour metabolism"[tiab]
  OR "metabolic reprogramming"[tiab] OR "metabolic rewiring"[tiab] OR "metabolic plasticity"[tiab]
  OR "energy metabolism"[tiab] OR bioenergetic*[tiab]
  OR "Metabolism"[MeSH Terms] OR "Metabolomics"[MeSH Terms]
)
AND
("2015/01/01"[Date - Publication] : "2026/05/05"[Date - Publication])
```

## 本地文件

推荐本地目录结构：

```text
literature_corpus/<release_id>/pmids.txt
literature_corpus/<release_id>/manifest.json
literature_evidence/<release_id>/
```

`manifest.json` 应记录获取日期、来源平台、检索式哈希、发表时间范围、影响因子来源和年份、纳入/排除数量、checksum，以及任何本地后处理脚本。若影响因子来自授权数据库，影响因子表及其派生整表应保留在 GitHub 仓库之外，除非相应许可证明确允许再分发。

## 分发边界

GitHub 仓库只包含检索式、筛选条件、schema 和重建说明；不包含下载得到的 PubMed/PMC 记录、摘要、全文、授权影响因子表或本地抽取的句级证据。需要文献证据 overlay 的用户，应在本地依据 PubMed、PMC/AWS Open Data、相关出版方以及影响因子数据提供方的条款自行重建。

文献层仅作为研究解释的证据 overlay，不会自动创建规范图谱事实。新颖、弱支持、不确定或冲突证据在独立复核前应保持探索性标记。

## 官方入口

- PubMed baseline/updatefiles 下载说明：https://pubmed.ncbi.nlm.nih.gov/download/
- PMC article datasets on AWS 说明：https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/
