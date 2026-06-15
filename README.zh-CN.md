# 土壤微生物与污染梯度可视化 Demo

语言： [English](README.md) | 中文

这个仓库是一个适合公开展示的 R portfolio demo，用 simulated / desensitized toy data 展示土壤微生物、污染梯度、环境因子、功能注释、多样性、网络、富集、PLS-PM、LEfSe-style 和 DESeq2-style 等可视化 workflow。

这个项目不是论文代码仓库，不复现论文图，也不发布真实研究数据、私有分析表或未公开科学结论。所有公开 demo 输出都来自 toy inputs，应理解为 workflow organization、visualization logic 和 reproducible figure generation 的展示，而不是实际环境微生物学发现。

## 项目重点

- 用统一的 shared toy data 支撑 11 个独立 numbered modules；
- 从 raw-like toy input tables 重新计算中间表、统计结果和图形；
- 展示微生物群落、污染物梯度、土壤理化因子、功能注释和机制模型相关图形；
- 明确区分 public-safe demo outputs 和真实研究结果；
- 提供项目级完整性检查、output manifest 和 lightweight smoke test；
- 保守处理 private / raw / original / large sequencing files，避免误提交敏感材料。

## 公开边界

本仓库中的数据和图形仅用于公开演示：

- 不包含真实样本测量值；
- 不包含 raw sequencing data；
- 不包含私有原始分析表；
- 不包含外部数据库再分发文件；
- 不包含论文敏感结果；
- 不复现已发表论文图；
- 不应被解释为实际 biological conclusion 或 ecological evidence。

KEGG / FAPROTAX / BacMet / sulfur-cycle 等相关示例均使用 simulated toy annotations 或脚本内生成的 toy mapping，只用于展示分析和作图结构。

## 仓库结构

```text
.
|-- README.md
|-- README.zh-CN.md
|-- docs/
|   |-- data_privacy.md
|   |-- demo_selection.md
|   |-- module_registry.csv
|   |-- output_manifest.csv
|   |-- project_overview.md
|   |-- public_release_checklist.md
|   |-- r_package_requirements.md
|   |-- reproducible_r_environment.md
|   `-- shared_toy_data_schema.md
|-- scripts/
|   |-- check_project_integrity.R
|   |-- create_shared_toy_data.R
|   |-- install_r_dependencies.R
|   |-- run_all_demos.R
|   |-- run_smoke_test.R
|   `-- write_output_manifest.R
|-- data/
|   `-- toy_shared/
`-- 01_... to 11_.../
    |-- README.md
    |-- docs/
    |-- data/toy/
    |-- scripts/run_demo.R
    |-- results/
    `-- figures/
```

每个 numbered module 都可以从自己的文件夹独立运行。模块默认读取仓库根目录下的 `data/toy_shared/`，然后把结果写入自己的 `results/` 和 `figures/`。

## 共享 Toy Data

从仓库根目录运行：

```bash
Rscript scripts/create_shared_toy_data.R
```

该脚本会生成：

- `data/toy_shared/sample_metadata.csv`
- `data/toy_shared/environmental_variables.csv`
- `data/toy_shared/taxonomy_table.csv`
- `data/toy_shared/abundance_table.csv`
- `data/toy_shared/functional_annotation_table.csv`

这些表是 simulated / desensitized toy data，用于支撑所有 demo modules。

## 模块概览

| Module | 主要展示内容 |
|---|---|
| `01_rf_correlation_heatmap` | Random forest feature screening 与环境相关性热图 |
| `02_microbe_env_network` | 微生物-环境因子 association network |
| `03_ternary_taxa_distribution` | 三元图展示 taxa 在不同污染背景组中的分布 |
| `04_faprotax_functional_profile` | FAPROTAX-style functional profile 与 group summary |
| `05_lefse_biomarker` | LEfSe-style biomarker visualization |
| `06_plspm_mechanism_model` | PLS-PM-style mechanism-oriented path model |
| `07_differential_volcano_heatmap` | DESeq2-style differential abundance volcano plots 与 heatmap |
| `08_kegg_enrichment` | KEGG-style enrichment bubble / bar plots |
| `09_sulfur_gene_contaminant_association` | sulfur-cycle gene 与污染因子 association visualization |
| `10_alpha_beta_diversity` | Alpha / beta diversity、PCoA、NMDS 和 distance plot |
| `11_vpa_mantel_partitioning` | VPA 与 Mantel-style environmental partitioning |

机器可读模块信息见：

- [docs/module_registry.csv](docs/module_registry.csv)
- [docs/output_manifest.csv](docs/output_manifest.csv)

## Demo 图形

**Random-forest correlation workflow**

![Random-forest correlation workflow](01_rf_correlation_heatmap/figures/rf_correlation_combined.png)

**DESeq2-style differential workflow**

![DESeq2-style differential workflow](07_differential_volcano_heatmap/figures/deseq2_volcano_combined.png)

**PLS-PM-style mechanism diagram**

![PLS-PM-style mechanism diagram](06_plspm_mechanism_model/figures/plspm_inner_path_model.png)

**KEGG-style enrichment workflow**

![KEGG-style enrichment workflow](08_kegg_enrichment/figures/kegg_pathway_bubble_plot.png)

**Beta-diversity ordination workflow**

![Beta-diversity ordination workflow](10_alpha_beta_diversity/figures/fig2_pcoa_bray_curtis.png)

这些图形均来自 toy data，只用于展示可复现作图流程。

## 快速运行

完整运行路径：

```bash
Rscript scripts/install_r_dependencies.R
Rscript scripts/create_shared_toy_data.R
Rscript scripts/run_all_demos.R
```

单独运行一个模块：

```bash
cd 01_rf_correlation_heatmap
Rscript scripts/run_demo.R
```

返回根目录后可以继续运行其他模块：

```bash
cd ..
cd 02_microbe_env_network
Rscript scripts/run_demo.R
```

## Lightweight Smoke Test

如果只想快速确认项目结构、toy data 和轻量模块是否可运行，可以从根目录运行：

```bash
Rscript scripts/run_smoke_test.R
```

该 smoke test 会：

- 重新生成 shared toy data；
- 运行不依赖 heavy Bioconductor 的 `01` 到 `04` 模块；
- 检查关键 CSV / PDF / PNG 是否生成且非空；
- 运行 project integrity checker；
- 刷新 `docs/output_manifest.csv`。

这个 smoke test 不替代 full run，只是为了让 CI 和本地检查更快发现结构性问题。

## 项目完整性检查

从仓库根目录运行：

```bash
Rscript scripts/check_project_integrity.R
Rscript scripts/write_output_manifest.R
```

完整性检查会确认：

- shared toy-data files 是否存在；
- `sample_id` 是否跨 metadata、environmental variables、abundance columns、functional annotation rows 一致；
- `feature_id` 是否跨 abundance、taxonomy、functional annotation tables 一致；
- 11 个 numbered modules 是否都有 README 和 `scripts/run_demo.R`；
- R scripts 是否包含本地绝对路径或 raw/private/original-data 依赖；
- declared outputs 是否存在。

## R 依赖

依赖说明见：

- [docs/r_package_requirements.md](docs/r_package_requirements.md)
- [docs/reproducible_r_environment.md](docs/reproducible_r_environment.md)

大多数模块使用 CRAN packages。`07_differential_volcano_heatmap` 和 `08_kegg_enrichment` 使用 Bioconductor packages，因此 full run 比 smoke test 更重。

## GitHub Actions

当前包含两个 workflow：

- `.github/workflows/smoke-test.yml`：轻量 smoke test；
- `.github/workflows/run-demos.yml`：完整安装依赖并运行 11 个 modules。

smoke workflow 用于快速结构检查，full workflow 用于完整复现检查。

## License And Reuse

代码、simulated toy data、generated demo results 和 generated demo figures 使用 MIT License，见 [LICENSE](LICENSE)。

仓库中的 demo outputs 仅用于流程展示，不应被解释为真实研究结果，也不应用于科学推断。

## Limitations

- toy data 只是为了让 workflow 可以运行和展示，不是真实测量值；
- 图形不代表实际污染场地、微生物类群、功能基因或环境机制；
- 这些模块展示的是 portfolio-oriented figure workflow，不替代正式研究中的实验设计、测序质控、统计验证和领域审查；
- 如果替换为真实数据，需要自行确认数据许可、隐私边界、统计假设和领域解释。
