# Ep_ISA_NEW 数据需求说明

Ep_ISA_NEW 的 `report()` 和 `train()` 调用的 exploring/plotting/preprocess 模块依赖外部参考数据文件。这些文件**不随包分发**——请从 Google Drive 下载后放入 `data/` 目录，或在代码中通过 `get_data_resource()` 的搜索路径放置。

`get_data_resource(filename)` 搜索顺序：
1. 从 `utils.py` 向上遍历目录树找 `data/<filename>`
2. 包内 `data/<filename>`
3. 当前工作目录 `./data/<filename>`

---

## 必需文件清单（果蝇 dm6 适配）

### 1. 训练 & 噪声估计（`train()` 路径）

| 文件名 | 用途 | 果蝇数据来源建议 |
|--------|------|------------------|
| `non_cCRE_non_blacklist_non_exon.bed` | 背景区域（噪声估计 + 负样本采样） | 从 dm6 ENCODE blacklist + non-cCRE 构建排除外显子后的 BED |

**列格式**：标准 BED3（chrom, start, end），无 header。dm6 染色体名（chr2L, chr2R, chr3L, chr3R, chrX, chr4）。

### 2. TF 表达推断（`scoring/infer_tf_expr.py`）

| 文件名 | 用途 | 果蝇数据来源建议 |
|--------|------|------------------|
| `hg38_TF_promoters_500bp.bed` → **请重命名为 `dm6_TF_promoters_500bp.bed`** | TF 启动子坐标 | 从 FlyBase gene model 提取 TF gene TSS ±250bp |

**列格式**：BED4（chrom, start, end, name）。`name` = TF gene symbol（大写）。注意代码中硬编码了文件名 `hg38_TF_promoters_500bp.bed`，你需要修改 `infer_tf_expr.py` 中的文件名，或直接用同名文件放果蝇数据。

### 3. TF 家族分析（`exploring/tf_family.py`）

| 文件名 | 用途 | 果蝇数据来源建议 |
|--------|------|------------------|
| `HTFs_with_JASPAR_Families.csv` | TF → DNA-binding domain (DBD) / 家族映射 | 从 FlyFactorSurvey 或 DNASE2 TF 分类整理 |

**列格式**：CSV，至少含 `gene_symbol` 列（大写 TF 名）、`DBD` 列（DNA-binding domain 类别名）、`JASPAR_Family` 列（家族名）。可额外含 `JASPAR_Class`。

### 4. TF 功能分类（`exploring/tf_function.py`）

| 文件名 | 用途 | 果蝇数据来源建议 |
|--------|------|------------------|
| `universal_stripe_factors.txt` | 泛表达 TF（USF）列表 | 从 FlyAtlas 跨组织表达一致性分析提取 |
| `pioneer_factors.txt` | 先锋因子列表（如 Zelda） | 文献整理：Zld, Dl, etc. |
| `context_only_tfs.txt` | 条件依赖型 TF | FlyAtlas 组织特异性高的 TF |
| `gtex.dispersionEstimates.tab` → **替换为果蝇版本** | 组织表达离散度（Gini 系数） | 从 FlyAtlas 提取各 TF 跨组织表达，计算 Gini |

**txt 文件格式**：每行一个 TF gene symbol（大写），无 header。
**tab 文件格式**：TSV，至少含 `gene_symbol` 列和 `gini` 列。

### 5. PPI 验证（`exploring/tf_pair_ppi.py`）

| 文件名 | 用途 | 果蝇数据来源建议 |
|--------|------|------------------|
| `TF_TF_I.txt` | 已知 TF-TF 物理交互 | DroID 或 BioGRID 果蝇 PPI |
| `TF_Cof_I.txt` | TF-辅因子交互矩阵 | 文献 / DroID 过滤 |
| `TF_binding_coop_cleaned.csv` | DNA 介导的 TF 共结合 | ChIP-seq 共结合分析或文献 |

**TF_TF_I.txt 格式**：TSV，含 `TF1`、`TF2` 列（大写 gene symbol）。
**TF_Cof_I.txt 格式**：TSV，首列 `TF`，后续列为各辅因子名，值为 0/1。
**TF_binding_coop_cleaned.csv 格式**：CSV，含 `prey`、`bait` 列（gene symbol），交互值列。

### 6. Motif GC 分析（`plotting/tf.py`）

| 文件名 | 用途 | 果蝇数据来源建议 |
|--------|------|------------------|
| `JASPAR2026_CORE_non-redundant_pfms_jaspar.txt` | JASPAR PFM（计算 motif GC%） | JASPAR 下载 insects 子集或直接用全量（果蝇 motif 在其中） |

**格式**：标准 JASPAR PFM 格式（`>ID\tName` + 4 行 ACGT counts）。

---

## 快速设置（Colab）

```python
# 方法 1：从 Google Drive 挂载
from google.colab import drive
drive.mount('/content/drive')
!ln -s /content/drive/MyDrive/Ep_ISA_NEW_data Ep_ISA_NEW/data

# 方法 2：直接创建 data/ 目录并复制
import os
os.makedirs('Ep_ISA_NEW/data', exist_ok=True)
# 然后把准备好的参考文件放入 Ep_ISA_NEW/data/
```

---

## 果蝇数据准备优先级

如果暂时无法准备全部文件，按以下优先级：

1. **必须**：`non_cCRE_non_blacklist_non_exon.bed`（否则 `train()` 崩溃）
2. **高优**：`HTFs_with_JASPAR_Families.csv`（否则 `report()` 中 family 分析崩溃）
3. **高优**：`JASPAR2026_CORE_non-redundant_pfms_jaspar.txt`（否则 GC 分析崩溃）
4. **中优**：PPI 三件套（否则 PPI 验证图缺失，但不影响核心 ISA 结果）
5. **低优**：`universal_stripe_factors.txt` / `pioneer_factors.txt` / `context_only_tfs.txt` / `gtex.dispersionEstimates.tab`（功能分类图缺失）

> **注意**：即使参考数据缺失，`load_finemo()` → `run_isa()` → `calc_coop_score()` 的核心 ISA pipeline **不受影响**。只有 `report()` 的部分图会因缺少参考数据而跳过或报错。
