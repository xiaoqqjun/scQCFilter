# scQCFilter v0.1.2

## 简介

**scQCFilter** 是一个功能强大的 R 包，用于单细胞 RNA 测序 (scRNA-seq) 数据的质量控制 (QC) 分析。提供一键式自动化 QC 分析、智能参数推荐和交互式 HTML 报告。

## 核心特性

✨ **完全自动化** - 自动检测和计算线粒体基因比例
✨ **5 种物种支持** - 自动识别人、小鼠、大鼠、猕猴、斑马鱼
✨ **三维分析框架** - 同时分析参数、分组、样本三个维度
✨ **智能参数推荐** - 自动提供质量诊断和改进建议
✨ **交互式报告** - 自动生成美观的 HTML 报告
✨ **多参数对比** - 支持多个参数组合的同时分析和比较

## 安装

```r
# 从源代码安装
install.packages("scQCFilter.tar.gz", repos = NULL, type = "source")

# install.packages("devtools")
devtools::install_github("xiaoqqjun/scQCFilter")

# 加载包
library(scQCFilter)
```

## 快速开始

### 最简洁的用法（推荐）

```r
library(scQCFilter)

# 一行代码完成 QC 分析！
result <- scQCFilter(scRNA_obj, organism = "human")

# 查看交互式报告
browseURL(result$report_file)
```

### 完整的用法示例

```r
# 准备你的 Seurat 对象
scRNA_obj <- readRDS("your_data.rds")

# 运行 QC 分析（所有参数可选）
result <- scQCFilter(
  scRNA_obj,
  
  # 质控参数
  percent_mt_max = 20,        # 线粒体百分比上限（%）
  nCount_min = 500,           # 最小 UMI 数
  nCount_max = 25000,         # 最大 UMI 数
  nFeature_min = 200,         # 最小基因数
  nFeature_max = 6000,        # 最大基因数
  
  # 分析选项
  organism = "mmu",           # 物种：human, mouse, rat, macaque, zebrafish
  groups.by = "cell_type",    # 按细胞类型分组分析（可选）
  sample.by = "orig.ident",   # 按样本分析（可选）
  
  # 输出选项
  return.filtered = TRUE,     # 返回过滤后的对象
  report.dir = "./qc_report", # 报告保存目录
  plot = TRUE,                # 生成图表
  verbose = TRUE              # 打印详细信息
)

# 查看报告
browseURL(result$report_file)

# 获取过滤后的数据
filtered_obj <- result$filtered_object

# 继续下游分析
filtered_obj <- NormalizeData(filtered_obj)
filtered_obj <- FindVariableFeatures(filtered_obj)
filtered_obj <- ScaleData(filtered_obj)
```

## 物种支持

| 物种 | 参数值 | 线粒体基因模式 |
|------|--------|-----------------|
| 人类 | `"human"`, `"homo"`, `"hs"` | `^MT-` |
| 小鼠 | `"mouse"`, `"mice"`, `"mmu"`, `"mm"` | `^mt-` |
| 大鼠 | `"rat"`, `"rno"` | `^mt-` |
| 猕猴 | `"macaque"`, `"mfa"` | `^MT-` |
| 斑马鱼 | `"zebrafish"`, `"dre"` | `^mt-` |

## 参数说明

### 质控参数（数值）

- **`percent_mt_max`** (默认: 20)  
  线粒体基因百分比上限。单位为 %。
  - 推荐值：10-30% 取决于细胞类型
  - 免疫细胞通常较低（5-15%）
  - 代谢活跃细胞可能较高（15-30%）

- **`nCount_min`** (默认: 500)  
  最小 UMI（Unique Molecular Identifier）数。
  - 过低的值可能包含低质量细胞
  - 推荐值：300-1000

- **`nCount_max`** (默认: Inf)  
  最大 UMI 数。
  - 用于排除可能的细胞双联体（doublets）
  - 推荐值：根据数据确定，通常 20000-50000

- **`nFeature_min`** (默认: 200)  
  最小基因数（检测到表达的基因数）。
  - 太低表示测序深度不足
  - 推荐值：200-500

- **`nFeature_max`** (默认: Inf)  
  最大基因数。
  - 通常不需要限制上限
  - 可选设置以排除异常细胞

### 分析参数（字符串）

- **`organism`** (默认: "human")  
  物种名称。自动识别线粒体基因模式。

- **`groups.by`** (默认: NULL)  
  按特定列分组分析。如 `"cell_type"` 会分别为每个细胞类型生成 QC 统计。

- **`sample.by`** (默认: "orig.ident")  
  按样本分析。通常是实验批次或样本 ID。

### 输出参数（逻辑值和路径）

- **`return.filtered`** (默认: FALSE)  
  是否返回过滤后的 Seurat 对象。

- **`report.dir`** (默认: "./qc_report")  
  生成 HTML 报告的保存目录。

- **`plot`** (默认: TRUE)  
  是否生成图表。

- **`verbose`** (默认: TRUE)  
  是否打印详细的处理信息。

## 结果对象

`scQCFilter()` 返回一个列表，包含以下元素：

```r
result$report_file              # HTML 报告的完整路径
result$filtered_object          # 过滤后的 Seurat 对象（如果 return.filtered=TRUE）
result$summary_table            # QC 摘要表格
result$detailed_statistics      # 详细的统计信息
result$diagnostics$suggestions  # 质量改进建议
result$diagnostics$quality_score # 总体质量分数
result$plots                    # 生成的所有图表列表
```

## 高级用法

### 多参数对比分析

```r
# 同时尝试多个参数组合，找到最优的 QC 阈值
result <- scQCFilter(
  scRNA_obj,
  percent_mt_max = c(10, 15, 20, 30),    # 多个值
  nCount_min = c(300, 500, 1000),        # 多个值
  nFeature_min = 200,
  organism = "human",
  groups.by = "cell_type",
  sample.by = "batch"
)

# 查看比较结果
head(result$summary_table)
```

### 按细胞类型的分层分析

```r
# 为不同的细胞类型使用不同的 QC 标准
result <- scQCFilter(
  scRNA_obj,
  percent_mt_max = 20,
  nCount_min = 500,
  nFeature_min = 200,
  groups.by = "cell_type",  # 按细胞类型分别分析
  organism = "human"
)

# 每个细胞类型会有独立的 QC 统计和建议
```

### 手动线粒体计算

```r
# 如果需要手动计算（通常不需要）
seurat_obj <- addMitoRatio(seurat_obj, organism = "human")

# 或检查线粒体统计
checkMitoStatus(seurat_obj)
```

## 输出示例

### 控制台输出

```
========== scQCFilter: Starting QC Analysis ==========
Step 0: Checking mitochondrial percentage...
  -> percent.mt not found in metadata.
  -> Running addMitoRatio() with organism...
  -> Successfully added percent.mt
Step 1: Parsing parameters...
Step 2: Extracting metadata...
Step 3-5: Performing QC Analysis...
Step 6: Generating HTML Report...

========== QC Analysis Complete ==========
Report saved to: ./qc_report/scQCFilter_report_2024-01-16_123456.html
```

### HTML 报告内容

- **QC 统计汇总** - 各参数的分布统计
- **交互式图表** - Plotly 交互式可视化
- **细胞过滤结果** - 保留/过滤的细胞数
- **分组分析** - 按细胞类型的详细统计
- **质量诊断** - 自动生成的改进建议
- **质量评分** - 0-100 的总体质量评分

## 常见问题

### Q: 我的数据没有线粒体基因？
**A:** scQCFilter 会自动处理。如果未找到线粒体基因，将返回 0。这通常表示：
- 基因名称格式不符合预期
- 使用了不同的物种注释
- 可以手动指定 `percent_mito` 列

### Q: 应该使用什么 QC 阈值？
**A:** 推荐从宽松的阈值开始，然后根据生物学背景逐步调整：
- 一般起点：`percent_mt_max=20, nCount_min=500, nFeature_min=200`
- 严格标准：`percent_mt_max=10, nCount_min=1000, nFeature_min=500`
- 检查报告中的质量诊断建议

### Q: 支持自定义的 Seurat 对象吗？
**A:** 是的！只要对象包含 `RNA` assay 和 `meta.data` 就可以。scQCFilter 会自动处理不同的数据格式。

### Q: 如何在论文中引用？
**A:** 详见 `CITATION.txt` 文件。

## 许可证

MIT License

## 联系方式

**作者:** Zhijun Feng (冯志军)  
**邮箱:** xiaoqqjun@sina.com  
**GitHub:** https://github.com/xiaoqqjun/scQCFilter  
**ORCID:** 0000-0003-1813-1669  
**WeChat:** 博士后的小酒馆  

## 致谢

感谢 Seurat 团队提供的优秀分析框架和 PercentageFeatureSet 函数。

---

**版本:** v0.1.2  
**发布日期:** 2024-01-16  
**状态:** 生产就绪 ✅
