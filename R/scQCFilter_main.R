# ============================================================================
# scQCFilter: Single-cell QC Filtering and Analysis Package
# ============================================================================
# 主函数 - 支持三维度分层分析 + 多参数对比 + 自动HTML报告生成


# ============================================================================
# 1. 主函数: scQCFilter
# ============================================================================

#' Single-cell Quality Control Filtering with Multi-dimensional Analysis
#'
#' @param seurat_obj Seurat object for filtering
#' @param percent_mt_max Maximum mitochondrial percentage (default: 20)
#' @param nCount_min Minimum UMI count (default: 500)
#' @param nCount_max Maximum UMI count (default: Inf)
#' @param nFeature_min Minimum gene count (default: 200)
#' @param nFeature_max Maximum gene count (default: Inf)
#' @param groups.by Column name for group analysis (default: NULL)
#' @param sample.by Column name for sample analysis (default: "orig.ident")
#' @param return.filtered Return filtered Seurat object (default: FALSE)
#' @param report.dir Directory to save HTML report (default: "./qc_report")
#' @param plot Generate plots (default: TRUE)
#' @param organism Species for automatic mitochondrial calculation ('human' or 'mouse'). Only used if percent_mito not in metadata.
#' @param verbose Print messages (default: TRUE)
#' @param show.qc.distribution Show pre-filtering QC distribution plots (default: TRUE)
#' @param qc.distribution.percentiles Percentile thresholds for distribution plots (default: c(2.5, 5, 95, 97.5))
#'
#' @return S3 object of class "scQCFilter" containing analysis results

scQCFilter <- function(
  seurat_obj,
  percent_mt_max = 20,
  nCount_min = 500,
  nCount_max = Inf,
  nFeature_min = 200,
  nFeature_max = Inf,
  groups.by = NULL,
  sample.by = "orig.ident",
  organism = "human",
  return.filtered = FALSE,
  report.dir = "./qc_report",
  plot = TRUE,
  verbose = TRUE,
  show.qc.distribution = TRUE,
  qc.distribution.percentiles = c(2.5, 5, 95, 97.5)
) {
  
  if (verbose) cat("========== scQCFilter: Starting QC Analysis ==========\n")
  
  # 步骤1: 解析参数（支持单值或向量）
  # =========================================================================
  # Auto-check and calculate mitochondrial percentage if not present
  # =========================================================================
  if (verbose) cat("Step 0: Checking mitochondrial percentage...\n")
  
  if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
    if (verbose) cat("  -> percent.mt not found in metadata.\n")
    if (verbose) cat("  -> Running addMitoRatio() with organism...\n")
    seurat_obj <- addMitoRatio(seurat_obj, organism = organism)
    if (verbose) cat("  -> Successfully added percent.mt\n")
  } else {
    if (verbose) cat("  -> percent.mt already exists in metadata. Skipping.\n")
  }
  
  if (verbose) cat("Step 1: Parsing parameters...\n")
  param_list <- .parse_params(
    percent_mt_max = percent_mt_max,
    nCount_min = nCount_min,
    nCount_max = nCount_max,
    nFeature_min = nFeature_min,
    nFeature_max = nFeature_max
  )
  
  n_params <- length(param_list)
  is_multi_params <- n_params > 1
  
  if (verbose) {
    if (is_multi_params) {
      cat(sprintf("  → %d parameter combinations detected\n", n_params))
    } else {
      cat("  → Single parameter set\n")
    }
  }
  
  # 步骤2: 提取元数据
  if (verbose) cat("Step 2: Extracting metadata...\n")
  metadata <- seurat_obj@meta.data
  n_cells_before <- nrow(metadata)
  
  if (verbose) {
    cat(sprintf("  → Total cells: %d\n", n_cells_before))
    if (!is.null(groups.by) && groups.by %in% colnames(metadata)) {
      cat(sprintf("  → Groups: %s\n", paste(unique(metadata[[groups.by]]), collapse=", ")))
    }
    if (!is.null(sample.by) && sample.by %in% colnames(metadata)) {
      cat(sprintf("  → Samples: %s\n", paste(unique(metadata[[sample.by]]), collapse=", ")))
    }
  }
  
  # 步骤3: 对每个参数组合执行分析
  if (verbose) cat("Step 3: Executing filtering for each parameter set...\n")
  
  all_results <- list()
  
  for (i in seq_along(param_list)) {
    params <- param_list[[i]]
    param_name <- names(param_list)[i]
    
    if (verbose) cat(sprintf("  [%d/%d] Processing %s\n", i, n_params, param_name))
    
    # 执行过滤
    keep_cells <- .apply_filters(metadata, params)
    
    # 维度1: Parameter分层统计
    by_param <- .analyze_by_parameter(metadata, keep_cells, params)
    
    # 维度2: Group分析
    by_group <- NULL
    if (!is.null(groups.by) && groups.by %in% colnames(metadata)) {
      by_group <- .analyze_by_group(metadata, keep_cells, params, groups.by)
    }
    
    # 维度3: Sample分析
    by_sample <- NULL
    if (!is.null(sample.by) && sample.by %in% colnames(metadata)) {
      by_sample <- .analyze_by_sample(metadata, keep_cells, params, sample.by)
    }
    
    # 维度4: Sample × Group
    sample_group <- NULL
    if (!is.null(groups.by) && !is.null(sample.by) && 
        groups.by %in% colnames(metadata) && sample.by %in% colnames(metadata)) {
      sample_group <- .analyze_sample_group(metadata, keep_cells, params, groups.by, sample.by)
    }
    
    # 整体统计
    n_cells_after <- sum(keep_cells)
    n_filtered <- n_cells_before - n_cells_after
    keep_rate <- n_cells_after / n_cells_before
    
    all_results[[param_name]] <- list(
      params = params,
      param_name = param_name,
      keep_cells = keep_cells,
      n_cells_before = n_cells_before,
      n_cells_after = n_cells_after,
      n_filtered = n_filtered,
      keep_rate = keep_rate,
      by_parameter = by_param,
      by_group = by_group,
      by_sample = by_sample,
      sample_group_composition = sample_group
    )
  }
  
  # 步骤4: 多参数对比
  comparison <- NULL
  if (is_multi_params) {
    if (verbose) cat("Step 4: Comparing multiple parameter sets...\n")
    comparison <- .compare_all_params(all_results, param_list)
  }
  
  # 步骤5: 参数推荐
  if (verbose) cat("Step 5: Recommending best parameters...\n")
  recommend <- .recommend_best_params(all_results, comparison, is_multi_params)
  
  # 步骤6: 诊断分析
  if (verbose) cat("Step 6: Running diagnostics...\n")
  diagnostics <- .run_diagnostics(all_results, groups.by, sample.by)
  
  # 步骤7: 获取第一个参数集的信息作为主要报告
  main_result <- all_results[[1]]

  # 步骤8: 生成QC分布图和统计（新功能）
  qc_distribution <- NULL
  if (show.qc.distribution) {
    if (verbose) cat("Step 7: Generating QC distribution plots...\n")
    qc_distribution <- .create_qc_distribution(
      metadata = metadata,
      groups.by = groups.by,
      params = main_result$params,
      percentiles = qc.distribution.percentiles,
      report_dir = report.dir,
      verbose = verbose
    )
  }

  # 步骤9: 生成HTML报告
  if (verbose) cat("Step 8: Generating HTML report...\n")
  
  if (!base::dir.exists(report.dir)) {
    base::dir.create(report.dir, recursive = TRUE)
  }
  
  report_file <- .generate_html_report(
    all_results = all_results,
    comparison = comparison,
    recommend = recommend,
    diagnostics = diagnostics,
    report_dir = report.dir,
    groups.by = groups.by,
    sample.by = sample.by,
    is_multi_params = is_multi_params,
    qc_distribution = qc_distribution,
    show.qc.distribution = show.qc.distribution
  )

  if (verbose) cat(sprintf("  → Report saved to: %s\n", report_file))

  # 步骤10: 生成图表
  if (plot) {
    if (verbose) cat("Step 9: Generating plots...\n")
    plots <- .create_visualizations(
      all_results = all_results,
      comparison = comparison,
      diagnostics = diagnostics,
      groups.by = groups.by,
      sample.by = sample.by,
      report_dir = report.dir,
      metadata = metadata
    )
  } else {
    plots <- NULL
  }

  # 步骤11: 组织返回值
  if (verbose) cat("Step 10: Organizing results...\n")
  
  # 获取过滤后的对象（如果需要）
  filtered_object <- NULL
  if (return.filtered) {
    filtered_object <- seurat_obj[, main_result$keep_cells]
  }
  
  # 创建返回对象
  output <- structure(
    list(
      summary = base::data.frame(
        metric = c("cells_before", "cells_after", "cells_filtered",
                  "filter_rate", "keep_rate"),
        value = c(n_cells_before, main_result$n_cells_after,
                 main_result$n_filtered, main_result$n_filtered/n_cells_before,
                 main_result$keep_rate)
      ),
      main_result = main_result,
      all_results = all_results,
      comparison = comparison,
      recommend = recommend,
      diagnostics = diagnostics,
      report_file = report_file,
      plots = plots,
      filtered_object = filtered_object,
      metadata = metadata,
      is_multi_params = is_multi_params,
      groups.by = groups.by,
      sample.by = sample.by,
      qc_distribution = qc_distribution
    ),
    class = "scQCFilter"
  )
  
  if (verbose) cat("========== QC Analysis Complete ==========\n\n")
  
  return(output)
}

# ============================================================================
# 2. 参数解析函数
# ============================================================================

.parse_params <- function(percent_mt_max, nCount_min, nCount_max,
                         nFeature_min, nFeature_max) {
  
  # 创建所有参数的组合
  param_grid <- base::expand.grid(
    percent_mt_max = percent_mt_max,
    nCount_min = nCount_min,
    nCount_max = nCount_max,
    nFeature_min = nFeature_min,
    nFeature_max = nFeature_max
  ) %>%
    base::data.frame()
  
  # 转换为list
  param_list <- list()
  for (i in seq_len(nrow(param_grid))) {
    param_name <- sprintf("Set_%d", i)
    param_list[[param_name]] <- as.list(param_grid[i, ])
  }
  
  return(param_list)
}

# ============================================================================
# 3. 过滤逻辑函数
# ============================================================================

.apply_filters <- function(metadata, params) {
  
  keep_cells <- (
    metadata$percent.mt <= params$percent_mt_max &
    metadata$nCount_RNA >= params$nCount_min &
    metadata$nCount_RNA <= params$nCount_max &
    metadata$nFeature_RNA >= params$nFeature_min &
    metadata$nFeature_RNA <= params$nFeature_max
  )
  
  return(keep_cells)
}

# ============================================================================
# 4. Parameter分层统计
# ============================================================================

.analyze_by_parameter <- function(metadata, keep_cells, params) {
  
  # 计算每个参数单独导致的过滤
  filter_mt <- metadata$percent.mt > params$percent_mt_max
  filter_count_min <- metadata$nCount_RNA < params$nCount_min
  filter_count_max <- metadata$nCount_RNA > params$nCount_max
  filter_count <- filter_count_min | filter_count_max
  filter_feature_min <- metadata$nFeature_RNA < params$nFeature_min
  filter_feature_max <- metadata$nFeature_RNA > params$nFeature_max
  filter_feature <- filter_feature_min | filter_feature_max
  
  # 计算重叠
  only_mt <- filter_mt & !filter_count & !filter_feature
  only_count <- filter_count & !filter_mt & !filter_feature
  only_feature <- filter_feature & !filter_mt & !filter_count
  multiple <- (filter_mt + filter_count + filter_feature) >= 2
  
  result <- base::data.frame(
    parameter = c("percent.mt", "nCount_RNA", "nFeature_RNA", "combined"),
    n_filtered = c(sum(filter_mt), sum(filter_count), sum(filter_feature), sum(multiple)),
    pct_of_total = c(sum(only_mt), sum(only_count), sum(only_feature), sum(multiple)) / sum(!keep_cells),
    severity_rank = c(NA, NA, NA, NA)
  )
  
  # 计算severity rank
  filtered_counts <- c(sum(filter_mt), sum(filter_count), sum(filter_feature))
  result$severity_rank[1:3] <- rank(-filtered_counts, ties.method = "min")
  
  return(result)
}

# ============================================================================
# 5. Group分析
# ============================================================================

.analyze_by_group <- function(metadata, keep_cells, params, groups.by) {
  
  groups <- unique(metadata[[groups.by]])
  
  result <- base::data.frame(
    group = character(),
    cells_before = numeric(),
    cells_after = numeric(),
    keep_rate = numeric(),
    mean_nCount_after = numeric(),
    mean_nFeature_after = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (group in groups) {
    idx <- metadata[[groups.by]] == group
    cells_before <- sum(idx)
    cells_after <- sum(idx & keep_cells)
    keep_rate <- cells_after / cells_before
    
    mean_count <- ifelse(cells_after > 0, 
                        mean(metadata$nCount_RNA[idx & keep_cells]), 
                        NA)
    mean_feature <- ifelse(cells_after > 0,
                          mean(metadata$nFeature_RNA[idx & keep_cells]),
                          NA)
    
    result <- rbind(result, base::data.frame(
      group = group,
      cells_before = cells_before,
      cells_after = cells_after,
      keep_rate = keep_rate,
      mean_nCount_after = mean_count,
      mean_nFeature_after = mean_feature
    ))
  }
  
  # 标记异常的group
  result$flag <- ifelse(
    result$keep_rate < mean(result$keep_rate) - sd(result$keep_rate),
    "CHECK",
    "OK"
  )
  
  return(result)
}

# ============================================================================
# 6. Sample分析
# ============================================================================

.analyze_by_sample <- function(metadata, keep_cells, params, sample.by) {
  
  samples <- unique(metadata[[sample.by]])
  
  result <- base::data.frame(
    sample = character(),
    cells_before = numeric(),
    cells_after = numeric(),
    keep_rate = numeric(),
    mean_nCount_after = numeric(),
    mean_nFeature_after = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (sample in samples) {
    idx <- metadata[[sample.by]] == sample
    cells_before <- sum(idx)
    cells_after <- sum(idx & keep_cells)
    keep_rate <- cells_after / cells_before
    
    mean_count <- ifelse(cells_after > 0,
                        mean(metadata$nCount_RNA[idx & keep_cells]),
                        NA)
    mean_feature <- ifelse(cells_after > 0,
                          mean(metadata$nFeature_RNA[idx & keep_cells]),
                          NA)
    
    result <- rbind(result, base::data.frame(
      sample = sample,
      cells_before = cells_before,
      cells_after = cells_after,
      keep_rate = keep_rate,
      mean_nCount_after = mean_count,
      mean_nFeature_after = mean_feature
    ))
  }
  
  # 标记异常的sample
  result$status <- ifelse(
    result$keep_rate < mean(result$keep_rate) - sd(result$keep_rate),
    "PROBLEM",
    "OK"
  )
  
  return(result)
}

# ============================================================================
# 7. Sample × Group二维分析
# ============================================================================

.analyze_sample_group <- function(metadata, keep_cells, params, groups.by, sample.by) {
  
  groups <- unique(metadata[[groups.by]])
  samples <- unique(metadata[[sample.by]])
  
  result <- base::data.frame(
    sample = character(),
    group = character(),
    keep_rate = numeric(),
    n_cells_after = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (sample in samples) {
    for (group in groups) {
      idx <- (metadata[[sample.by]] == sample) & (metadata[[groups.by]] == group)
      cells_before <- sum(idx)
      cells_after <- sum(idx & keep_cells)
      
      keep_rate <- ifelse(cells_before > 0, cells_after / cells_before, NA)
      
      result <- rbind(result, base::data.frame(
        sample = sample,
        group = group,
        keep_rate = keep_rate,
        n_cells_after = cells_after
      ))
    }
  }
  
  return(result)
}

# ============================================================================
# 8. 多参数对比
# ============================================================================

.compare_all_params <- function(all_results, param_list) {
  
  comparison <- base::data.frame(
    param_set = names(all_results),
    cells_before = numeric(length(all_results)),
    cells_after = numeric(length(all_results)),
    cells_filtered = numeric(length(all_results)),
    keep_rate = numeric(length(all_results)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(all_results)) {
    result <- all_results[[i]]
    comparison$cells_before[i] <- result$n_cells_before
    comparison$cells_after[i] <- result$n_cells_after
    comparison$cells_filtered[i] <- result$n_filtered
    comparison$keep_rate[i] <- result$keep_rate
  }
  
  # 添加参数信息
  for (col in names(param_list[[1]])) {
    comparison[[col]] <- sapply(param_list, function(x) x[[col]])
  }
  
  return(comparison)
}

# ============================================================================
# 9. 参数推荐
# ============================================================================

.recommend_best_params <- function(all_results, comparison, is_multi_params) {
  
  if (!is_multi_params) {
    # 单参数情况，只返回使用的参数
    params <- all_results[[1]]$params
    return(list(
      best_set = names(all_results)[1],
      best_params = params,
      recommendation = "Using provided parameters"
    ))
  }
  
  # 多参数情况：选择平衡方案
  # 评分标准：保留细胞数 > 80% & 过滤严格度适中
  
  scores <- base::data.frame(
    param_set = comparison$param_set,
    keep_rate = comparison$keep_rate,
    score = NA
  )
  
  # 计算评分：偏好keep_rate在85%-95%之间
  scores$score <- ifelse(
    scores$keep_rate >= 0.85 & scores$keep_rate <= 0.95,
    100 - abs(scores$keep_rate - 0.90) * 200,
    ifelse(scores$keep_rate > 0.95, 50, 30)
  )
  
  best_idx <- which.max(scores$score)
  best_set <- scores$param_set[best_idx]
  
  return(list(
    best_balanced = best_set,
    scores = scores,
    all_recommendations = comparison %>%
      dplyr::arrange(desc(score)) %>%
      dplyr::mutate(rank = row_number())
  ))
}

# ============================================================================
# 10. 诊断分析
# ============================================================================

.run_diagnostics <- function(all_results, groups.by, sample.by) {
  
  main_result <- all_results[[1]]
  
  problems <- list()
  suggestions <- character()
  
  # 检查group问题
  if (!is.null(main_result$by_group)) {
    by_group <- main_result$by_group
    problem_groups <- by_group$group[by_group$flag == "CHECK"]
    if (length(problem_groups) > 0) {
      problems$problematic_groups <- problem_groups
      suggestions <- c(suggestions,
        paste("Group(s) with low keep rate:", paste(problem_groups, collapse=", "))
      )
    }
  }
  
  # 检查sample问题
  if (!is.null(main_result$by_sample)) {
    by_sample <- main_result$by_sample
    problem_samples <- by_sample$sample[by_sample$status == "PROBLEM"]
    if (length(problem_samples) > 0) {
      problems$problematic_samples <- problem_samples
      suggestions <- c(suggestions,
        paste("Sample(s) with low keep rate:", paste(problem_samples, collapse=", "))
      )
    }
  }
  
  # 总体质量评估
  keep_rate <- main_result$keep_rate
  overall_quality <- ifelse(keep_rate > 0.95, "Excellent",
                           ifelse(keep_rate > 0.90, "Good",
                                 ifelse(keep_rate > 0.80, "Acceptable",
                                       "Poor")))
  
  suggestions <- c(suggestions, 
                  paste("Overall quality:", overall_quality))
  
  return(list(
    problems = problems,
    overall_quality = overall_quality,
    keep_rate = keep_rate,
    suggestions = suggestions
  ))
}

# ============================================================================
# 10.5 QC分布图生成函数 (v1.0.0新增)
# ============================================================================

.create_qc_distribution <- function(metadata, groups.by, params, percentiles,
                                     report_dir, verbose = TRUE) {

  # 检查是否有groups.by且有效
  has_groups <- !is.null(groups.by) && groups.by %in% colnames(metadata)

  # 计算百分位统计
  qc_stats <- .calculate_qc_percentiles(
    metadata = metadata,
    groups.by = if(has_groups) groups.by else NULL,
    percentiles = percentiles
  )

  # 生成图表
  plots <- list()

  if (has_groups) {
    # 分组模式
    if (verbose) cat("  → Generating grouped QC distribution plots...\n")
    plots$by_group <- .create_grouped_qc_violin(
      metadata = metadata,
      groups.by = groups.by,
      percentiles = percentiles,
      params = params,
      report_dir = report_dir
    )
  }

  # 整体模式（始终生成）
  if (verbose) cat("  → Generating overall QC distribution plots...\n")
  plots$overall <- .create_overall_qc_violin(
    metadata = metadata,
    percentiles = percentiles,
    params = params,
    report_dir = report_dir,
    has_groups = has_groups
  )

  return(list(
    stats = qc_stats,
    plots = plots,
    has_groups = has_groups,
    groups.by = if(has_groups) groups.by else NULL
  ))
}

# ============================================================================
# 10.6 计算QC指标百分位统计
# ============================================================================

.calculate_qc_percentiles <- function(metadata, groups.by, percentiles) {

  # 定义QC指标
  qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt")

  # 整体统计
  overall_stats <- list()

  for (metric in qc_metrics) {
    if (metric %in% colnames(metadata)) {
      values <- metadata[[metric]]
      overall_stats[[metric]] <- list(
        values = values,
        quantiles = stats::quantile(values, probs = percentiles / 100, na.rm = TRUE),
        mean = mean(values, na.rm = TRUE),
        median = stats::median(values, na.rm = TRUE),
        sd = sd(values, na.rm = TRUE)
      )
    }
  }

  # 分组统计
  group_stats <- NULL
  if (!is.null(groups.by) && groups.by %in% colnames(metadata)) {
    group_stats <- list()
    groups <- unique(metadata[[groups.by]])

    for (metric in qc_metrics) {
      if (metric %in% colnames(metadata)) {
        group_stats[[metric]] <- list()
        for (grp in groups) {
          idx <- metadata[[groups.by]] == grp
          values <- metadata[[metric]][idx]
          group_stats[[metric]][[grp]] <- list(
            values = values,
            quantiles = stats::quantile(values, probs = percentiles / 100, na.rm = TRUE),
            mean = mean(values, na.rm = TRUE),
            median = stats::median(values, na.rm = TRUE),
            n_cells = sum(idx)
          )
        }
      }
    }
  }

  # 计算各区间的细胞数量
  cell_counts <- .calculate_cell_counts(metadata, percentiles, qc_metrics)

  return(list(
    overall = overall_stats,
    by_group = group_stats,
    cell_counts = cell_counts,
    percentiles = percentiles
  ))
}

# ============================================================================
# 10.7 计算各百分位区间的细胞数量（v2.1.1优化）
# ============================================================================

.calculate_cell_counts <- function(metadata, percentiles, qc_metrics) {

  results <- list()
  n_cells_total <- nrow(metadata)

  # 为下限表格准备数据
  low_data <- list()

  for (metric in c("nFeature_RNA", "nCount_RNA")) {
    if (metric %in% colnames(metadata)) {
      values <- metadata[[metric]]

      # 计算百分位阈值
      q_2.5 <- stats::quantile(values, 0.025, na.rm = TRUE)
      q_5 <- stats::quantile(values, 0.05, na.rm = TRUE)

      # 计算累计移除细胞数
      cells_2.5 <- sum(values < q_2.5, na.rm = TRUE)
      cells_5 <- sum(values < q_5, na.rm = TRUE)

      # 计算Delta
      delta_value <- q_5 - q_2.5
      delta_value_pct <- if (q_2.5 > 0) (delta_value / q_2.5) * 100 else 0

      delta_cells <- cells_5 - cells_2.5
      delta_cells_pct <- if (cells_2.5 > 0) (delta_cells / cells_2.5) * 100 else 0

      # 细胞数占总数的百分比
      cells_2.5_pct <- (cells_2.5 / n_cells_total) * 100
      cells_5_pct <- (cells_5 / n_cells_total) * 100

      metric_name <- if (metric == "nFeature_RNA") "nFeature" else "nCount"

      low_data[[metric_name]] <- list(
        range_2.5_5 = sprintf("%.0f ~ %.0f", q_2.5, q_5),
        delta_value = sprintf("+%.0f (%.1f%%)", delta_value, delta_value_pct),
        cells_remove = sprintf("%.0f (%.1f%%) / %.0f (%.1f%%)",
                               cells_2.5, cells_2.5_pct,
                               cells_5, cells_5_pct),
        delta_cells = sprintf("+%.0f (%.1f%%)", delta_cells, delta_cells_pct),
        # 原始数据用于建议
        raw = list(
          q_2.5 = q_2.5, q_5 = q_5,
          cells_2.5 = cells_2.5, cells_5 = cells_5,
          delta_value = delta_value, delta_value_pct = delta_value_pct,
          delta_cells = delta_cells, delta_cells_pct = delta_cells_pct
        )
      )
    }
  }

  results$low_threshold <- low_data

  # 为上限表格准备数据
  high_data <- list()

  for (metric in c("nFeature_RNA", "nCount_RNA")) {
    if (metric %in% colnames(metadata)) {
      values <- metadata[[metric]]

      # 计算百分位阈值
      q_95 <- stats::quantile(values, 0.95, na.rm = TRUE)
      q_97.5 <- stats::quantile(values, 0.975, na.rm = TRUE)

      # 计算累计移除细胞数
      cells_95 <- sum(values > q_95, na.rm = TRUE)
      cells_97.5 <- sum(values > q_97.5, na.rm = TRUE)

      # 计算Delta
      delta_value <- q_97.5 - q_95
      delta_value_pct <- if (q_95 > 0) (delta_value / q_95) * 100 else 0

      delta_cells <- cells_95 - cells_97.5  # 从95%到97.5%减少的细胞
      delta_cells_pct <- if (cells_95 > 0) (delta_cells / cells_95) * 100 else 0

      # 细胞数占总数的百分比
      cells_95_pct <- (cells_95 / n_cells_total) * 100
      cells_97.5_pct <- (cells_97.5 / n_cells_total) * 100

      metric_name <- if (metric == "nFeature_RNA") "nFeature" else "nCount"

      high_data[[metric_name]] <- list(
        range_95_97.5 = sprintf("%.0f ~ %.0f", q_95, q_97.5),
        delta_value = sprintf("+%.0f (%.1f%%)", delta_value, delta_value_pct),
        cells_remove = sprintf("%.0f (%.1f%%) / %.0f (%.1f%%)",
                               cells_95, cells_95_pct,
                               cells_97.5, cells_97.5_pct),
        delta_cells = sprintf("-%.0f (%.1f%%)", delta_cells, delta_cells_pct),
        # 原始数据用于建议
        raw = list(
          q_95 = q_95, q_97.5 = q_97.5,
          cells_95 = cells_95, cells_97.5 = cells_97.5,
          delta_value = delta_value, delta_value_pct = delta_value_pct,
          delta_cells = delta_cells, delta_cells_pct = delta_cells_pct
        )
      )
    }
  }

  results$high_threshold <- high_data

  # 线粒体数据（用于文字说明）
  if ("percent.mt" %in% colnames(metadata)) {
    values <- metadata$percent.mt
    mt_summary <- list(
      gt_5 = list(count = sum(values > 5, na.rm = TRUE), pct = (sum(values > 5, na.rm = TRUE) / n_cells_total) * 100),
      gt_10 = list(count = sum(values > 10, na.rm = TRUE), pct = (sum(values > 10, na.rm = TRUE) / n_cells_total) * 100),
      gt_15 = list(count = sum(values > 15, na.rm = TRUE), pct = (sum(values > 15, na.rm = TRUE) / n_cells_total) * 100),
      gt_20 = list(count = sum(values > 20, na.rm = TRUE), pct = (sum(values > 20, na.rm = TRUE) / n_cells_total) * 100)
    )
    results$mt_summary <- mt_summary
  }

  results$n_cells_total <- n_cells_total

  return(results)
}

# ============================================================================
# 10.8 生成整体QC小提琴图（v2.1.1优化）
# ============================================================================

.create_overall_qc_violin <- function(metadata, percentiles, params, report_dir, has_groups) {

  # 分别创建三个独立的图
  plots <- list()

  # 1. nCount_RNA
  if ("nCount_RNA" %in% colnames(metadata)) {
    values <- metadata$nCount_RNA
    q_5 <- stats::quantile(values, 0.05, na.rm = TRUE)
    q_95 <- stats::quantile(values, 0.95, na.rm = TRUE)
    y_range <- range(values, na.rm = TRUE)
    y_expansion <- (y_range[2] - y_range[1]) * 0.1

    p1 <- ggplot2::ggplot(data.frame(value = values), ggplot2::aes(x = 1, y = value)) +
      ggplot2::geom_violin(fill = "#667eea", alpha = 0.6, scale = "width", width = 0.7) +
      ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8,
                            fill = "white", color = "#333333") +
      ggplot2::geom_hline(yintercept = q_5, linetype = "dashed", color = "red", linewidth = 0.4) +
      ggplot2::geom_hline(yintercept = q_95, linetype = "dashed", color = "red", linewidth = 0.4) +
      ggplot2::annotate("text", x = 1.6, y = q_5, label = sprintf("%.0f", q_5),
                        color = "red", size = 3.2, hjust = 0, fontface = "bold") +
      ggplot2::annotate("text", x = 1.6, y = q_95, label = sprintf("%.0f", q_95),
                        color = "red", size = 3.2, hjust = 0, fontface = "bold") +
      ggplot2::coord_flip(ylim = c(y_range[1] - y_expansion, y_range[2] + y_expansion * 2)) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
        plot.margin = ggplot2::margin(5, 50, 5, 5)
      ) +
      ggplot2::labs(title = "nCount (UMI)") +
      ggplot2::xlim(0.5, 2.8)

    plots$nCount <- p1
  }

  # 2. nFeature_RNA
  if ("nFeature_RNA" %in% colnames(metadata)) {
    values <- metadata$nFeature_RNA
    q_5 <- stats::quantile(values, 0.05, na.rm = TRUE)
    q_95 <- stats::quantile(values, 0.95, na.rm = TRUE)
    y_range <- range(values, na.rm = TRUE)
    y_expansion <- (y_range[2] - y_range[1]) * 0.1

    p2 <- ggplot2::ggplot(data.frame(value = values), ggplot2::aes(x = 1, y = value)) +
      ggplot2::geom_violin(fill = "#667eea", alpha = 0.6, scale = "width", width = 0.7) +
      ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8,
                            fill = "white", color = "#333333") +
      ggplot2::geom_hline(yintercept = q_5, linetype = "dashed", color = "red", linewidth = 0.4) +
      ggplot2::geom_hline(yintercept = q_95, linetype = "dashed", color = "red", linewidth = 0.4) +
      ggplot2::annotate("text", x = 1.6, y = q_5, label = sprintf("%.0f", q_5),
                        color = "red", size = 3.2, hjust = 0, fontface = "bold") +
      ggplot2::annotate("text", x = 1.6, y = q_95, label = sprintf("%.0f", q_95),
                        color = "red", size = 3.2, hjust = 0, fontface = "bold") +
      ggplot2::coord_flip(ylim = c(y_range[1] - y_expansion, y_range[2] + y_expansion * 2)) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
        plot.margin = ggplot2::margin(5, 50, 5, 5)
      ) +
      ggplot2::labs(title = "nFeature (Genes)") +
      ggplot2::xlim(0.5, 2.8)

    plots$nFeature <- p2
  }

  # 3. percent.mt - y轴从0开始，无负坐标
  if ("percent.mt" %in% colnames(metadata)) {
    values <- metadata$percent.mt
    y_max <- max(values, na.rm = TRUE)
    y_expansion <- y_max * 0.15

    p3 <- ggplot2::ggplot(data.frame(value = values), ggplot2::aes(x = 1, y = value)) +
      ggplot2::geom_violin(fill = "#667eea", alpha = 0.6, scale = "width", width = 0.7) +
      ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8,
                            fill = "white", color = "#333333") +
      # 关键：expand = FALSE 防止自动扩展到负值
      ggplot2::coord_flip(ylim = c(0, y_max + y_expansion), expand = FALSE) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
        plot.margin = ggplot2::margin(5, 10, 5, 5)
      ) +
      ggplot2::labs(title = "Mitochondrial %%") +
      ggplot2::xlim(0.5, 1.5)

    plots$percent.mt <- p3
  }

  # 组合三个图
  p_combined <- patchwork::wrap_plots(plots, ncol = 3) +
    patchwork::plot_annotation(
      title = "Pre-filtering QC Distribution",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12)
      )
    )

  # 保存图表 - 压缩尺寸
  plot_file <- base::file.path(report_dir, "00_qc_distribution_overall.png")
  ggplot2::ggsave(plot_file, p_combined, width = 9, height = 2.5, dpi = 300)

  return(list(plot = p_combined, file = plot_file))
}

# ============================================================================
# 10.9 生成分组QC小提琴图（v2.1.0优化）
# ============================================================================

.create_grouped_qc_violin <- function(metadata, groups.by, percentiles, params, report_dir) {

  # 准备数据
  plot_data <- data.frame(
    cell_id = rownames(metadata),
    group = as.character(metadata[[groups.by]]),
    nCount_RNA = metadata$nCount_RNA,
    nFeature_RNA = metadata$nFeature_RNA,
    percent.mt = metadata$percent.mt
  )

  # 分别创建三个独立的图，以精确控制Y轴
  plots <- list()

  # 通用主题
  common_theme <- ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 11),
      legend.position = "none"
    )

  # 1. nCount_RNA
  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = nCount_RNA, fill = group)) +
    ggplot2::geom_violin(alpha = 0.6, scale = "width", width = 0.8) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8,
                          fill = "white", color = "#333333") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(title = "nCount (UMI)") +
    common_theme

  plots$nCount <- p1

  # 2. nFeature_RNA
  p2 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = nFeature_RNA, fill = group)) +
    ggplot2::geom_violin(alpha = 0.6, scale = "width", width = 0.8) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8,
                          fill = "white", color = "#333333") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(title = "nFeature (Genes)") +
    common_theme

  plots$nFeature <- p2

  # 3. percent.mt - Y轴从0开始，无负值
  y_max <- max(plot_data$percent.mt, na.rm = TRUE)
  y_expansion <- y_max * 0.15

  p3 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = percent.mt, fill = group)) +
    ggplot2::geom_violin(alpha = 0.6, scale = "width", width = 0.8) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8,
                          fill = "white", color = "#333333") +
    # 关键：使用 coord_cartesian 的 expand = FALSE 防止Y轴扩展到负值
    ggplot2::coord_cartesian(ylim = c(0, y_max + y_expansion), expand = FALSE) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(title = "Mitochondrial %") +
    common_theme

  plots$percent.mt <- p3

  # 组合三个图 - 按照 nCount, nFeature, MT 的顺序
  p_combined <- patchwork::wrap_plots(plots$nCount, plots$nFeature, plots$percent.mt, ncol = 3) +
    patchwork::plot_annotation(
      title = paste0("QC Distribution by ", groups.by),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 12)
      )
    )

  # 保存图表
  plot_file <- base::file.path(report_dir, "00_qc_distribution_by_group.png")
  ggplot2::ggsave(plot_file, p_combined, width = 9, height = 3.5, dpi = 300)

  return(list(plot = p_combined, file = plot_file))
}

# ============================================================================
# 11. HTML报告生成
# ============================================================================

.generate_html_report <- function(all_results, comparison, recommend, diagnostics,
                                 report_dir, groups.by, sample.by, is_multi_params,
                                 qc_distribution = NULL, show.qc.distribution = TRUE) {

  report_file <- base::file.path(report_dir, "scQCFilter_report.html")

  # 获取主结果
  main_result <- all_results[[1]]
  
  # 开始构建HTML
  html_content <- '
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>scQCFilter Quality Control Report</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            min-height: 100vh;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            overflow: hidden;
        }
        
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px 30px;
            text-align: center;
        }
        
        .header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
        }
        
        .header p {
            font-size: 1.1em;
            opacity: 0.9;
        }
        
        .content {
            padding: 40px 30px;
        }
        
        .section {
            margin-bottom: 40px;
        }
        
        .section h2 {
            color: #667eea;
            border-bottom: 3px solid #667eea;
            padding-bottom: 15px;
            margin-bottom: 20px;
            font-size: 1.8em;
        }
        
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }
        
        .summary-card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }
        
        .summary-card .label {
            font-size: 0.9em;
            opacity: 0.9;
            margin-bottom: 10px;
        }
        
        .summary-card .value {
            font-size: 2em;
            font-weight: bold;
        }
        
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
        }
        
        thead {
            background: #f8f9fa;
            border-bottom: 2px solid #dee2e6;
        }
        
        th {
            padding: 15px;
            text-align: left;
            font-weight: 600;
            color: #333;
        }
        
        td {
            padding: 12px 15px;
            border-bottom: 1px solid #dee2e6;
        }
        
        tbody tr:hover {
            background: #f8f9fa;
        }
        
        .badge {
            display: inline-block;
            padding: 5px 12px;
            border-radius: 20px;
            font-size: 0.85em;
            font-weight: 600;
        }
        
        .badge-ok {
            background: #d4edda;
            color: #155724;
        }
        
        .badge-check {
            background: #fff3cd;
            color: #856404;
        }
        
        .badge-problem {
            background: #f8d7da;
            color: #721c24;
        }
        
        .badge-good {
            background: #d1ecf1;
            color: #0c5460;
        }
        
        .parameters {
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            margin: 20px 0;
            border-left: 4px solid #667eea;
        }
        
        .parameters h3 {
            color: #667eea;
            margin-bottom: 15px;
        }
        
        .parameters-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
        }
        
        .parameter-item {
            background: white;
            padding: 12px;
            border-radius: 6px;
            border-left: 3px solid #667eea;
        }
        
        .parameter-item .name {
            font-size: 0.9em;
            color: #666;
            margin-bottom: 5px;
        }
        
        .parameter-item .value {
            font-size: 1.2em;
            font-weight: bold;
            color: #667eea;
        }
        
        .diagnostics {
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #dc3545;
        }
        
        .diagnostics h3 {
            color: #dc3545;
            margin-bottom: 15px;
        }
        
        .suggestion {
            background: white;
            padding: 10px 15px;
            margin: 10px 0;
            border-radius: 6px;
            border-left: 3px solid #dc3545;
        }
        
        .suggestion.info {
            border-left-color: #0066cc;
        }
        
        .suggestion.warning {
            border-left-color: #ff9800;
        }
        
        .footer {
            background: #f8f9fa;
            padding: 20px 30px;
            text-align: center;
            color: #666;
            border-top: 1px solid #dee2e6;
            font-size: 0.9em;
        }
        
        .figure {
            text-align: center;
            margin: 30px 0;
        }
        
        .figure img {
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }
        
        .figure-caption {
            font-size: 0.9em;
            color: #666;
            margin-top: 10px;
            font-style: italic;
        }
        
        .tabs {
            display: flex;
            flex-wrap: wrap;
            border-bottom: 2px solid #dee2e6;
            margin-bottom: 20px;
        }
        
        .tab {
            padding: 12px 20px;
            cursor: pointer;
            border: none;
            background: transparent;
            color: #666;
            font-weight: 500;
            border-bottom: 3px solid transparent;
            transition: all 0.3s;
        }
        
        .tab.active {
            color: #667eea;
            border-bottom-color: #667eea;
        }
        
        .tab-content {
            display: none;
        }
        
        .tab-content.active {
            display: block;
        }
    </style>
</head>
<body>
    <div class="container">
'
  
  # Header
  html_content <- paste0(html_content, '
        <div class="header">
            <h1>📊 scQCFilter Quality Control Report</h1>
            <p>Single-cell RNA-seq Quality Control and Analysis</p>
        </div>
        
        <div class="content">
')
  
  # Summary Section
  html_content <- paste0(html_content, sprintf('
            <div class="section">
                <h2>Summary Statistics</h2>
                <div class="summary-grid">
                    <div class="summary-card">
                        <div class="label">Total Cells Before</div>
                        <div class="value">%d</div>
                    </div>
                    <div class="summary-card">
                        <div class="label">Cells Retained</div>
                        <div class="value">%d</div>
                    </div>
                    <div class="summary-card">
                        <div class="label">Cells Filtered</div>
                        <div class="value">%d</div>
                    </div>
                    <div class="summary-card">
                        <div class="label">Keep Rate</div>
                        <div class="value">%.1f%%</div>
                    </div>
                </div>
            </div>
',
    main_result$n_cells_before,
    main_result$n_cells_after,
    main_result$n_filtered,
    main_result$keep_rate * 100
  ))

  # QC Distribution Section (v2.1.0 优化)
  if (show.qc.distribution && !is.null(qc_distribution)) {
    html_content <- paste0(html_content, '
            <div class="section">
                <h2>📊 Pre-filtering QC Distribution</h2>
                <p style="color: #666; margin-bottom: 15px;">
                    Distribution of QC metrics before filtering. Red dashed lines indicate recommended threshold references.
                </p>
')

    # 整体分布图
    overall_plot_file <- qc_distribution$plots$overall$file
    if (!is.null(overall_plot_file) && base::file.exists(overall_plot_file)) {
      html_content <- paste0(html_content, '
                <div class="figure">
                    <img src="00_qc_distribution_overall.png" alt="QC Distribution Overview" style="max-width: 100%%;">
                </div>
')
    }

    # 分组分布图（如果有）
    if (qc_distribution$has_groups && !is.null(qc_distribution$plots$by_group)) {
      group_plot_file <- qc_distribution$plots$by_group$file
      if (!is.null(group_plot_file) && base::file.exists(group_plot_file)) {
        html_content <- paste0(html_content, sprintf('
                <h3 style="color: #667eea; margin-top: 20px;">Distribution by %s</h3>
                <div class="figure">
                    <img src="00_qc_distribution_by_group.png" alt="QC Distribution by Group" style="max-width: 100%%;">
                </div>
',
          qc_distribution$groups.by
        ))
      }
    }

    # 统计汇总表 - 两个表格横向排列
    if (!is.null(qc_distribution$stats$cell_counts)) {
      html_content <- paste0(html_content, '
                <h3 style="color: #667eea; margin-top: 20px;">Cell Count Statistics</h3>
                <div style="display: flex; flex-wrap: wrap; gap: 30px; margin-top: 15px;">
')

      # 表格1: nCount 和 nFeature（横向合并）
      if (!is.null(qc_distribution$stats$cell_counts$nCount_RNA) &&
          !is.null(qc_distribution$stats$cell_counts$nFeature_RNA)) {

        ncount_data <- qc_distribution$stats$cell_counts$nCount_RNA
        nfeature_data <- qc_distribution$stats$cell_counts$nFeature_RNA

        html_content <- paste0(html_content, '
                  <div style="flex: 1; min-width: 400px;">
                    <table style="width: 100%%;">
                        <thead>
                            <tr>
                                <th>Range</th>
                                <th>nCount Threshold</th>
                                <th>nFeature Threshold</th>
                                <th>Cells</th>
                                <th>%%</th>
                            </tr>
                        </thead>
                        <tbody>
')

        for (i in seq_len(nrow(ncount_data))) {
          row_class <- ""
          if (i == 1 || i == 5) row_class <- 'style="background-color: #fff3cd;"'

          html_content <- paste0(html_content, sprintf('
                            <tr %s>
                                <td>%s</td>
                                <td>%s</td>
                                <td>%s</td>
                                <td>%d / %d</td>
                                <td>%.1f%% / %.1f%%</td>
                            </tr>
',
            row_class,
            ncount_data$range[i],
            ncount_data$threshold[i],
            nfeature_data$threshold[i],
            ncount_data$n_cells[i],
            nfeature_data$n_cells[i],
            ncount_data$percentage[i],
            nfeature_data$percentage[i]
          ))
        }

        html_content <- paste0(html_content, '
                        </tbody>
                    </table>
                  </div>
')
      }

      # 表格2: percent.mt
      if (!is.null(qc_distribution$stats$cell_counts$percent.mt)) {
        mt_data <- qc_distribution$stats$cell_counts$percent.mt

        html_content <- paste0(html_content, '
                  <div style="flex: 0 0 280px;">
                    <table style="width: 100%%;">
                        <thead>
                            <tr>
                                <th>MT %%</th>
                                <th>Cells</th>
                                <th>%%</th>
                            </tr>
                        </thead>
                        <tbody>
')

        for (i in seq_len(nrow(mt_data))) {
          row_class <- ""
          if (i == 4) row_class <- 'style="background-color: #f8d7da;"'

          html_content <- paste0(html_content, sprintf('
                            <tr %s>
                                <td>%s</td>
                                <td>%d</td>
                                <td>%.1f%%</td>
                            </tr>
',
            row_class,
            mt_data$range[i],
            mt_data$n_cells[i],
            mt_data$percentage[i]
          ))
        }

        html_content <- paste0(html_content, '
                        </tbody>
                    </table>
                  </div>
')
      }

      html_content <- paste0(html_content, '
                </div>
')
    }

    # 如果没有groups.by，显示提示
    if (!qc_distribution$has_groups) {
      html_content <- paste0(html_content, '
                <div class="suggestion info" style="margin-top: 20px;">
                    <strong>💡 Tip:</strong> Consider using the <code>groups.by</code> parameter to analyze QC distribution
                    by cell type or experimental condition. Different groups may have distinct QC profiles.
                    <br><br>
                    Example: <code>scQCFilter(seurat_obj, groups.by = "cell_type", ...)</code>
                </div>
')
    }

    html_content <- paste0(html_content, '
            </div>
')
  }

  # Parameters Section
  html_content <- paste0(html_content, '
            <div class="section">
                <h2>Filtering Parameters</h2>
                <div class="parameters">
                    <h3>Applied QC Thresholds</h3>
                    <div class="parameters-grid">
')
  
  params <- main_result$params
  html_content <- paste0(html_content, sprintf('
                        <div class="parameter-item">
                            <div class="name">Mitochondrial %% (max)</div>
                            <div class="value">%.1f%%</div>
                        </div>
                        <div class="parameter-item">
                            <div class="name">UMI Count (min)</div>
                            <div class="value">%d</div>
                        </div>
                        <div class="parameter-item">
                            <div class="name">UMI Count (max)</div>
                            <div class="value">%s</div>
                        </div>
                        <div class="parameter-item">
                            <div class="name">Gene Count (min)</div>
                            <div class="value">%d</div>
                        </div>
                        <div class="parameter-item">
                            <div class="name">Gene Count (max)</div>
                            <div class="value">%s</div>
                        </div>
',
    params$percent_mt_max,
    params$nCount_min,
    ifelse(is.infinite(params$nCount_max), "Inf", params$nCount_max),
    params$nFeature_min,
    ifelse(is.infinite(params$nFeature_max), "Inf", params$nFeature_max)
  ))
  
  html_content <- paste0(html_content, '
                    </div>
                </div>
            </div>
')
  
  # Parameter分层表格
  if (!is.null(main_result$by_parameter)) {
    html_content <- paste0(html_content, '
            <div class="section">
                <h2>Filtering Contribution Analysis</h2>
                <p>Shows how much each parameter contributed to cell filtering:</p>
                <table>
                    <thead>
                        <tr>
                            <th>Parameter</th>
                            <th>Cells Filtered</th>
                            <th>% of Total Filtered</th>
                            <th>Severity Rank</th>
                        </tr>
                    </thead>
                    <tbody>
')
    
    by_param <- main_result$by_parameter
    for (i in seq_len(nrow(by_param))) {
      rank_text <- ifelse(is.na(by_param$severity_rank[i]), "-", by_param$severity_rank[i])
      html_content <- paste0(html_content, sprintf('
                        <tr>
                            <td>%s</td>
                            <td>%d</td>
                            <td>%.1f%%</td>
                            <td>%s</td>
                        </tr>
',
        by_param$parameter[i],
        by_param$n_filtered[i],
        by_param$pct_of_total[i] * 100,
        rank_text
      ))
    }
    
    html_content <- paste0(html_content, '
                    </tbody>
                </table>
            </div>
')
  }
  
  # Group分析表格
  if (!is.null(main_result$by_group)) {
    html_content <- paste0(html_content, sprintf('
            <div class="section">
                <h2>Analysis by %s</h2>
                <table>
                    <thead>
                        <tr>
                            <th>%s</th>
                            <th>Cells Before</th>
                            <th>Cells After</th>
                            <th>Keep Rate</th>
                            <th>Status</th>
                        </tr>
                    </thead>
                    <tbody>
',
    gsub("_", " ", groups.by),
    gsub("_", " ", groups.by)
  ))
    
    by_group <- main_result$by_group
    for (i in seq_len(nrow(by_group))) {
      badge_class <- ifelse(by_group$flag[i] == "OK", "badge-ok", "badge-check")
      html_content <- paste0(html_content, sprintf('
                        <tr>
                            <td>%s</td>
                            <td>%d</td>
                            <td>%d</td>
                            <td>%.1f%%</td>
                            <td><span class="badge %s">%s</span></td>
                        </tr>
',
        by_group$group[i],
        by_group$cells_before[i],
        by_group$cells_after[i],
        by_group$keep_rate[i] * 100,
        badge_class,
        by_group$flag[i]
      ))
    }
    
    html_content <- paste0(html_content, '
                    </tbody>
                </table>
            </div>
')
  }
  
  # Sample分析表格
  if (!is.null(main_result$by_sample)) {
    html_content <- paste0(html_content, sprintf('
            <div class="section">
                <h2>Analysis by %s</h2>
                <table>
                    <thead>
                        <tr>
                            <th>%s</th>
                            <th>Cells Before</th>
                            <th>Cells After</th>
                            <th>Keep Rate</th>
                            <th>Status</th>
                        </tr>
                    </thead>
                    <tbody>
',
    gsub("_", " ", sample.by),
    gsub("_", " ", sample.by)
  ))
    
    by_sample <- main_result$by_sample
    for (i in seq_len(nrow(by_sample))) {
      badge_class <- ifelse(by_sample$status[i] == "OK", "badge-ok", "badge-problem")
      html_content <- paste0(html_content, sprintf('
                        <tr>
                            <td>%s</td>
                            <td>%d</td>
                            <td>%d</td>
                            <td>%.1f%%</td>
                            <td><span class="badge %s">%s</span></td>
                        </tr>
',
        by_sample$sample[i],
        by_sample$cells_before[i],
        by_sample$cells_after[i],
        by_sample$keep_rate[i] * 100,
        badge_class,
        by_sample$status[i]
      ))
    }
    
    html_content <- paste0(html_content, '
                    </tbody>
                </table>
            </div>
')
  }
  
  # 多参数对比表格（如果有）
  if (!is.null(comparison)) {
    html_content <- paste0(html_content, '
            <div class="section">
                <h2>Multi-Parameter Comparison</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Parameter Set</th>
                            <th>mt.max</th>
                            <th>count.min</th>
                            <th>feature.min</th>
                            <th>Cells After</th>
                            <th>Keep Rate</th>
                        </tr>
                    </thead>
                    <tbody>
')
    
    for (i in seq_len(nrow(comparison))) {
      html_content <- paste0(html_content, sprintf('
                        <tr>
                            <td>%s</td>
                            <td>%.1f</td>
                            <td>%d</td>
                            <td>%d</td>
                            <td>%d</td>
                            <td>%.1f%%</td>
                        </tr>
',
        comparison$param_set[i],
        comparison$percent_mt_max[i],
        comparison$nCount_min[i],
        comparison$nFeature_min[i],
        comparison$cells_after[i],
        comparison$keep_rate[i] * 100
      ))
    }
    
    html_content <- paste0(html_content, '
                    </tbody>
                </table>
            </div>
')
  }
  
  # 诊断和建议
  html_content <- paste0(html_content, sprintf('
            <div class="section">
                <h2>Diagnosis & Recommendations</h2>
                <div class="diagnostics">
                    <h3>Overall Quality: <span class="badge badge-good">%s</span></h3>
                    <div class="suggestion info">
                        <strong>Keep Rate:</strong> %.1f%% of cells retained
                    </div>
',
    diagnostics$overall_quality,
    diagnostics$keep_rate * 100
  ))
    
    if (length(diagnostics$suggestions) > 0) {
      for (sugg in diagnostics$suggestions) {
        html_content <- paste0(html_content, sprintf('
                    <div class="suggestion warning">
                        %s
                    </div>
',
          sugg
        ))
      }
    }
    
    html_content <- paste0(html_content, '
                </div>
            </div>
')
  
  # 关闭
  html_content <- paste0(html_content, '
        </div>
        
        <div class="footer">
            <p>Generated by <strong>scQCFilter</strong> | Single-cell Quality Control Tool</p>
            <p>Report generated on ', base::Sys.time(), '</p>
        </div>
    </div>
</body>
</html>
')
  
  # 写入文件
  base::writeLines(html_content, report_file)
  
  return(report_file)
}

# ============================================================================
# 12. 可视化函数
# ============================================================================

.create_visualizations <- function(all_results, comparison, diagnostics,
                                  groups.by, sample.by, report_dir, metadata) {
  
  main_result <- all_results[[1]]
  
  plots <- list()
  
  # 1. Parameter分层饼图
  if (!is.null(main_result$by_parameter)) {
    by_param <- main_result$by_parameter
    by_param_plot <- by_param %>%
      dplyr::filter(parameter != "combined") %>%
      ggplot2::ggplot(ggplot2::aes(x = "", y = n_filtered, fill = parameter)) +
      ggplot2::geom_bar(stat = "identity", width = 1) +
      ggplot2::coord_polar("y", start = 0) +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "right") +
      ggtitle("Parameter Contribution to Filtering")
    
    ggplot2::ggsave(base::file.path(report_dir, "01_parameter_contribution.png"), 
           by_param_plot, width = 8, height = 6, dpi = 300)
    plots$parameter_contribution <- by_param_plot
  }
  
  # 2. Group Keep Rate柱状图
  if (!is.null(main_result$by_group)) {
    group_plot <- main_result$by_group %>%
      ggplot2::ggplot(ggplot2::aes(x = reorder(group, -keep_rate), y = keep_rate, fill = flag)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = c("OK" = "#28a745", "CHECK" = "#ffc107")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(x = "Group", y = "Keep Rate", title = paste("Keep Rate by", groups.by))
    
    ggplot2::ggsave(base::file.path(report_dir, "02_group_analysis.png"),
           group_plot, width = 8, height = 6, dpi = 300)
    plots$group_analysis <- group_plot
  }
  
  # 3. Sample Keep Rate柱状图
  if (!is.null(main_result$by_sample)) {
    sample_plot <- main_result$by_sample %>%
      ggplot2::ggplot(ggplot2::aes(x = reorder(sample, -keep_rate), y = keep_rate, fill = status)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = c("OK" = "#28a745", "PROBLEM" = "#dc3545")) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(x = "Sample", y = "Keep Rate", title = paste("Keep Rate by", sample.by))
    
    ggplot2::ggsave(base::file.path(report_dir, "03_sample_analysis.png"),
           sample_plot, width = 8, height = 6, dpi = 300)
    plots$sample_analysis <- sample_plot
  }
  
  return(plots)
}

# ============================================================================
# 13. S3 方法：print
# ============================================================================

#' @export
print.scQCFilter <- function(x, ...) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════\n")
  cat("                    scQCFilter Results                  \n")
  cat("═══════════════════════════════════════════════════════\n\n")
  
  # 摘要统计
  cat("Summary Statistics:\n")
  cat("─────────────────────────────────────────────────────\n")
  cat(sprintf("  Total cells before filtering: %d\n", x$main_result$n_cells_before))
  cat(sprintf("  Total cells after filtering:  %d\n", x$main_result$n_cells_after))
  cat(sprintf("  Cells filtered out:           %d\n", x$main_result$n_filtered))
  cat(sprintf("  Keep rate:                   %.1f%%\n", x$main_result$keep_rate * 100))
  cat("\n")
  
  # 参数
  cat("Filtering Parameters:\n")
  cat("─────────────────────────────────────────────────────\n")
  params <- x$main_result$params
  cat(sprintf("  Mitochondrial %% (max):       %.1f%%\n", params$percent_mt_max))
  cat(sprintf("  UMI count (min):             %d\n", params$nCount_min))
  cat(sprintf("  Gene count (min):            %d\n", params$nFeature_min))
  cat("\n")
  
  # 报告文件
  cat("Report:\n")
  cat("─────────────────────────────────────────────────────\n")
  cat(sprintf("  HTML report saved to: %s\n", x$report_file))
  cat("\n")
  
  cat("═══════════════════════════════════════════════════════\n\n")
  
}

# ============================================================================
# 导出函数
# ============================================================================


# ============================================================================
# 14. scQCExplore - QC分布探索函数（v2.1.0新增）
# ============================================================================

#' Explore QC Distribution Before Filtering
#'
#' @param seurat_obj Seurat object to explore
#' @param groups.by Column name for group analysis (default: NULL)
#' @param organism Species for automatic mitochondrial calculation (default: "human")
#' @param report.dir Directory to save exploration report (default: "./qc_explore")
#' @param verbose Print messages (default: TRUE)
#'
#' @return S3 object of class "scQCExplore" containing exploration results
#' @export

scQCExplore <- function(
  seurat_obj,
  groups.by = NULL,
  organism = "human",
  report.dir = "./qc_explore",
  verbose = TRUE
) {

  if (verbose) cat("========== scQCExplore: Starting QC Exploration ==========\n")

  # Step 1: 检查并计算线粒体比例
  if (verbose) cat("Step 1: Checking mitochondrial percentage...\n")

  if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
    if (verbose) cat("  -> percent.mt not found in metadata.\n")
    if (verbose) cat("  -> Running addMitoRatio() with organism...\n")
    seurat_obj <- addMitoRatio(seurat_obj, organism = organism)
    if (verbose) cat("  -> Successfully added percent.mt\n")
  } else {
    if (verbose) cat("  -> percent.mt already exists in metadata. Skipping.\n")
  }

  # Step 2: 提取元数据
  if (verbose) cat("Step 2: Extracting metadata...\n")
  metadata <- seurat_obj@meta.data
  n_cells <- nrow(metadata)

  if (verbose) {
    cat(sprintf("  -> Total cells: %d\n", n_cells))
    if (!is.null(groups.by) && groups.by %in% colnames(metadata)) {
      cat(sprintf("  -> Groups: %s\n", paste(unique(metadata[[groups.by]]), collapse=", ")))
    }
  }

  # Step 3: 创建报告目录
  if (!base::dir.exists(report.dir)) {
    base::dir.create(report.dir, recursive = TRUE)
  }

  # Step 4: 计算百分位统计
  if (verbose) cat("Step 3: Calculating percentile statistics...\n")
  qc_stats <- .calculate_qc_percentiles(
    metadata = metadata,
    groups.by = if(!is.null(groups.by) && groups.by %in% colnames(metadata)) groups.by else NULL,
    percentiles = c(2.5, 5, 95, 97.5)
  )

  # Step 5: 计算细胞计数
  cell_counts <- .calculate_cell_counts(
    metadata = metadata,
    percentiles = c(2.5, 5, 95, 97.5),
    qc_metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt")
  )

  # Step 6: 生成图表
  if (verbose) cat("Step 4: Generating distribution plots...\n")

  # 检查是否有有效的groups.by
  has_groups <- !is.null(groups.by) && groups.by %in% colnames(metadata)

  plots <- list()

  # 整体分布图
  plots$overall <- .create_overall_qc_violin(
    metadata = metadata,
    percentiles = c(2.5, 5, 95, 97.5),
    params = list(percent_mt_max = 20, nCount_min = 500, nFeature_min = 200),
    report_dir = report.dir,
    has_groups = has_groups
  )

  # 分组分布图（如果有groups.by）
  if (has_groups) {
    plots$by_group <- .create_grouped_qc_violin(
      metadata = metadata,
      groups.by = groups.by,
      percentiles = c(2.5, 5, 95, 97.5),
      params = list(percent_mt_max = 20, nCount_min = 500, nFeature_min = 200),
      report_dir = report.dir
    )
  }

  # Step 7: 生成阈值建议
  if (verbose) cat("Step 5: Generating threshold suggestions...\n")
  suggestions <- .generate_threshold_suggestions(qc_stats, cell_counts)

  # Step 8: 生成HTML报告
  if (verbose) cat("Step 6: Generating exploration report...\n")
  report_file <- .generate_explore_html(
    qc_stats = qc_stats,
    cell_counts = cell_counts,
    suggestions = suggestions,
    report_dir = report.dir,
    groups.by = if(has_groups) groups.by else NULL,
    n_cells = n_cells
  )

  if (verbose) cat(sprintf("  -> Report saved to: %s\n", report_file))

  # 组织返回值
  output <- structure(
    list(
      n_cells = n_cells,
      stats = qc_stats,
      cell_counts = cell_counts,
      suggestions = suggestions,
      plots = plots,
      report_file = report_file,
      groups.by = if(has_groups) groups.by else NULL
    ),
    class = "scQCExplore"
  )

  if (verbose) cat("========== QC Exploration Complete ==========\n\n")

  return(output)
}


# ============================================================================
# 14.1 生成阈值建议（v2.1.1优化）
# ============================================================================

.generate_threshold_suggestions <- function(qc_stats, cell_counts) {

  suggestions <- list()

  # 辅助函数：将值取整到最近的50的倍数（范围150-1000）- 用于nFeature下限
  round_to_50 <- function(x) {
    rounded <- round(x / 50) * 50
    rounded <- pmax(150, pmin(1000, rounded))
    return(rounded)
  }

  # 辅助函数：将值取整到最近的100的倍数 - 用于nFeature上限
  round_to_100 <- function(x) {
    rounded <- round(x / 100) * 100
    rounded <- pmax(100, rounded)
    return(rounded)
  }

  # 辅助函数：将值取整到最近的1000的倍数 - 用于nCount
  round_to_1000 <- function(x) {
    rounded <- round(x / 1000) * 1000
    rounded <- pmax(1000, rounded)
    return(rounded)
  }

  # 辅助函数：生成范围字符串
  make_range <- function(low, high) {
    return(sprintf("%d-%d", low, high))
  }

  # 判断是否超过阈值
  # Delta Value > 30% 或 Delta Cells > 120%
  check_threshold_exceeded <- function(delta_value_pct, delta_cells_pct) {
    return(delta_value_pct > 30 || delta_cells_pct > 120)
  }

  # ========== 下限建议 ==========
  # 检查nFeature和nCount的变化幅度
  nfeature_exceeded <- FALSE
  ncount_exceeded <- FALSE

  if (!is.null(cell_counts$low_threshold$nFeature)) {
    raw <- cell_counts$low_threshold$nFeature$raw
    nfeature_exceeded <- check_threshold_exceeded(raw$delta_value_pct, raw$delta_cells_pct)
  }

  if (!is.null(cell_counts$low_threshold$nCount)) {
    raw <- cell_counts$low_threshold$nCount$raw
    ncount_exceeded <- check_threshold_exceeded(raw$delta_value_pct, raw$delta_cells_pct)
  }

  # nFeature 下限建议
  # 注意：下限阈值中，moderate（宽松）值更小，conservative（严格）值更大
  # 范围格式应该是 "小-大"，即 "moderate-conservative"
  if (!is.null(cell_counts$low_threshold$nFeature)) {
    raw <- cell_counts$low_threshold$nFeature$raw

    if (!nfeature_exceeded) {
      # 变化平缓：宽松用2.5%，严格用5%
      moderate_low <- round_to_50(raw$q_2.5)
      conservative_low <- round_to_50(raw$q_5)
    } else {
      # 变化剧烈：宽松用2.5%，严格缩减到30%水平
      moderate_low <- round_to_50(raw$q_2.5)
      # 计算缩减后的值
      reduced_value <- raw$q_2.5 + (raw$q_5 - raw$q_2.5) * 0.3
      conservative_low <- round_to_50(reduced_value)
    }

    # 范围格式：小-大（moderate-conservative）
    suggestions$nFeature$min <- make_range(moderate_low, conservative_low)
  }

  # nCount 下限建议（取整到50的倍数，与nFeature一致）
  if (!is.null(cell_counts$low_threshold$nCount)) {
    raw <- cell_counts$low_threshold$nCount$raw

    if (!ncount_exceeded) {
      # 变化平缓：宽松用2.5%，严格用5%
      moderate_low <- round_to_50(raw$q_2.5)
      conservative_low <- round_to_50(raw$q_5)
    } else {
      # 变化剧烈：宽松用2.5%，严格缩减到30%水平
      moderate_low <- round_to_50(raw$q_2.5)
      reduced_value <- raw$q_2.5 + (raw$q_5 - raw$q_2.5) * 0.3
      conservative_low <- round_to_50(reduced_value)
    }

    # 范围格式：小-大（moderate-conservative）
    suggestions$nCount$min <- make_range(moderate_low, conservative_low)
  }

  # ========== 上限建议 ==========
  # 检查nFeature和nCount的变化幅度（上限）
  nfeature_high_exceeded <- FALSE
  ncount_high_exceeded <- FALSE

  if (!is.null(cell_counts$high_threshold$nFeature)) {
    raw <- cell_counts$high_threshold$nFeature$raw
    nfeature_high_exceeded <- check_threshold_exceeded(raw$delta_value_pct, raw$delta_cells_pct)
  }

  if (!is.null(cell_counts$high_threshold$nCount)) {
    raw <- cell_counts$high_threshold$nCount$raw
    ncount_high_exceeded <- check_threshold_exceeded(raw$delta_value_pct, raw$delta_cells_pct)
  }

  # nFeature 上限建议
  if (!is.null(cell_counts$high_threshold$nFeature)) {
    raw <- cell_counts$high_threshold$nFeature$raw

    if (!nfeature_high_exceeded) {
      # 变化平缓：宽松用97.5%，严格用95%
      moderate_high <- round_to_100(raw$q_97.5)
      conservative_high <- round_to_100(raw$q_95)
    } else {
      # 变化剧烈：宽松用97.5%，严格缩减到30%水平
      moderate_high <- round_to_100(raw$q_97.5)
      reduced_value <- raw$q_95 + (raw$q_97.5 - raw$q_95) * 0.3
      conservative_high <- round_to_100(reduced_value)
    }

    suggestions$nFeature$max <- make_range(conservative_high, moderate_high)
  }

  # nCount 上限建议
  if (!is.null(cell_counts$high_threshold$nCount)) {
    raw <- cell_counts$high_threshold$nCount$raw

    if (!ncount_high_exceeded) {
      # 变化平缓：宽松用97.5%，严格用95%
      moderate_high <- round_to_1000(raw$q_97.5)
      conservative_high <- round_to_1000(raw$q_95)
    } else {
      # 变化剧烈：宽松用97.5%，严格缩减到30%水平
      moderate_high <- round_to_1000(raw$q_97.5)
      reduced_value <- raw$q_95 + (raw$q_97.5 - raw$q_95) * 0.3
      conservative_high <- round_to_1000(reduced_value)
    }

    suggestions$nCount$max <- make_range(conservative_high, moderate_high)
  }

  # percent.mt 建议：保守10%，适中20%
  suggestions$percent_mt <- list(
    conservative = 10,
    moderate = 20
  )

  return(suggestions)
}


# ============================================================================
# 14.2 生成探索HTML报告
# ============================================================================

.generate_explore_html <- function(qc_stats, cell_counts, suggestions,
                                    report_dir, groups.by, n_cells) {

  report_file <- base::file.path(report_dir, "scQCExplore_report.html")

  html_content <- '
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>scQCExplore Report</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #28a745 0%, #20c997 100%);
            padding: 20px;
            min-height: 100vh;
        }
        .container {
            max-width: 1000px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            overflow: hidden;
        }
        .header {
            background: linear-gradient(135deg, #28a745 0%, #20c997 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }
        .header h1 { font-size: 2em; margin-bottom: 10px; }
        .content { padding: 30px; }
        .section { margin-bottom: 30px; }
        .section h2 {
            color: #28a745;
            border-bottom: 3px solid #28a745;
            padding-bottom: 10px;
            margin-bottom: 15px;
        }
        .figure { text-align: center; margin: 20px 0; }
        .figure img { max-width: 100%; border-radius: 8px; }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
            background: white;
        }
        th { background: #f8f9fa; padding: 12px; text-align: left; font-weight: 600; }
        td { padding: 10px 12px; border-bottom: 1px solid #dee2e6; }
        .suggestion-box {
            background: #f8f9fa;
            border-left: 4px solid #28a745;
            padding: 15px;
            margin: 15px 0;
            border-radius: 6px;
        }
        .suggestion-box h4 { color: #28a745; margin-bottom: 10px; }
        code {
            background: #e9ecef;
            padding: 2px 6px;
            border-radius: 4px;
            font-family: Consolas, monospace;
        }
        .footer {
            background: #f8f9fa;
            padding: 15px 30px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>📊 scQCExplore Report</h1>
            <p>QC Distribution Exploration & Threshold Suggestions</p>
        </div>
        <div class="content">
'

  # 添加基本信息（合并为一行）
  html_content <- paste0(html_content, sprintf('
            <div class="section">
                <h2>Overview</h2>
                <p><strong>Total Cells:</strong> %d%s</p>
            </div>
',
    n_cells,
    if(!is.null(groups.by)) sprintf(" | <strong>Grouped by:</strong> %s", groups.by) else ""
  ))

  # 添加分布图
  html_content <- paste0(html_content, '
            <div class="section">
                <h2>QC Distribution</h2>
                <div class="figure">
                    <img src="00_qc_distribution_overall.png" alt="QC Distribution">
                </div>
')

  # 如果有分组图
  if (!is.null(groups.by) && base::file.exists(base::file.path(report_dir, "00_qc_distribution_by_group.png"))) {
    html_content <- paste0(html_content, sprintf('
                <h3 style="color: #28a745; margin-top: 20px;">Distribution by %s</h3>
                <div class="figure">
                    <img src="00_qc_distribution_by_group.png" alt="QC Distribution by Group">
                </div>
', groups.by))
  }

  html_content <- paste0(html_content, '
            </div>
')

  # 添加统计表格 - 下限表格
  html_content <- paste0(html_content, '
            <div class="section">
                <h2>Low Threshold Assessment (min)</h2>
                <table style="width: 100%%; max-width: 800px;">
                    <thead>
                        <tr>
                            <th>Items</th>
                            <th>2.5%%~5%%</th>
                            <th>Delta Value</th>
                            <th>Cells to Remove</th>
                            <th>Delta Cells</th>
                        </tr>
                    </thead>
                    <tbody>
')

  # nFeature 行
  if (!is.null(cell_counts$low_threshold$nFeature)) {
    d <- cell_counts$low_threshold$nFeature
    html_content <- paste0(html_content, sprintf('
                        <tr>
                            <td>nFeature</td>
                            <td>%s</td>
                            <td>%s</td>
                            <td>%s</td>
                            <td>%s</td>
                        </tr>
',
      d$range_2.5_5,
      d$delta_value,
      d$cells_remove,
      d$delta_cells
    ))
  }

  # nCount 行
  if (!is.null(cell_counts$low_threshold$nCount)) {
    d <- cell_counts$low_threshold$nCount
    html_content <- paste0(html_content, sprintf('
                        <tr>
                            <td>nCount</td>
                            <td>%s</td>
                            <td>%s</td>
                            <td>%s</td>
                            <td>%s</td>
                        </tr>
',
      d$range_2.5_5,
      d$delta_value,
      d$cells_remove,
      d$delta_cells
    ))
  }

  html_content <- paste0(html_content, '
                    </tbody>
                </table>
            </div>
')

  # 上限表格
  html_content <- paste0(html_content, '
            <div class="section">
                <h2>High Threshold Assessment (max)</h2>
                <table style="width: 100%%; max-width: 800px;">
                    <thead>
                        <tr>
                            <th>Items</th>
                            <th>95%%~97.5%%</th>
                            <th>Delta Value</th>
                            <th>Cells to Remove</th>
                            <th>Delta Cells</th>
                        </tr>
                    </thead>
                    <tbody>
')

  # nFeature 行
  if (!is.null(cell_counts$high_threshold$nFeature)) {
    d <- cell_counts$high_threshold$nFeature
    html_content <- paste0(html_content, sprintf('
                        <tr>
                            <td>nFeature</td>
                            <td>%s</td>
                            <td>%s</td>
                            <td>%s</td>
                            <td>%s</td>
                        </tr>
',
      d$range_95_97.5,
      d$delta_value,
      d$cells_remove,
      d$delta_cells
    ))
  }

  # nCount 行
  if (!is.null(cell_counts$high_threshold$nCount)) {
    d <- cell_counts$high_threshold$nCount
    html_content <- paste0(html_content, sprintf('
                        <tr>
                            <td>nCount</td>
                            <td>%s</td>
                            <td>%s</td>
                            <td>%s</td>
                            <td>%s</td>
                        </tr>
',
      d$range_95_97.5,
      d$delta_value,
      d$cells_remove,
      d$delta_cells
    ))
  }

  html_content <- paste0(html_content, '
                    </tbody>
                </table>
            </div>
')

  # 线粒体文字说明
  if (!is.null(cell_counts$mt_summary)) {
    mt <- cell_counts$mt_summary
    html_content <- paste0(html_content, sprintf('
            <div class="section">
                <h2>Mitochondrial Percentage</h2>
                <p style="color: #333; line-height: 1.8;">
                    <strong>Cells with high mitochondrial content:</strong><br>
                    &gt;5%%: %d (%.1f%%) | &gt;10%%: %d (%.1f%%) | &gt;15%%: %d (%.1f%%) | &gt;20%%: %d (%.1f%%)
                </p>
            </div>
',
      mt$gt_5$count, mt$gt_5$pct,
      mt$gt_10$count, mt$gt_10$pct,
      mt$gt_15$count, mt$gt_15$pct,
      mt$gt_20$count, mt$gt_20$pct
    ))
  }

  # 添加阈值建议 - 简洁格式
  html_content <- paste0(html_content, '
            <div class="section">
                <h2>Threshold Suggestions</h2>
')

  # nFeature 建议
  if (!is.null(suggestions$nFeature)) {
    min_range <- suggestions$nFeature$min
    max_range <- if (!is.null(suggestions$nFeature$max)) suggestions$nFeature$max else "Inf"
    html_content <- paste0(html_content, sprintf('
                <div class="suggestion-box">
                    <h4>nFeature (Genes)</h4>
                    <p><strong>Moderate:</strong> <code>nFeature_min = %s</code></p>
                    <p><strong>Conservative:</strong> <code>nFeature_min = %s, nFeature_max = %s</code></p>
                </div>
',
      min_range,
      strsplit(min_range, "-")[[1]][1],  # 取范围的下限
      max_range
    ))
  }

  # nCount 建议
  if (!is.null(suggestions$nCount)) {
    min_range <- suggestions$nCount$min
    max_range <- if (!is.null(suggestions$nCount$max)) suggestions$nCount$max else "Inf"
    html_content <- paste0(html_content, sprintf('
                <div class="suggestion-box">
                    <h4>nCount (UMI)</h4>
                    <p><strong>Moderate:</strong> <code>nCount_min = %s</code></p>
                    <p><strong>Conservative:</strong> <code>nCount_min = %s, nCount_max = %s</code></p>
                </div>
',
      min_range,
      strsplit(min_range, "-")[[1]][1],  # 取范围的下限
      max_range
    ))
  }

  # percent.mt 建议：保守10%，适中20%
  if (!is.null(suggestions$percent_mt)) {
    html_content <- paste0(html_content, sprintf('
                <div class="suggestion-box">
                    <h4>Mitochondrial %%</h4>
                    <p><strong>Moderate:</strong> <code>percent_mt_max = %d</code></p>
                    <p><strong>Conservative:</strong> <code>percent_mt_max = %d</code></p>
                </div>
',
      suggestions$percent_mt$moderate,
      suggestions$percent_mt$conservative
    ))
  }

  # 使用示例 + 双语提醒
  # 提取建议值用于示例代码
  nfeature_min <- if (!is.null(suggestions$nFeature$min)) {
    parts <- strsplit(suggestions$nFeature$min, "-")[[1]]
    parts[length(parts)]  # 取范围的上限作为示例
  } else "200"

  ncount_min <- if (!is.null(suggestions$nCount$min)) {
    parts <- strsplit(suggestions$nCount$min, "-")[[1]]
    parts[length(parts)]
  } else "500"

  html_content <- paste0(html_content, sprintf('
                <div class="suggestion-box" style="border-left-color: #667eea;">
                    <h4 style="color: #667eea;">Example Usage</h4>
                    <pre style="background: #f1f3f5; padding: 10px; border-radius: 4px; margin-top: 10px; overflow-x: auto;"><code>result <- scQCFilter(
  seurat_obj,
  nCount_min = %s,
  nFeature_min = %s,
  percent_mt_max = 20,
  groups.by = "your_group_column"
)</code></pre>
                    <div style="margin-top: 15px; padding: 10px; background: #fff3cd; border-radius: 4px; border-left: 3px solid #ffc107;">
                        <strong>⚠️ Reminder / 提醒:</strong><br>
                        <span style="color: #856404;">The above suggestions are AI-generated based on data distribution.
                        Please verify and adjust the QC thresholds according to your specific biological context.</span><br>
                        <span style="color: #856404;">以上建议是基于数据分布的智能化推荐，请根据您的具体生物学背景验证并调整质控阈值。</span><br>
                        <span style="color: #28a745; font-weight: bold; margin-top: 8px; display: inline-block;">Thank you for using scQCFilter! 感谢您使用 scQCFilter！</span>
                    </div>
                </div>
            </div>
',
    ncount_min,
    nfeature_min
  ))

  # 结束HTML
  html_content <- paste0(html_content, '
        </div>
        <div class="footer">
            <p>Generated by <strong>scQCExplore</strong> v2.1.0</p>
            <p>Report generated on ', base::Sys.time(), '</p>
        </div>
    </div>
</body>
</html>
')

  # 写入文件
  base::writeLines(html_content, report_file)

  return(report_file)
}


# ============================================================================
# 14.3 S3 方法：print.scQCExplore
# ============================================================================

#' @export
print.scQCExplore <- function(x, ...) {

  cat("\n")
  cat("═══════════════════════════════════════════════════════\n")
  cat("                    scQCExplore Results                \n")
  cat("═══════════════════════════════════════════════════════\n\n")

  cat("Data Overview:\n")
  cat("─────────────────────────────────────────────────────\n")
  cat(sprintf("  Total cells: %d\n", x$n_cells))
  if (!is.null(x$groups.by)) {
    cat(sprintf("  Grouped by: %s\n", x$groups.by))
  }
  cat("\n")

  cat("Suggested Thresholds (Moderate):\n")
  cat("─────────────────────────────────────────────────────\n")
  if (!is.null(x$suggestions$nCount)) {
    cat(sprintf("  nCount_min: %d\n", x$suggestions$nCount$moderate$min))
  }
  if (!is.null(x$suggestions$nFeature)) {
    cat(sprintf("  nFeature_min: %d\n", x$suggestions$nFeature$moderate$min))
  }
  if (!is.null(x$suggestions$percent_mt)) {
    cat(sprintf("  percent_mt_max: %d\n", x$suggestions$percent_mt$moderate$max))
  }
  cat("\n")

  cat("Report:\n")
  cat("─────────────────────────────────────────────────────\n")
  cat(sprintf("  HTML report: %s\n", x$report_file))
  cat("\n")

  cat("═══════════════════════════════════════════════════════\n\n")

}


# ============================================================================
# 15. scQCReport - 生成组合PDF报告（v2.1.1新增）
# ============================================================================

#' Generate Combined PDF Report with scQCExplore and scQCFilter
#'
#' @param seurat_obj Seurat object to analyze
#' @param groups.by Column name for group analysis (default: NULL)
#' @param sample.by Column name for sample analysis (default: "orig.ident")
#' @param organism Species for mitochondrial calculation (default: "human")
#' @param percent_mt_max Maximum mitochondrial percentage (default: 20)
#' @param nCount_min Minimum UMI count (default: 500)
#' @param nCount_max Maximum UMI count (default: Inf)
#' @param nFeature_min Minimum gene count (default: 200)
#' @param nFeature_max Maximum gene count (default: Inf)
#' @param report.dir Directory to save reports (default: "./qc_report")
#' @param use.suggestions Use scQCExplore suggestions for scQCFilter (default: TRUE)
#' @param verbose Print messages (default: TRUE)
#'
#' @return S3 object of class "scQCReport" containing both results and PDF path
#' @export

scQCReport <- function(
  seurat_obj,
  groups.by = NULL,
  sample.by = "orig.ident",
  organism = "human",
  percent_mt_max = 20,
  nCount_min = 500,
  nCount_max = Inf,
  nFeature_min = 200,
  nFeature_max = Inf,
  report.dir = "./qc_report",
  use.suggestions = TRUE,
  verbose = TRUE
) {

  if (verbose) cat("========== scQCReport: Generating Combined PDF Report ==========\n")

  # Step 1: 创建报告目录
  if (!base::dir.exists(report.dir)) {
    base::dir.create(report.dir, recursive = TRUE)
  }

  # Step 2: 运行 scQCExplore
  if (verbose) cat("Step 1: Running scQCExplore...\n")
  explore_result <- scQCExplore(
    seurat_obj = seurat_obj,
    groups.by = groups.by,
    organism = organism,
    report.dir = report.dir,
    verbose = verbose
  )

  # Step 3: 确定过滤参数
  if (use.suggestions && !is.null(explore_result$suggestions)) {
    # 使用建议值（取范围的上限作为宽松值）
    if (!is.null(explore_result$suggestions$nCount$min)) {
      parts <- strsplit(explore_result$suggestions$nCount$min, "-")[[1]]
      nCount_min <- as.numeric(parts[length(parts)])
    }
    if (!is.null(explore_result$suggestions$nFeature$min)) {
      parts <- strsplit(explore_result$suggestions$nFeature$min, "-")[[1]]
      nFeature_min <- as.numeric(parts[length(parts)])
    }
    if (!is.null(explore_result$suggestions$percent_mt$moderate)) {
      percent_mt_max <- explore_result$suggestions$percent_mt$moderate
    }
    if (verbose) cat(sprintf("  -> Using suggested parameters: nCount_min=%d, nFeature_min=%d, percent_mt_max=%d\n",
                             nCount_min, nFeature_min, percent_mt_max))
  }

  # Step 4: 运行 scQCFilter
  if (verbose) cat("Step 2: Running scQCFilter...\n")
  filter_result <- scQCFilter(
    seurat_obj = seurat_obj,
    percent_mt_max = percent_mt_max,
    nCount_min = nCount_min,
    nCount_max = nCount_max,
    nFeature_min = nFeature_min,
    nFeature_max = nFeature_max,
    groups.by = groups.by,
    sample.by = sample.by,
    organism = organism,
    return.filtered = FALSE,
    report.dir = report.dir,
    plot = TRUE,
    show.qc.distribution = FALSE,
    verbose = verbose
  )

  # Step 5: 生成组合PDF
  if (verbose) cat("Step 3: Generating combined PDF report...\n")
  pdf_file <- .generate_combined_pdf(
    explore_result = explore_result,
    filter_result = filter_result,
    report.dir = report.dir,
    groups.by = groups.by,
    params = list(
      percent_mt_max = percent_mt_max,
      nCount_min = nCount_min,
      nCount_max = nCount_max,
      nFeature_min = nFeature_min,
      nFeature_max = nFeature_max
    )
  )

  if (verbose) cat(sprintf("  -> PDF saved to: %s\n", pdf_file))
  if (verbose) cat("========== Combined PDF Report Complete ==========\n\n")

  # 组织返回值
  output <- structure(
    list(
      explore = explore_result,
      filter = filter_result,
      pdf_file = pdf_file,
      parameters = list(
        percent_mt_max = percent_mt_max,
        nCount_min = nCount_min,
        nCount_max = nCount_max,
        nFeature_min = nFeature_min,
        nFeature_max = nFeature_max,
        groups.by = groups.by,
        sample.by = sample.by,
        use.suggestions = use.suggestions
      )
    ),
    class = "scQCReport"
  )

  return(output)
}


# ============================================================================
# 15.1 生成组合PDF内部函数
# ============================================================================

.generate_combined_pdf <- function(explore_result, filter_result, report.dir,
                                   groups.by, params) {

  pdf_file <- base::file.path(report.dir, "scQCReport_combined.pdf")

  # 使用grDevices生成PDF
  grDevices::pdf(pdf_file, width = 11, height = 8.5)

  # ========== 第一页: scQCExplore ==========
  .draw_explore_page(explore_result, groups.by)

  # ========== 第二页: scQCFilter ==========
  .draw_filter_page(filter_result, params, groups.by)

  grDevices::dev.off()

  return(pdf_file)
}


# ============================================================================
# 15.2 绘制Explore页面
# ============================================================================

.draw_explore_page <- function(explore_result, groups.by) {

  # 设置页面布局
  graphics::par(mar = c(0.5, 0.5, 1, 0.5), bg = "white")

  # 标题区域
  graphics::plot.new()
  graphics::rect(0, 0.85, 1, 1, col = "#28a745", border = NA)
  graphics::text(0.5, 0.925, "scQCExplore Report", cex = 2, font = 2, col = "white")
  graphics::text(0.5, 0.87, "QC Distribution Exploration & Threshold Suggestions", cex = 1, col = "white")

  # 基本信息区域
  graphics::rect(0.02, 0.78, 0.98, 0.84, col = "#f8f9fa", border = NA, lty = 0)
  info_text <- sprintf("Total Cells: %d", explore_result$n_cells)
  if (!is.null(groups.by)) {
    info_text <- paste0(info_text, sprintf("  |  Grouped by: %s", groups.by))
  }
  graphics::text(0.5, 0.81, info_text, cex = 1.1, font = 1)

  # 加载并显示分布图
  plot_file <- base::file.path(dirname(explore_result$report_file), "00_qc_distribution_overall.png")
  if (base::file.exists(plot_file)) {
    img <- png::readPNG(plot_file)
    # 显示图片在中央区域
    graphics::par(fig = c(0.05, 0.95, 0.35, 0.76), new = TRUE, mar = c(0, 0, 0.5, 0))
    graphics::plot.new()
    graphics::rasterImage(img, 0, 0, 1, 1)
  }

  # 恢复主布局
  graphics::par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0.5, 0.5, 1, 0.5))

  # 统计表格区域
  graphics::rect(0.02, 0.02, 0.98, 0.32, col = "#f8f9fa", border = NA)

  # 表格标题
  graphics::text(0.08, 0.29, "Low Threshold Assessment (min)", cex = 0.9, font = 2, adj = 0)

  # 绘制表格
  tbl_x <- 0.04
  tbl_y <- 0.26
  row_h <- 0.035
  col_widths <- c(0.12, 0.16, 0.15, 0.18, 0.15)

  # 表头
  headers <- c("Items", "2.5%~5%", "Delta Value", "Cells to Remove", "Delta Cells")
  graphics::rect(tbl_x, tbl_y, tbl_x + sum(col_widths), tbl_y + row_h, col = "#e9ecef", border = "gray80")
  x_pos <- tbl_x
  for (i in seq_along(headers)) {
    graphics::text(x_pos + col_widths[i]/2, tbl_y + row_h/2, headers[i], cex = 0.75, font = 2)
    x_pos <- x_pos + col_widths[i]
  }

  # 数据行
  tbl_y <- tbl_y - row_h
  cc <- explore_result$cell_counts

  # nFeature行
  if (!is.null(cc$low_threshold$nFeature)) {
    d <- cc$low_threshold$nFeature
    graphics::rect(tbl_x, tbl_y, tbl_x + sum(col_widths), tbl_y + row_h, col = "white", border = "gray80")
    x_pos <- tbl_x
    values <- c("nFeature", d$range_2.5_5, d$delta_value, d$cells_remove, d$delta_cells)
    for (i in seq_along(values)) {
      graphics::text(x_pos + col_widths[i]/2, tbl_y + row_h/2, values[i], cex = 0.7)
      x_pos <- x_pos + col_widths[i]
    }
  }

  # nCount行
  tbl_y <- tbl_y - row_h
  if (!is.null(cc$low_threshold$nCount)) {
    d <- cc$low_threshold$nCount
    graphics::rect(tbl_x, tbl_y, tbl_x + sum(col_widths), tbl_y + row_h, col = "white", border = "gray80")
    x_pos <- tbl_x
    values <- c("nCount", d$range_2.5_5, d$delta_value, d$cells_remove, d$delta_cells)
    for (i in seq_along(values)) {
      graphics::text(x_pos + col_widths[i]/2, tbl_y + row_h/2, values[i], cex = 0.7)
      x_pos <- x_pos + col_widths[i]
    }
  }

  # 阈值建议
  sug_x <- 0.55
  sug_y <- 0.26
  graphics::text(sug_x, sug_y, "Threshold Suggestions", cex = 0.9, font = 2, adj = 0)

  suggestions <- explore_result$suggestions
  sug_y <- sug_y - 0.04

  # nFeature建议
  if (!is.null(suggestions$nFeature$min)) {
    parts <- strsplit(suggestions$nFeature$min, "-")[[1]]
    conservative <- parts[1]
    moderate <- parts[2]
    graphics::text(sug_x, sug_y, sprintf("nFeature: Moderate = %s | Conservative = %s", moderate, conservative),
                   cex = 0.7, adj = 0, col = "#333")
    sug_y <- sug_y - 0.035
  }

  # nCount建议
  if (!is.null(suggestions$nCount$min)) {
    parts <- strsplit(suggestions$nCount$min, "-")[[1]]
    conservative <- parts[1]
    moderate <- parts[2]
    graphics::text(sug_x, sug_y, sprintf("nCount: Moderate = %s | Conservative = %s", moderate, conservative),
                   cex = 0.7, adj = 0, col = "#333")
    sug_y <- sug_y - 0.035
  }

  # percent.mt建议
  if (!is.null(suggestions$percent_mt)) {
    graphics::text(sug_x, sug_y, sprintf("percent.mt: Moderate = %d%% | Conservative = %d%%",
                                         suggestions$percent_mt$moderate, suggestions$percent_mt$conservative),
                   cex = 0.7, adj = 0, col = "#333")
    sug_y <- sug_y - 0.035
  }

  # 线粒体统计
  if (!is.null(cc$mt_summary)) {
    mt <- cc$mt_summary
    sug_y <- sug_y - 0.02
    graphics::text(sug_x, sug_y, "Mitochondrial Content:", cex = 0.75, font = 2, adj = 0)
    sug_y <- sug_y - 0.03
    mt_text <- sprintf(">5%%: %d (%.1f%%)  |  >10%%: %d (%.1f%%)  |  >15%%: %d (%.1f%%)  |  >20%%: %d (%.1f%%)",
                       mt$gt_5$count, mt$gt_5$pct, mt$gt_10$count, mt$gt_10$pct,
                       mt$gt_15$count, mt$gt_15$pct, mt$gt_20$count, mt$gt_20$pct)
    graphics::text(sug_x, sug_y, mt_text, cex = 0.65, adj = 0, col = "#666")
  }
}


# ============================================================================
# 15.3 绘制Filter页面
# ============================================================================

.draw_filter_page <- function(filter_result, params, groups.by) {

  # 获取主结果
  main_result <- filter_result$detailed_statistics[[1]]

  # 设置页面布局
  graphics::par(mar = c(0.5, 0.5, 1, 0.5), bg = "white")

  # 标题区域
  graphics::plot.new()
  graphics::rect(0, 0.85, 1, 1, col = "#667eea", border = NA)
  graphics::text(0.5, 0.925, "scQCFilter Report", cex = 2, font = 2, col = "white")
  graphics::text(0.5, 0.87, "Quality Control Filtering Results", cex = 1, col = "white")

  # 统计卡片区域
  card_y <- 0.72
  card_h <- 0.1
  card_widths <- c(0.22, 0.22, 0.22, 0.22)
  card_x <- 0.03

  cards <- list(
    list(label = "Total Cells", value = main_result$n_cells_before),
    list(label = "Cells Retained", value = main_result$n_cells_after),
    list(label = "Cells Filtered", value = main_result$n_filtered),
    list(label = "Keep Rate", value = sprintf("%.1f%%", main_result$keep_rate * 100))
  )

  for (i in seq_along(cards)) {
    graphics::rect(card_x, card_y, card_x + card_widths[i], card_y + card_h,
                   col = "#667eea", border = NA, radius = 0.02)
    graphics::text(card_x + card_widths[i]/2, card_y + card_h - 0.025,
                   cards[[i]]$label, cex = 0.7, col = "white")
    graphics::text(card_x + card_widths[i]/2, card_y + card_h/2 - 0.01,
                   cards[[i]]$value, cex = 1.3, font = 2, col = "white")
    card_x <- card_x + card_widths[i] + 0.02
  }

  # 过滤参数区域
  graphics::rect(0.02, 0.56, 0.48, 0.70, col = "#f8f9fa", border = NA)
  graphics::text(0.04, 0.68, "Filtering Parameters", cex = 0.9, font = 2, adj = 0)

  param_text <- sprintf(
    "percent_mt_max: %s\nnCount: %s ~ %s\nnFeature: %s ~ %s",
    params$percent_mt_max,
    if(is.finite(params$nCount_min)) params$nCount_min else "-Inf",
    if(is.finite(params$nCount_max)) params$nCount_max else "Inf",
    if(is.finite(params$nFeature_min)) params$nFeature_min else "-Inf",
    if(is.finite(params$nFeature_max)) params$nFeature_max else "Inf"
  )
  graphics::text(0.04, 0.60, param_text, cex = 0.75, adj = c(0, 1), col = "#333")

  # 诊断建议区域
  if (!is.null(filter_result$diagnostics)) {
    graphics::rect(0.50, 0.56, 0.98, 0.70, col = "#fff3cd", border = NA)
    graphics::text(0.52, 0.68, "Diagnostics", cex = 0.9, font = 2, adj = 0, col = "#856404")

    diag_y <- 0.65
    if (!is.null(filter_result$diagnostics$suggestions)) {
      for (i in seq_along(filter_result$diagnostics$suggestions)) {
        if (i <= 3) {  # 最多显示3条建议
          graphics::text(0.52, diag_y, paste0("- ", filter_result$diagnostics$suggestions[i]),
                        cex = 0.65, adj = c(0, 1), col = "#856404")
          diag_y <- diag_y - 0.03
        }
      }
    }
  }

  # 加载并显示过滤后的分布图
  filter_plot_file <- base::file.path(dirname(filter_result$report_file), "filtering_violin_plot.png")
  if (!base::file.exists(filter_plot_file)) {
    filter_plot_file <- base::file.path(dirname(filter_result$report_file), "00_qc_distribution_overall.png")
  }

  if (base::file.exists(filter_plot_file)) {
    img <- png::readPNG(filter_plot_file)
    graphics::par(fig = c(0.05, 0.95, 0.08, 0.54), new = TRUE, mar = c(0, 0, 0.5, 0))
    graphics::plot.new()
    graphics::rasterImage(img, 0, 0, 1, 1)
  }

  # 恢复主布局
  graphics::par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0.5, 0.5, 1, 0.5))

  # 页脚
  graphics::rect(0, 0, 1, 0.05, col = "#f8f9fa", border = NA)
  graphics::text(0.5, 0.025,
                 sprintf("Generated by scQCFilter v2.1.1 | %s", base::Sys.time()),
                 cex = 0.7, col = "#666")
}


# ============================================================================
# 15.4 S3 方法：print.scQCReport
# ============================================================================

#' @export
print.scQCReport <- function(x, ...) {

  cat("\n")
  cat("═══════════════════════════════════════════════════════\n")
  cat("                    scQCReport Results                \n")
  cat("═══════════════════════════════════════════════════════\n\n")

  cat("Combined PDF Report:\n")
  cat("─────────────────────────────────────────────────────\n")
  cat(sprintf("  PDF file: %s\n", x$pdf_file))
  cat("\n")

  cat("Parameters Used:\n")
  cat("─────────────────────────────────────────────────────\n")
  cat(sprintf("  percent_mt_max: %s\n", x$parameters$percent_mt_max))
  cat(sprintf("  nCount_min: %s\n", x$parameters$nCount_min))
  cat(sprintf("  nCount_max: %s\n", if(is.finite(x$parameters$nCount_max)) x$parameters$nCount_max else "Inf"))
  cat(sprintf("  nFeature_min: %s\n", x$parameters$nFeature_min))
  cat(sprintf("  nFeature_max: %s\n", if(is.finite(x$parameters$nFeature_max)) x$parameters$nFeature_max else "Inf"))
  if (!is.null(x$parameters$groups.by)) {
    cat(sprintf("  groups.by: %s\n", x$parameters$groups.by))
  }
  cat("\n")

  cat("═══════════════════════════════════════════════════════\n\n")

}
