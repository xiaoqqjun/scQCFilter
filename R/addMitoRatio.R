#' Calculate Mitochondrial Gene Percentage Using Seurat's PercentageFeatureSet
#'
#' Simple wrapper around Seurat::PercentageFeatureSet()
#' Automatically detects organism and applies correct pattern
#'
#' @param seurat_obj A Seurat object
#' @param organism Species: "human", "mouse", "rat", "macaque", "zebrafish"
#'                 (default: "human")
#'
#' @return Seurat object with percent.mt added to metadata
#'
#' @details
#' Uses Seurat::PercentageFeatureSet() with species-specific patterns:
#' - human: "^MT-"
#' - mouse: "^mt-"
#' - rat: "^mt-"
#' - macaque: "^MT-"
#' - zebrafish: "^mt-"
#'
#' @examples
#' \dontrun{
#' seurat_obj <- addMitoRatio(seurat_obj, organism = "human")
#' }
#'
#' @export
addMitoRatio <- function(seurat_obj, organism = "human") {

  # 标准化物种名称
  organism <- tolower(trimws(organism))

  # 物种映射 - 正则表达式模式
  organism_patterns <- list(
    human = "^MT-",
    homo = "^MT-",
    homo_sapiens = "^MT-",
    hs = "^MT-",
    mouse = "^mt-",
    mice = "^mt-",
    mmu = "^mt-",
    mus_musculus = "^mt-",
    mm = "^mt-",
    rat = "^mt-",
    rno = "^mt-",
    rattus_norvegicus = "^mt-",
    macaque = "^MT-",
    mfa = "^MT-",
    macaca_fascicularis = "^MT-",
    zebrafish = "^mt-",
    dre = "^mt-",
    danio_rerio = "^mt-"
  )

  # 查找对应的模式
  pattern <- organism_patterns[[organism]]
  
  if (is.null(pattern)) {
    stop(sprintf(
      "Unsupported organism: '%s'. Supported: human, mouse, rat, macaque, zebrafish",
      organism
    ))
  }

  # 使用 Seurat 的 PercentageFeatureSet 函数
  seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, 
                                              pattern = pattern,
                                              col.name = "percent.mt")

  return(seurat_obj)
}


#' Check Mitochondrial Status
#'
#' Quick check of mitochondrial percentage statistics
#'
#' @param seurat_obj A Seurat object
#'
#' @examples
#' \dontrun{
#' checkMitoStatus(seurat_obj)
#' }
#'
#' @export
checkMitoStatus <- function(seurat_obj) {
  if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {
    cat("Mitochondrial percentage statistics:\n")
    cat(sprintf("  Mean:   %.2f%%\n", mean(seurat_obj@meta.data$percent.mt)))
    cat(sprintf("  Median: %.2f%%\n", stats::median(seurat_obj@meta.data$percent.mt)))
    cat(sprintf("  Min:    %.2f%%\n", min(seurat_obj@meta.data$percent.mt)))
    cat(sprintf("  Max:    %.2f%%\n", max(seurat_obj@meta.data$percent.mt)))
    cat(sprintf("  SD:     %.2f%%\n", stats::sd(seurat_obj@meta.data$percent.mt)))
    invisible(NULL)
  } else {
    warning("'percent.mt' column not found in metadata")
    invisible(NULL)
  }
}
