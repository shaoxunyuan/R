SE_DEseq2 <- function(SE, assayname = "Count", groupname = "group") {
  suppressPackageStartupMessages({
    library(DESeq2)
    library(SummarizedExperiment)
  })
    
  # 验证SE对象类型
  if (!is(SE, "SummarizedExperiment")) {
    stop("输入对象必须是SummarizedExperiment类型")
  }
  
  # 验证assay是否存在
  if (!assayname %in% assayNames(SE)) {
    stop("指定的assay不存在: ", assayname)
  }
  
  # 验证分组列是否存在
  if (!groupname %in% colnames(colData(SE))) {
    stop("分组列不存在于colData中: ", groupname)
  }

  # 提取Count矩阵并强制转为整数 (DESeq2要求整数)
  countData <- assay(SE, assayname)
  if (!all(countData == floor(countData))) {
    message("检测到非整数值，执行ceiling转换...")
    countData <- ceiling(countData)
  }

  # 提取样本信息并动态构建分组公式
  sample_info <- colData(SE)
  design_formula <- reformulate(groupname)  # 动态适配分组列名
    
  # --------------------------
  # DESeq2差异分析
  # --------------------------
  # 创建DESeq2对象
  dds <- DESeqDataSetFromMatrix(countData = countData,colData = sample_info,design = design_formula)
  
  # 运行分析并抑制冗余信息
  suppressMessages({dds <- DESeq(dds, quiet = TRUE)})
  
  # --------------------------
  # 结果处理
  # --------------------------
  # 获取分组信息并校验
    groups <- unique(sample_info$group)
    
    num_groups <- length(groups)
    
    # 存储所有差异分析结果的列表
    all_results <- list()
    # 进行所有可能的两两比较
    for (i in 1:(num_groups - 1)) {
      for (j in (i + 1):num_groups) {
        contrast <- c("group", groups[i], groups[j])
        res <- results(dds, contrast = contrast)
        # 找出log2FoldChange为NA的行索引
        na_rows <- which(is.na(res$pvalue))
        # 将这些行的指定列值进行修改
        res[na_rows, "log2FoldChange"] <- 0
        res[na_rows, "lfcSE"] <- 0
        res[na_rows, "stat"] <- 0
        res[na_rows, "pvalue"] <- 1
        res[na_rows, "padj"] <- 1
        res = res[order(res$pvalue,decreasing = F),]
        comparison_name <- paste0(groups[i], "_vs_", groups[j])
        all_results[[comparison_name]] <- res
      }
    }
 
# 添加差异分析结果到新对象的metadata中
metadata(SE)$DEresults <- all_results
 
# 返回新对象（原SE保持不变）
return(SE)
}