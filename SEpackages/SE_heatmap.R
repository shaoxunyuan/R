SE_heatmap <- function(SE, assayname = "TPM", group_col = NULL, normalization = "none", genes_of_interest = NULL, use_raster = TRUE) {
	#SE_heatmap(SE, assayname = "TPM", group_col = NULL, normalization = "scale", genes_of_interest = c("AAGAB","AARS"), use_raster = TRUE)  
    library(SummarizedExperiment)  
    library(ComplexHeatmap)  
    library(circlize)  
    library(RColorBrewer)  # 确保可以使用调色板  
    # 检查输入有效性  
    if (!inherits(SE, "SummarizedExperiment")) {  
        stop("Input SE must be a SummarizedExperiment object.")  
    }  
    
    # 提取表达数据  
    exp_data <- assay(SE, assayname)  
    
    # 检查基因是否存在于数据中  
    if (!is.null(genes_of_interest)) {  
        missing_genes <- setdiff(genes_of_interest, rownames(exp_data))  
        if (length(missing_genes) > 0) {  
            warning(paste("The following genes are not found in the dataset:", paste(missing_genes, collapse = ", ")))  
            genes_of_interest <- intersect(genes_of_interest, rownames(exp_data))  
        }  
        exp_data <- exp_data[genes_of_interest, , drop = FALSE]  
    }  
    
    # 数据归一化  
    if (normalization == "scale") {  
        exp_data <- t(apply(exp_data, 1, scale))  # z-score 标准化  
    } else if (normalization == "log") {  
        exp_data <- log1p(exp_data)  # log1p 转换  
    }  
    
    # 获取分组信息  
    if (!is.null(group_col) && group_col %in% names(colData(SE))) {  
        group_info <- colData(SE)[[group_col]]  
        unique_groups <- unique(group_info)  
        group_colors <- brewer.pal(length(unique_groups), "Set1")[1:length(unique_groups)]  
        names(group_colors) <- unique_groups  # 添加名称以匹配组  
        
        # 创建注释数据框  
        ha <- HeatmapAnnotation(df = data.frame(Group = factor(group_info)),   
                                 col = list(Group = group_colors))  
    } else {  
        ha <- NULL  # 没有分组信息时不添加注释  
    }  
    
    # 关闭消息提示  
    ht_opt$message <- FALSE  

    # 绘制热图  
    Heatmap(exp_data,   
        name = "Expression",   
        top_annotation = ha,  # 添加顶部的分组注释  
        show_column_names = TRUE,  
        cluster_rows = TRUE,  # 行聚类  
        cluster_columns = TRUE,  # 列聚类  
        show_row_names = TRUE,  # 显示行名  
        show_column_dend = TRUE,  # 显示列聚类树  
        show_row_dend = TRUE,  # 显示行聚类树  
        heatmap_legend_param = list(  
          title = "Expression",   
          title_position = "topcenter",   
          direction = "vertical",  # 垂直方向，让图例上下显示
          legend_width = unit(10, "cm"),  # 调整图例宽度  
          title_gp = gpar(fontsize = 10, fontface = "bold"),  # 图例标题字体设置  
          labels_gp = gpar(fontsize = 8)  # 标签字体设置  
        ),  
        use_raster = use_raster,  # 使用光栅化的选项  
        column_title = "samples",  
        row_title = "features")
}  

