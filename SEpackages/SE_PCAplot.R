SE_PCAplot = function(SE, assayname = "TPM", groupname = "group", outlier_threshold = 2, scale = TRUE){  
  message("可以通过 SE_subset <- SE[rownames(SE) %in% feature, ]筛选一个SE子集")  
  message("可以通过 SE_subset <- SE[,  colnames(SE) %in% sample]筛选一个SE子集")  
  library(ggplot2)  
  library(plyr)  
  library(gridExtra)  
  library(cluster)  

  SCvalue = function(data){ # 计算轮廓系数  
    num_clusters <- length(unique(data$group)) # 使用已知的组数作为聚类数目,  
    kmeans_result <- kmeans(data[, c("PC1", "PC2")], centers = num_clusters, nstart = 25)  # nstart 增加可以提高结果的稳定性  
    data$cluster <- as.factor(kmeans_result$cluster) # 将kmeans的结果加到pca_data_filtered中  
    distance_matrix <- dist(data[, c("PC1", "PC2")])  
    silhouette_result <- silhouette(as.numeric(data$cluster), distance_matrix)  
    mean_silhouette_width <- round(mean(silhouette_result[, "sil_width"]),3)  
    mean_silhouette_width  
  }  

  expdata = assay(SE, assayname) 
  row_variances <- apply(expdata, 1, var, na.rm = TRUE)  
  constant_rows <- which(row_variances == 0)  

  # 检查是否存在方差为 0 的特征
  if (length(constant_rows) > 0) {  
    print(paste0("Delete ", length(constant_rows), " common express features"))  
    expdata <- expdata[-constant_rows, ]  
  } else {  
    print("No common express features to delete.")  
  }  
  num_feature <- nrow(expdata) 

  sample_info = colData(SE)  
  pca_result <- prcomp(as.data.frame(t(expdata)), scale. = scale) 
  pca_data <- as.data.frame(pca_result$x)  

  missing_samples <- setdiff(rownames(pca_data), rownames(sample_info))  
  if (length(missing_samples) > 0) {  
    print(paste("Warning: Missing sample information for:", paste(missing_samples, collapse = ", ")))  

    pca_data <- pca_data[!(rownames(pca_data) %in% missing_samples), ]  
  }  

  pca_data$group <- mapvalues(rownames(pca_data), rownames(sample_info), sample_info[,groupname], warn_missing = FALSE) # 使用groupname参数  
 
  # 离群值过滤
  mean_pc1 <- mean(pca_data$PC1)  
  sd_pc1 <- sd(pca_data$PC1)  
  mean_pc2 <- mean(pca_data$PC2)  
  sd_pc2 <- sd(pca_data$PC2)  
  pca_data_filter <- pca_data[(pca_data$PC1 > (mean_pc1 - outlier_threshold * sd_pc1)) &(pca_data$PC1 < (mean_pc1 + outlier_threshold * sd_pc1)) &  
                                 (pca_data$PC2 > (mean_pc2 - outlier_threshold * sd_pc2)) & (pca_data$PC2 < (mean_pc2 + outlier_threshold * sd_pc2)),]  

  group_counts <- table(pca_data$group)  
  caption_text1 <- paste("feature:", num_feature,";","sample:", paste(names(group_counts), ":", group_counts, collapse = ", "))  

  group_counts_filtered <- table(pca_data_filter$group)  
  caption_text2 <- paste("feature:", num_feature,";","sample:", paste(names(group_counts_filtered), ":", group_counts_filtered, collapse = ", "))  

  plot_pca <- function(data, title, caption) {  
    ggplot(data, aes(x = PC1, y = PC2, color = group)) +  
      geom_point(size = 3) +  
      labs(title = paste0("PCAplot: ", SCvalue(data)), x = "PCA1", y = "PCA2", caption = caption) +  
      theme_minimal() +  
      stat_ellipse(type = "norm", level = 0.95, linetype = 2) +  
      theme(axis.title.x = element_text(size = 12),  
            axis.title.y = element_text(size = 12),  
            axis.text = element_text(size = 12),  
            plot.title = element_text(size = 12),  
            panel.grid.major = element_blank(),  
            panel.grid.minor = element_blank(),  
            strip.text = element_text(size = 12),  
            panel.border = element_rect(color = "gray", fill = NA, size = 1),  
            legend.position = c(0.9,0.2), legend.text = element_text(size = 12),  
            legend.title = element_text(size = 12))  
  }  

  pca_plot1 <- plot_pca(pca_data, "PCAplot", caption_text1)  
  pca_plot2 <- plot_pca(pca_data_filter, "PCAplot Filtered", caption_text2)  

  grid.arrange(pca_plot1, pca_plot2, nrow=1)  
  invisible(list(pca_plot1 = pca_plot1, pca_plot2 = pca_plot2)) 
}