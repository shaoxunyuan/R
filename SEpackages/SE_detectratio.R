SE_detectratio <- function(SE, assayname = "TPM") {  
    # 加载必要的库  
    library(SummarizedExperiment)  
    
    # 检查输入的有效性  
    if (!inherits(SE, "SummarizedExperiment")) {  
        stop("Input SE must be a SummarizedExperiment object.")  
    }  
    if (!assayname %in% names(assays(SE))) {  
        stop(paste("Assay", assayname, "not found in SE."))  
    }  
    
    # 获取 feature_info 和数据矩阵  
    feature_info <- rowData(SE)  
    expdata <- assay(SE, assayname)  
    
    # 计算每一行不为0的个数  
    detectsample_counts <- rowSums(expdata != 0)  
    
    # 更新 feature_info 中的 detectsample 列  
    feature_info$detectsample <- detectsample_counts  
    
    # 计算 detectratio  
    total_samples <- nrow(colData(SE))  
    feature_info$detectratio <- round(feature_info$detectsample / total_samples, 4)   
    
    # 使用更新后的 feature_info 更新 SE 的 rowData  
    rowData(SE) <- feature_info  
    
    # 返回SE  
    return(SE)  
}
