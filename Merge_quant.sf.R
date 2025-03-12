# 读取命令行参数
# 该脚本将单个数据集的quant.sf结果整合为linear和circular两个水平的SummarizedExperiment格式
# linear包含isoform和gene水平，circular包含isoform，bsj和gene水平
args <- commandArgs(trailingOnly = TRUE)

# 获取传入的数据集名字如PRJNA429023
datasetname <- args[1]

options(warn = -1)
library(data.table)
library(plyr)
library(dplyr)
library(SummarizedExperiment)  

source("/data2/shaoxun/reference/functions.R") 
BloodCircR_circRNA = fread("/data2/shaoxun/BloodCircleR/dataset/circRNA_results/BloodCircR_circRNA.txt",data.table = F)

SEmake = function(data,sampletable){
# 样本信息
    sample_info = sampletable
    rownames(sample_info) = sample_info$BioSample
	overlapsample = intersect(unique(data$BioSample),unique(sample_info$BioSample))
	sample_info = sample_info[sample_info$BioSample %in% overlapsample,]
# 表达矩阵
	data = data[data$BioSample %in% overlapsample,]
    TPMmatrix <- dcast(data = data,formula = marker ~ BioSample,value.var = "TPM")
    rownames(TPMmatrix) = TPMmatrix$marker;TPMmatrix=TPMmatrix[,-1];TPMmatrix[is.na(TPMmatrix)] = 0
    Countmatrix <- dcast(data = data,formula = marker ~ BioSample,value.var = "NumReads")
    rownames(Countmatrix) = Countmatrix$marker;Countmatrix=Countmatrix[,-1];Countmatrix[is.na(Countmatrix)] = 0
# 特征信息
    feature_info <- data[,setdiff(names(data),c("Length","EffectiveLength", "NumReads", "TPM", "log2TPM", "BioSample"))];
    feature_info = distinct(feature_info)
    rownames(feature_info) = feature_info$marker
# 确保矩阵的列和sample_info的行顺序一致  
    TPMmatrix = TPMmatrix[,rownames(sample_info)]  
    Countmatrix = Countmatrix[,rownames(sample_info)]  
# 确保矩阵的行和feature_info的行一致  
    TPMmatrix = TPMmatrix[rownames(feature_info),]  
    Countmatrix = Countmatrix[rownames(feature_info),]
# 生成和保存结果
    se <- SummarizedExperiment(  
      assays = list(TPM = as.matrix(TPMmatrix), Count = as.matrix(Countmatrix)), # 将两个矩阵添加到 assays 槽中  
      rowData = feature_info,  
      colData = sample_info,  
      metadata = list() # 创建一个空的 metadata 列表  
    )  
    return(se)
}

Function_quant.sf_raw = function(datasetname){
    print(datasetname)
    sampletable_file = paste0("/data2/shaoxun/BloodCircleR/sampletable/", datasetname, "_sampletable.txt")  
    sampletable = fread(sampletable_file, data.table = F)  
    sampletable$BioSample = gsub("\\(.*","",sampletable$BioSample)
    datadir = "/data2/shaoxun/BloodCircleR/dataset/"
    quant.sf.list <- list()
      for (sampleindex in 1:length(sampletable$BioSample)) {
      print(sampletable$BioSample[sampleindex])
      quant_sf_path <- paste0(datadir,datasetname,"/",sampletable$BioSample[sampleindex],"/quant/profile_results/quant.sf")
      if (file.exists(quant_sf_path)) {
          quant.sf <- fread(quant_sf_path, data.table = F)
          quant.sf$BioSample <- sampletable$BioSample[sampleindex]
          quant.sf$datasetname <- datasetname
          quant.sf.list[[sampleindex]] <- quant.sf
      } else { 
          print("quant.sf 文件不存在，不执行操作...")}    
      }
    quant.sf.df <- do.call(rbind, quant.sf.list)
    rownames(quant.sf.df) <- NULL
    saveRDS(quant.sf.df,file=paste0(datadir,"/circRNA_quant/rawdata/quant.sf.",datasetname,".allsample.raw.rds"))
}

Function_quant.sf_filter = function(datasetname){
  record.list = list()
  datadir = "/data2/shaoxun/BloodCircleR/dataset/"
  message("初步过滤，数据集为:", datasetname)
  # 读取数据
  quant.sf.df = readRDS(paste0(datadir, "circRNA_quant/rawdata/quant.sf.", datasetname, ".allsample.raw.rds"))
  # 筛选线性和环状 RNA
    quant.linear = quant.sf.df[grep("ENST", quant.sf.df$Name),]
    quant.circular = quant.sf.df[grep("chr", quant.sf.df$Name),]
    quant.linear$gene_name = mapvalues(quant.linear$Name,gtf_transcript$transcript_id,gtf_transcript$gene_name,warn_missing = FALSE)
    quant.linear$gene_biotype = mapvalues(quant.linear$Name,gtf_transcript$transcript_id,gtf_transcript$gene_biotype,warn_missing = FALSE)
    quant.linear$BSJ_ID = NA
    quant.circular$gene_name = mapvalues(quant.circular$Name,BloodCircR_circRNA$isoformID,BloodCircR_circRNA$host_gene,warn_missing = FALSE)
    quant.circular$gene_biotype = mapvalues(quant.circular$gene_name,gtf_transcript$gene_name,gtf_transcript$gene_biotype,warn_missing = FALSE)
    quant.circular$BSJ_ID = mapvalues(quant.circular$Name,BloodCircR_circRNA$isoformID,BloodCircR_circRNA$BSJ_ID,warn_missing = FALSE)
    quant.linear$type = "linear"
    quant.circular$type = "cicular"
    quant.sf <- rbind(quant.linear, quant.circular)
    record.list[["Raw"]] = dplyr::count(quant.sf,type) # 保留原始
    quant.sf_filter <- quant.sf[!quant.sf$gene_name %in% c(HBgene$Approved.symbol, RPgene$Approved.symbol), ]
    quant.sf_filter = quant.sf_filter[quant.sf_filter$gene_biotype %in% c("protein_coding", "lincRNA"),]
    record.list[["TypeFilter"]] = dplyr::count(quant.sf_filter,type) # 过滤后
    quant.sf_filter_adj <- TPMrecalculate(quant.sf_filter)
    quant.sf_filter_adj <- quant.sf_filter_adj[, !names(quant.sf_filter_adj) %in% c("TPM")]
    names(quant.sf_filter_adj)[ncol(quant.sf_filter_adj)] <- "TPM"
    quant.sf_filter_adj$TPM = round(quant.sf_filter_adj$TPM, 4)
    # 保存结果
    saveRDS(quant.sf_filter_adj, file = paste0(datadir, "circRNA_quant/rawdata/quant.sf.", datasetname, ".allsample.TypeFilter.rds"))
    # 保存过滤过程
    saveRDS(record.list,file = paste0(datadir, "circRNA_quant/rawdata/quant.sf.", datasetname, ".record.list.rds")) 
    message("过滤完成，数据集为:", datasetname)            
}

# 读取数据，一次即可
datadir <- "/data2/shaoxun/BloodCircleR/dataset/"  
sampletable_file <- paste0("/data2/shaoxun/BloodCircleR/sampletable/", datasetname, "_sampletable.txt")  
sampletable <- fread(sampletable_file, data.table = FALSE)  
quant_sf_file <- paste0(datadir, "/circRNA_quant/rawdata/quant.sf.", datasetname, ".allsample.TypeFilter.rds")  
quant.sf_filter_adj <- readRDS(quant_sf_file) 
	
Function_quant.sf_linear <- function(datasetname) {  
  library(plyr) # 确保加载 plyr 包  
  library(dplyr) # 确保加载 dplyr 包  
  datadir <- "/data2/shaoxun/BloodCircleR/dataset/"  
  quant.linear_isoform <- quant.sf_filter_adj[grep("ENST", quant.sf_filter_adj$Name), ]  
  quant.linear_isoform$marker <- quant.linear_isoform$Name 
  quant.linear_gene <- quant.linear_isoform %>% dplyr::group_by(BioSample, gene_name) %>% dplyr::summarise(NumReads = sum(NumReads),TPM = sum(TPM) )  
  quant.linear_gene$marker <- quant.linear_gene$gene_name  
  quant.linear_gene$gene_biotype <- mapvalues(quant.linear_gene$marker, gtf_gene$gene_name, gtf_gene$gene_biotype, warn_missing = FALSE)  
  quant.linear_isoform$TPM <- round(quant.linear_isoform$TPM, 4)  
  quant.linear_gene$TPM <- round(quant.linear_gene$TPM, 4)  
  filepath <- paste0(datadir, "/circRNA_quant/quant.sf.", datasetname, ".linear.isoform.se.rds")  
  saveRDS(SEmake(quant.linear_isoform,sampletable), file = filepath)  
  filepath <- paste0(datadir, "/circRNA_quant/quant.sf.", datasetname, ".linear.gene.se.rds")  
  saveRDS(SEmake(quant.linear_gene,sampletable), file = filepath)  
}

Function_quant.sf_circular = function(datasetname){
    datadir <- "/data2/shaoxun/BloodCircleR/dataset/"  
message("cirular数据加载:", datasetname)  
    quant.circular_isoform <- quant.sf_filter_adj[grep("chr", quant.sf_filter_adj$Name), ]  
    quant.circular_isoform$marker = mapvalues(quant.circular_isoform$Name,BloodCircR_circRNA$isoformID,BloodCircR_circRNA$BloodCircR_ID,warn_missing = FALSE)# 赋值一次
message("cirular-isoform处理:",datasetname)  
    quant.circular_isoform_count <- dplyr::count(quant.circular_isoform, marker) %>% arrange(desc(n))  
    quant.circular_isoform$detectsample <- mapvalues(quant.circular_isoform$marker, quant.circular_isoform_count$marker, quant.circular_isoform_count$n, warn_missing = FALSE)   
    quant.circular_isoform$detectratio <- round(as.numeric(quant.circular_isoform$detectsample) / nrow(sampletable), 2)  
    filepath <- paste0(datadir, "/circRNA_quant/quant.sf.", datasetname, ".circular.isoform.se.rds")  
    saveRDS( SEmake(quant.circular_isoform,sampletable),file = filepath)  
message("cirular-isoform处理完成:",datasetname)   
message("cirular-bsj处理:",datasetname)  
    quant.circular_bsj <- quant.circular_isoform %>% dplyr::group_by(BioSample, BSJ_ID) %>% dplyr::summarise(NumReads = sum(NumReads), TPM = sum(TPM)) %>% as.data.frame()  
    quant.circular_bsj$marker = quant.circular_bsj$BSJ_ID
    quant.circular_bsj_count <- dplyr::count(quant.circular_bsj, marker) %>% arrange(desc(n))  
    quant.circular_bsj$detectsample <- mapvalues(quant.circular_bsj$marker, quant.circular_bsj_count$marker, quant.circular_bsj_count$n, warn_missing = FALSE)   
    quant.circular_bsj$detectratio <- round(as.numeric(quant.circular_bsj$detectsample) / nrow(sampletable), 2)  
    quant.circular_bsj$TPM <- round(quant.circular_bsj$TPM, 4)  
    quant.circular_bsj$log2TPM <- round(log2(quant.circular_bsj$TPM+1), 4)  
    quant.circular_bsj$gene_name <- mapvalues(quant.circular_bsj$marker, BloodCircR_circRNA$BSJ_ID, BloodCircR_circRNA$host_gene, warn_missing = FALSE) 
    quant.circular_bsj$gene_biotype <- mapvalues(quant.circular_bsj$gene_name, gtf_gene$gene_name, gtf_gene$gene_biotype, warn_missing = FALSE) 
    filepath <- paste0(datadir, "/circRNA_quant/quant.sf.", datasetname, ".circular.bsj.se.rds")  
    saveRDS( SEmake(quant.circular_bsj,sampletable),file = filepath) 
message("cirular-bsj处理完成:",datasetname)   
message("cirular-gene处理:",datasetname)   
    quant.circular_gene <- quant.circular_isoform %>% dplyr::group_by(BioSample, gene_name) %>% dplyr::summarise(NumReads = sum(NumReads), TPM = sum(TPM)) %>% as.data.frame()  
    quant.circular_gene$marker = quant.circular_gene$gene_name
    quant.circular_gene_count <- dplyr::count(quant.circular_gene, marker) %>% arrange(desc(n))  
    quant.circular_gene$detectsample <- mapvalues(quant.circular_gene$marker, quant.circular_gene_count$marker, quant.circular_gene_count$n, warn_missing = FALSE)   
    quant.circular_gene$detectratio <- round(as.numeric(quant.circular_gene$detectsample) / nrow(sampletable), 2)  
    quant.circular_gene$TPM <- round(quant.circular_gene$TPM, 4)  
    quant.circular_gene$log2TPM <- round(log2(quant.circular_gene$TPM+1), 4)  
    quant.circular_gene$gene_biotype <- mapvalues(quant.circular_gene$gene_name, gtf_gene$gene_name, gtf_gene$gene_biotype, warn_missing = FALSE) 
    filepath <- paste0(datadir, "/circRNA_quant/quant.sf.", datasetname, ".circular.gene.se.rds")  
    saveRDS( SEmake(quant.circular_gene,sampletable),file = filepath) 
message("cirular-gene处理完成:",datasetname)   
}

Function_quant.sf_filter(datasetname)
Function_quant.sf_linear(datasetname)
Function_quant.sf_circular(datasetname)
