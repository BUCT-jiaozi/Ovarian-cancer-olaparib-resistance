library(Seurat)
library(harmony)
library(clustree)
library(ggpubr)
library(dplyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(future)
library(tidyverse)
library(gridExtra)
library(ggridges)
library(ggExtra)
library(data.table)
library(readxl)
library(SeuratDisk)
rm(list = ls());gc()
set.seed(101)
plan("multisession", workers = 2) 
options(future.globals.maxSize = 10000 * 1024^2) # set 50G RAM
setwd('/home/Rstudio/yhb/Data/OC_data/kura_new_fq/')
dir <- "./"

### 分开读取单个样本创建seurat对象
if(T){
  samples <- dir(path="./kb_allresult/counts_filte",
                 # pattern="^CID", # 遍历字符开头所有文件夹下的文件
                 pattern="^SRR", 
                 full.names=T, 
                 recursive=T, include.dirs=T)# 递归搜索
  
  

  # 读取seurat对象文件，meta.data信息为空
  samples_objects <- lapply(samples, function(s) {
    matrix_dir = s
    
    # 设置文件的路径
    h5seurat_path <- paste0(matrix_dir, "/adata.h5seurat")
    
    # 读取
    seurat_object <- LoadH5Seurat(file = h5seurat_path)
  })
  names(samples_objects) <- c("T320", "T160", "T80", "T40", "T20", "T10", "T5", "T2.5", "T1", "C") # samples的顺序
  
  rownames(samples_objects$T320)
  colnames(samples_objects$T320)
  
  # samples_objects <- lapply(seq_along(samples_raw_data), function(i) {
  #   
  #   seurat_object <- CreateSeuratObject(counts = samples_raw_data[[i]], project = names(samples_raw_data)[i], min.cells = 5)
  #   # seurat_object <- CreateSeuratObject(counts = samples_raw_data[[i]], project = names(samples_raw_data)[i], names.delim = "_", names.field = 2)
  # }) # 基因会自动去重
  # names(samples_objects) <- names(samples_raw_data)
  new_order <-  c("C", "T1", "T2.5", "T5", "T10", "T20", "T40", "T80", "T160", "T320") ## 正确排列样本顺序
  samples_objects <- samples_objects[match(new_order, names(samples_objects))]
  
  dim(samples_objects$T320)
}


# 基因名转换，修改细胞信息
if(F){
  
  source(file = "../code/Functions/RenameGenesSeurat.R")
  gene_use <- read.table('/home/Rstudio/yhb/sup_files/human_genes.tsv',sep='\t',header=F,stringsAsFactor=F) # 比对生成的结果，39546基因
  gene_use <- gene_use[gene_use$V2 != "", ] # 去除空值 25705
  gene_use <- gene_use[!duplicated(gene_use$V2), ] # 对 V2 列进行去重 25668
  
  for (i in seq_along(samples_objects)){
    #colnames(samples_objects[[i]]) <- paste0(names(samples_objects)[i], "_", colnames(samples_objects[[i]])) # 细胞名加样本信息
    #samples_objects[[i]]$orig.ident <- names(samples_objects)[i] # 添加组别信息
    # 基因转换
    samples_objects[[i]] <- RenameGenesSeurat(obj=samples_objects[[i]], newnames=gene_use$V2, gene.use=gene_use$V1, de.assay="RNA")

  }
  
  # 修改细胞名，添加样本信息
  for (i in seq_along(samples_objects)){
    samples_objects[[i]] <- RenameCells(samples_objects[[i]], new.names=paste0(names(samples_objects)[i], "_", colnames(samples_objects[[i]])))
    samples_objects[[i]]$orig.ident <- names(samples_objects)[i] # 添加组别信息
  }
}



### 单样本顺序质控
if(T){
  # 计算线粒体含量
  for (i in seq_along(samples_objects)){
    mito_genes=rownames(samples_objects[[i]])[grep("^MT-", rownames(samples_objects[[i]]),ignore.case = T)] 
    ribosomal_genes=rownames(samples_objects[[i]])[grep("^RP[SL]", rownames(samples_objects[[i]]),ignore.case = T)]
    
    samples_objects[[i]] <- PercentageFeatureSet(samples_objects[[i]], 
                                                 #pattern = "^MT-",
                                                 features = mito_genes,
                                                 col.name = "percent.mt")
    samples_objects[[i]] <- PercentageFeatureSet(samples_objects[[i]], 
                                                 features = ribosomal_genes,
                                                 col.name = "percent.ribo")
  }
  # names(samples_objects) <- names(samples_raw_data)
  
  dir.create('kb_fitle_R')
  setwd("./kb_fitle_R")
  dir.create('./1.QualityControl')
  pdf(file = "./1.QualityControl/1.1count.feature.mt.pdf", width = 8, height = 6)
  for (i in seq_along(samples_objects)){
    #质控小提琴图
    print(VlnPlot(samples_objects[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0, ncol = 2, same.y.lims=F) +
            scale_y_continuous(breaks=seq(0, 100, 10)) +
            NoLegend()) ## 数据已经过滤过
    plot1 <- FeatureScatter(samples_objects[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
    ploT1 <- FeatureScatter(samples_objects[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    print(plot1 + ploT1)
    
    print(ggdensity(samples_objects[[i]]@meta.data, x = "nCount_RNA", title = names(samples_objects[i])))
    print(ggdensity(samples_objects[[i]]@meta.data, x = "nFeature_RNA", title = names(samples_objects[i])))
    print(ggdensity(samples_objects[[i]]@meta.data, x = "percent.mt", title = names(samples_objects[i])))
  }
  dev.off()
}

### 质控，初步标准化分群
if(T){
  # 质控
  samples_objects <- lapply(seq_along(samples_objects), function(i) {
    samples_objects[[i]] <- subset(x = samples_objects[[i]],
                                   subset = nFeature_RNA > 500 & percent.mt < 15 & nCount_RNA < 20000 & nFeature_RNA < 5000 & percent.ribo < 30)
  })
  names(samples_objects) <- c("C", "T1", "T2.5", "T5", "T10", "T20", "T40", "T80", "T160", "T320")
  dim(samples_objects$T160)
  
  samples_objects <- lapply(seq_along(samples_objects), function(i) {
    samples_objects[[i]] <- NormalizeData(object = samples_objects[[i]], normalization.method= "LogNormalize", verbose = FALSE)
    samples_objects[[i]] <- FindVariableFeatures(object = samples_objects[[i]],
                                                 selection.method = "vst", 
                                                 nfeatures = 2000, verbose = FALSE)
    samples_objects[[i]] <- ScaleData(object = samples_objects[[i]], vars.to.regress = c("nCount_RNA", "percent.mt", "percent.ribo"))
  })
  names(samples_objects) <- c("C", "T1", "T2.5", "T5", "T10", "T20", "T40", "T80", "T160", "T320")
  
  set.resolutions <- seq(0.1, 1, by = 0.1)
  
  samples_objects <- lapply(seq_along(samples_objects), function(i) {
    samples_objects[[i]] <- RunPCA(samples_objects[[i]], npcs = 100, verbose = FALSE)
    #ElbowPlot(object = samples_objects[[i]], ndims = 100)
    samples_objects[[i]] <- FindNeighbors(object = samples_objects[[i]], dims = 1:50, verbose = FALSE)
    samples_objects[[i]] <- FindClusters(object = samples_objects[[i]], resolution = set.resolutions, verbose = FALSE) 
    #clustree(samples_objects[[i]])
    samples_objects[[i]] <- RunUMAP(samples_objects[[i]], dims = 1:50)
  })
  names(samples_objects) <- c("C", "T1", "T2.5", "T5", "T10", "T20", "T40", "T80", "T160", "T320")
  
  pdf(file = "./1.QualityControl/1.2clusters.pdf", width = 8, height = 6)
  for (i in seq_along(samples_objects)){
    print(ElbowPlot(object = samples_objects[[i]], ndims = 100))
    print(clustree(samples_objects[[i]]))
    
    samples_objects.res <- sapply(set.resolutions, function(x){
      p <- DimPlot(object = samples_objects[[i]], reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
      print(p)
    })
  }
  dev.off()
  saveRDS(samples_objects, file = "1.0data.rds")
}

samples_objects <- readRDS("./1.0data.rds")


### 单样本去除双细胞
if(T){
  #### 去除双细胞
  # 方式1
  # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
  # remotes::install_github('lzmcboy/DoubletFinder_204_fix')  ##Seurat是V4
  library(DoubletFinder) # Require cleanup of low-quality cells in advance
  source(file = "../code/Functions/doubletDetect.R")
  
  samples_double <- list()
  samples_double[[1]] <- doubletDetect(Seurat.object = samples_objects[[1]], PCs = 1:50, doublet.rate = 0.132, annotation = "RNA_snn_res.0.3", sct = T) #12373  1839
  samples_double[[2]] <- doubletDetect(Seurat.object = samples_objects[[2]], PCs = 1:50, doublet.rate = 0.061, annotation = "RNA_snn_res.0.3", sct = T) #2294
  samples_double[[3]] <- doubletDetect(Seurat.object = samples_objects[[3]], PCs = 1:50, doublet.rate = 0.031, annotation = "RNA_snn_res.0.4", sct = T) #1544
  samples_double[[4]] <- doubletDetect(Seurat.object = samples_objects[[4]], PCs = 1:50, doublet.rate = 0.031, annotation = "RNA_snn_res.0.3", sct = T) #716
  samples_double[[5]] <- doubletDetect(Seurat.object = samples_objects[[5]], PCs = 1:50, doublet.rate = 0.061, annotation = "RNA_snn_res.0.3", sct = T) #1974
  samples_double[[6]] <- doubletDetect(Seurat.object = samples_objects[[6]], PCs = 1:50, doublet.rate = 0.041, annotation = "RNA_snn_res.0.4", sct = T) #1452
  samples_double[[7]] <- doubletDetect(Seurat.object = samples_objects[[7]], PCs = 1:50, doublet.rate = 0.061, annotation = "RNA_snn_res.0.2", sct = T) #2052
  samples_double[[8]] <- doubletDetect(Seurat.object = samples_objects[[8]], PCs = 1:50, doublet.rate = 0.061, annotation = "RNA_snn_res.0.1", sct = T) #2183
  samples_double[[9]] <- doubletDetect(Seurat.object = samples_objects[[9]], PCs = 1:50, doublet.rate = 0.091, annotation = "RNA_snn_res.0.3", sct = T) #2611
  samples_double[[10]] <- doubletDetect(Seurat.object = samples_objects[[10]], PCs = 1:50, doublet.rate = 0.031, annotation = "RNA_snn_res.0.7", sct = T) #1281
  names(samples_double) <- c("C", "T1", "T2.5", "T5", "T10", "T20", "T40", "T80", "T160", "T320")
  # saveRDS(samples_double, "1.1data_double.rds")
  
  
  samples_singlet <- list()
  samples_singlet <- lapply(seq_along(samples_double), function(i) {
    samples_singlet[[i]] <- subset(x = samples_double[[i]],subset = Doublet == "Singlet")
  })
  names(samples_singlet) <- c("C", "T1", "T2.5", "T5", "T10", "T20", "T40", "T80", "T160", "T320")
  saveRDS(samples_singlet, file = "1.2data_singlecell.rds")
  
  
  table(samples_double[[1]]$Doublet)
  # Doublet Singlet 
  # 94    1745  112 1184
  table(samples_double[[2]]$Doublet)
  # Doublet Singlet 
  # 83    2211  73 1914
  table(samples_double[[3]]$Doublet)
  # Doublet Singlet 
  # 48    1496  14 763 
  table(samples_double[[4]]$Doublet)
  # Doublet Singlet 
  # 5     711  12  621 
  table(samples_double[[5]]$Doublet)
  # Doublet Singlet 
  # 85    1889  110 2466
  table(samples_double[[6]]$Doublet)
  # Doublet Singlet 
  # 42    1410  33 1092
  table(samples_double[[7]]$Doublet)
  # Doublet Singlet 
  # 88    1964  161 3599 
  table(samples_double[[8]]$Doublet)
  # Doublet Singlet 
  # 89    2094  112 2518 
  table(samples_double[[9]]$Doublet)
  # Doublet Singlet 
  # 143    2468  110 1897 
  table(samples_double[[10]]$Doublet)
  # Doublet Singlet 
  # 25    1256  21 1088
  
  # pdf("1.QC/1.4doublet.cell.pdf")
  # DimPlot(object = samples_double[[10]], reduction = 'umap', group.by = "Doublet")
  # dev.off()
  
  
}


### 单样本的 infercnv
if(T){
  samples_singlet <- readRDS("./1.2data_singlecell.rds")
  
  # umap定制颜色
  if(F){
    # 定义颜色顺序列表
    all_cluster_colors <- list(
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041"),   # C
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041", "3" = "#501d8a"),  # T1
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041"),   # T2
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041"),   # T3
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041", "3" = "#501d8a"),  # T4
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041", "3" = "#501d8a"),  # T5
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041", "3" = "#501d8a"),  # T6
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041", "3" = "#501d8a"),  # T7
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041"),   # T8
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041")    # T9
    )
    all_annotation_map <- list(
      C     = "RNA_snn_res.0.3",
      T1    = "RNA_snn_res.0.3",
      T2.5  = "RNA_snn_res.0.4",
      T5    = "RNA_snn_res.0.3",
      T10   = "RNA_snn_res.0.3",
      T20   = "RNA_snn_res.0.4",
      T40   = "RNA_snn_res.0.2",
      T80   = "RNA_snn_res.0.1",
      T160  = "RNA_snn_res.0.3",
      T320  = "RNA_snn_res.0.7"
    )
    # 遍历每个 Seurat 对象并绘制图
    plots <- list()
    for (i in 1:length(samples_singlet)) {
      # 选择相应的 Seurat 对象
      seurat_object <- samples_singlet[[i]]
      
      # 从 cluster_colors 中选择相应的颜色
      colors <- all_cluster_colors[[i]]
      annotation_column <- all_annotation_map[[i]]
      
      seurat_object$orig.ident <- seurat_object[[annotation_column]]

      
      data.merge_umap_sample_plot <- CellDimPlot(seurat_object, group.by = "orig.ident", reduction = "UMAP", theme_use = "theme_blank", #theme_blank theme_void theme_linedraw theme_light
                  #title = "Kuramochi",
                  palcolor = colors,
                  #theme_use = ggplot2::theme_classic, 
                  #cells.highlight = colnames(pancreas_sub)[pancreas_sub$SubCellType == "Epsilon"],
                  #label = TRUE, label_insitu = TRUE, label_repel = TRUE, label_segment_color = "red",
                  legend.position  = "right", legend.direction = "vertical",
                  theme_args = list(plot.title = element_text(hjust = 0.5, # 标题居中（关键）
                                                              size  = 14, face  = "bold"), legend.title = element_text(size = 12), 
                                    legend.key.height = unit(0.6, "cm"),  # 图例间距
                                    legend.spacing.y  = unit(4, "pt")))+
        scale_color_manual(
          name =  names(samples_singlet)[i],
          values = colors,
          guide = guide_legend(
            override.aes = list(size = 3)  # 👈 图例圆点大小
          ))
      # 将每个图保存到列表中
      plots[[i]] <- data.merge_umap_sample_plot
    }
    
    # 查看某个图，假设查看第一个
    print(plots[[1]])
    
    # 如果需要将所有图保存为图片
    dir.create("./2.Cluster")
    pdf("./2.Cluster/2.3all_clusters.pdf", width = 5, height = 3)
    for (i in 1:length(plots)) {
      #ggsave(paste0("UMAP_Plot_", i, ".png"), plot = plots[[i]], width = 8, height = 6)
      print(plots[[i]])
    }
    dev.off()
    pdf("./2.Cluster/2.6all_clusters.pdf", width = 12, height = 4)
    wrap_plots(plots,ncol = 5)
    dev.off()
  }
  
  
  ##Predict the CNV of potential cancer cells based on the annotated normal cells
  library(Seurat)
  library(infercnv) #Version:1.6
  #data.merge <- readRDS(file = "./1.8data.merge.pro.rds")
  
  # 对应 annotation 的映射
  annotation_map <- list(
    T1    = "RNA_snn_res.0.3",
    T2.5  = "RNA_snn_res.0.4",
    T5    = "RNA_snn_res.0.3",
    T10   = "RNA_snn_res.0.3",
    T20   = "RNA_snn_res.0.4",
    T40   = "RNA_snn_res.0.2",
    T80   = "RNA_snn_res.0.1",
    T160  = "RNA_snn_res.0.3",
    T320  = "RNA_snn_res.0.7"
  )
  
  # 参考样本
  ref_obj <- samples_singlet[["C"]]
  ref_cells <- colnames(ref_obj)
  ref_label <- rep("C", length(ref_cells))
  
  dir.create('3.CNV')
  setwd("./3.CNV")
  
  # 循环其它样本
  for (s in names(annotation_map)) {
    message("Running inferCNV for sample: ", s)
    
    obj <- samples_singlet[[s]]
    
    # counts 矩阵：拼接参考和肿瘤
    counts_matrix <- cbind(GetAssayData(ref_obj, slot = "counts"),
                           GetAssayData(obj, slot = "counts"))
    
    # annoFile: 参考 + 当前样本
    anno_ref <- data.frame(cell.id = ref_cells, celltype = ref_label)
    anno_tumor <- data.frame(cell.id = colnames(obj),
                             celltype = obj@meta.data[[annotation_map[[s]]]])
    annoFile <- rbind(anno_ref, anno_tumor)
    rownames(annoFile) <- annoFile$cell.id
    
    # 保存到临时文件
    write.table(annoFile, file = paste0("anno_", s, ".txt"),
                col.names = FALSE, row.names = F, quote = FALSE, sep = "\t")
    
    # 创建 inferCNV 对象
    infercnv_obj <- CreateInfercnvObject(
      raw_counts_matrix = counts_matrix,
      annotations_file = paste0("anno_", s, ".txt"),
      delim = "\t",
      gene_order_file = "../../code/sup_files/hg38_gencode_v27.txt",
      ref_group_names = "C"
    )
    
    # 运行 infercnv
    out_dir <- paste0("infercnv_", s)
    dir.create(out_dir, showWarnings = FALSE)
    
    infercnv_obj <- infercnv::run(infercnv_obj,
                                  cutoff = 0.1,
                                  out_dir = out_dir,
                                  cluster_by_groups = TRUE,
                                  denoise = F,
                                  HMM = T, # 基于HNN预测cnv，速度慢
                                  #analysis_mode = "cells", # samples|subclusters|cells 图像过滤或HMM预测的分组级别。默认值：samples
                                  num_threads = 8,
                                  write_expr_matrix = T, # 生成绘图时是否写入包含矩阵内容的文本文件
                                  cluster_references = TRUE # 在注释中对引用进行聚类，不显示树状图，默认T
    )
    
    saveRDS(infercnv_obj, file = paste0("infercnv_", s, ".rds"))
  }
}






