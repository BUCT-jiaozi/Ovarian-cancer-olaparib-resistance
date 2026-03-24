library(Seurat)

rm(list = ls());gc()
set.seed(101)
plan("multisession", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)
setwd('/home/Rstudio/yhb/Data/ovsaho_sc/')

raw.data <- read.csv('./ovsaho_all_C_to_T40.csv', header = TRUE, row.names = 1)

### 分离每一个样本，
if(T){
  dir.create("./data")
  # 定义前缀
  prefixes <- c('C_', 'T5_', 'T10_', 'T20_', 'T40_')
  
  # 创建一个空列表用于存储分割后的数据框
  dfs <- list()
  
  # 按前缀筛选列，并保存到不同数据框中
  for (prefix in prefixes) {
    # 使用 grep() 函数根据前缀筛选列
    selected_columns <- grep(paste0("^", prefix), colnames(raw.data), value = TRUE)
    
    # 将筛选后的数据框保存到列表中
    dfs[[prefix]] <- raw.data[, selected_columns, drop = FALSE]
    
    # 保存每个数据框为单独的 CSV 文件
    write.csv(dfs[[prefix]], file = paste0("./data/bc_", prefix, "data.csv"), row.names = TRUE)
  }
}


### 分开读取单个样本创建seurat对象
if(T){
  samples <- dir(path="./data",
                 # pattern="^CID", # 遍历字符开头所有文件夹下的文件
                 pattern="^bc", 
                 full.names=T, 
                 recursive=T, include.dirs=T)# 递归搜索
  
  samples_raw_data <- lapply(samples, function(s) {
    
    matrix_dir = s
    matrix_path <- paste0(matrix_dir)
    mat <- read.csv(file = matrix_path, header = T, row.names = 1)
    
    mat
  })
  # names(samples_raw_data) <- c(dir("~/projects/breast/GSE176078/CID",pattern="CID."))
  names(samples_raw_data) <- c("C", "T10", "T20", "T40", "T5") # samples的顺序
  
  samples_objects <- lapply(seq_along(samples_raw_data), function(i) {
    seurat_object <- CreateSeuratObject(counts = samples_raw_data[[i]], project = names(samples_raw_data)[i])
    # seurat_object <- CreateSeuratObject(counts = samples_raw_data[[i]], project = names(samples_raw_data)[i], names.delim = "_", names.field = 2)
  })
  names(samples_objects) <- names(samples_raw_data)
  new_order <-  c("C", "T5", "T10", "T20", "T40") ## 正确排列样本顺序
  samples_objects <- samples_objects[match(new_order, names(samples_objects))]
  dim(samples_objects$T5)
}

### 单样本顺序质控
if(T){
  # 计算线粒体含量
  for (i in seq_along(samples_objects)){
    samples_objects[[i]] <- PercentageFeatureSet(samples_objects[[i]], 
                                                 pattern = "^MT-",
                                                 col.name = "percent.mt")
    samples_objects[[i]] <- PercentageFeatureSet(samples_objects[[i]], 
                                                 pattern = "^RP[SL]",
                                                 col.name = "percent.ribo")
  }
  # names(samples_objects) <- names(samples_raw_data)
  
  # dir.create('1.QualityControl')
  pdf(file = "1.QualityControl/1.1count.feature.mt.pdf", width = 8, height = 6)
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
                                   subset = nFeature_RNA > 300 & percent.mt < 18 & nCount_RNA < 8500 & nFeature_RNA < 3500 & percent.ribo < 48 & percent.ribo > 8)
  })
  names(samples_objects) <- c("C", "T5", "T10", "T20", "T40")
  dim(samples_objects$T5)
  
  samples_objects <- lapply(seq_along(samples_objects), function(i) {
    samples_objects[[i]] <- NormalizeData(object = samples_objects[[i]], normalization.method= "LogNormalize", verbose = FALSE)
    samples_objects[[i]] <- FindVariableFeatures(object = samples_objects[[i]],
                                                 selection.method = "vst", 
                                                 nfeatures = 2000, verbose = FALSE)
    samples_objects[[i]] <- ScaleData(object = samples_objects[[i]], vars.to.regress = c("nCount_RNA", "percent.mt", "percent.ribo"))
  })
  names(samples_objects) <- c("C", "T5", "T10", "T20", "T40")
  
  set.resolutions <- seq(0.05, 0.3, by = 0.05)
  
  samples_objects <- lapply(seq_along(samples_objects), function(i) {
    samples_objects[[i]] <- RunPCA(samples_objects[[i]], npcs = 100, verbose = FALSE)
    #ElbowPlot(object = samples_objects[[i]], ndims = 100)
    samples_objects[[i]] <- FindNeighbors(object = samples_objects[[i]], dims = 1:30, verbose = FALSE)
    samples_objects[[i]] <- FindClusters(object = samples_objects[[i]], resolution = set.resolutions, verbose = FALSE) 
    #clustree(samples_objects[[i]])
    samples_objects[[i]] <- RunUMAP(samples_objects[[i]], dims = 1:30)
  })
  names(samples_objects) <- c("C", "T5", "T10", "T20", "T40")
  
  pdf(file = "1.QualityControl/1.2.1single_clusters.pdf", width = 5, height = 3)
  for (i in seq_along(samples_objects)){
    print(ElbowPlot(object = samples_objects[[i]], ndims = 100))
    print(clustree(samples_objects[[i]]))
    
    samples_objects.res <- sapply(set.resolutions, function(x){
      p <- DimPlot(object = samples_objects[[i]], reduction = 'umap',label = TRUE, group.by = paste0("RNA_snn_res.", x))
      print(p)
    })
  }
  dev.off()
  saveRDS(samples_objects, file = "1.2.1sample_data.rds")
}

samples_singlet <- readRDS("./1.2.1sample_data.rds")



### 单样本的 infercnv
if(T){
  #samples_singlet <- readRDS("./1.2data_singlecell.rds")
  
  # umap定制颜色
  if(F){
    # 定义颜色顺序列表
    all_cluster_colors <- list(
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041", "3" = "#501d8a"),   # C
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041", "3" = "#501d8a"),  # T5
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041"),   # T10
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041", "3" = "#501d8a"),   # T20
      c("0" = "#aa3474", "1" = "#ee8c7d", "2" = "#1c8041")  # T40

    )
    all_annotation_map <- list(
      C     = "RNA_snn_res.0.15",
      T5    = "RNA_snn_res.0.2",
      T10   = "RNA_snn_res.0.25",
      T20   = "RNA_snn_res.0.2",
      T40   = "RNA_snn_res.0.2"

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
    # dir.create("./2.Cluster")
    # pdf("./2.Cluster/2.3all_clusters.pdf", width = 5, height = 3)
    # for (i in 1:length(plots)) {
    #   #ggsave(paste0("UMAP_Plot_", i, ".png"), plot = plots[[i]], width = 8, height = 6)
    #   print(plots[[i]])
    # }
    # dev.off()
    pdf("./2.Cluster/2.5all_clusters.pdf", width = 12, height = 2)
    wrap_plots(plots,ncol = 5)
    dev.off()
  }
 




