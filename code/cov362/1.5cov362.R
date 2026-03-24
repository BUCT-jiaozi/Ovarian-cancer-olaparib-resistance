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
setwd('/home/Rstudio/yhb/Data/cov362_sc/')
dir <- "./"
#devtools::install_version("ggplot2", "3.5.2", force = TRUE)

### 分开读取单个样本创建seurat对象
if(T){
  samples <- dir(path="./outputs/matrix",
                 # pattern="^CID", # 遍历字符开头所有文件夹下的文件
                 pattern="^SRR", 
                 full.names=T, 
                 recursive=T, include.dirs=T)# 递归搜索
  
  htos <- dir(path="./d4citeseq",
                 # pattern="^CID", # 遍历字符开头所有文件夹下的文件
                 pattern="^SRR", 
                 full.names=T, 
                 recursive=T, include.dirs=T)# 递归搜索
  
  
  # 读取seurat对象文件，meta.data信息为空
  samples_raw_data <- lapply(samples, function(s) {
    matrix_dir = s
    
    # 读取counts
    mat_counts <- Read10X(data.dir = paste0(matrix_dir))
  })
  samples
  #names(samples_raw_data) <- c("SRR26822983", "SRR26822984", "SRR26822985", "SRR26822986", "SRR26822987") # samples的顺序
  # hto信息
  hto_data <- lapply(htos, function(s) {
    hto_dir = s
    
    hto_path <- paste0(hto_dir, "/umi_count")
    # 读取hto
    mat_htos <- Read10X(data.dir = hto_path, gene.column=1)
  })
  htos
  #names(hto_data) <- c("SRR26822983", "SRR26822984", "SRR26822985", "SRR26822986", "SRR26822987")
  
  # 统一hto细胞名称格式
  hto_data <- lapply(seq_along(hto_data), function(i) {
    # colnames(samples_raw_data[[i]]) <- paste0(names(samples_raw_data)[i], "_", colnames(samples_raw_data[[i]]))
    colnames(hto_data[[i]]) <- paste0(colnames(hto_data[[i]]), "-1")
    return(hto_data[[i]])  # 返回修改后的结果
  })

}


# 对于每个样本和对应的 HTO 数据进行操作合并
samples_objects <- lapply(seq_along(samples_raw_data), function(i) {
  
  # 获取对应的 RNA 数据和 HTO 数据
  rna_data <- samples_raw_data[[i]]
  hto_mat <- hto_data[[1]] # 只用第一个有哈希信息
  
  # Select cell barcodes detected by both RNA and HTO
  joint.bcs <- intersect(colnames(rna_data), colnames(hto_mat))
  
  # Subset RNA and HTO counts by joint cell barcodes
  pbmc.umis <- rna_data[, joint.bcs]
  pbmc.htos <- as.matrix(hto_mat[, joint.bcs])
  
  # Setup Seurat object
  pbmc.hashtag <- CreateSeuratObject(counts = Matrix::Matrix(as.matrix(pbmc.umis), sparse = TRUE))
  
  # Normalize RNA data with log normalization
  pbmc.hashtag <- NormalizeData(pbmc.hashtag)
  pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
  pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))
  
  # Add HTO data as a new assay independent from RNA
  pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
  
  # Normalize HTO data using CLR transformation
  pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")
  
  # Demux cells based on HTO
  pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)
  
  # Assign classification results
  table(pbmc.hashtag$HTO_classification.global)
  table(pbmc.hashtag$HTO_maxID) # 药物处理组别信息
  
  
  # Create violin plot for UMI counts
  Idents(pbmc.hashtag) <- "HTO_classification.global"
  VlnPlot(pbmc.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
  
  # 给 singlet 数据添加对应的 orig.ident 信息
  pbmc.hashtag@meta.data$orig.ident <- recode(pbmc.hashtag@meta.data$HTO_maxID,
                                              "HTO-6-GGTTGCCAGATGTCA" = "C",
                                              "HTO-7-TGTCTTTCCTGCCAG" = "T5",
                                              "HTO-8-CTCCTCTGCAATTAC" = "T10",
                                              "HTO-9-CAGTAGTCACGGTCA" = "T20",
                                              "HTO-10-ATTGACCCGCGTTAG" = "T40")
  
  # 返回处理后的 pbmc.hashtag 对象
  return(pbmc.hashtag)
})

# 修改细胞名，添加样本信息
for (i in seq_along(samples_objects)){
    samples_objects[[i]] <- RenameCells(samples_objects[[i]], new.names=paste0(samples_objects[[i]]$orig.ident, "_", colnames(samples_objects[[i]]), "_s", i))
    # samples_objects[[i]]$orig.ident <- names(samples_objects)[i] # 添加组别信息
}

names(samples_objects) <- c("SRR26822983", "SRR26822984", "SRR26822985", "SRR26822986", "SRR26822987")
saveRDS(samples_objects, file = "1.0data_allcell.rds")


samples_singlet <- list()
samples_singlet <- lapply(seq_along(samples_objects), function(i) {
  samples_singlet[[i]] <- subset(x = samples_objects[[i]],subset = HTO_classification.global == "Singlet") ## Extract singlets (single cells)
})

names(samples_singlet) <- c("SRR26822983", "SRR26822984", "SRR26822985", "SRR26822986", "SRR26822987")

# 先导出一次细胞信息用于速率分析loom文件改名
for (i in seq_along(samples_singlet)){
  write.csv(Cells(samples_singlet[[i]]), file = paste0(names(samples_singlet)[i],"_cellID_obs.csv"), row.names = FALSE)
}


saveRDS(samples_singlet, file = "1.1data_singlecell.rds")

samples_singlet <- readRDS('./1.1data_singlecell.rds')

data.merge <- merge(samples_singlet[[1]], y = samples_singlet[2:length(samples_singlet)], project = "scRNA")


# 计算线粒体含量
mito_genes <- rownames(data.merge)[grep("^MT-", rownames(data.merge),ignore.case = T)] 
ribosomal_genes <- rownames(data.merge)[grep("^RP[SL]", rownames(data.merge),ignore.case = T)]
data.merge <- PercentageFeatureSet(data.merge, 
                                             #pattern = "^MT-",
                                             features = mito_genes,
                                             col.name = "percent.mt")
data.merge <- PercentageFeatureSet(data.merge, 
                                             features = ribosomal_genes,
                                             col.name = "percent.ribo")

dir.create('./1.QualityControl')
pdf(file = "./1.QualityControl/1.1count.feature.mt.pdf", width = 8, height = 6)
print(VlnPlot(data.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0, ncol = 2, same.y.lims=F) +
        scale_y_continuous(breaks=seq(0, 100, 10)) +
        NoLegend())
plot1 <- FeatureScatter(data.merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
ploT1 <- FeatureScatter(data.merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + ploT1)

print(ggdensity(data.merge@meta.data, x = "nCount_RNA", title = "cov362"))
print(ggdensity(data.merge@meta.data, x = "nFeature_RNA", title = "cov362"))
print(ggdensity(data.merge@meta.data, x = "percent.mt", title = "cov362"))
dev.off()



### 质控，初步标准化分群
# 质控
data.merge <- subset(x = data.merge,subset = nFeature_RNA > 450 & percent.mt < 15 & nCount_RNA < 8500 & nFeature_RNA < 4000 & percent.ribo < 40 & percent.ribo > 5)
dim(data.merge)
  
data.merge <- NormalizeData(object = data.merge, normalization.method= "LogNormalize", verbose = FALSE)
data.merge <- FindVariableFeatures(object = data.merge,selection.method = "vst", 
                                                 nfeatures = 2000, verbose = FALSE)
data.merge <- ScaleData(object = data.merge, vars.to.regress = c("nCount_RNA", "percent.mt", "percent.ribo"))

set.resolutions <- seq(0.1, 0.3, by = 0.05)

data.merge <- RunPCA(data.merge, npcs = 100, verbose = FALSE)
ElbowPlot(object = data.merge, ndims = 100)
data.merge <- FindNeighbors(object = data.merge, dims = 1:30, verbose = FALSE)
data.merge <- FindClusters(object = data.merge, resolution = set.resolutions, verbose = FALSE) 
clustree(data.merge)
data.merge <- RunUMAP(data.merge, dims = 1:30)

data.merge@meta.data$orig.ident <- factor(data.merge@meta.data$orig.ident, 
                                          levels = c("C", "T5", "T10", "T20", "T40"))

DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "RNA_snn_res.0.2")
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "orig.ident")


### 修改cluster编号顺序
if(F){
  Idents(data.merge) <- data.merge$RNA_snn_res.0.2
  # rename clusters to put in order of more resistant
  data.merge <- RenameIdents(object = data.merge,
                             "0" = "c1",
                             "1" = "c2",
                             "2" = "c0",
                             "3" = "c0",
                             "4" = "c3")
  # assign new order
  data.merge@active.ident <- factor(x = data.merge@active.ident,
                                    levels = c("c0", "c1", "c2", "c3"))
  
  data.merge@meta.data$cellType_low <- data.merge@active.ident
  Idents(data.merge) <- data.merge@meta.data$cellType_low
  data.merge@meta.data$cellType_low <- factor(data.merge@meta.data$cellType_low, 
                                        levels = c("c0", "c1", "c2", "c3"))
  
  DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "cellType_low")
}


source(file = "../OC_data/kura_new_fq/code/Functions/ratio.plot.R")

pdf(file = "./1.QualityControl/1.4clusters.pdf", width = 8, height = 6)
VlnPlot(data.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0, ncol = 2, same.y.lims=F) +
  scale_y_continuous(breaks=seq(0, 100, 10)) +
  NoLegend() ## 数据已经过滤过
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "RNA_snn_res.0.2")
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "orig.ident")
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "cellType_low")
ratio.plot(seurat.object = data.merge, id.vars1 = "orig.ident", id.vars2 = "cellType_low", angle = 60)
dev.off()
saveRDS(data.merge,"./1.5data.merge.rds")

pdf(file = "./2.Cluster/2.0clusters_ratio.pdf", width = 4, height = 4)
ratio.plot(seurat.object = data.merge, id.vars1 = "cellType_low", id.vars2 = "orig.ident", angle = 0, 
           color.len = c("#8DA3A6", "#7EB875", "#947D99", "#57A2AC"))#c("#8DA3A6", "#7EB875", "#57A2AC", "#D1B9C5", "#947D99", "#6D4D6E")
dev.off()



data.merge <- readRDS("./1.5data.merge.rds")
dir.create("2.Cluster")
#umap展示
if(F){

  # devtools::install_github("zhanghao-njmu/SCP")
  library(SCP)
  pdf("2.Cluster/2.1sample_umap.pdf", width = 5.5, height = 5)
  CellDimPlot(data.merge, group.by = "orig.ident", reduction = "UMAP", theme_use = "theme_blank", #theme_blank theme_void theme_linedraw theme_light
              title = "COV362",
              palcolor = c("#7EB875", "#B997C7", "#824D99", "#57A2AC", "#4E78C4"),
              #theme_use = ggplot2::theme_classic, 
              #cells.highlight = colnames(pancreas_sub)[pancreas_sub$SubCellType == "Epsilon"],
              #label = TRUE, label_insitu = TRUE, label_repel = TRUE, label_segment_color = "red",
              legend.position  = "right", legend.direction = "vertical",
              theme_args = list(plot.title = element_text(hjust = 0.5, # 标题居中
                  size  = 14, face  = "bold"), legend.title = element_text(size = 12), 
                legend.key.height = unit(0.6, "cm"),  # 图例间距
                legend.spacing.y  = unit(4, "pt"))
              )+
    scale_color_manual(
      name = "orig.ident:",
      values = c("#7EB875", "#B997C7", "#824D99", "#57A2AC", "#4E78C4"),
      guide = guide_legend(
        override.aes = list(size = 4)  # 👈 图例圆点大小
      ))
  dev.off()
  
  
  data.merge$Clusters <- data.merge$cellType_low
  pdf("2.Cluster/2.2cluster_umap.pdf", width = 5.5, height = 5)
  CellDimPlot(data.merge, group.by = "Clusters", reduction = "UMAP", theme_use = "theme_blank", #theme_blank theme_void theme_linedraw theme_light
              title = "COV362",
              palcolor = c("#57A2AC", "#4E78C4", "#824D99", "#D0B541"),
              #theme_use = ggplot2::theme_classic, 
              #cells.highlight = colnames(pancreas_sub)[pancreas_sub$SubCellType == "Epsilon"],
              #label = TRUE, label_insitu = TRUE, label_repel = TRUE, label_segment_color = "red",
              legend.position  = "right", legend.direction = "vertical",
              theme_args = list(plot.title = element_text(hjust = 0.5, # 标题居中
                                                          size  = 14, face  = "bold"), legend.title = element_text(size = 12), 
                                legend.key.height = unit(0.6, "cm"),  # 图例间距
                                legend.spacing.y  = unit(4, "pt")))+
    scale_color_manual(
      name = "Clusters:",
      values = c("#57A2AC", "#4E78C4", "#824D99", "#D0B541"),
      guide = guide_legend(
        override.aes = list(size = 4)  # 👈 图例圆点大小
      ))
  dev.off()
  
  FeatureDimPlot(
    srt = data.merge.pro,
    features = c("SOX17", "WT1", "IFI6", "IFI27"),
    reduction = "UMAP",
    theme_use = "theme_blank",ncol = 4
  )
  
  FeatureDimPlot(
    srt = data.merge,
    features = c("SOX17", "WT1", "IFI6", "SLC1A3"),
    compare_features = TRUE,
    label = TRUE,
    label_insitu = TRUE,
    reduction = "UMAP",
    theme_use = "theme_blank",pt.size = 2
  )
  
  # 气泡图
  ht <- GroupHeatmap(
    srt = data.merge,
    features = c(
      "SOX17", "WT1", "IFI6", "IFI27", "CD44", "KLF4", "SGK1", "CITED2",
      "PARP1", "RAD23A", "CDK1", "ID1", "RPL17", "DAB2","SLC1A3"
    ),
    group.by = c("cellType_low"),
    heatmap_palette = "YlOrRd",
    # cell_annotation = c("Phase", "G2M_score", "Cdh2"),
    cell_annotation_palette = c("Dark2", "Paired", "Paired"),
    #show_row_names = TRUE, 
    row_names_side = "left",
    add_dot = TRUE, 
    add_reticle = TRUE
  )
  print(ht$plot)
  
  # 小提琴
  FeatureStatPlot(
    srt = data.merge,
    group.by = "cellType_low",
    #bg.by = "celltype",
    stat.by = c("SOX17", "WT1", "IFI6", "IFI27", "CD44", "KLF4", "SGK1", "CITED2",
                "PARP1", "RAD23A", "CDK1", "ID1", "RPL17", "DAB2", "SLC1A3"),
    add_box = TRUE,
  )
}

library(scales)
dir.create("d5velocyto")
# 获得每个细胞的UMAP或TSNE坐标，使用 Embeddings函数
write.csv(Embeddings(data.merge, reduction = "umap"), file = "./d5velocyto/cell_embeddings.csv")
# 获取每个细胞的barcode
write.csv(Cells(data.merge), file = "./d5velocyto/cellID_obs.csv", row.names = FALSE)
# 提取每个细胞的样本信息
write.csv(data.merge@meta.data[, 'orig.ident', drop = FALSE], file = "./d5velocyto/cell_sample.csv")
# 提取每个细胞的celltype信息
write.csv(data.merge@meta.data[, 'Clusters', drop = FALSE], file = "./d5velocyto/cell_celltype.csv")
# 获取celltype的颜色信息
hue_pal()(length(levels(data.merge$Clusters)))
# 获取样本的颜色信息
hue_pal()(length(levels(data.merge$orig.ident)))


library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(cowplot)

new_metadata = data.merge@meta.data
labels = as.character(unique(new_metadata$orig.ident))
list_means = list()
for(label in labels) {
  cells = rownames(new_metadata[new_metadata$orig.ident == label, ])
  avg_cells = rowMeans(as.matrix(data.merge@assays$RNA@data[, cells]))
  list_means[[label]] <- avg_cells
}
data.merge_avg_bulk = do.call(cbind, list_means)
data.merge_avg_bulk = data.merge_avg_bulk[, c(3,4,5,1,2)] # 调整顺序

# try with variable genes
var_genes <- VariableFeatures(data.merge)
cormM = cor((data.merge_avg_bulk[var_genes, ]), method="spearman")


library(dendextend)
dend = as.dendrogram(hclust(dist(cormM)), type = "average") # can be used as cluster_rows and columns arg instead of T
dend <- click_rotate(dend, continue = TRUE)

h2 = Heatmap(cormM, show_column_names = T, show_row_dend = T, 
             show_column_dend = F, show_row_names = T, 
             name = "Corr", row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 12),
             col = colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(50),
             cluster_rows = dend, row_dend_reorder = F, cluster_columns = dend,
             column_dend_reorder = F,
             heatmap_width = unit(7, "cm"), heatmap_height = unit(6, "cm"),
             heatmap_legend_param = list(title_gp = gpar(fontsize = 10),
                                         title = "Correlation",
                                         labels_gp = gpar(fontsize = 10),
                                         legend_height = unit(2, "cm")))

pdf("2.Cluster/2.3Correlation.pdf", width = 4, height = 3)
draw(h2,padding = unit(c(1, 1, 6, 1), "mm"))  # 下、左、上、右
grid::grid.text(
  "COV362",
  y = unit(1, "npc") - unit(5, "mm"),
  gp = gpar(fontsize = 14, fontface = "bold")
)
dev.off()



# Head to check
levels(data.merge)



# 计算差异基因
if(F){
  cluster.all.markers <- FindAllMarkers(data.merge, only.pos = TRUE, group.by = "cellType_low", test.use = "MAST", latent.vars = "orig.ident")
  idents <- as.character(levels(data.merge))
  cluster.sig.markers <- cluster.all.markers[which(cluster.all.markers$p_val_adj<0.05),]
  saveFormat <- lapply(idents, function(x){
    index <- which(cluster.sig.markers$cluster == x)
    DEGs <- cluster.sig.markers[index,]
    DEGs.up <- DEGs %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))
    DEGs.down <- DEGs %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
    DEGs <- rbind(DEGs.up, DEGs.down)
    return(DEGs)
  })
  library(openxlsx)
  write.xlsx(saveFormat, file = "2.Cluster/AnnotateCellType/celltype.all.DEGs.xlsx", sheetName = idents, rowNames = F)
  # saveRDS(cellType.all.markers, file = "2.Cluster/AnnotateCellType/cellType.all.DEGs.rds")
  
 
}



new_metadata <- data.merge@meta.data
new_metadata$cellType_low <- factor(new_metadata$cellType_low, 
                                    levels = c("c0", "c1", "c2", "c3"))
#labels = as.character(unique(new_metadata$cellType_low))
labels = as.character(levels(new_metadata$cellType_low))
list_means = list()
for(label in labels) {
  cells = rownames(new_metadata[new_metadata$cellType_low == label, ])
  avg_cells = rowMeans(as.matrix(data.merge@assays$RNA@data[, cells]))
  list_means[[label]] <- avg_cells
}
sc_avg_data = do.call(cbind, list_means)
#sc_avg_data = sc_avg_data[, c(2,3,1)] # 调整顺序

genes_to_show = c("SOX17", "KRT8", "PAX8", "PARP1", "S100A4",
                  "CD24", "CDKN2A", "FN1", "CD44", "LBH",
                  "RAD21", "RAD23A", "CCNB1", "CCNB2", "CCND1", "SMAD3", "TOP2A", 
                  "IFI27", "IFI6", "IFI35", "PARP14", "CYP1B1", "SERPINE1")
# , "PARD6B", 
# "VAMP8", "TM4SF1", "S100A6", "DPYSL3", "CLIC1", "DAZAP2", "CRABP2", "ANXA4", "FN1", # "KRT18", 
# "TGFBI", "SIX1", "VCAM1", "CD44", "CCDC80", "SERPINE1", "CDC42EP3", "AHNAK", "DAB2", "SMAD3",
# "GDF15", "ID1", "CYP1B1", "RPL36A", "NQO1", "PRR11", "DDIT4", "GPX4", "TOP2A"

avg_expression = DotPlot(data.merge, features = genes_to_show)
avg_expression = avg_expression$data
avg_expression$features.plot <- factor(avg_expression$features.plot, levels = (rownames(avg_expression)))


#avg_expression$id <- factor(avg_expression$id, levels = c(2,1,0))
avg_expression$id <- factor(avg_expression$id, levels = c("c0", "c1", "c2", "c3"))

# 自定义颜色映射
id_colors <- c("#57A2AC", "#4E78C4", "#824D99", "#D0B541")

# 根据你的要求对基因进行分组颜色映射
gene_colors <- c(rep("#57A2AC", 5), rep("#4E78C4", 5), rep("#824D99", 7),
                 rep("#D0B541", 6))


# 绘制图形
pdf("2.Cluster/2.4cluster_markers.pdf", width = 7, height = 3.5)
ggplot(avg_expression,
       aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled)) + 
  theme_bw() +
  geom_point(shape = 21) +
  scale_fill_gradient2(low = "#053061",
                       mid = "#eeeeee",
                       high = "#67001F",
                       breaks = c(-1, 0, 1)) +
  scale_color_manual(values = id_colors) +  # 应用id颜色映射
  ylab("") +
  xlab("") +
  labs(size = "Percent of cell", fill = "Scaled expression") +
  theme(axis.text = element_text(color = "black"), 
        text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, 
                                   color = gene_colors), # 设置x轴字体颜色
        # axis.title.x = element_text(color = "#8DA3A6", size = 14), # 设置x轴标签字体颜色
        # axis.title.y = element_text(color = "#DD4E70", size = 14), # 设置y轴标签字体颜色
        axis.text.y = element_text(size = 12, color = id_colors), # 设置y轴字体颜色
        legend.key.size = unit(0.4, 'cm'), legend.position = "top")
dev.off()


### 富集分析

## 根据细胞类型重新计算差异基因
#### Cell type specific gene
# Idents(data.merge) <- data.merge$cellType_low
# idents <- as.character(levels(data.merge))
cellType.all.markers <- FindAllMarkers(data.merge, 
                                       group.by = "cellType_low", 
                                       logfc.threshold = 0, 
                                       min.pct = 0.2, 
                                       test.use = "MAST", 
                                       latent.vars = "orig.ident")

saveRDS(cellType.all.markers, file = "2.Cluster/AnnotateCellType/celltype.all.DEGs.enrich.rds")

cellType.all.DEGs <- readRDS("2.Cluster/AnnotateCellType/celltype.all.DEGs.enrich.rds")
Tumor.DEGs <- cellType.all.DEGs[cellType.all.DEGs$cluster=="R5",]
# geneList <- Tumor.DEGs$avg_log2FC
# names(geneList) <- Tumor.DEGs$gene # 5439 genes
geneList <- data.frame(
  gene = Tumor.DEGs$gene,
  avg_log2FC = Tumor.DEGs$avg_log2FC
)