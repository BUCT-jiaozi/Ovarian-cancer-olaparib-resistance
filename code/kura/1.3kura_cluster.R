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

setwd('/home/Rstudio/yhb/Data/OC_data/kura_new_fq/kb_fitle_R')
dir <- "./"
source(file = "../code/Functions/ratio.plot.R")

samples_singlet <- readRDS('./1.2data_singlecell.rds')

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

dir.create('./5.new_result')
pdf(file = "./5.new_result/1.1count.feature.mt.pdf", width = 8, height = 6)
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



## 质控，初步标准化分群
# 质控
data.merge <- subset(x = data.merge,subset = nFeature_RNA > 500 & percent.mt < 15 & nCount_RNA < 20000 & nFeature_RNA < 5000 & percent.ribo < 30)
dim(data.merge)

data.merge <- NormalizeData(object = data.merge, normalization.method= "LogNormalize", verbose = FALSE)
data.merge <- FindVariableFeatures(object = data.merge,selection.method = "vst",
                                   nfeatures = 3000, verbose = FALSE)
data.merge <- ScaleData(object = data.merge, vars.to.regress = c("nCount_RNA", "percent.mt", "percent.ribo"))

set.resolutions <- seq(0.1, 0.45, by = 0.05)

data.merge <- RunPCA(data.merge, npcs = 100, verbose = FALSE)
ElbowPlot(object = data.merge, ndims = 100)
data.merge <- FindNeighbors(object = data.merge, dims = 1:40, verbose = FALSE)
data.merge <- FindClusters(object = data.merge, resolution = set.resolutions, verbose = FALSE) 
data.merge <- FindClusters(object = data.merge, resolution = 1.2, verbose = FALSE) 
clustree(data.merge)
data.merge <- RunUMAP(data.merge, dims = 1:40)

data.merge@meta.data$orig.ident <- factor(data.merge@meta.data$orig.ident, 
                                          levels = c("C", "T1", "T2.5", "T5", "T10", "T20", "T40", "T80", "T160", "T320"))

DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "RNA_snn_res.1.2")
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "orig.ident")


### 修改cluster编号顺序
if(F){
  Idents(data.merge) <- data.merge$RNA_snn_res.1.2
  # rename clusters to put in order of more resistant
  data.merge <- RenameIdents(object = data.merge,
                             "0" = "R3","1" = "R1","2" = "R4","3" = "R5","4" = "R2","5" = "R3",
                             "6" = "C","7" = "R3","8" = "R4","9" = "R5","10" = "R4","11" = "R1",
                             "12" = "R2","13" = "R5","14" = "R4","15" = "R1","16" = "R5","17" = "R2",
                             "18" = "R4","19" = "R1","20" = "R1","21" = "R2"
                             )
  # assign new order
  data.merge@active.ident <- factor(x = data.merge@active.ident,
                                    levels = c("C", "R1", "R2", "R3", "R4", "R5"))
  
  data.merge@meta.data$cellType_low <- data.merge@active.ident
  Idents(data.merge) <- data.merge@meta.data$cellType_low
  data.merge@meta.data$cellType_low <- factor(data.merge@meta.data$cellType_low, 
                                              levels = c("C", "R1", "R2", "R3", "R4", "R5"))
  
  DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "cellType_low")
}


source(file = "../code/Functions/ratio.plot.R")

pdf(file = "./5.new_result/1.4clusters.pdf", width = 8, height = 6)
VlnPlot(data.merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0, ncol = 2, same.y.lims=F) +
  scale_y_continuous(breaks=seq(0, 100, 10)) +
  NoLegend() ## 数据已经过滤过
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "RNA_snn_res.1.2")
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "orig.ident")
DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "cellType_low")
#ratio.plot(seurat.object = data.merge, id.vars1 = "orig.ident", id.vars2 = "cellType_low", angle = 0, color.len = c("#8DA3A6", "#B997C7","#824D99","#4E78C4","#57A2AC","#7EB875","#D0B541", "#E67F33", "#DD4E70", "#CE2220"))
ratio.plot(seurat.object = data.merge, id.vars1 = "cellType_low", id.vars2 = "orig.ident", angle = 0, color.len = c("#8DA3A6", "#7EB875", "#57A2AC", "#D1B9C5", "#947D99", "#6D4D6E"))
dev.off()


saveRDS(data.merge,"./1.5data.merge.pro.rds")

data.merge <- readRDS("./1.5data.merge.pro.rds")

dir.create("2.Cluster")

pdf(file = "./2.Cluster/2.0clusters_ratio.pdf", width = 6, height = 4)
ratio.plot(seurat.object = data.merge, id.vars1 = "cellType_low", id.vars2 = "orig.ident", angle = 0, color.len = c("#8DA3A6", "#7EB875", "#57A2AC", "#D1B9C5", "#947D99", "#6D4D6E"))
#ratio.plot(seurat.object = data.merge, id.vars1 = "orig.ident", id.vars2 = "cellType_low", angle = 0, color.len = c("#8DA3A6", "#B997C7","#824D99","#4E78C4","#57A2AC","#7EB875","#D0B541", "#E67F33", "#DD4E70", "#CE2220"))
dev.off()

#umap展示
if(F){
  ######################
  # Plot Umap Sample   #
  ######################
  # devtools::install_github("zhanghao-njmu/SCP")
  library(SCP)
  pdf("2.Cluster/2.1sample_umap.pdf", width = 5.5, height = 5)
  CellDimPlot(data.merge, group.by = "orig.ident", reduction = "UMAP", theme_use = "theme_blank", #theme_blank theme_void theme_linedraw theme_light
              title = "Kuramochi",
              palcolor = c("#8DA3A6", "#B997C7","#824D99","#4E78C4","#57A2AC","#7EB875","#D0B541", "#E67F33", "#DD4E70", "#CE2220"),
              #theme_use = ggplot2::theme_classic, 
              #cells.highlight = colnames(pancreas_sub)[pancreas_sub$SubCellType == "Epsilon"],
              #label = TRUE, label_insitu = TRUE, label_repel = TRUE, label_segment_color = "red",
              legend.position  = "right", legend.direction = "vertical",
              theme_args = list(plot.title = element_text(hjust = 0.5, # 标题居中（关键）
                                                          size  = 14, face  = "bold"), legend.title = element_text(size = 12),
                                legend.key.height = unit(0.6, "cm"),  # 图例间距
                                legend.spacing.y  = unit(4, "pt")))+
    scale_color_manual(
      name = "orig.ident:",
      values = c("#8DA3A6", "#B997C7","#824D99","#4E78C4","#57A2AC","#7EB875","#D0B541", "#E67F33", "#DD4E70", "#CE2220"),
      guide = guide_legend(
        override.aes = list(size = 4)  # 👈 图例圆点大小
      ))
  dev.off()
  
  
  data.merge$Clusters <- data.merge$cellType_low
  pdf("2.Cluster/2.2cluster_umap.pdf", width = 5.5, height = 5)
  CellDimPlot(data.merge, group.by = "Clusters", reduction = "UMAP", theme_use = "theme_blank", #theme_blank theme_void theme_linedraw theme_light
              title = "Kuramochi",
              palcolor = c("#8DA3A6","#B997C7","#57A2AC", "#7EB875", "#4E78C4", "#DD4E70"),
              #theme_use = ggplot2::theme_classic, 
              #cells.highlight = colnames(pancreas_sub)[pancreas_sub$SubCellType == "Epsilon"],
              #label = TRUE, label_insitu = TRUE, label_repel = TRUE, label_segment_color = "red",
              legend.position  = "right", legend.direction = "vertical",
              theme_args = list(plot.title = element_text(hjust = 0.5, # 标题居中（关键）
                                                          size  = 14, face  = "bold"), legend.title = element_text(size = 12), 
                                legend.key.height = unit(0.6, "cm"),  # 图例间距
                                legend.spacing.y  = unit(4, "pt"))
  )
  dev.off()
  
  FeatureDimPlot(
    srt = RNA,
    features = c("PARP1", "RAD51C", "FANCI", "RAD51"),
    reduction = "UMAP",
    theme_use = "theme_blank",ncol = 4
  )
  
  FeatureDimPlot(
    srt = data.merge.pro,
    features = c("SOX17", "WT1", "IFI6", "IFI27"),
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
      "SOX17", "RAD51", "RAD51C", "RAD21", "RAD23B", "FANCA","FANCB", "FANCC", "FANCG", "FANCD2", "FANCI",
      "PARP1", "RAD23A", "RPA1", "RPA2", "RPA3", "CCNB1", "CCND2", "BRCA1", "SLC1A3"
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
    stat.by = c("SOX17", "RAD51", "RAD51C", "RAD21", "RAD23B", "FANCA","FANCB", "FANCC", "FANCG", "FANCD2", "FANCI",
                "PARP1", "RAD23A", "RPA1", "RPA2", "RPA3", "CCNB1", "CCND2", "BRCA1", "SLC1A3"),
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
#data.merge_avg_bulk = data.merge_avg_bulk[, c(3,4,5,1,2)] # 调整顺序

# try with variable genes
var_genes <- VariableFeatures(data.merge)
cormM = cor((data.merge_avg_bulk[var_genes, ]), method="spearman")


library(dendextend)
dend = as.dendrogram(hclust(dist(cormM)), type = "average") # can be used as cluster_rows and columns arg instead of T
#dend <- click_rotate(dend, continue = TRUE)
desired_order <- c("C", "T1", "T2.5", "T5", "T10", "T20", "T40", "T80", "T160", "T320")
dend <- rotate(dend, order = desired_order)

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
  "Kuramochi",
  y = unit(1, "npc") - unit(3, "mm"),
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
  
  
  ## 火山图
  ################ 亚群差异分析
  table(data.merge$cellType_low)

  # sce.markers <- bind_rows(R1_DEGs, R2_DEGs, R3_DEGs, R4_DEGs, R5_DEGs)
  sce.markers <- FindAllMarkers(data.merge, only.pos = F,return.thresh = 0.01,logfc.threshold = 0.3,min.pct = 0.2) #logfc.threshold 调整中间留白大小
  head(sce.markers)

  ## 添加显著性 红色 p_val_adj < 0.01,  黑色 p_val_adj >= 0.01
  sce.markers$sig <- if_else(sce.markers$p_val_adj < 0.01, "P_adj < 0.01", "P_adj >= 0.01")
  sce.markers$reg <- if_else(sce.markers$avg_log2FC > 0, "up", "down")
  table(sce.markers$sig)
  table(sce.markers$cluster, sce.markers$sig)

  ## 挑选 红色 p_val_adj < 0.01 中log2FC的top10
  # top10 <- sce.markers %>%
  #   filter(sig =="P_adj < 0.01") %>%
  #   group_by(cluster,reg) %>%
  #   top_n(abs(avg_log2FC),n = 5) %>%
  #   ungroup()
  
  if(F){
    up_genes_manual <- list(
      C  = c("SOX17", "WT1", "KRT8", "CD24", "CDKN2A"),
      R1 = c("S100A4", "PAX8", "IFI27", "IFI6", "CCND1"),
      R2 = c("COL3A1", "SERPINE2", "LBH", "HMGA2", "HES1"),
      R3 = c("TM4SF1", "S100A6", "DPYSL3", "DAZAP2", "FN1"),
      R4 = c("SIX1", "VCAM1", "CD44", "CCDC80", "AHNAK"),
      R5 = c("GDF15", "CYP1B1", "RPL36A", "DDIT4", "TOP2A"))
    down_genes <- sce.markers %>%
      filter(sig == "P_adj < 0.01", reg == "down") %>%
      group_by(cluster) %>%
      slice_max(order_by = abs(avg_log2FC), n = 5) %>%
      ungroup()
    up_genes <- sce.markers %>%
      filter(sig == "P_adj < 0.01", reg == "up") %>%
      rowwise() %>%
      filter(gene %in% up_genes_manual[[cluster]]) %>%
      ungroup()
    top10 <- bind_rows(down_genes, up_genes)
  }
  
  head(top10)

  ## 绘图整体
  # 绘制每个Cluster 的散点火山图：
  library(ggrepel)

  ## 灰色背景柱子
  #根据图p中log2FC区间确定背景柱长度：
  top_log2FC <- sce.markers %>%
    group_by(cluster) %>%       # 按cluster分组
    slice_max(avg_log2FC, n = 1) %>%# 在每个分组中选择log2FC最大的值
    ungroup()                   # 取消分组

  down_log2FC <- sce.markers %>%
    group_by(cluster) %>%       # 按cluster分组
    slice_min(avg_log2FC, n = 1) %>%# 在每个分组中选择log2FC最大的值
    ungroup()                   # 取消分组

  dfbar <- data.frame(cluster=unique(sce.markers$cluster), up=top_log2FC$avg_log2FC, down=down_log2FC$avg_log2FC)

  ## 绘制背景柱：这里注意要将背景画在底部，不然放在点图后面会遮住点
  p <- ggplot() +
    geom_col(data = dfbar, mapping = aes(x = cluster,y = up), fill = "#efefef") +
    geom_col(data = dfbar,mapping = aes(x = cluster,y = down), fill = "#efefef") +
    geom_jitter(data = sce.markers, aes(x = cluster, y = avg_log2FC, color = sig),
                size = 0.6, width =0.4) +
    geom_jitter(data = top10, aes(x = cluster, y = avg_log2FC, color = sig),
                size = 1, width =0.4) + ## top10的点，大小突出一下
    scale_color_manual(name=NULL, values = c("#C42F40","gray")) + ## 点的颜色调整 CE2220
    geom_text_repel(data=top10, aes(x=cluster,y=avg_log2FC,label=gene), force = 1.2,
                    arrow = arrow(length = unit(0.008, "npc"),type = "open", ends = "last") ) ## 添加文字标签

  p

  ## 添加cluster方框
  ## 方框的高度为前面 log2FC的阈值的两倍，好看一点可以*0.8
  # 添加X轴的cluster色块标签：
  dfcol <- data.frame(x= unique(sce.markers$cluster), y=0, label=unique(sce.markers$cluster),
                      labelcol = c("white",rep("black",2),"white",rep("black",2))
  )
  dfcol
  mycol <- c("#8DA3A6","#B997C7","#57A2AC", "#7EB875", "#4E78C4", "#DD4E70")

  p2 <- p +
    geom_tile(data = dfcol, aes(x=x,y=y), height=0.3 * 2 * 0.8, color = "black", fill = mycol,  show.legend = F) +
    geom_text(data=dfcol, aes(x=x,y=y,label=label), size =6, color = dfcol$labelcol) ## 添加方框中的cluster标签

  p2

  ## 修改其他主题
  p3 <- p2 +
    labs(x="Clusters",y="average logFC") +
    theme_minimal()+
    theme( axis.title = element_text(size = 13, color = "black",face = "bold"),
           axis.line.y = element_line(color = "black",size = 1.2),
           axis.line.x = element_blank(),
           axis.text.x = element_blank(),
           panel.grid = element_blank(),
           legend.position = "top",
           legend.direction = "vertical",
           legend.justification = c(1,0),
           legend.text = element_text(size = 15)
    ) +
    guides(color=guide_legend(override.aes = list(size=4.5))) ## 修改图例中图标的大小

  p3

  ## 保存
  ggsave(filename = "./2.Cluster/2.5cluster_heatmap.pdf",width = 8, height = 6,plot = p3)

}



new_metadata <- data.merge@meta.data
new_metadata$cellType_low <- factor(new_metadata$cellType_low, 
                                                levels = c("C", "R1", "R2", "R3", "R4", "R5"))
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

genes_to_show = c("SOX17", "WT1", "KRT8", "S100A16", "CD24", "PARP1", "CDKN2A", 
                  "S100A4", "PAX8", "IFI27", "IFI6", "BCAM", "CCND1", 
                  "COL3A1", "SERPINE2", "LBH", "HMGA2", "HES1", "SLIT2", "PARD6B", "RAD21", "RAD23B",
                  "VAMP8", "TM4SF1", "S100A6", "DPYSL3", "CLIC1", "DAZAP2", "CRABP2", "ANXA4", "FN1", # "KRT18", 
                  "SIX1", "VCAM1", "CD44", "CCDC80", "CDC42EP3", "AHNAK", "DAB2", "SMAD3",
                  "GDF15", "ID1", "CYP1B1", "EPRS1", "NQO1", "PRR11", "DDIT4", "GPX4", "TOP2A")


avg_expression = DotPlot(data.merge, features = genes_to_show)
avg_expression = avg_expression$data
avg_expression$features.plot <- factor(avg_expression$features.plot, levels = (rownames(avg_expression)))


#avg_expression$id <- factor(avg_expression$id, levels = c(2,1,0))
avg_expression$id <- factor(avg_expression$id, levels = c("C", "R1", "R2", "R3", "R4", "R5"))

# 自定义颜色映射
id_colors <- c("#8DA3A6","#B997C7","#57A2AC", "#7EB875", "#4E78C4", "#DD4E70")

# 根据你的要求对基因进行分组颜色映射
gene_colors <- c(rep("#8DA3A6", 7), rep("#B997C7", 6), rep("#57A2AC", 9),
                 rep("#7EB875", 9), rep("#4E78C4", 8), rep("#DD4E70", 9))


# 绘制图形
pdf("2.Cluster/2.4cluster_markers.pdf", width = 10, height = 4)
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

save_plot("cov362_markers_bubble.pdf", cov362_markers_bubble, base_height = 3, base_width = 9.5)




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
source(file = "../code/Functions/clusterProfiler.enricher.R")
dir.create('./2.Cluster/AnnotateCellType/GSEA')
pdf("2.Cluster/AnnotateCellType/GSEA/MsigDB_H.pdf")
gProfiler.res_H <- cluterProfiler.enricher(gene = geneList$gene, geneType = "SYMBOL", db.type = "MsigDB", MsigDB.category = list(H = c("All")), 
                                           title ="MsigDB_H", saveDir = paste0(getwd(), "/2.Cluster/AnnotateCellType/GSEA"))
dev.off()

gProfiler.res_C2 <- cluterProfiler.enricher(gene = geneList$gene, geneType = "SYMBOL", db.type = "MsigDB", MsigDB.category = list(C2 = c("REACTOME")), 
                                            title ="MsigDB_C2", saveDir = paste0(getwd(), "/2.Cluster/AnnotateCellType/GSEA"))

gProfiler.res_GO <- cluterProfiler.enricher(gene = geneList$gene, geneType = "SYMBOL", db.type = "GO", GO.ont = "ALL", simplify = F, 
                                            title ="GO", saveDir = paste0(getwd(), "/2.Cluster/AnnotateCellType/GSEA"))

gProfiler.res_kegg <- cluterProfiler.enricher(gene = geneList$gene, geneType = "SYMBOL", db.type = "KEGG",
                                              title ="kegg", saveDir = paste0(getwd(), "/2.Cluster/AnnotateCellType/GSEA"))



