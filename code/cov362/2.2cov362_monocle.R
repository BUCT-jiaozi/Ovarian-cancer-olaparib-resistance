library(Seurat)
library(ggplot2)
library(tidyverse)
rm(list = ls());gc()
set.seed(101)
plan("multisession", workers = 2) 
options(future.globals.maxSize = 10000 * 1024^2) # set 50G RAM
setwd('/home/Rstudio/yhb/Data/cov362_sc/')

dir <- "./"
dir.create("4.monocle")
setwd("./4.monocle/")
# data.merge <- readRDS("../1.5data.merge.pro.rds")
# DimPlot(data.merge, reduction = "umap", label = T, group.by = 'cellType_low')

# data.merge$orig.ident <- data.merge$Patient_ID

dir.create("./plot_out")


RNA <- readRDS("../1.5data.merge.rds")
RNA$Neu_cluster_annotation <- factor(RNA$cellType_low,
                                     levels=c("c0","c1", "c2", "c3"))
####### trajactory analysis ######
library(monocle)
data <- as(as.matrix(RNA@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = RNA@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
# Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())

monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
## 
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(monocle_cds),num_cells_expressed >= 10))

# Neu_DEG_genes <- differentialGeneTest(monocle_cds[expressed_genes,], fullModelFormulaStr = '~Neu_cluster_annotation', cores = 12)
# write.csv(Neu_DEG_genes,"./Neu_monocle_DEG2.csv",quote=F)

Neu_DEG_genes <- read.csv("./Neu_monocle_DEG2.csv",row.names = 1)
ordering_genes <- row.names(subset(Neu_DEG_genes, qval < 0.001))
length(ordering_genes)
gene_id <- row.names(Neu_DEG_genes)[order(Neu_DEG_genes$qval)][1:1500] # 1500 1000 1200 2300

##
monocle_cds <- setOrderingFilter(monocle_cds, gene_id)
monocle_cds <- reduceDimension(
  monocle_cds,
  max_components = 2,
  norm_method = "vstExprs", #c("log", "vstExprs", "none")
  #num_dim = 3, 
  #residualModelFormulaStr = "~orig.ident", # 指定在聚类之前从数据减去的效果，orig.ident即去除样本影响，有利于减少分支数
  reduction_method = 'DDRTree') #c("DDRTree", "ICA", "tSNE", "SimplePPT", "L1-graph", "SGL-tree"),


monocle_cds <- orderCells(monocle_cds)

ClusterName_color_panel <- c("c0" = "#1c8041", "c1" = "#fddf8b", "c2" = "#c77cff",
                             "c3" = "#4E78C4") #, "R4"="#c77cff", "R5"= "#DD4E70"

library(cowplot)

# 选择根节点，重新排序
#monocle_cds <- orderCells(monocle_cds, root_state = 2) # 

#pdf(paste0(dir,"plot_out/2.2Neu_Pseudotime_new.pdf"),width = 6,height = 4)
plot_cell_trajectory(monocle_cds,color_by = "State")
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",show_branch_points=F)+theme_cowplot() + theme(legend.position = "top")#+ scale_color_viridis(option = "plasma")
plot_cell_trajectory(monocle_cds, color_by = "Neu_cluster_annotation",cell_size = 1,show_branch_points=F)+theme_cowplot()+
  scale_color_manual(name = "", values = ClusterName_color_panel)+guides(color = guide_legend(override.aes = list(size=4)))+
  theme(axis.text.y = element_text(size=12),axis.text.x = element_text(size=12),axis.title.y = element_text(size=12),axis.title.x = element_text(size=12),axis.ticks.length = unit(0.2,"cm"))
#dev.off()

pdf(paste0(dir,"plot_out/4.1Pseudotime_cell.pdf"),width = 4.5,height = 3.5)
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime",show_branch_points=F)+theme_cowplot() + theme(legend.position = "top")#+ scale_color_viridis(option = "plasma")
dev.off()

pdf(paste0(dir,"plot_out/4.2Pseudotime_cluster.pdf"),width = 6,height = 4)
plot_cell_trajectory(monocle_cds, color_by = "Neu_cluster_annotation",cell_size = 1,show_branch_points=F)+theme_cowplot()+
  scale_color_manual(name = "", values = ClusterName_color_panel)+guides(color = guide_legend(override.aes = list(size=4)))+
  theme(axis.text.y = element_text(size=12),axis.text.x = element_text(size=12),axis.title.y = element_text(size=12),axis.title.x = element_text(size=12),axis.ticks.length = unit(0.2,"cm"))
dev.off()

saveRDS(monocle_cds,"./Neu_monocle_new.rds")


# 展示基因在拟时序过程的表达变化
plot_genes <- c("KRT8","RAD23A","CD44") # "SOX17","WT1","ID1","SMAD3","RAD21", "RAD23B", "DDIT4"
cds_subset <- monocle_cds[plot_genes,] 


source("../../OC_data/kura_new_fq/code/Functions/Monocle_plot_gene.R")
monocle_theme_opts <- function () {
  theme(strip.background = element_rect(colour = "white", 
                                        fill = "white")) + theme(panel.border = element_blank()) + 
    theme(axis.line.x = element_line(size = 0.25, color = "black")) + 
    theme(axis.line.y = element_line(size = 0.25, color = "black")) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill = "white")) + 
    theme(legend.key = element_blank())
}
pdf(paste0(dir,"plot_out/4.3pseudo_CXCR4_new.pdf"),width = 5.5,height = 2)
My_plot_gene_pseudotime(cds_subset, color_by = "Neu_cluster_annotation",ncol=3)
dev.off()


# 展示在拟时序过程中每种细胞类型比例的峰峦图
### density
plotdf=pData(monocle_cds)
library(ggridges)
ggplot(plotdf, aes(x=Pseudotime,y=Neu_cluster_annotation,fill=Neu_cluster_annotation))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank())+
  scale_fill_manual(
    name = "Cluster",
    values = c("#57A2AC", "#4E78C4", "#824D99", "#D0B541")
    )
ggsave("plot_out/4.4monocle_cluster_density.pdf",width = 6,height = 2)




monocle_cds <- readRDS("./Neu_monocle_new.rds")
# 会返回每个基因的显著性，显著的基因就是那些随不同branch变化的基因
#### branch heatmap
# BEAM_res <- BEAM(monocle_cds, branch_point = 1, cores = 1)
BEAM_res <- BEAM(monocle_cds, branch_point = 1, cores = 7, progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
saveRDS(BEAM_res, file = "BEAM_res.rds")

# 热图
tmp1=plot_genes_branched_heatmap(monocle_cds[row.names(subset(BEAM_res,
                                                              qval < 1e-6)),],
                                 branch_point = 1,
                                 num_clusters = 3, # 这些基因被分成几个group
                                 cores = 1,
                                 use_gene_short_name = T,
                                 show_rownames = F,
                                 return_heatmap = T)
pdf(paste0(dir,"plot_out/4.5.1monocle_brach_heatmap.pdf"),width=5,height=7)
tmp1$ph_res
dev.off()
## GO-BP analysis 
gene_group <- tmp1$annotation_row
gene_group$gene <- rownames(gene_group)
library(clusterProfiler)
library(org.Hs.eg.db)
allcluster_go <- data.frame()
for(i in unique(gene_group$Cluster)){
  samll_gene_group= filter(gene_group,gene_group$Cluster==i)
  df_name = bitr(samll_gene_group$gene,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene=unique(df_name$ENTREZID),
                 OrgDb=org.Hs.eg.db,
                 keyType="ENTREZID",
                 ont="BP",
                 pAdjustMethod="BH",
                 pvalueCutoff=0.05,
                 qvalueCutoff=0.05,
                 readable=T)
  go_res=go@result
  if(dim(go_res)[1] !=0){
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
write.csv(allcluster_go,"./plot_out/4.5.1monocle_brach_BP.csv")
saveRDS(allcluster_go, "./allcluster_go.rds")


Neu_DEG_genes <- read.csv("./Neu_monocle_DEG2.csv",row.names = 1)
ordering_genes <- row.names(subset(Neu_DEG_genes, qval < 0.001))
length(ordering_genes)
gene_id <- row.names(Neu_DEG_genes)[order(Neu_DEG_genes$qval)][1:1500] # 1500 1000 1200 2300


###### pseudo-time heatmap
heat_cds <- monocle_cds[gene_id ,]
tmp2=plot_pseudotime_heatmap(heat_cds,
                             num_clusters = 3,
                             cores = 1,
                             show_rownames = F,
                             return_heatmap = T)
pdf(paste0(dir,"plot_out/4.5.2monocle_pseudotime_heatmap.pdf"),width=5,height=6)
print(tmp2)
dev.off()
## GO-BP analysis 
gene_group <- data.frame(Cluster = factor(cutree(tmp2$tree_row,3)))
gene_group$gene <- rownames(gene_group)
head(gene_group)
library(clusterProfiler)
library(org.Hs.eg.db)
allcluster_go <- data.frame()
for(i in unique(gene_group$Cluster)){
  samll_gene_group= filter(gene_group,gene_group$Cluster==i)
  df_name = bitr(samll_gene_group$gene,fromType="SYMBOL",toType=c("ENTREZID"),OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene=unique(df_name$ENTREZID),
                 OrgDb=org.Hs.eg.db,
                 keyType="ENTREZID",
                 ont="BP",
                 pAdjustMethod="BH",
                 pvalueCutoff=0.05,
                 qvalueCutoff=0.05,
                 readable=T)
  go_res=go@result
  if(dim(go_res)[1] !=0){
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
write.csv(allcluster_go,"./plot_out/4.5.2monocle_pseudotime_BP.csv")



