library(devtools)
library(DESeq2)
library(magrittr)
library(gmodels)
library(pheatmap)
library(EnhancedVolcano)
library(MASS) 
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(ReactomePA)
library(ggridges)
library(cowplot)
library(patchwork)
library(org.Hs.eg.db)  #人 hsa
library(ggpubr)
library(ggthemes)


###data division into groups###
rm(list = ls());gc()
setwd('C:/Data/work_Data/A2780_RNAseq/20250531A2780_GSE153867/')
raw_data <- read.table("./all_exon_A2780.txt",header = T,skip = 0,check.names = F) ##Need to be modified

data <- as.data.frame(raw_data[,6:ncol(raw_data)]) #Contains Length
rownames(data) <- raw_data[,1]



###Separate analysis###

resistant_ids <- paste0("SRR134899", 81:88)
#wt_ids <- paste0("SRR134899", 89:96)


#Create database
ID <- colnames(data)[2:ncol(data)]
database <- as.data.frame(ID)
#database$Sample_Type <- ifelse(grepl("Con", database$ID), "Con", "Ola+Ida") ##Need to be modified Olaparib-resistant
#database$Sample_Type <- ifelse(database$ID %in% resistant_ids, "Olaparib_resistant", "WT")
database$Sample_Type <- ifelse(database$ID %in% resistant_ids, "OR", "WT")
write.csv(database, "./database.csv")

kb <- data$Length/1000
countData <- as.matrix(data[,2:ncol(data)])
#rownames(countData) <- data[,1]
rownames(countData) <- rownames(data)



rpk <- countData / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)#
write.csv(tpm, "tpm.csv")

colData <- data.frame(sample = colnames(countData), condition = database$Sample_Type)
dds <- DESeqDataSetFromMatrix(countData, colData, design= ~ condition)
dds <- dds[rowSums(counts(dds)) > 1, ] 
dds <- DESeq(dds)##标准化结果
res <- results(dds,contrast = c("condition", "OR", "WT"))  ###
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)

gene_df <- bitr(resdata$Row.names, fromType = "ENSEMBL", toType = c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db) # 修改
resdata_symbol <- merge(resdata,gene_df,by.x="Row.names",by.y="ENSEMBL")
write.csv(resdata_symbol,"resdata_symbol.csv")


res_de <- subset(resdata_symbol, resdata_symbol$padj < 0.05, select=c( 'log2FoldChange', 'pvalue','padj','SYMBOL','ENTREZID'))
res_de_up <- subset(res_de, res_de$log2FoldChange >= 1)
res_de_dw <- subset(res_de, res_de$log2FoldChange <= (-1)*1)

res_de_up_id <- as.vector(res_de_up$SYMBOL)
res_up <- resdata_symbol[resdata_symbol$SYMBOL%in%res_de_up_id,]
write.csv(res_up,"diff_up.csv")
res_de_dw_id <- as.vector(res_de_dw$SYMBOL)
res_dw <- resdata_symbol[resdata_symbol$SYMBOL%in%res_de_dw_id,]
write.csv(res_dw,"diff_dw.csv")
res_de_top <- c(res_de_up_id, res_de_dw_id)
res_diff <- resdata_symbol[resdata_symbol$SYMBOL%in%res_de_top,]
write.csv(res_diff,"diff_gens.csv",row.names = F)

top10gene <- res_up %>% top_n(n = 10,wt = log2FoldChange) %>% pull(SYMBOL)  
bottom10gene <- res_dw %>% top_n(n = -10,wt = log2FoldChange) %>% pull(SYMBOL)
gene20 <- c(top10gene, bottom10gene)


### plot1
label_mi <- resdata_symbol$SYMBOL
P <- EnhancedVolcano(resdata_symbol,lab = label_mi,x = 'log2FoldChange',y = 'padj',labSize = 4,title = "",subtitle = "",pCutoff = 0.05,FCcutoff = 1,xlim = c(-7, 7),ylim = c(0,15),border = 'full',colAlpha = 1,cutoffLineType = 'twodash',   xlab = bquote(~Log[2]~ 'fold change'),ylab = bquote(~-Log[10]~adjusted~italic(P)))
plot1 <- P + ggplot2::coord_cartesian(xlim=c(-14, 14),ylim = c(0,300))+ggplot2::scale_x_continuous(breaks=seq(-12,12,2))+ggplot2::scale_y_continuous(breaks=seq(0,325,50))

pdf(file = 'Volcanomap.pdf',width = 12,height = 10)
print(plot1)
dev.off()

dir.create("./plot_out")

## plot2
# 另一个火山图

gene_vol <- c("CCNB1", "RAD51", "RPA3",  
              "FANCI", "BRCA1", "BRCA2", "PARP1", "ID1",
              "FOXR1", "TRPC4", "VCAN", "CD4", "OSM", "MAPK4")
dw_gene <- c("FOXR1", "TRPC4", "VCAN", "CD4", "OSM", "MAPK4")
up_gene <- c("CCNB1", "RAD51", "RPA3","FANCI", "BRCA1", "BRCA2") #, "PARP1", "ID1"

for (i in 1:nrow(resdata_symbol)) {
  if ( !(is.na(resdata_symbol[i,'log2FoldChange'])) & resdata_symbol[i,'log2FoldChange'] > 1 & !(is.na(resdata_symbol[i,'pvalue'])) & resdata_symbol[i,'pvalue'] < 0.05) resdata_symbol[i,'select_change'] <- 'Up'
  else if ( !(is.na(resdata_symbol[i,'log2FoldChange'])) & resdata_symbol[i,'log2FoldChange'] < -1 & !(is.na(resdata_symbol[i,'pvalue'])) & resdata_symbol[i,'pvalue'] < 0.05) resdata_symbol[i,'select_change'] <- 'Down'
  else resdata_symbol[i,'select_change'] <- 'None' 
  resdata_symbol[i,'select'] <- paste(resdata_symbol[i,'select_change'])
}
P <- ggplot(resdata_symbol,aes(log2FoldChange, -log(padj, 10))) +
  geom_point(aes(color = select_change), alpha = 0.9) +
  #scale_color_manual(values = c('blue2','gray30', 'red2')) +
  #scale_color_manual(values = c("steelblue","gray","brown")) +
  scale_color_manual(values = c("#4E78C4","gray","#CE2220")) +
  scale_x_continuous (limits = c (-4, 4), breaks = seq (-4, 4, 2))+
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent',size = 1.5)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent',size = 1.5)) +
  theme(legend.title = element_blank(), legend.key = element_rect(fill = 'transparent', color = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
  theme(legend.position = c(0.92,0.9)) + # Up Down标签位置
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.6) + 
  geom_hline(yintercept = -log(0.25, 10), color = 'gray', size = 0.6) +
  labs(x = 'log2(FoldChange)', y = '-log10(padj)') +
  theme(text = element_text(size = 16)) 
P
plot01 <- P + ggplot2::coord_cartesian(xlim=c(-16, 16),ylim = c(0,300))+ggplot2::scale_x_continuous(breaks=seq(-16,16,2))+
  ggplot2::scale_y_continuous(breaks=seq(0,350,50))
plot01




### plot3
nrDEG_limma_voom <- resdata_symbol
nrDEG_limma_voom$regulate <- ifelse(nrDEG_limma_voom$padj > 0.05, "unchanged",
                                    ifelse(nrDEG_limma_voom$log2FoldChange > 1, "up-regulated",
                                           ifelse(nrDEG_limma_voom$log2FoldChange < -1, "down-regulated", "unchanged")))

nrDEG_limma_voom <- nrDEG_limma_voom[!is.na(nrDEG_limma_voom$padj),]
for(i in 1:nrow(nrDEG_limma_voom)) {
  if(nrDEG_limma_voom$padj[i] == 0) {
    nrDEG_limma_voom$padj[i] <- 1e-300
  }
}

plot03 <- ggplot(nrDEG_limma_voom,aes(x=log2FoldChange, y=-log10(padj),color = regulate)) + 
  geom_point() + 
  theme_bw()+
  scale_color_manual(values = c("#4E78C4","gray","#CE2220")) +
  geom_hline(yintercept = -log10(0.05),linetype=2,size=0.6) +
  geom_vline(xintercept = c(-1,1),linetype=2,size=0.5)+
  theme(legend.title = element_text(size=20, face = "bold"))+ ##设置图例大小
  theme(legend.text=element_text(size=15))+
  theme(legend.key.size = unit(1.6, 'lines'))+
  theme(axis.text=element_text(vjust=1,size=15))+
  theme(axis.title = element_text(size = 17))+
  geom_text_repel(data = nrDEG_limma_voom[(abs(nrDEG_limma_voom$log2FoldChange) > 1.5 & nrDEG_limma_voom$padj < 1e-260) |
                                            (nrDEG_limma_voom$log2FoldChange > 6.02 & nrDEG_limma_voom$padj < 0.03 ) |
                                            (nrDEG_limma_voom$log2FoldChange < -4 & nrDEG_limma_voom$padj < 0.03 ), ],
                  aes(label = SYMBOL),
                  size = 3,
                  show.legend = F)+
  scale_x_continuous(limits = c(-16,16))+
  scale_y_continuous(limits = c(-1,300))

pdf(file = 'Volcanomap02.pdf',width = 8,height = 6)
print(plot03)
dev.off()


##GSEA分析
#geneList三部曲
#1.获取基因logFC
geneList <- resdata_symbol$log2FoldChange
#2.命名
names(geneList) <- resdata_symbol$ENTREZID
#3.排序很重要
geneList <- sort(geneList, decreasing = TRUE)
#geneList <- geneList[!duplicated(names(geneList))]

## msigdb数据库下载人类的C5注释集
# library(msigdbr)
# m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
#   dplyr::select(gs_name, entrez_gene)
# head(m_t2g)
# gsea_res <- GSEA(genelist, TERM2GENE = m_t2g, minGSSize = 10, maxGSSize = 500,pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = 456)

kk_gsea <- gseKEGG(geneList, organism = "hsa", 
                   minGSSize = 10, maxGSSize = 10000, 
                   pvalueCutoff = 0.5,verbose = FALSE)
#kk_gsea <- gseKEGG(geneList, organism = "mmu",pvalueCutoff = 1.0)
saveRDS(kk_gsea,"kk_gsea.rds")
kk_gsea <- readRDS("./kk_gsea.rds")
kk <- kk_gsea

head(kk)

#按照enrichment score从高到低排序，便于查看富集通路
sortkk <- kk[order(kk$enrichmentScore, decreasing = T),]
head(sortkk)

#write.csv(sortkk, "pathway_gsea_output2.csv", row.names = FALSE)
kk@result <- setReadable(kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame() # 转换基因ID信息为symbol
# kk <- setReadable(kk,OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", drop = FALSE) # 

gene_conversion <- AnnotationDbi::select(org.Hs.eg.db, keys = names(kk@geneList),keytype = "ENTREZID",columns = c("ENTREZID", "SYMBOL")) # 
kk@gene2Symbol <- setNames(gene_conversion$SYMBOL, gene_conversion$ENTREZID)

GOresult <- kk@result
# GOresult$Description <- sub(" - Mus musculus \\(house mouse\\)", "", GOresult$Description) #去除通路多余的描述
# GOresult$Description <- gsub("/", "_", GOresult$Description)

write.csv(GOresult, "gsea_pathway.csv",row.names = FALSE)


### 绘制所有通路的KEGG通路或者前25条：
dir.create('gseaplot')  #创建文件夹储存结果
for (i in 1:length(GOresult$Description)) {
  plot1 <- gseaplot(kk, geneSetID = i, title = GOresult$Description[i])
  pdf(paste0('gseaplot/', GOresult$Description[i], ".pdf"), height = 10, width = 10)
  print(plot1)
  dev.off()
}
for (i in 1:length(GOresult$Description)) {
  # 替换 Description 中的 / 为下划线 _
  safe_title <- gsub("/", "_", GOresult$Description[i])
  # 创建 GSEA 图
  plot1 <- gseaplot(kk, geneSetID = i, title = GOresult$Description[i])
  # 保存为 PDF 文件，使用 safe_title 避免 '/' 导致的路径问题
  pdf(paste0('gseaplot/', safe_title, ".pdf"), height = 10, width = 10)
  print(plot1)
  dev.off()
}

dir.create('gseaplot2')  #创建文件夹储存结果
for (i in 1:length(GOresult$Description)) {
  plot1 <- gseaplot2(kk, geneSetID = i, title = GOresult$Description[i],base_size = 14, pvalue_table = F)
  pdf(paste0('gseaplot2/', GOresult$Description[i], ".pdf"), height = 8, width = 8)
  print(plot1)
  dev.off()
}
for (i in 1:length(GOresult$Description)) {
  # 替换 Description 中的 / 为下划线 _
  safe_title <- gsub("/", "_", GOresult$Description[i])
  # 创建 GSEA 图
  plot1 <- gseaplot2(kk, geneSetID = i, title = GOresult$Description[i],base_size = 14, pvalue_table = F)
  # 保存为 PDF 文件，使用 safe_title 避免 '/' 导致的路径问题
  pdf(paste0('gseaplot2/', safe_title, ".pdf"), height = 8, width = 8)
  print(plot1)
  dev.off()
}


dir.create('gseaplot')  #创建文件夹储存结果
for (i in 1:min(25, length(GOresult$Description))) {
  plot1 <- gseaplot(kk, geneSetID = i, title = GOresult$Description[i])
  pdf(paste0('gseaplot/', GOresult$Description[i], ".pdf"), height = 10, width = 10)
  print(plot1)
  dev.off()
}
dir.create('gseaplot2')  #创建文件夹储存结果
for (i in 1:min(25, length(GOresult$Description))) {
  plot1 <- gseaplot2(kk, geneSetID = i, title = GOresult$Description[i],base_size = 14, pvalue_table = F)
  pdf(paste0('gseaplot2/', GOresult$Description[i], ".pdf"), height = 8, width = 8)
  print(plot1)
  dev.off()
}



### 单独指定通路绘制及多条通路绘制到一起：
gseaplot2(kk,c("hsa03430","hsa03410","hsa03460","hsa05204", "hsa04110"),
          base_size = 18,
          rel_heights = c(1, 0.2),subplots = 1:2,
          pvalue_table = FALSE,
          ES_geom = "line",
          color = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"))




gene_kk <- union(res_up$ENTREZID, res_dw$ENTREZID)

##GO分析
GO_BP <- enrichGO(gene = gene_kk,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "ENTREZID",  #设定读取的gene ID类型
                  ont = "BP",  #(ont包括 Biological Process,Cellular Component,Mollecular Function三部分）
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,  #设定p值阈值
                  qvalueCutoff = 0.05  #设定q值阈值
)
# trans ID
GO_BP@result <- setReadable(GO_BP, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame() # 转换基因ID信息为symbol
# kk <- setReadable(kk,OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", drop = FALSE) # 直接转换基因会变少

# gene_conversion <- AnnotationDbi::select(org.Hs.eg.db, keys = names(GO_BP@geneList),keytype = "ENTREZID",columns = c("ENTREZID", "SYMBOL")) # AnnotationDbi不会去重
# GO_BP@gene2Symbol <- setNames(gene_conversion$SYMBOL, gene_conversion$ENTREZID)
saveRDS(GO_BP, "./GO_BP.rds")
write.csv(GO_BP@result, "GO_BP.csv")


GO_CC <- enrichGO(gene = gene_kk,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "ENTREZID",  #设定读取的gene ID类型
                  ont = "CC",  #(ont包括 Biological Process,Cellular Component,Mollecular Function三部分）
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,  #设定p值阈值
                  qvalueCutoff = 0.05  #设定q值阈值
)
# trans ID
GO_CC@result <- setReadable(GO_CC, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame() # 转换基因ID信息为symbol
# kk <- setReadable(kk,OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", drop = FALSE) # 直接转换基因会变少

# gene_conversion <- AnnotationDbi::select(org.Hs.eg.db, keys = names(GO_CC@geneList),keytype = "ENTREZID",columns = c("ENTREZID", "SYMBOL")) # AnnotationDbi不会去重
# GO_CC@gene2Symbol <- setNames(gene_conversion$SYMBOL, gene_conversion$ENTREZID)
saveRDS(GO_CC, "./GO_CC.rds")
write.csv(GO_CC@result, "GO_CC.csv")


GO_MF <- enrichGO(gene = gene_kk,
                  OrgDb = "org.Hs.eg.db",
                  keyType = "ENTREZID",  #设定读取的gene ID类型
                  ont = "MF",  #(ont包括 Biological Process,Cellular Component,Mollecular Function三部分）
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,  #设定p值阈值
                  qvalueCutoff = 0.05  #设定q值阈值
)
# trans ID
GO_MF@result <- setReadable(GO_MF, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame() # 转换基因ID信息为symbol
# kk <- setReadable(kk,OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", drop = FALSE) # 直接转换基因会变少

# gene_conversion <- AnnotationDbi::select(org.Hs.eg.db, keys = names(GO_MF@geneList),keytype = "ENTREZID",columns = c("ENTREZID", "SYMBOL")) # AnnotationDbi不会去重
# GO_MF@gene2Symbol <- setNames(gene_conversion$SYMBOL, gene_conversion$ENTREZID)
saveRDS(GO_MF, "./GO_MF.rds")
write.csv(GO_MF@result, "GO_MF.csv")


# 限制基因数200个，也可以不限制
kk_KEGG <- enrichKEGG(gene = gene_kk, organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff =0.05)
# trans ID
kk_KEGG@result <- setReadable(kk_KEGG, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID") %>% as.data.frame() # 转换基因ID信息为symbol
# kk <- setReadable(kk,OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", drop = FALSE) # 直接转换基因会变少

# gene_conversion <- AnnotationDbi::select(org.Hs.eg.db, keys = names(kk_KEGG@geneList),keytype = "ENTREZID",columns = c("ENTREZID", "SYMBOL")) # AnnotationDbi不会去重
# kk_KEGG@gene2Symbol <- setNames(gene_conversion$SYMBOL, gene_conversion$ENTREZID)
saveRDS(kk_KEGG, "./kk_KEGG.rds")
write.csv(kk_KEGG@result, "KEGG.csv")

KEGG <- read.csv("KEGG.csv")
KEGG$Ratio <- KEGG$Count/as.numeric(sub(".*/", "", KEGG$GeneRatio[1]))

# KEGG$Description <- sub(" - Mus musculus \\(house mouse\\)", "", KEGG$Description) #去除通路多余的描述
KEGG$Description <- gsub("/", "_", KEGG$Description)

write.csv(KEGG, "KEGG_new.csv",row.names = FALSE)

p7 <- ggplot(KEGG[1:20,], aes(x=Ratio, y=Description, size = Count, color=-log10(pvalue)))+
  geom_point(alpha=0.5)+
  theme_bw()+
  scale_size(range = c(2, 9), name="Gene Number")+
  theme(
    axis.text.x = element_text(size = 15, color = "black", face = "bold"),
    axis.text.y = element_text(size = 15, color = "black", face = "bold"),
    legend.title = element_text(size = 15, color = "black", face = "bold"),
    legend.text = element_text(size = 13, face = "bold"),
    axis.title.x = element_text(size = 18, color = "black", face = "bold")
  )+
  labs(x="GeneRatio", y="")+
  scale_colour_gradientn(colours = c("#132B43","#56B1F7"),
                         guide = guide_colourbar(reverse = TRUE)) +
  plot_annotation(
    title = "KEGG Enrichment ScatterPlot",
    theme = theme(plot.title = element_text(size = 23, hjust = 0.6,face = "bold"))
  )

pdf(file = 'KEGG_qipao.pdf',width = 8,height = 6)
print(p7)
dev.off()


# 限制200个

sortup <- res_up[order(res_up$log2FoldChange, decreasing = T),]
sortdw <- res_dw[order(res_dw$log2FoldChange, decreasing = F),]
res_up200 <- sortup[1:200, ] #限制富集分析使用的基因数
res_dw200 <- sortdw[1:200, ]
kk_up <- enrichKEGG(gene = res_up200$ENTREZID,organism = 'hsa',pAdjustMethod = "BH",pvalueCutoff = 1)
kk_dw <- enrichKEGG(gene = res_dw200$ENTREZID,organism = 'hsa',pAdjustMethod = "BH",pvalueCutoff = 1)
# options(clusterProfiler.download.method = "wininet")
kk_diff <- enrichKEGG(gene = res_diff$ENTREZID,organism = 'hsa',pAdjustMethod = "BH",pvalueCutoff = 1)
# kk_up <- enrichKEGG(gene = res_up$ENTREZID,organism = 'rno',pAdjustMethod = "BH",pvalueCutoff = 1)
# kk_dw <- enrichKEGG(gene = res_dw$ENTREZID,organism = 'rno',pAdjustMethod = "BH",pvalueCutoff = 1)

enrich_kegg_up <- as.data.frame(kk_up@result)
# enrich_kegg_up$Description <- sub(" - Mus musculus \\(house mouse\\)", "", enrich_kegg_up$Description) #去除通路多余的描述
enrich_kegg_up$Description <- gsub("/", "_", enrich_kegg_up$Description)
write.csv(enrich_kegg_up,"pathway_kegg_up.csv",row.names = F)
enrich_kegg_dw <- as.data.frame(kk_dw@result)
#enrich_kegg_dw$Description <- sub(" - Mus musculus \\(house mouse\\)", "", enrich_kegg_dw$Description) #去除通路多余的描述
enrich_kegg_dw$Description <- gsub("/", "_", enrich_kegg_dw$Description)
write.csv(enrich_kegg_dw,"pathway_kegg_dw.csv",row.names = F)

# 条形图
plot7 <- barplot(kk_dw,title="Enrichment KEGG_dw",showCategory=20)
pdf(file = 'kegg_bar_dw.pdf',width = 10,height = 8)
print(plot7)
dev.off()
plot8 <- barplot(kk_up,title="Enrichment KEGG_up",showCategory=20)
pdf(file = 'kegg_bar_up.pdf',width = 10,height = 8)
print(plot8)
dev.off()

# 点图展示，使用一个函数
erich2plot <- function(data4plot){
  library(ggplot2)
  data4plot <- data4plot[order(data4plot$qvalue,decreasing = F)[1:20],]  #取qvalue最小的20个通路
  data4plot$BgRatio <- 
    apply(data4plot,1,function(x){
      as.numeric(strsplit(x[3],'/')[[1]][1])
    })/apply(data4plot,1,function(x){
      as.numeric(strsplit(x[4],'/')[[1]][1])
    })  #计算提供基因富集到通路的比例
  
  p <- ggplot(data4plot,aes(BgRatio,Description))  #画图BgRatio为富集到的比例，Description为通路名称
  p <- p + geom_point()
  
  pbubble <- p + geom_point(aes(size = Count,color= -1*log10(qvalue)))  #增加维度Count为hit到的数量，color表示其显著性
  
  pr <- pbubble + scale_colour_gradient(low = "blue",high = "red") +
    labs(color=expression(-log[10](qvalue)),size="observed.gene.count",
         x = "Richfactor",y = "term.description",title = "Enrichment Process")
  
  pr <- pr + theme_bw()
  pr  #颜色表示通路富集的显著性，大小表示通路富集到的基因数量，横坐标表示富集到的基因占通路的比例
}






## 筛选之后出图
if(F){
  library(stringr)
  library(DOSE)
  dir.create("GO_KEGG_result")
  setwd("./GO_KEGG_result/")
  goBP <- read.csv("../GO_BP.csv")
  # GO:0006261 DNA-templated DNA replication
  # GO:0044786 cell cycle DNA replication
  # GO:0051446 positive regulation of meiotic cell cycle
  # GO:0000727 double-strand break repair via break-induced replication
  # GO:1902969 GO:0006270 GO:0030174 GO:0006268 GO:0033260 GO:0090329 GO:0006271 GO:0006261 复制相关
  
  gobp_id <- c("GO:0006261", "GO:0044786", "GO:0051446", "GO:0000727", "GO:1902969", "GO:0033260", "GO:0006268")
  goBP <- goBP[goBP$ID %in% gobp_id, ]
  
  goBP$ONTOLOGY <- rep("BP")
  goBP$Ratio <- goBP$Count/as.numeric(sub(".*/", "", goBP$GeneRatio[1]))
  goBP <- goBP[order(goBP$Ratio,decreasing = T),]
  goBP$Description <- factor(goBP$Description,levels = rev(goBP$Description))
  goBP <- mutate(goBP, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))# 计算富集倍数，如果使用其画图，需要重新排序通路顺序
  
  
  goCC <- read.csv("../GO_CC.csv")
  gocc_id <- c("GO:0031261", "GO:0098636", "GO:0030426", "GO:0005923", "GO:0030427")
  goCC <- goCC[goCC$ID %in% gocc_id, ]
  
  goCC$ONTOLOGY <- rep("CC")
  goCC$Ratio <- goCC$Count/as.numeric(sub(".*/", "", goCC$GeneRatio[1]))
  goCC <- goCC[order(goCC$Ratio,decreasing = T),]
  goCC$Description <- factor(goCC$Description,levels = rev(goCC$Description))
  goCC <- mutate(goCC, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))# 计算富集倍数，如果使用其画图，需要重新排序通路顺序
  
  
  goMF <- read.csv("../GO_MF.csv")
  gomf_id <- c("GO:0001216", "GO:0098631", "GO:0019838", "GO:0098632", "GO:0008083", "GO:0070851") # 生长因子相关 GO:0008083 GO:0070851
  goMF <- goMF[goMF$ID %in% gomf_id, ]
  
  goMF$ONTOLOGY <- rep("MF")
  goMF$Ratio <- goMF$Count/as.numeric(sub(".*/", "", goMF$GeneRatio[1]))
  goMF <- goMF[order(goMF$Ratio,decreasing = T),]
  goMF$Description <- factor(goMF$Description,levels = rev(goMF$Description))
  goMF <- mutate(goMF, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))# 计算富集倍数，如果使用其画图，需要重新排序通路顺序
  
  
  KEGG <- read.csv("../KEGG_new.csv")
  kegg_id <- c("hsa03030", "hsa04514", "hsa04115", "hsa04015", "hsa04974", "hsa04060", "hsa04520", "hsa04151", "hsa04024", "hsa04810") # 后3条重新加的 hsa00330
  KEGG <- KEGG[KEGG$ID %in% kegg_id, ]
  
  KEGG$ONTOLOGY <- rep("KEGG")
  #KEGG$Ratio <- KEGG$Count/as.numeric(sub(".*/", "", KEGG$GeneRatio[1]))
  KEGG <- KEGG[order(KEGG$Ratio,decreasing = T),]
  KEGG$Description <- factor(KEGG$Description,levels = rev(KEGG$Description))
  KEGG <- mutate(KEGG, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))# 计算富集倍数，如果使用其画图，需要重新排序通路顺序
}


### 可视化方法1
if(F){
  # 创建颜色梯度映射
  color_gradient_BP <- colorRampPalette(c("#f4b183", "#FCEDE3"))
  color_gradient_CC <- colorRampPalette(c("#a9d18e", "#EBF4E5"))
  color_gradient_MF <- colorRampPalette(c("#9dc3e6", "#E9F1F9"))
  
  
  p1 <- ggplot(data = goBP, aes(x = Description, y = Ratio, fill = pvalue))+ 
    geom_bar(stat = "identity",width = 0.9)+ 
    coord_flip()+
    theme_bw()+
    labs(x = "",y = "",fill = "Pvalue \n\nBP")+ # 设置坐标轴标题及标题
    scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+
    theme(axis.text.y = element_text(size = 15, color = "black", face = "bold"), 
          legend.title = element_text(size = 14, color = "black", face="bold"), 
          legend.text = element_text(size = 10, color = "black", face = "bold"),
          legend.key.size = unit(1, 'lines'),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.ticks.length.x = unit(0, "pt"), 
          axis.text.x = element_blank(),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 1.5, 
                               unit = "pt"))+
    scale_fill_gradientn(colors =  color_gradient_BP(10),
                         guide = guide_colourbar(reverse = TRUE))+
    scale_y_continuous(limits = c(0, 0.03)) # 调整横坐标
  print(p1)
  
  p2 <- ggplot(data = goCC, aes(x = Description, y = Ratio, fill = pvalue))+ 
    geom_bar(stat = "identity",width = 0.9)+ 
    coord_flip()+
    theme_bw()+ 
    labs(x = "",y = "",fill = "CC")+ # 设置坐标轴标题及标题
    scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+
    theme(axis.text.y = element_text(size = 15, color = "black", face = "bold"), # 坐标轴标签大小
          legend.title = element_text(size = 14, color = "black", face="bold"), # 图例标题大小
          legend.text = element_text(size = 10, color = "black", face = "bold"),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          legend.key.size = unit(1, 'lines'),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.length.x = unit(0, "pt"),
          plot.margin = margin(t = 3, r = 0, b = 0, l = 1.5, 
                               unit = "pt"))+
    scale_fill_gradientn(colors =  color_gradient_CC(10),
                         guide = guide_colourbar(reverse = TRUE))+
    scale_y_continuous(limits = c(0, 0.03)) # 调整横坐标
  print(p2)
  
  
  p3 <- ggplot(data = goMF, aes(x = Description, y = Ratio, fill = pvalue))+ 
    geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
    coord_flip()+
    theme_bw()+ 
    labs(x = "",y = "GeneRatio",fill = "MF")+ # 设置坐标轴标题及标题
    scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+
    theme(axis.title.x = element_text(size = 20, color = "black", face = "bold"),
          axis.text = element_text(size = 15, color = "black", face = "bold"), # 坐标轴标签大小
          legend.title = element_text(size = 14, color = "black",face="bold"), # 图例标题大小
          legend.text = element_text(size = 10, color = "black", face = "bold"),
          legend.key.size = unit(1, 'lines'),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(t = 3, r = 0, b = 0, l = 1.5, 
                               unit = "pt"))+
    scale_fill_gradientn(colors =  color_gradient_MF(10),
                         guide = guide_colourbar(reverse = TRUE))+
    scale_y_continuous(limits = c(0, 0.03)) # 调整横坐标
  print(p3)
  
  combined_plot <- wrap_plots(p1, p2, p3, ncol = 1, nrow = 3)
  print(combined_plot)
  combined_plot +
    plot_annotation(
      title = "Barplot of Enriched GO Terms",
      theme = theme(plot.title = element_text(size = 23, hjust = 0.5,face = "bold"))
    )
  
}


### 可视化方法2
if(F){
  GO_data <- rbind(goBP, goCC, goMF)
  
  ### 条形图
  #ontology_colors <- c("BP" = "#e55709", "CC" = "#1c8041", "MF" = "#501d8a")
  ontology_colors <- c("BP" = "#4E78C4", "CC" = "#7EB875", "MF" = "#824D99")
  
  pdf(file = 'GO_p5.pdf',width = 12,height = 8)
  ggplot(GO_data, aes(x = reorder(Description, Count), y = Count, fill = ONTOLOGY)) +
    geom_bar(stat = "identity") + 
    coord_flip() +  # 旋转为横向柱状图
    facet_grid(ONTOLOGY ~ ., scales = "free") +  # 根据 BP/CC/MF 分类
    scale_fill_manual(values = ontology_colors) +  # 指定颜色
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12)
    ) +
    xlab("GO Term") +
    ylab("Gene Count") +
    ggtitle("GO Enrichment Analysis")
  dev.off()
  
  ### 气泡图
  pdf(file = 'GO_p6.pdf',width = 12,height = 8)
  ggplot(GO_data, aes(x = Ratio, y = reorder(Description, Ratio), size = Count, color = ONTOLOGY)) +
    geom_point(alpha = 0.8) +  # 透明度 0.8
    scale_color_manual(values = ontology_colors) +  # 指定分类颜色
    facet_grid(ONTOLOGY ~ ., scales = "free") +  # 分面 BP/CC/MF
    theme_minimal() +
    theme(
      text = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      strip.text = element_text(size = 12),  # 调整分面标签字体
      legend.text = element_text(size = 12)  # 调整图例字体
    ) +
    labs(x = "Gene Ratio", y = "GO Term", title = "GO Enrichment Bubble Plot") +
    scale_size(range = c(3, 10))  # 调整气泡大小范围
  dev.off()
}

### 可视化方法3
if(F){
  library(tidyverse)
  library(ggh4x)  #ggnewscale_0.5.1 ggh4x_0.3.0
  library(ggfun)
  library(ggnewscale)
  library(grid)
  library(clusterProfiler)
  #remotes::install_version("ggnewscale", "0.5.1")
  
  plot_df <- bind_rows(goBP, goCC, goMF) %>%
    dplyr::mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c("BP", "CC", "MF")), ordered = T)) %>%
    dplyr::arrange(ONTOLOGY, desc(Count)) %>% # 按照 ONTOLOGY 和 Count 列对数据进行排序
    # dplyr::mutate(Description = str_remove(Description, pattern = ",.*")) %>% # 删除逗号后所有字符
    dplyr::mutate(Description = factor(Description, levels = rev(Description), ordered = T))
  
  ####----Plot----####
  
  ### Count
  plot <- plot_df %>%
    ggplot() + 
    # geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "KEGG"),
    #            aes(x = Count, y = interaction(Description, ONTOLOGY), fill = p.adjust, size = Count), shape = 21) + 
    # scale_fill_gradient(low = "#a1d99b", high = "#238b45", name = "KEGG p.adjust") + 
    # ggnewscale::new_scale_fill() + 
    geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "MF"),
               aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) + 
    scale_fill_gradient(low = "#b7c1e4", high = "#4E78C4", name = "MF pvalue") + 
    ggnewscale::new_scale_fill() + 
    geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "CC"),
               aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) + 
    scale_fill_gradient(low = "#bbd2bc", high = "#7EB875", name = "CC pvalue") + 
    ggnewscale::new_scale_fill() +
    geom_point(data = plot_df %>% dplyr::filter(ONTOLOGY == "BP"),
               aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) + 
    scale_fill_gradient(low = "#c8b3d6", high = "#824D99",  name = "BP pvalue") + 
    guides(y = "axis_nested",
           y.sec = guide_axis_manual(breaks = 1:18, # 通路数
                                     labels = plot_df %>%
                                       group_by(ONTOLOGY) %>%
                                       mutate(Description = rev(Description)) %>%
                                       pull(Description) # 提取反转后的名称
           )) +  # 反转顺序,检查通路展示顺序是否有误
    ggtitle(label = "GO annotation") + 
    labs(x = "Count", y = "Description") + 
    scale_size(range = c(3, 7),
               guide = guide_legend(override.aes = list(fill = "#000000"))) + 
    theme_bw() + 
    theme(
      ggh4x.axis.nestline.y = element_line(size = 3, color = c("#4E78C4", "#7EB875", "#824D99")), # "#74c476",
      ggh4x.axis.nesttext.y = element_text(colour = c("#4E78C4", "#7EB875", "#824D99")), # "#74c476", 
      legend.background = element_roundrect(color = "#969696"),
      panel.border = element_rect(size = 0.5),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
      axis.text = element_text(color = "#000000", size = 11),
      #axis.text.y = element_text(color = rep(c("#225ea8", "#fc4e2a", "#88419d"), each = 5)), # "#41ae76", 
      axis.text.y = element_text(color = c(rep("#4E78C4",6), rep("#7EB875",5), rep("#824D99",7))),
      axis.text.y.left = element_blank(),
      axis.ticks.length.y.left = unit(10, "pt"),
      axis.ticks.y.left = element_line(color = NA),
      axis.title = element_text(color = "#000000", size = 15),
      plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
    ) + 
    coord_cartesian(clip = "off") + 
    annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                           gp = gpar(col = "#969696", lwd = 1.5)),
                      xmin = unit(3, "native"),
                      xmax = unit(15, "native"),
                      ymin = unit(40.85, "native"),
                      ymax = unit(42.25, "native"))
  
  plot
  pdf(file = '../plot_out/p3_3_GO_p1.pdf',width = 8,height = 8)
  print(plot)
  dev.off()
  # ggsave(filename = "GO_KEGG.pdf",
  #        plot = plot,
  #        height = 11,
  #        width = 12.5)
  
  ## 条形图
  plot1 <- plot_df %>%
    ggplot() + 
    # MF 条形图
    geom_col(data = plot_df %>% filter(ONTOLOGY == "MF"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue), width = 0.7) + 
    scale_fill_gradient(low = "#b7c1e4", high = "#4E78C4", name = "MF pvalue") + 
    ggnewscale::new_scale_fill() +
    # CC 条形图
    geom_col(data = plot_df %>% filter(ONTOLOGY == "CC"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue), width = 0.7) + 
    scale_fill_gradient(low = "#bbd2bc", high = "#7EB875", name = "CC pvalue") + 
    ggnewscale::new_scale_fill() +
    # BP 条形图
    geom_col(data = plot_df %>% filter(ONTOLOGY == "BP"),
             aes(x = Count, y = interaction(Description, ONTOLOGY), fill = pvalue), width = 0.7) + 
    scale_fill_gradient(low = "#c8b3d6", high = "#824D99", name = "BP pvalue") + 
    
    guides(y = "axis_nested",
           y.sec = guide_axis_manual(breaks = 1:18,
                                     labels = plot_df %>%
                                       group_by(ONTOLOGY) %>%
                                       mutate(Description = rev(Description)) %>%
                                       pull(Description) # 提取反转后的名称
           )) +  # 反转顺序,检查通路展示顺序是否有误
    ggtitle(label = "GO annotation") + 
    labs(x = "Count", y = "Description") + 
    scale_size(range = c(3, 7),
               guide = guide_legend(override.aes = list(fill = "#000000"))) + 
    theme_bw() + 
    theme(
      ggh4x.axis.nestline.y = element_line(size = 3, color = c("#4E78C4", "#7EB875", "#824D99")), # "#74c476",
      ggh4x.axis.nesttext.y = element_text(colour = c("#4E78C4", "#7EB875", "#824D99")), # "#74c476", 
      legend.background = element_roundrect(color = "#969696"),
      panel.border = element_rect(size = 0.5),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
      axis.text = element_text(color = "#000000", size = 11),
      #axis.text.y = element_text(color = rep(c("#225ea8", "#fc4e2a", "#88419d"), each = 5)), # "#41ae76", 
      axis.text.y = element_text(color = c(rep("#4E78C4",6), rep("#7EB875",5), rep("#824D99",7))),
      axis.text.y.left = element_blank(),
      axis.ticks.length.y.left = unit(10, "pt"),
      axis.ticks.y.left = element_line(color = NA),
      axis.title = element_text(color = "#000000", size = 15),
      plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
    ) + 
    coord_cartesian(clip = "off") + 
    annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                           gp = gpar(col = "#969696", lwd = 1.5)),
                      xmin = unit(3, "native"),
                      xmax = unit(15, "native"),
                      ymin = unit(40.85, "native"),
                      ymax = unit(42.25, "native"))
  
  plot1
  pdf(file = '../plot_out/p3_3_GO_p2.pdf',width = 8,height = 6)
  print(plot1)
  dev.off()
}


### KEGG
if(F){
  ### Ratio
  plot5 <- KEGG %>%
    ggplot() + 
    geom_point(data = KEGG %>% dplyr::filter(ONTOLOGY == "KEGG"),
               aes(x = Ratio, y = interaction(Description, ONTOLOGY), fill = pvalue, size = Count), shape = 21) + 
    scale_fill_gradient(low = "#de8c92", high = "#ce2220",  name = "pvalue") + 
    guides(y = "axis_nested",
           y.sec = guide_axis_manual(breaks = 1:10, ## 通路条数
                                     labels = KEGG %>%
                                       group_by(ONTOLOGY) %>%
                                       mutate(Description = rev(Description)) %>%
                                       pull(Description) # 提取反转后的名称
           )) +  # 反转顺序,检查通路展示顺序是否有误
    ggtitle(label = "KEGG annotation") + 
    labs(x = "Ratio", y = "Description") + 
    scale_size(range = c(3, 7),
               guide = guide_legend(override.aes = list(fill = "#000000"))) + 
    theme_bw() + 
    theme(
      ggh4x.axis.nestline.y = element_blank(), # 去除左侧线
      ggh4x.axis.nesttext.y = element_text(size = 0, colour = "#ce2220"),  
      legend.background = element_roundrect(color = "#969696"),
      panel.border = element_rect(size = 0.5),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
      axis.text = element_text(color = "#000000", size = 11),
      # axis.text.y = element_text(color = rep(c("#225ea8", "#fc4e2a", "#88419d"), each = 5)), # "#41ae76", 
      axis.text.y = element_text(color = "#000000"), # "#88419d"
      axis.text.y.left = element_blank(),
      axis.ticks.length.y.left = unit(10, "pt"),
      axis.ticks.y.left = element_line(color = NA),
      axis.title = element_text(color = "#000000", size = 15),
      plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
    ) + 
    coord_cartesian(clip = "off") + 
    annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                           gp = gpar(col = "#969696", lwd = 1.5)),
                      #xmin = unit(3, "native"),xmax = unit(15, "native"),ymin = unit(40.85, "native"),ymax = unit(42.25, "native")
                      xmin = 3, xmax = 15, ymin = 40.85, ymax = 42.25)
  
  plot5
  pdf(file = '../plot_out/p3_4_KEGG_dot.pdf',width = 8,height = 4.5)
  print(plot5)
  dev.off()
  # ggsave(filename = "GO_KEGG.pdf",
  #        plot = plot,
  #        height = 11,
  #        width = 12.5)
  
  plot6 <- KEGG %>%
    ggplot() + 
    geom_col(data = KEGG %>% filter(ONTOLOGY == "KEGG"),
             aes(x = Ratio, y = interaction(Description, ONTOLOGY), fill = pvalue), width = 0.7) + 
    scale_fill_gradient(low = "#CCC8E5", high = "#8c6bb1", name = "pvalue") + 
    
    guides(y = "axis_nested",
           y.sec = guide_axis_manual(breaks = 1:10,
                                     labels = KEGG %>%
                                       group_by(ONTOLOGY) %>%
                                       mutate(Description = rev(Description)) %>%
                                       pull(Description) # 提取反转后的名称
           )) +  # 反转顺序,检查通路展示顺序是否有误
    ggtitle(label = "KEGG annotation") + 
    labs(x = "Ratio", y = "Description") + 
    scale_size(range = c(3, 7),
               guide = guide_legend(override.aes = list(fill = "#000000"))) + 
    theme_bw() + 
    theme(
      ggh4x.axis.nestline.y = element_blank(), # 去除左侧线
      ggh4x.axis.nesttext.y = element_text(size = 0, colour = "#9e9ac8"),  
      legend.background = element_roundrect(color = "#969696"),
      panel.border = element_rect(size = 0.5),
      plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "cm"),
      axis.text = element_text(color = "#000000", size = 11),
      # axis.text.y = element_text(color = rep(c("#225ea8", "#fc4e2a", "#88419d"), each = 5)), # "#41ae76", 
      axis.text.y = element_text(color = "#000000"), # "#88419d"
      axis.text.y.left = element_blank(),
      axis.ticks.length.y.left = unit(10, "pt"),
      axis.ticks.y.left = element_line(color = NA),
      axis.title = element_text(color = "#000000", size = 15),
      plot.title = element_text(color = "#000000", size = 20, hjust = 0.5)
    ) + 
    coord_cartesian(clip = "off") + 
    annotation_custom(grob = roundrectGrob(r = unit(0.2, "snpc"),
                                           gp = gpar(col = "#969696", lwd = 1.5)),
                      xmin = unit(3, "native"),
                      xmax = unit(15, "native"),
                      ymin = unit(40.85, "native"),
                      ymax = unit(42.25, "native"))
  
  plot6
  pdf(file = 'KEGG_col.pdf',width = 8,height = 6)
  print(plot2)
  dev.off()
}




setwd("../")
#tpm <- tpm[, -c(1:8, (ncol(tpm)-2):ncol(tpm))] # wangyue_data
tpm_diff <- tpm[rownames(tpm)%in%res_diff$Row.names,]
write.csv(tpm_diff, "tpm_diff.csv")
tpm_up <- tpm[rownames(tpm)%in%res_up$Row.names,]
tpm_dw <- tpm[rownames(tpm)%in%res_dw$Row.names,]

pca.info <- fast.prcomp(tpm_diff)
head(pca.info$rotation) #显示PCA计算结果
pca.data <- data.frame(sample = rownames(pca.info$rotation),Type = database$Sample_Type,pca.info$rotation) #加入标签信息
plot2 <- ggscatter(pca.data,x = "PC1",y = "PC2",color = "Type") + theme_base() #绘图
pdf(file = 'pc_1.pdf',width = 8,height = 6)
print(plot2)
dev.off()

# 另一个PCA图：
ggscatter(pca.data,x = "PC1",y = "PC2",color = "Type") + 
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"),
    legend.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 13),
    axis.title = element_text(size = 18, color = "black"))

# 计算95%置信，适用于多样本查看组间差异，带置信圈的PCA图
#pcaData <- plotPCA(pca.data, intgroup=c("condition", "type"), returnData=T) 
percentVar <- round(100*attr(pca.data, "percentVar"))
plot3 <- ggplot(pca.data, aes(PC1, PC2, color=Type)) + 
  geom_point(size=3) +
  scale_color_manual(values=c("#CE2220","#4E78C4"))+ # 
  ggtitle("PCA") + 
  xlab(paste0("PC1")) + 
  ylab(paste0("PC2"))+
  theme_bw()+stat_ellipse(aes(color = Type), level = 0.95, show.legend = FALSE) ##计算95%置信
pdf(file = './plot_out/p3_1B_pc_2.pdf',width = 7,height = 6)
print(plot3)
dev.off()




df <- database
rownames(df) <- df[,1]
df[,1] <- NULL
tpm_diff1 <- merge(tpm_diff,gene_df,by.x="row.names",by.y="ENSEMBL")
#tpm_diff1 <- tpm_diff1[, -8] # wangyue_data
#tpm_diff1 <- tpm_diff1[, -1] # wangyue_data
tpm_diff1 <- filter(tpm_diff1,!duplicated(tpm_diff1$SYMBOL))
rownames(tpm_diff1) <- tpm_diff1$SYMBOL
tpm_diff1 <- tpm_diff1[,-c(1,ncol(tpm_diff1)-1,ncol(tpm_diff1))]
#tpm_diff1 <- tpm_diff1[, -c(1:6, (ncol(tpm_diff1)-2):ncol(tpm_diff1))] # wangyue_data
tpm_diff1 <- t(tpm_diff1)
database_cut <- subset(database,select=c(1,Sample_Type))  #
tpm_out <- merge(tpm_diff1,database_cut,by.x="row.names",by.y="ID")
tpm_out <- as.data.frame(tpm_out)
write.csv(tpm_out,"tpm_out.csv")



## 5.1计算矩阵中大于一的元素比例，画热图

calc_non_zero_ratio <- function(x) {
  return(sum(x > 1 ) / length(x))
}
tpm_out$Sample_Type <- ifelse(tpm_out$Sample_Type=="WT",0,1)
non_zero_ratios <- apply(tpm_out, 2, calc_non_zero_ratio)
selected_columns <- which(non_zero_ratios > 1/2)
filtered_mat <- tpm_out[, selected_columns]  # 选择 non_zero_ratio > 1/2 的列
b <- database_cut$Sample_Type
filtered_mat$condition <- b  # 将条件添加到数据框中
tpm_out_b <- filtered_mat
tpm_out_heat <- tpm_out_b
rownames(tpm_out_heat) <- tpm_out_heat[,1]  # 将第一列作为行名
tpm_out_heat <- tpm_out_heat[,-1]  # 删除第一列
tpm_out_heat <- tpm_out_heat[,-ncol(tpm_out_heat)]  # 删除最后一列
tpm_out_heat <- t(tpm_out_heat)  # 转置数据框
a <- rownames(tpm_out_heat)  # 获取行名
b <- res_up$SYMBOL  # 获取显著上调基因的SYMBOL
c <- res_dw$SYMBOL  # 获取显著下调基因的SYMBOL
filter_up <- intersect(a,b)  # 获取在两个数据框中共有的基因
filter_dw <- intersect(a,c)
res_up_path <- res_up[res_up$SYMBOL%in%filter_up,]  # 筛选显著上调的基因
res_up_path <- subset(res_up_path,select=c("SYMBOL","ENTREZID","log2FoldChange"))
res_dw_path <- res_dw[res_dw$SYMBOL%in%filter_dw,]  # 筛选显著下调的基因
res_dw_path <- subset(res_dw_path,select=c("SYMBOL","ENTREZID","log2FoldChange"))
col <- colorRampPalette(c("white", "#4CAC8E"))(70)
ann_colors <- list(condition = c(N="cyan2",P="darkmagenta"),gender=c(Male="cyan2",Female="darkmagenta"),age=col)#配色分组
coul <- colorRampPalette(brewer.pal(11, "PiYG"))(25)
new_df <- subset(df, select = c(Sample_Type))  #设置热图只展示Sample_Type的信息
plot4 <- pheatmap(tpm_diff, cluster_row=T,cluster_cols=F,cutree_cols= NA,treeheight_row = 0, cutree_rows = NA, scale="row",annotation_col=new_df,annotation_colors = ann_colors,color = coul,show_rownames = F,show_colnames = F,border_color = NA)
plot04 <- pheatmap(tpm_diff, cluster_row=T, scale="row",annotation_col=new_df,annotation_colors = ann_colors,color = coul,show_rownames = F,show_colnames = F,border_color = NA)

pdf(file = 'pheatmap.pdf',width = 8,height = 6)
print(plot4)
dev.off()
pdf(file = 'pheatmap01.pdf',width = 8,height = 6)
print(plot04)
dev.off()
#tpm_out_sc <- tpm_out[c("Row.names", top10gene, bottom10gene, "Sample_Type")]


### 定制热图
gene <- c("CCNB1", "RAD51", "RPA3",  
          "FANCA", "FANCB", "FANCD2", "FANCG", "FANCI", "BRCA1", "BRCA2", "PARP1","ID1",
          "FOXR1", "TRPC4", "VCAN", "CD4", "OSM", "MAPK4",
          "CDKN1B", "CDKN1C", "CASP16P", "PAPPA", "USP2", "FRG2", "BCL2L10"
) #"CDK1", "CDK4", "ZEB1", "BCL2", "CDH2", "DLAT", "XRCC4", "NAT10", "PI3K", "MYC", "FANCL"


gene_list <- gene_df[gene_df$SYMBOL %in% gene, ]

tpm_selected <- merge(tpm, gene_list, by.x = "row.names", by.y = "ENSEMBL")
rownames(tpm_selected) <- tpm_selected$SYMBOL
tpm_selected <-  tpm_selected %>% 
  dplyr::select(-c("Row.names","SYMBOL","ENTREZID"))
# colnames(tpm_selected) <- c("Con 1", "Con 2", "Con 3", 
#                             "Ola+Ida 1", "Ola+Ida 2", "Ola+Ida 3")

ann_colors <- list(
  Type = c("WT"="#4E78C4", "OR"="#CE2220")
)

# 修改样本顺序
if(F){
  cols_to_front <- c("SRR13489989", "SRR13489990", "SRR13489991", "SRR13489992",
                     "SRR13489993", "SRR13489994", "SRR13489995", "SRR13489996")
  tpm_selected <- tpm_selected[, c(cols_to_front, setdiff(names(tpm_selected), cols_to_front))]
  
}

pheatmap(tpm_selected,
         scale = "row",
         cluster_cols = F,
         show_colnames = T,
         annotation_colors = ann_colors,
         border = T,
         border_color = "gray",
         fontsize_row = 15,
         fontsize = 17)

#annotation_col_data = data.frame(Sample_Type = factor(rep(c("A1", "A2"), c(3,3))))
annotation_col_data <- database
annotation_col_data <- annotation_col_data %>% rename(Type = Sample_Type) # 修改列名
row.names(annotation_col_data) <- annotation_col_data[,1]
annotation_col_data[, 1] <- NULL
plot5 <- pheatmap(tpm_selected,
                  #color = colorRampPalette(c("#A1C8E3", "#FFE5E0", "#C42F40"))(60), #设置渐变色及等级 c("navy", "white", "firebrick3")
                  color = colorRampPalette(c("#4E78C4", "#ffffff", "#CE2220"))(60), #设置渐变色及等级 c("navy", "white", "firebrick3")
                  annotation_col = annotation_col_data,
                  scale = "row",
                  cluster_cols = F,
                  cluster_rows = T,
                  show_colnames = T,
                  angle_col = "45", #列名45°
                  annotation_colors = ann_colors,
                  border = T,
                  border_color = "gray",
                  fontsize_row = 14,
                  fontsize_col = 14,
                  fontsize = 12)

pdf(file = './plot_out/p3_8_pheatmap02.pdf',width = 12,height = 8)
print(plot5)
dev.off()


### 用resdata_symbol画热图
tpm_selected <- res_diff
tpm_selected <- filter(tpm_selected,!duplicated(tpm_selected$SYMBOL))
rownames(tpm_selected) <- tpm_selected$SYMBOL
tpm_selected <- tpm_selected %>%
  dplyr::select(-c(1:7, "SYMBOL", "ENTREZID"))






library(ComplexHeatmap)
library(tidyverse)
library(pheatmap)
library(corrplot)
library(circlize)
library(GGally)
library(linkET)
library(ggraph)
library(ggpubr)
library(ggside)
library(ggcor)
library(colorRamp2)

### 定制热图
gene <- c("CCNB1", "RAD51", "RPA3",  
          "FANCA", "FANCB", "FANCD2", "FANCG", "FANCI", "BRCA1", "BRCA2", "PARP1","ID1",
          "FOXR1", "TRPC4", "VCAN", "CD4", "OSM", "MAPK4",
          "CDKN1B", "CDKN1C", "CASP16P", "PAPPA", "USP2", "FRG2", "BCL2L10"
) #"CDK1", "CDK4", "ZEB1", "BCL2", "CDH2", "DLAT", "XRCC4", "NAT10", "PI3K", "MYC", "FANCL"


gene_list <- gene_df[gene_df$SYMBOL %in% gene, ]

tpm_selected <- merge(tpm, gene_list, by.x = "row.names", by.y = "ENSEMBL")
#tpm_selected <- tpm_selected[match(gene, tpm_selected$SYMBOL), , drop = FALSE]
# tpm_selected <- merge(countData, gene_list, by.x = "row.names", by.y = "ENSEMBL") ## tpm换为countData
rownames(tpm_selected) <- tpm_selected$SYMBOL
tpm_selected <-  tpm_selected %>% 
  dplyr::select(-c("Row.names","SYMBOL","ENTREZID"))
# colnames(tpm_selected) <- c("Con 1", "Con 2", "Con 3", 
#                             "Ola+Ida 1", "Ola+Ida 2", "Ola+Ida 3")
if(F){
  # common_genes <- intersect(colnames(lihc.tumor), ferroptosis_markers$Symbol)
  # hippo <- c("YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4",
  #            "SMAD1", "SMAD2", "SMAD3", "SMAD4")
  # cor_data <- cor(lihc.tumor[, common_genes], lihc.tumor[, hippo], method = 'spearman') %>%
  #   as.data.frame() %>%
  #   dplyr::filter(if_any(everything(), ~ abs(.) > 0.4))
}


# 计算相关性矩阵
cor_data <- cor(t(tpm_selected), t(tpm_selected), method = 'spearman') %>%
  as.data.frame() %>%
  dplyr::filter(if_any(everything(), ~ abs(.) > 0.4))

#p <- corr.test(t(tpm_selected), t(tpm_selected), method="spearman")

#my_palette  <-  colorRampPalette(c( "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582"))(50)
gene <- c('BRCA2','FANCD2','RAD51', "PARP1") #"PARP1", "FANCI"
lihc.tumor <- t(tpm_selected)


# 应该使用 correlate 来计算相关系数和显著性水平。注意需要将前两列名称改为 spec 和 env。
mantel.data <- ggcor::correlate(lihc.tumor[, gene], 
                                lihc.tumor[, rownames(cor_data)], 
                                method = "spearman", cor.test = TRUE) %>%
  ggcor::as_cor_tbl() %>%
  rename_with(~c("spec", "env"), .cols = 1:2) %>%
  mutate(
    size = cut(r, breaks = c(-Inf, 0.5, 0.7, Inf),
               labels = c("<0.5", "0.5~0.7", ">=0.7"),
               right = FALSE),
    p.value = cut(p.value, right = FALSE,
                  breaks = c(-Inf,  0.05, Inf),
                  labels = c("<0.05", ">=0.05"))
  )
my_palette  <-  colorRampPalette(c("#9dc6ff",  "#ffffff",  "#ff9999"))(30)

#"#B997C7","#824D99"

quickcor(lihc.tumor[, rownames(cor_data)], type = "upper") + 
  geom_square() + 
  add_link(mantel.data, mapping = aes(colour = p.value, size = size)) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  # remove_axis("x") + # 去除上方的基因名
  scale_fill_gradientn(colours = my_palette) +
  scale_colour_manual(values = c("#A1C9F0", "#dcb0f2"))

pdf(file = "./relation_heatmap.pdf", width = 14, height = 10)
quickcor(lihc.tumor[, rownames(cor_data)], type = "upper") + 
  geom_circle2() +  #geom_square geom_star geom_shade geom_pie2 geom_mark geom_link2 geom_hc_rect geom_elipse2 geom_cross geom_confbox geom_colour geom_circle2
  add_link(mantel.data, mapping = aes(colour = p.value, size = size)) +
  scale_size_manual(values = c(0.4, 0.8, 1.2)) +
  scale_fill_gradientn(colours = my_palette) +
  scale_colour_manual(values = c("#bfcfff", "#dcb0f2"))+
  theme(
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black")
  )
dev.off()
confbox
# R 自带的 pairs 绘图函数来绘制相关性图。
panel.raters<-function(x,y,corr=NULL,...){
  if(!is.null(corr))
    return()
  plot.xy(xy.coords(x,y),
          type="p",#点
          pch=21,#点的形状
          cex=1.5,#点的大小
          bg="#4B8643",#点的填充色
          col="#4B8643",#点的边颜色
          lwd=0.5,
          ...)
  abline(lm(y~x),lwd=2,col="#f4b183")#画拟合线
  box(col="black",lwd=2)#黑色粗边框
}

#上三角相关系数热图
panel.fill.cor <- function(x,y,corr=NULL,...)
{
  #计算相关系数，可以换成"kendall"或"pearson", "spearman"
  corr<-round(cor(x,y,use="pairwise",method="spearman"),2)
  pal<-my_palette#颜色范围
  ncol<-length(pal)
  col.ind<-as.numeric(cut(corr,breaks=seq(from=-1,to=1,length.out=ncol+1),
                          include.lowest=TRUE))
  
  #画背景
  par(new=TRUE)
  plot(0,type='n',xlim=c(-1,1),ylim=c(-1,1),axes=FALSE,asp=1)
  usr<-par("usr")
  rect(usr[1],usr[3],usr[2],usr[4],col=pal[col.ind],
       border=NA)
  
  #相关系数值
  text(0,0,labels=corr,cex=2.5,
       col=ifelse(corr>0,"black","white"))
  box(col="black")#边框颜色
}
#对角线文本
textPanel <- function(x=0.5,y=0.5,txt,cex,font){
  text(x,y,txt,cex=cex,font=font)
  box(col="black",lwd=2)
}

# 该函数不适合展示太多的基因。
par(bg = "#fdfdfd")
gene2 <- c('FANCD2','BRCA2','PARP1','RAD51', 'FANCI','CDKN1B', 'USP2', 'CASP16P')
pdf(file = "./relation_heatmap1.pdf", width = 8, height = 6)
pairs(lihc.tumor[,gene2], 
      gap = .5,                      # 小图之间的空隙
      text.panel = textPanel,        # 对角线文本
      lower.panel = panel.raters,    # 左下角
      upper.panel = panel.fill.cor)  # 右上角
dev.off()

