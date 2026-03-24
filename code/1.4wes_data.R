rm(list = ls());gc()
setwd("/home/Rstudio/yhb/Data/WES/")
library(future)
set.seed(101)
plan("multisession", workers = 2) 
options(future.globals.maxSize = 10000 * 1024^2) # set RAM

library(maftools)
options(stringsAsFactors = F)
## annovar
annovar.laml <- annovarToMaf(annovar = "./d8annotation/annovar/annovar_merge.vcf", 
                             refBuild = 'hg38',
                             tsbCol = 'Tumor_Sample_Barcode', 
                             table = 'refGene',
                             MAFobj = T)
## gatk
# gatk.laml = read.maf(maf = 'gatk/gatk4.1.4.0_merge.maf')

# library(data.table)
# tmp=fread('./7.annotation/funcatator/funcatator_merge.maf')
# gatk.laml = read.maf(maf = tmp)

## vep
vep.laml = read.maf(maf = './d8annotation/vep/vep_merge.maf')
## for vep.laml
library(stringr)
vep.laml@data$Protein_Change = paste0("p.",
                                      str_sub(vep.laml@data$Amino_acids,1,1),
                                      vep.laml@data$Protein_position,
                                      str_sub(vep.laml@data$Amino_acids,3,3))

## save Rdata
#save(annovar.laml, gatk.laml, vep.laml, file = 'laml.Rdata')
save(annovar.laml, vep.laml, file = 'laml.Rdata')

## Summary
laml=annovar.laml
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml)
getGeneSummary(laml)
getFields(laml)

# laml=gatk.laml
# unique(laml@data$Tumor_Sample_Barcode)
# getSampleSummary(laml)
# getGeneSummary(laml)
# getFields(laml)

laml=vep.laml
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml)
getGeneSummary(laml)
getFields(laml)


load(file = 'laml.Rdata')
#laml=c(annovar.laml,gatk.laml,vep.laml)
laml=c(annovar.laml,vep.laml)

## mafsummary
#anno=c('annovar','gatk','vep')
anno=c('annovar','vep')
for (i in 1:2) {
  #i = 1
  #png(paste0('plotmafSummary_', anno[i], '.png'),res = 150,width = 1080,height = 1080)
  pdf(paste0('plotmafSummary_', anno[i], '.pdf'),width = 8,height = 6)
  plotmafSummary(maf = laml[[i]],
                 rmOutlier = TRUE,
                 showBarcodes = T,
                 textSize = 0.4,
                 addStat = 'median',
                 dashboard = TRUE,
                 titvRaw = FALSE)
  dev.off()
}

## oncoplot_top30
for (i in 1:2) {
  #i = 1
  #png(paste0('oncoplot_top30_', anno[i], '.png'),res = 150,width = 1080,height = 1080)
  pdf(paste0('oncoplot_top30_', anno[i], '.pdf'),width = 3,height = 8)
  oncoplot(maf = laml[[i]],
           top = 30,
           fontSize = 0.5,
           sampleOrder = laml[[i]]@clinical.data$Tumor_Sample_Barcode,
           showTumorSampleBarcodes = T)
  dev.off()
}


## lollipopPlot for SP140
gene='MAPK11'
protein=c("AAChange.refGene","Protein_Change","Protein_Change")
for (i in 1:2) {
  #i=3
  #png(paste0(gene,'_', anno[i], '.png'),res = 150,width = 1080,height = 1080)
  pdf(paste0(gene,'_', anno[i], '.pdf'),width = 8,height = 6)
  maftools::lollipopPlot(maf = laml[[i]],
                         gene = gene,
                         AACol = protein[i],
                         labelPos = 'all')
  dev.off()
}

## lollipopPlot for TP53
gene='BRAP'
protein=c("AAChange.refGene","Protein_Change","Protein_Change")
for (i in 1:2) {
  #i=3
  #png(paste0(gene,'_', anno[i], '.png'),res = 150,width = 1080,height = 1080)
  pdf(paste0(gene,'_', anno[i], '.pdf'),width = 8,height = 6)
  maftools::lollipopPlot(maf = laml[[i]],
                         gene = gene,
                         AACol = protein[i],
                         labelPos = 'all')
  dev.off()
}




# options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# devtools::install_github('raerose01/deconstructSigs')
# BiocManager::install('BSgenome', force = TRUE)
# BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggrepel)
## https://github.com/raerose01/deconstructSigs

## 读入数据
maf <- read.table('./d8annotation/vep/vep_merge.maf',header = T,sep = '\t',quote = "")
maf[1:5,1:5]

dir.create("./1.4result")
### 补充
if(F){
  DLD_WES <- maf
  
  mut_genes <- unique(DLD_WES$Hugo_Symbol)
  
  # ==== 1. 计算突变频率 ====
  DLD_WES <- DLD_WES %>%
    mutate(
      VAF_tumor = ifelse(t_depth > 0, t_alt_count / t_depth, NA),
      VAF_normal = ifelse(n_depth > 0, n_alt_count / n_depth, NA),
      delta_VAF = VAF_tumor - VAF_normal
    )
  #耐药细胞相对于对照细胞中的等位基因频率变化
  
  # ==== 1. 突变类型筛选 ====
  mutation_types <- c("Frame_Shift_Ins", "Frame_Shift_Del", 
                      "Missense_Mutation", "Nonsense_Mutation", "Splice_Site")
  
  DLD_WES <- DLD_WES %>%
    filter(Variant_Classification %in% mutation_types) # 437到317个
  
  
  # ==== 4. Fisher 精确检验 ====
  DLD_WES <- DLD_WES %>%
    rowwise() %>%
    mutate(
      p_value = fisher.test(matrix(
        c(t_alt_count, t_ref_count, n_alt_count, n_ref_count),
        nrow = 2
      ))$p.value
    ) %>%
    ungroup() %>%
    mutate(
      log10_p = -log10(p_value),
      MutType = case_when(
        p_value < 0.05 & delta_VAF > 0 ~ "Treat enrich",
        p_value < 0.05 & delta_VAF < 0 ~ "C enrich",
        TRUE ~ "None sign"
      )
    )
  
  # ==== 火山图 ====
  # ==== 5. 标注显著突变（可选）====
  gene_id <- c("ERCC6L", "RAD51AP2", "PPARGC1A", "KLF14")
  top_labels <- subset(DLD_WES, Hugo_Symbol %in% gene_id)
  
  # top_labels <- DLD_WES %>%
  #   arrange(p_value) %>%
  #   filter(p_value < 0.01 & abs(delta_VAF) > 0.2) %>%
  #   head(10)
  
  plot07 <- ggplot(DLD_WES, aes(x = delta_VAF, y = log10_p, color = MutType)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_text_repel(
      data = top_labels,
      aes(label = Hugo_Symbol),
      size = 3, max.overlaps = 10
    ) +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
    scale_color_manual(values = c(
      "Resistant-enriched" = "#E74C3C",
      "WildType-enriched" = "#2E86C1",
      "Not significant" = "grey70"
    )) +
    labs(
      title = "Significant Mutations",
      x = expression(Delta~VAF~"(Treat - C)"),
      y = expression(-log[10](pvalue))
    ) +
    theme_minimal(base_size = 12)
  
  pdf("1.4result/p1_VAF_Volcano_Plot.pdf", width = 6, height = 4)
  print(plot07)
  dev.off()
  
  plot07 <- ggplot(DLD_WES, aes(x = delta_VAF, y = log10_p, color = MutType)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_text_repel(data = top_labels,   # 添加带边框的标签,上调基因
                     # geom_text_repel(data = subset(resdata_symbol, SYMBOL %in% gene_vol),   # 基因不带框
                     aes(label = Hugo_Symbol),
                     seed = 233,                           # 设置随机数种子，用于确定标签位置
                     size = 3.5,                           # 设置标签的字体大小
                     color = 'black',                      # 设置标签的字体颜色
                     min.segment.length = 0,               # 始终为标签添加指引线段
                     #min.segment.length = Inf, # 去除引线
                     force = 2,                            # 设置标签重叠时的排斥力
                     force_pull = 3,                       # 设置标签与数据点之间的吸引力
                     box.padding = 0.4,                    # 设置标签周围的填充量
                     max.overlaps = Inf,                   # 设置排斥重叠过多标签的阈值为无穷大，保持始终显示所有标签
                     # 两侧排列基因
                     segment.linetype = 2,                 # 设置线段类型为虚线
                     segment.color = 'black',              # 设置线段颜色 'black'
                     segment.alpha = 0.8,                  # 设置线段不透明度
                     nudge_x = 0,   # 调整标签x轴起始位置
                     nudge_y = 80, # 调整标签y轴起始位置
                     direction = "y",                      # 按y轴调整标签位置方向，若想水平对齐则为x
                     hjust = 1                         # 对齐标签：0右对齐，1左对齐，0.5居中
    ) +
    geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
    scale_color_manual(values = c(
      "Treat enrich" = "#E74C3C",
      "C enrich" = "#2E86C1",
      "None sign" = "grey70"
    )) +
    labs(
      title = "Significant mutation",
      x = expression(Delta~VAF~"(Treat - C)"),
      y = expression(-log[10](pvalue))
    ) +
    theme_minimal(base_size = 12)
  # #=====新的基因列表加和====
  # genes <- unique(c(genes, top_labels$Hugo_Symbol))
  # 
  # wes_sub <- DLD_WES %>%
  #   filter(Hugo_Symbol %in% genes)
  # 
  # # ==== 2. 计算每个基因、每种突变类型的平均 delta_VAF ====
  # wes_summary <- wes_sub %>%
  #   group_by(Hugo_Symbol, Variant_Classification) %>%
  #   summarise(mean_delta_VAF = mean(delta_VAF, na.rm = TRUE), .groups = "drop")
  # 
  # # ==== 3. 确保所有关注基因都出现（未突变补0） ====
  # all_combinations <- expand.grid(
  #   Hugo_Symbol = genes,
  #   Variant_Classification = mutation_types,
  #   stringsAsFactors = FALSE
  # )
  # 
  # wes_summary <- all_combinations %>%
  #   left_join(wes_summary, by = c("Hugo_Symbol", "Variant_Classification")) %>%
  #   mutate(mean_delta_VAF = ifelse(is.na(mean_delta_VAF), 0, mean_delta_VAF))
  # 
  # # ==== ★★ 新增：筛除没有突变的基因 ====
  # wes_summary <- wes_summary %>%
  #   group_by(Hugo_Symbol) %>%
  #   mutate(total_delta = sum(mean_delta_VAF)) %>%
  #   ungroup() %>%
  #   filter(total_delta > 0) %>%     # ← ← 只保留有突变的基因
  #   select(-total_delta)
  # 
  # # ==== 4. 计算每个基因的总突变频率（用于排序） ====
  # wes_summary <- wes_summary %>%
  #   group_by(Hugo_Symbol) %>%
  #   mutate(total_VAF = sum(mean_delta_VAF, na.rm = TRUE)) %>%
  #   ungroup()
  # 
  # # ==== 5. 按总突变频率从大到小排序 ====
  # wes_summary$Hugo_Symbol <- factor(
  #   wes_summary$Hugo_Symbol,
  #   levels = wes_summary %>%
  #     distinct(Hugo_Symbol, total_VAF) %>%
  #     arrange(total_VAF) %>%
  #     pull(Hugo_Symbol)
  # )
  # 
  # # ==== 6. 定义配色 ====
  # mutation_colors <- c(
  #   "Frame_Shift_Ins"  = "#FF7F00BB",  # 半透明橙
  #   "Frame_Shift_Del"  = "#984EA3BB",  # 半透明紫
  #   "Missense_Mutation" = "#1B9E77BB", # 半透明青绿
  #   "Nonsense_Mutation" = "#B0C4DEBB", # 淡蓝
  #   "Splice_Site"       = "#FFC0CBBB"  # 淡粉
  # )
  # 
  # # ==== 7. 绘制纵向堆叠条形图 ====
  # plot08 <- ggplot(wes_summary, aes(y = Hugo_Symbol, x = mean_delta_VAF, fill = Variant_Classification)) +
  #   geom_bar(stat = "identity", position = "stack", width = 0.6, color = "black", linewidth = 0.3) +
  #   geom_text(
  #     aes(label = ifelse(mean_delta_VAF > 0, round(mean_delta_VAF, 2), "")),
  #     position = position_stack(vjust = 0.5), size = 3.5, color = "black"
  #   ) +
  #   scale_fill_manual(values = mutation_colors) +
  #   theme_minimal(base_size = 14) +
  #   labs(
  #     y = "Genes",
  #     x = expression(Delta~VAF~"(Resistance - Wild_Type)"),
  #     fill = "Mutation Type",
  #     title = "Mutation Profile"
  #   ) +
  #   theme(
  #     axis.text.y = element_text(face = "bold"),
  #     panel.grid.minor = element_blank(),
  #     panel.grid.major.x = element_blank(),  # 去除纵线
  #     panel.grid.major.y = element_blank(),  # 去除横线
  #     plot.title = element_text(face = "bold", hjust = 0.5)
  #   )
  # pdf("07-Mutation Profile.pdf", width = 10, height = 7)
  # print(plot08)
  # dev.off()
  # 
  # WES_output <- DLD_WES[DLD_WES$Hugo_Symbol%in%(wes_summary$Hugo_Symbol),]
  # WES_output$Exon_Number <- NULL
  # write.csv(WES_output,"all_WES.csv")
}

## 构建肿瘤突变数据框，需要5列信息: sample.ID,chr,pos,ref,alt 
sample.mut.ref <- data.frame(Sample=maf[,16], 
                             chr = maf[,5],
                             pos = maf[,6],
                             ref = maf[,11],
                             alt = maf[,13])

sample.mut.ref[1:5,1:5]

sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref,
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)

#这样就拿到了每个样本的三联核苷酸突变类型分布了，如前 5 行前 5 列
sigs.input[1:3,1:3]


sigs.output <- whichSignatures(tumor.ref = sigs.input,
                               signatures.ref = signatures.cosmic, 
                               sample.id = 'SRR26823039',
                               contexts.needed = TRUE)


# Plot output
pdf("./sigs.pdf", width = 8, height = 6)
plotSignatures(sigs.output)
dev.off()

pdf("./sigs_pie.pdf", width = 8, height = 6)
makePie(sigs.output)
dev.off()


# pheatmap
df = data.frame()
for (i in rownames(sigs.input)) {
  sigs.output <- whichSignatures(tumor.ref = sigs.input,
                                 signatures.ref = signatures.cosmic, 
                                 sample.id = i,
                                 contexts.needed = TRUE)
  df = rbind(df,sigs.output$weights)
}
df = df[ , apply(df, 2, function(x){sum(x>0)})>0]

pdf("./sigs_heatmap.pdf", width = 8, height = 4)
pheatmap::pheatmap(df,cluster_cols = F,cluster_rows = F,fontsize = 16)
dev.off()


### 肿瘤异质性
rm(list = ls())
library(dplyr)
library(stringr)
library(tidyr)
# 读入数据
laml = read.maf('./d8annotation/vep/vep_merge.maf')
laml@data = laml@data[!grepl('^MT-', laml@data$Hugo_Symbol),] 
# 增加一列t_vaf，即肿瘤样本中突变位点的覆盖深度t_alt_count占测序覆盖深度t_depth的比值
laml@data$t_vaf = (laml@data$t_alt_count/laml@data$t_depth)
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml) 
getGeneSummary(laml) 
getFields(laml)


mut = laml@data[laml@data$t_alt_count >= 5 &
                  laml@data$t_vaf >= 0.05, c("Hugo_Symbol",
                                             "Chromosome",
                                             "Start_Position",
                                             "Tumor_Sample_Barcode",
                                             "t_vaf")]

mut$patient = substr(mut$Tumor_Sample_Barcode, 1, 9)

mut[1:6,1:6]

pid = unique(mut$patient)


## 红色是snv位点
lapply(pid , function(p){
  p = 'SRR268230'
  print(p)
  mat = unique(mut[mut$patient == p,c("Tumor_Sample_Barcode",'Hugo_Symbol')]) 
  mat$tmp = 1
  # 长变扁
  mat = spread(mat,Tumor_Sample_Barcode,tmp,fill = 0)
  class(mat)
  mat = as.data.frame(mat)
  rownames(mat) = mat$Hugo_Symbol
  mat = mat[ , -1, drop = FALSE]
  #dat = mat[order(mat[,1],mat[,2],mat[,3],mat[,4]),]
  dat = mat[order(mat[,1]),]
  #png(paste0('overlap_', p, '.png'),width = 300, height = 1500)
  pdf(paste0('overlap_', p, '.pdf'),width = 3,height = 16)
  pheatmap::pheatmap(mat = mat, cluster_cols = F, cluster_rows = T, show_rownames = T, legend = F,fontsize_row = 10,fontsize_col = 12)
  dev.off()
})


# 肿瘤样本中的异质性：变异等位基因频率（Variant Allele Frequecy） t_vaf
laml = read.maf('./d8annotation/vep/vep_merge.maf')
laml@data = laml@data[!grepl('^MT-',laml@data$Hugo_Symbol),] 
# 增加一列 t_vaf，即肿瘤样本中突变位点的覆盖深度t_alt_count占测序覆盖深度t_depth的比值
laml@data$t_vaf = (laml@data$t_alt_count/laml@data$t_depth)
heter = inferHeterogeneity(maf = laml, top = 3, vafCol = 't_vaf') # 3个样本
heter$clusterData
heter$clusterMeans

plotClusters(clusters = heter,
             showCNvars = T,
             tsb = 'SRR26823039')

for (i in unique(laml@data$Tumor_Sample_Barcode)) {
  # i = 'SRR26823039'
  #png(paste0('vaf_', i, '.png'),width = 500, height = 330)
  pdf(paste0('vaf_', i, '.pdf'),width = 8,height = 6)
  plotClusters(clusters = heter,
               tsb = i,
               showCNvars = T)
  dev.off()
}


### KEGG注释
laml = read.maf('./d8annotation/vep/vep_merge.maf')
laml@data=laml@data[!grepl('^MT-',laml@data$Hugo_Symbol),] 
# 增加一列t_vaf，即肿瘤样本中突变位点的覆盖深度t_alt_count占测序覆盖深度t_depth的比值
laml@data$t_vaf = (laml@data$t_alt_count/laml@data$t_depth)
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml) 
getGeneSummary(laml) 
getFields(laml)

mut = laml@data[laml@data$t_alt_count >= 5 &
                  laml@data$t_vaf >= 0.05, c("Hugo_Symbol",
                                             "Chromosome",
                                             "Start_Position",
                                             "Tumor_Sample_Barcode",
                                             "t_vaf")]

mut$patient = substr(mut$Tumor_Sample_Barcode, 1, 9)
mut[1:6,1:6]

# 得到每一个病人的 3 个样本（列）和突变基因（行）的矩阵
pid = unique(mut$patient)

all_snv = lapply(pid , function(p){
  # p='SRR268230'
  print(p)
  mat=unique(mut[mut$patient %in% p,c("Tumor_Sample_Barcode",'Hugo_Symbol')]) 
  mat$tmp = 1
  # 长变扁
  mat = spread(mat,Tumor_Sample_Barcode,tmp,fill = 0)
  class(mat)
  mat = as.data.frame(mat)
  rownames(mat) = mat$Hugo_Symbol
  mat=mat[,-1]
  #dat = mat[order(mat[,1],mat[,2],mat[,3],mat[,4]),]
  dat = mat[order(mat[,1]),]
  return(dat)
})

# 突变分类：样本出现次数
trunk_gene = unlist(sapply(all_snv, function(x) rownames(x[rowSums(x) == 3,])))
#branch_gene = unlist(sapply(all_snv, function(x) rownames(x[rowSums(x) == 3|2,])))
branch_gene = unlist(sapply(all_snv, function(x) rownames(x[rowSums(x) == 2,])))
private_gene = unlist(sapply(all_snv, function(x) rownames(x[rowSums(x) == 1,])))
trunk_gene
branch_gene
private_gene


library(org.Hs.eg.db)
library(clusterProfiler)
kegg_SYMBOL_hsa <- function(genes){ 
  gene.df <- bitr(genes, fromType = "SYMBOL",
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)
  head(gene.df) 
  diff.kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                        organism     = 'hsa',
                        pvalueCutoff = 0.99,
                        qvalueCutoff = 0.99
  )
  return(setReadable(diff.kk, OrgDb = org.Hs.eg.db,keyType = 'ENTREZID'))
}

# 对于trunk_gene
trunk_kk=kegg_SYMBOL_hsa(trunk_gene)
trunk_df=trunk_kk@result
write.csv(trunk_df,file = 'trunk_kegg.csv')
png(paste0('trunk_kegg', '.png'),width = 1080,height = 540)
barplot(trunk_kk,font.size = 20)
dev.off()

# 对于branch_gene
branch_kk=kegg_SYMBOL_hsa(branch_gene)
branch_df=branch_kk@result
write.csv(branch_df,file = 'branch_kegg.csv')
png(paste0('branch_kegg', '.png'),width = 1080,height = 540)
barplot(branch_kk,font.size = 20)
dev.off()

# 对于private_gene
private_kk=kegg_SYMBOL_hsa(private_gene)
private_df=private_kk@result
write.csv(private_df,file = 'private_kegg.csv')
png(paste0('private_kegg', '.png'),width = 1080,height = 540)
barplot(private_kk,font.size = 20)
dev.off()



### ABSOLUTE可以直接从体细胞突变和拷贝数变异的结果中推断肿瘤纯度和倍性
library(ABSOLUTE)
dir.create("./absolute_re")

RunAbsolute(
  seg.dat.fn = "./d9cnv/gatk/segments/SRR26823039.cr.igv.seg",
  sigma.p = 0,
  max.sigma.h = 0.2,
  min.ploidy = 0.5,
  max.ploidy = 8,
  primary.disease = "OC",
  platform = "Illumina_WES",  #Illumina_WES SNP_6.0
  results.dir = "./absolute_re/SRR26823039",
  sample.name = "SRR26823039",
  max.as.seg.count = 1500,
  max.neg.genome = 0.005,
  max.non.clonal = 1,
  copy_num_type = "total",
  maf.fn = "./d8annotation/vep/absolute_maf/SRR26823039_absolute.maf",
  #maf.fn = "./d8annotation/vep/SRR26823039_vep.maf",
  min.mut.af = 0.05)

## 批量处理
config = read.table('config3.txt')
for (i in config$V1) {
  # i = "SRR26823039"
  segmets = paste0("./d9cnv/gatk/segments/", i, ".cr.igv.seg")
  maf = paste0("./d8annotation/vep/absolute_maf/", i, "_absolute.maf")
  RunAbsolute(
    seg.dat.fn = segmets, 
    sigma.p = 0,
    max.sigma.h = 0.2,
    min.ploidy = 0.5, 
    max.ploidy = 8, 
    primary.disease = "OC", 
    platform = "Illumina_WES", 
    results.dir = paste0("./absolute_re/", i), 
    sample.name = i, 
    max.as.seg.count = 1500, 
    max.neg.genome = 0.005,
    max.non.clonal = 1, 
    copy_num_type = "total",
    maf.fn = maf,
    min.mut.af = 0.05
  )
}


### pyclone推断肿瘤细胞的克隆组成
rm(list = ls())
options(stringsAsFactors = F)
case1_loci = read.table("./d10pyclone/allSRR_pyclone_analysis/tables/loci.tsv",
                        header = T)
# 获取clusters 的分组信息
clusters_list = unique(case1_loci[, c(1, 3)])
rownames(clusters_list) = clusters_list[, 1]
cluster_id = data.frame(cluster_id = as.character(clusters_list$cluster_id))
rownames(cluster_id) = clusters_list[, 1]
# 获取同一个突变位点在不同样本中的cellular_prevalence，然后画热图可视化
library(tidyr)
cellular_prevalence = spread(case1_loci[, c(1, 2, 4)], key = sample_id, value = cellular_prevalence)
rownames(cellular_prevalence) = cellular_prevalence[, 1]
cellular_prevalence = cellular_prevalence[,-1]
sampe_id = colnames(cellular_prevalence)
cellular_prevalence = as.data.frame(t(apply(cellular_prevalence, 1, as.numeric)))
colnames(cellular_prevalence) = sampe_id

cluster_id$cluster_id <- as.character(cluster_id$cluster_id)
ann_colors <- list(
  cluster_id = c("0"="#D4477D","1"="#7DD06F", "2"="#844081", "3"="#688EC1", "4"="#C17E73", "5"="#484125", "6"="#6CD3A7", "7"="#597873","8"="#7B6FD0", "9"="#CF4A31", "10"="#D0CD47",
                 "11"="#722A2D", "12"="#CBC594", "13"="#D19EC4", "14"="#5A7E36")
)
pdf("pycolne_heatmap.pdf",width = 4,height = 6)
pheatmap::pheatmap(
  cellular_prevalence,
  annotation_row = cluster_id,
  annotation_colors = ann_colors,
  show_rownames = F,
  cluster_cols = F,
  cluster_rows = T,
  clustering_method = 'median',
  angle_col = 45)
dev.off()

# 获取同一个突变位点在不同样本中的variant_allele_frequency，也就是vaf，同样可视化，为了聚类，采用了不同的聚类方法
library(tidyr)
variant_allele_frequency = spread(case1_loci[, c(1, 2, 6)], key = sample_id, value = variant_allele_frequency)
rownames(variant_allele_frequency) = variant_allele_frequency[, 1]
variant_allele_frequency = variant_allele_frequency[,-1]
sampe_id = colnames(variant_allele_frequency)
variant_allele_frequency = as.data.frame(t(apply(variant_allele_frequency, 1, as.numeric)))
colnames(variant_allele_frequency) = sampe_id

cluster_id$cluster_id <- as.character(cluster_id$cluster_id)
ann_colors <- list(
  cluster_id = c("0"="#D4477D","1"="#7DD06F", "2"="#844081", "3"="#688EC1", "4"="#C17E73", "5"="#484125", "6"="#6CD3A7", "7"="#597873","8"="#7B6FD0", "9"="#CF4A31", "10"="#D0CD47",
                 "11"="#722A2D", "12"="#CBC594", "13"="#D19EC4", "14"="#5A7E36")
)
pdf("pycolne_heatmap_vaf.pdf",width = 4,height = 6)
pheatmap::pheatmap(
  cellular_prevalence,
  annotation_row = cluster_id,
  annotation_colors = ann_colors,
  show_rownames = F,
  cluster_cols = F,
  clustering_method = 'average',
  angle_col = 45)
dev.off()

cor(case1_loci[,c(4,6)])
#cellular_prevalence variant_allele_frequency
#cellular_prevalence                1.0000000                0.9374382
#variant_allele_frequency           0.9374382                1.0000000



### 肿瘤进化分析
library(timescape)
options(stringsAsFactors = F)
example("timescape")
browseVignettes("timescape") 
library(plotly)
library(htmlwidgets)
library(webshot)

tree_edges = read.table("d10pyclone/allSRR_pyclone_analysis/tree.txt")
colnames(tree_edges) = c("source","target")
# clonal prevalences
cellfreq = read.table("d10pyclone/allSRR_pyclone_analysis/cellfreq.txt")
colnames(cellfreq) = 0:(length(cellfreq)-1)
sample_id = read.table("d10pyclone/allSRR_pyclone_analysis/sample_id")
cellfreq$timepoint = sample_id[ , 1]
library(tidyr)
clonal_prev = gather(cellfreq, key="clone_id", value = "clonal_prev", -timepoint)
clonal_prev = clonal_prev[order(clonal_prev$timepoint),]
# targeted mutations
# mutations <- read.csv(system.file("extdata", "AML_mutations.csv", package = "timescape"))
p = timescape(clonal_prev = clonal_prev, tree_edges = tree_edges,height=260)
p
saveWidget(p, "A2780_timescape.html")


# for (i in 1:6) {
#   # i = 1
#   # tree 
#   tree_edges = read.table(paste0("9.pyclone/case", i, "_pyclone_analysis/tree.txt"))
#   colnames(tree_edges) = c("source","target")
#   # clonal prevalences
#   cellfreq = read.table(paste0("9.pyclone/case", i, "_pyclone_analysis/cellfreq.txt"))
#   colnames(cellfreq) = 0:(length(cellfreq)-1)
#   sample_id = read.table(paste0("9.pyclone/case", i, "_pyclone_analysis/sample_id"))
#   cellfreq$timepoint = sample_id[ , 1]
#   library(tidyr)
#   clonal_prev = gather(cellfreq, key="clone_id", value = "clonal_prev", -timepoint)
#   clonal_prev = clonal_prev[order(clonal_prev$timepoint),]
#   # targeted mutations
#   # mutations <- read.csv(system.file("extdata", "AML_mutations.csv", package = "timescape"))
#   p = timescape(clonal_prev = clonal_prev, tree_edges = tree_edges,height=260)
#   saveWidget(p, paste0("case", i,"_timescape", ".html"))
# }




### 拷贝数变异 CNVkit
library(data.table)
library(ComplexHeatmap)
library(circlize)

# 读取多个样本的cnr
files <- list.files("./d9cnv/cnvkit/", pattern = "*.cnr", full.names = TRUE)

data_list <- lapply(files, fread)

# 给每个样本加上名字
names(data_list) <- gsub("_bqsr.cnr", "", basename(files))

# 1) 给每个样本挑出 log2 列
dt_list <- Map(function(dt, nm) {
  dt2 <- dt[, .(chromosome, start, end, log2)]
  setnames(dt2, "log2", nm)
  dt2
}, data_list, names(data_list))

# 2) 合并
merged <- Reduce(function(x,y) merge(x,y, by=c("chromosome","start","end")), dt_list)

# 3) 做矩阵（行=样本，列=bin）
# mat <- as.matrix(t(merged[, -c(1:3), with=FALSE]))
# rownames(mat) <- names(data_list)
# colnames(mat) <- paste0(merged$chromosome, ":", merged$start)


# 定义颜色
#col_fun = colorRamp2(c(-1,0,1), c("blue","white","red"))
col_fun = colorRamp2(c(-1,0,1), c("#000A87", "#ffffff", "#900512"))

# 染色体断点
chrom <- merged$chromosome
split_factor <- factor(chrom, levels = as.character(1:22))

# 绘制热图
# 只保留 chr1-22
merged2 <- merged[chromosome %in% paste0("chr", 1:22)]

# 提取矩阵：行=样本，列=bin
mat <- as.matrix(t(merged2[, -c(1:3), with=FALSE]))
#rownames(mat) <- names(data_list)
rownames(mat) <- c("T320", "T40", "T10")
colnames(mat) <- paste0(merged2$chromosome, ":", merged2$start)

# 定义颜色
#col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# 染色体分割因子
chrom <- merged2$chromosome
split_factor <- factor(chrom, levels = paste0("chr", 1:22))
levels(split_factor) <- sub("^chr", "", levels(split_factor)) # 去掉chr前缀

#sample_split <- c("SRR26823039", "SRR26823040", "SRR26823041")
sample_split <- factor(rownames(mat), levels = c("T10", "T40", "T320"))

# 画热图
Heatmap(mat,
        name = "CNA ratio",
        col = col_fun,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = TRUE, show_column_names = FALSE,
        column_split = split_factor, border = TRUE)

ht <- Heatmap(mat, show_column_names = F, show_row_dend = F, 
             show_column_dend = F, show_row_names = F, 
             name = "CNA ratio",
             row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 8),
             column_title_side = c("bottom"),
             col = col_fun,
             # col = colorRamp2(seq(-1, 1, length.out = 30), 
             #                  colorRampPalette(rev(brewer.pal(9, "RdBu")))(30)),
             cluster_rows = F, cluster_columns = F, border = "black",
             column_split = split_factor, row_split = sample_split, 
             #gap = unit(3, "mm"), # 染色体之间的间隔
             row_gap = unit(1, "mm"), column_gap = unit(1, "mm"), 
             #heatmap_width = unit(12, "cm"), heatmap_height = unit(10, "cm"),
             heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                         title = "CNA ratio",
                                         labels_gp = gpar(fontsize = 8),
                                         legend_height = unit(2, "cm")))
ht
#Cairo::CairoPDF("./p6_heat_cnvratio.pdf", width = 20, height = 10)
pdf("./1.4result/p2_heat_cnvratio.pdf", width = 10.5, height = 1.5)
draw(ht)
dev.off()

### 拷贝数变异 facets
if(F){
  install.packages('devtools')
  devtools::install_github("mskcc/facets") #, build_vignettes = TRUE
  devtools::install_github("mskcc/pctGCdata")
}
options(stringsAsFactors = F)
library("facets")
# check if .Random.seed exists
seedexists <- exists(".Random.seed")
# save seed
if(seedexists) 
  oldSeed <- .Random.seed
# Alway use the same random seed
set.seed(0xfade)
# read the data
if(F){
  library(vcfR)
  vcf_file='./d6gatk/merge.vcf'
  ### 直接读取群体gvcf文件即可
  vcf <- read.vcfR(vcf_file, verbose = FALSE)
  save(vcf,file = 'input_vcf.Rdata')
}
load(file = 'input_vcf.Rdata')
vcf@fix[1:4,1:4]
vcf@gt[1:14,1:4]
colnames(vcf@gt)
library(stringr)
tmp_f <- function(x){
  #x=vcf@gt[,2]
  gt=str_split(x,':',simplify = T)[,1]
  ad=str_split(str_split(x,':',simplify = T)[,2],',',simplify = T)[,2]
  dp=str_split(x,':',simplify = T)[,3]
  return(data.frame(gt=gt,dp=as.numeric(dp),ad=as.numeric(ad)))
}
colnames(vcf@gt)
tms=colnames(vcf@gt)[!grepl('SRR26823042',colnames(vcf@gt))][-1] # 去除对照正常样本
tmp <- lapply(tms, function(tm){
  #tm=tms[1]
  print(tm)
  #nm=paste0(strsplit(tm,'_')[[1]][1],'_germline')
  nm="SRR26823042"
  print(nm)
  if(nm %in% colnames(vcf@gt)){
    snpm=cbind(vcf@fix[,1:2],
               tmp_f(vcf@gt[,nm]),
               tmp_f(vcf@gt[,tm]))
    ## only keep tumor or normal have mutations
    kp=snpm[,3] %in% c('1/1','0/1') | snpm[,6] %in% c('1/1','0/1')  
    table(kp)
    snpm=snpm[kp,]
    ## remove those show ./. positions
    kp=snpm[,3] == './.' | snpm[,6]== './.'  
    print(table(!kp))
    snpm=snpm[!kp,c(1,2,4,5,7,8)]
    rcmat=snpm
    rcmat$POS=as.numeric(rcmat$POS)
    dim(rcmat)
    rcmat=na.omit(rcmat)
    colnames(rcmat)=c("Chromosome", "Position","NOR.DP","NOR.RD","TUM.DP","TUM.RD")
    rcmat[,1]=gsub('chr','',rcmat$Chrom)
    ## fit segmentation tree
    xx = preProcSample(rcmat)
    ## estimate allele specific copy numbers
    oo=procSample(xx,cval=150)
    oo$dipLogR
    ## EM fit version 1
    fit=emcncf(oo)
    tmp=fit$cncf
    write.table(tmp,file = paste0('facets_cnv_',tm,'.seg.txt'),
                row.names = F,col.names = F,quote = F)
    head(fit$cncf)
    fit$purity
    fit$ploidy
    #png(paste0('facets_cnv_',tm,'.png'),res=150,width = 880,height = 880)
    pdf(paste0('facets_cnv_',tm,'.pdf'),width = 8,height = 6)
    plotSample(x=oo,emfit=fit, sname=tm)
    dev.off()
    if(F){
      fit2=emcncf2(oo)
      head(fit2$cncf)
      fit2$purity
      fit2$ploidy
      #png(paste0('facets_cnv2_',tm,'.png'),res=150,width = 880,height = 880)
      pdf(paste0('facets_cnv2_',tm,'.pdf'),width = 8,height = 6)
      plotSample(x=oo,emfit=fit2, sname=tm)
      dev.off()
      
    }
    return(c(fit$dipLogR,fit$ploidy,fit$purity))
  }
  
})
tmp <- do.call(rbind,tmp)
rownames(tmp)=tms
colnames(tmp)=c('dipLogR', 'ploidy', 'purity')
write.csv(tmp,'ploidy_and_purity.csv')
