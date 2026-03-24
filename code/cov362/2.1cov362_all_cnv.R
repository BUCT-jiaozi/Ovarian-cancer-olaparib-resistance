
##Predict the CNV of potential cancer cells based on the annotated normal cells
setwd("/home/Rstudio/yhb/Data/cov362_sc")
library(Seurat)
data.merge <- readRDS(file = "./1.5data.merge.rds")


if(F){
  index <- which(data.merge@meta.data$orig.ident %in% c("C")) #
  Normal.cell <- rownames(data.merge@meta.data)[index] #
  Normal.cellType <- as.character(data.merge@meta.data$orig.ident[index])
}


##Need to guess the cell population of the copy number
dir.create('3.CNV')
setwd("./3.CNV")


#DimPlot(object = data.merge, reduction = 'umap',label = TRUE, group.by = "cellType_low")
#possible.cancer.cell <- rownames(data.merge@meta.data)[which(data.merge@meta.data$seurat_clusters %in% c("4"))] # 3564 cells
# possible.cancer.cell <- rownames(data.merge@meta.data)[which(data.merge@meta.data$cellType_low %in% c("c0", "c1", "c2", "c3"))] #

if(F){
  possible.cancer.cell <- rownames(data.merge@meta.data)[data.merge@meta.data$cellType_low %in% c("c0", "c1", "c2", "c3") &
                                                           data.merge@meta.data$orig.ident %in% c("T5", "T10", "T20", "T40")] #

}

DefaultAssay(data.merge) <- "RNA"
cancer.count <- GetAssayData(data.merge, slot = "counts")
counts_matrix <- cancer.count[, c(Normal.cell, possible.cancer.cell)] #12373 17244
# cell.id <- gsub("-", "\\.", colnames(counts_matrix)) # 将-转化为.
# colnames(counts_matrix) <- cell.id
# write.table(counts_matrix, file = "countmatrix.txt", sep = "\t", row.names = T, col.names = T, quote = F)
cell.id <- colnames(counts_matrix)

##Building cell type annotation file
index <- match(possible.cancer.cell, rownames(data.merge@meta.data))
cellLabel1 <- as.character(data.merge$cellType_low[index])
# cellLabel1 <- paste0("Sample_", cellLabel1)
annoFile <- data.frame(cell.id = cell.id, celltype = c(Normal.cellType, cellLabel1))
# write.table(annoFile, file = "annoFile.txt", col.names = F, row.names = F, quote = F, sep = "\t")
row.names(annoFile) <- annoFile$cell.id
annoFile <- annoFile[, -1, drop = FALSE]

write.table(annoFile, file = "anno_all.txt",col.names = FALSE, row.names = T, quote = FALSE, sep = "\t")

###231 Server
#setwd("./3.CNV")
library(infercnv) #Version:1.6
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= counts_matrix,
                                    annotations_file= annoFile,
                                    delim="\t",
                                    gene_order_file="../../OC_data/kura_new_fq/code/sup_files/hg38_gencode_v27.txt",
                                    ref_group_names = "C") # c("CD8+ T cell", "NK/NKT cell")
# infercnv_obj = infercnv::run(infercnv_obj,
#                              cutoff=0.1,# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
#                              out_dir= "./",
#                              cluster_by_groups=T, # 按照样本分组
#                              denoise=TRUE,
#                              HMM=TRUE, # 基于HNN预测cnv，速度慢
#                              #analysis_mode = "cells", # samples|subclusters|cells 图像过滤或HMM预测的分组级别。默认值：samples
#                              num_threads = 6)
# saveRDS(infercnv_obj,"infercnv_obj.rds")
# dir.create('./infercnv')
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir= "./",
                             cluster_by_groups=T,
                             denoise=F,
                             HMM=T,
                             num_threads = 6,
                             write_expr_matrix = T, # 生成绘图时是否写入包含矩阵内容的文本文件
                             cluster_references = TRUE # 在注释中对引用进行聚类，不显示树状图，默认T
                             # k_obs_groups = 8,
)
saveRDS(infercnv_obj,"infercnv_obj.rds")
# final_object <- readRDS("./run.final.infercnv_obj")




### cnv热图重新画
if(T){
  # annotation state for all
  cell_ann_all = read.table("anno_all.txt")
  
  cell_ann_list = list(cell_ann_all)
  # remove C from annotation
  cell_ann_list = lapply(cell_ann_list, function(x) x = x[x$V2 != "C", ])
  # get gene order
  gene_ordering = read.table("../../OC_data/kura_new_fq/code/sup_files/hg38_gencode_v27.txt")
  # read matrices
  infercnv_all = readRDS("./run.final.infercnv_obj")
  
  
  infercnv_obj_list = list(infercnv_all@expr.data)
  # shared genes in all
  all_genes = lapply(infercnv_obj_list, rownames)
  genes_shared = Reduce(intersect, all_genes)
  
  # subset infercnv with only treated cells and common genes in all
  for(i in 1:length(infercnv_obj_list)){
    cells = as.character(cell_ann_list[[i]]$V1) # 去除了 C组 细胞
    infercnv_obj_list[[i]] <- infercnv_obj_list[[i]][genes_shared, cells] # 保留交集基因
  }

  
  merge_exp <- function(exp_matrix_list) {

    
    merged_data = Reduce(function(x, y) transform(merge(x, y, by=0, all=TRUE), 
                                                  row.names=Row.names, Row.names=NULL),
                         exp_matrix_list)
    # Remove eventual NAs and replace by 0
    merged_data[is.na(merged_data)] <- 0
    return(merged_data)
  }
  
  # sample X % of cells in each group
  sample_cells <- function(exp_cnv_matrix, fraction = 0.3) {
    # sample cells to plot cnvs
    n_cells = round(length(colnames(exp_cnv_matrix)) * fraction)  # 取总细胞数的 30%
    sample_cells = sample(colnames(exp_cnv_matrix), n_cells)      # 随机抽样
    return(exp_cnv_matrix[, sample_cells])                        # 返回抽样后的矩阵
  }
  

  mapping <- c(
    "C" = "C",
    "T5"   = "T5",
    "T10" = "T10",
    "T20"   = "T20",
    "T40"  = "T40"
  )
  
  prepare_data_plot <- function(infercnv_obj_list_sample_ind, cell_ann_list_ind,
                                genes_shared, color_cluster_ind, color_sample_ind) {
    
    # 数据抽样
    # prepare data for individual samples to plot, the thing ready to plot?
    infercnv_obj_list_sample = lapply(list(infercnv_obj_list_sample_ind),
                                      sample_cells, fraction = 0.06)  
    # infercnv_obj_list_sample = lapply(list(infercnv_obj_list[[1]]),
    #                                   sample_cells, fraction = 0.9) # 示例运行
    
    # merge list of dataframes into one single dataframe
    infercnv_obj_sample_merged = merge_exp(infercnv_obj_list_sample)
    infercnv_obj_sample_merged = infercnv_obj_sample_merged[genes_shared, ]
    
    cell_ann_all = cell_ann_list_ind
    #cell_ann_all = cell_ann_list[[1]] # 示例运行
    rownames(cell_ann_all) <- as.character(cell_ann_all$V1)
    print(length(cell_ann_all))
    cell_ann_states = cell_ann_all[colnames(infercnv_obj_sample_merged), ]
    print(dim(infercnv_obj_sample_merged))
    
    # 细胞名包含分组信息
    #sample_state_ann = stringr::str_split(as.character(cell_ann_states$V2), pattern = "_")
    sample_state_ann = stringr::str_split(as.character(cell_ann_states$V1), pattern = "_")
    # extracting samples
    sample_ann = factor(sapply(sample_state_ann, "[[", 1)) # 提取样本名
    
    sample_ann <- factor(mapping[as.character(sample_ann)], levels = unique(mapping))
    
    # extracting clusters
    #state_ann = factor(sapply(sample_state_ann, "[[", 2)) # 提取 cluster/状态
    state_ann = factor(as.character(cell_ann_states$V2)) # 提取 cluster/状态
    
    # chr annotation 染色体注释
    rownames(gene_ordering) <- gene_ordering$V1
    gene_order_shared = gene_ordering[genes_shared, ]
    chr_ann = factor(gene_order_shared$V2, levels = c(paste("chr",1:22, sep = "")))
    levels(chr_ann) <- sub("^chr", "", levels(chr_ann)) # 去掉chr前缀
    
    # 细胞名称定义，细胞名_cluster
    # add a number in the end to make unique names
    rand_number = c(1:length(colnames(infercnv_obj_sample_merged)))
    # create new cell names that will be sorted after dendrogram! this helps to keep the
    # sample and state order according to the time points.
    new_col_names_for_sort = paste(as.character(sample_ann), as.character(state_ann), 
                                   1:length(colnames(infercnv_obj_sample_merged)), sep = "_")
    # add the new cell names
    colnames(infercnv_obj_sample_merged) <- new_col_names_for_sort
    
    # 直接调整细胞类型顺序
    # 获取列名并提取第二部分（state_ann）
    state_info <- sapply(strsplit(colnames(infercnv_obj_sample_merged), "_"), `[`, 2)
    sample_info <- sapply(strsplit(colnames(infercnv_obj_sample_merged), "_"), `[`, 1)
    
    # 按照 state_ann 的顺序重新排序列
    #infercnv_obj_sample_merged_new <- infercnv_obj_sample_merged[, order(factor(state_info, levels = levels(state_ann)))]
    # 按照两个因素排序：先按state，再按sample
    order_idx <- order(
      factor(state_info, levels = levels(state_ann)),
      factor(sample_info, levels = levels(sample_ann))
    )
    infercnv_obj_sample_merged_new <- infercnv_obj_sample_merged[, order_idx]
    
    state_ann_new <- sapply(strsplit(colnames(infercnv_obj_sample_merged_new), "_"), `[`, 2)
    sample_ann_new <- sapply(strsplit(colnames(infercnv_obj_sample_merged_new), "_"), `[`, 1)
    state_ann_new <- factor(as.character(state_ann_new))
    sample_ann_new <- factor(as.character(sample_ann_new), levels = unique(mapping))
    
    
    # 样本和分群颜色注释
    # sample state annotation
    sample_state_ann_heat = rowAnnotation(df = data.frame(sample = sample_ann_new, state = state_ann_new), 
                                          simple_anno_size = unit(2, "mm"),
                                          col = list(state = color_cluster_ind,
                                                     sample = color_sample_ind),
                                          show_annotation_name = F,
                                          annotation_legend_param = list(
                                            state  = list(title = "Cluster",title_gp = gpar(fontsize = 10),labels_gp = gpar(fontsize = 10)),
                                            sample = list(title = "Sample",title_gp = gpar(fontsize = 10),labels_gp = gpar(fontsize = 10))))
    
    # plot heatmap
    col_fun = colorRamp2(c(0.8, 1, 1.2), c("#000A87", "#ffffff", "#900512")) # this works well
    # max(infercnv_obj_list[[1]])
    # min(infercnv_obj_list[[1]])
    
    # 聚类信息
    library(dendextend)
    #dend = as.dendrogram(hclust(dist(t(infercnv_obj_sample_merged))), type = "average") # can be used as cluster_rows and columns arg instead of T
    # this does the nice trick to sort cells based on their names that I redefined by sample and state!
    # this orders strings based on numbers! and keeps the clustering constraints 数字大小排序
    #dend <- dendextend::rotate(dend, stringr::str_sort(labels(dend), numeric = T))
    
    
    heatmap_cnv = Heatmap(t(infercnv_obj_sample_merged_new), cluster_rows = F, row_dend_reorder = F,
                          cluster_columns = F,
                          show_row_names = F, show_column_names = F, show_row_dend = F,
                          use_raster = T, left_annotation = sample_state_ann_heat,
                          col = col_fun, column_split = chr_ann, column_gap = unit(1, "mm"),row_split = state_ann_new,row_gap = unit(0.5, "mm"),
                          row_title_rot = 0, row_title_gp = gpar(fontsize = 12),
                          border = TRUE, column_title_side = "bottom", 
                          column_title_gp = gpar(fontsize = 12),
                          row_names_gp = gpar(fontsize = 10),
                          column_names_gp = gpar(fontsize = 10),
                          heatmap_legend_param = list(title_gp = gpar(fontsize = 10),
                                                      title = "Expression",
                                                      labels_gp = gpar(fontsize = 10),
                                                      legend_height = unit(2, "cm")))
    return(heatmap_cnv)
  }
  # new_colors = c("#8DA3A6","#B997C7","#824D99","#4E78C4","#57A2AC","#7EB875",
  #                "#D0B541", "#E67F33", "#CE2220", "#80271E")
  
  cluster_colors = list(c("c0" = "#8DA3A6", "c1" = "#57A2AC", "c2" = "#7EB875", "c3" = "#4E78C4")) #, "R5" = "#421510"
  
  
  sample_colors = list(c("T5" = "#B997C7", "T10" = "#824D99", "T20" = "#D0B541", "T40" = "#E67F33"))

  
  #infercnv_obj_list_sample = lapply(infercnv_obj_list, sample_cells, fraction = 0.3)
  # infercnv_obj_list_sample = infercnv_obj_list
  
  
  list_heatmaps = list()
  
  for(i in 1:length(infercnv_obj_list)) {
    heat_cnv = prepare_data_plot(infercnv_obj_list[[i]], 
                                 cell_ann_list[[i]],
                                 genes_shared, 
                                 cluster_colors[[i]],
                                 sample_colors[[i]]
    )
    list_heatmaps[[i]] <- heat_cnv
  }
  
  
  pdf("./1.0all_cnv_heatmaps.pdf", width = 10, height = 3.5)
  draw(
    list_heatmaps[[1]],
    padding = unit(c(1, 1, 6, 1), "mm")  # 下、左、上、右
  )
  
  grid::grid.text(
    "COV362",
    y = unit(1, "npc") - unit(3, "mm"),
    gp = gpar(fontsize = 14, fontface = "bold")
  )
  dev.off()
}




