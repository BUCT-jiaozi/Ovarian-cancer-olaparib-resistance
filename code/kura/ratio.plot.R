ratio.plot <- function(seurat.object, id.vars1 = "orig.ident", id.vars2 = "seurat_clusters", angle = 45, 
                       color.len = c("#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#E5B350", "#F47D2B", "#83AD00", "#00C0BA", "#7B6FD0", "#C06CAB")){
  
  color.len <- color.len[1:length(unique(seurat.object@meta.data[,id.vars1]))]
  # a <- sort(table(seurat.object@meta.data[,id.vars2]), decreasing=T)
  a <- table(seurat.object@meta.data[,id.vars2])
  cluster.number <- data.frame(Cluster = names(a), number = as.numeric(a))
  # cluster.number$Cluster <- factor(cluster.number$Cluster, cluster.number$Cluster) # 横坐标按细胞数量大小排列
  
  patient.Cluster.ratio <- table(seurat.object@meta.data[,c(id.vars1, id.vars2)])
  patient.Cluster.ratio <- as.data.frame(patient.Cluster.ratio)
  cell.num <- tapply(patient.Cluster.ratio$Freq, patient.Cluster.ratio[, id.vars2], sum)
  cell.num <- rep(cell.num, each = length(unique(seurat.object@meta.data[,id.vars1]))) #each指定每个重复多少次
  patient.Cluster.ratio$ratio <- patient.Cluster.ratio$Freq/cell.num
  patient.Cluster.ratio[, id.vars2] <- factor(patient.Cluster.ratio[, id.vars2], levels = cluster.number$Cluster)
  #xlabs <- paste0(cluster.number$Cluster, " (n=", cluster.number$number, " cells)") # 加细胞数量
  xlabs <- cluster.number$Cluster # 加细胞数量
  library(ggplot2)
  colnames(patient.Cluster.ratio)[2] <- "Type2"
  colnames(patient.Cluster.ratio)[1] <- "Type1"
  p <- ggplot(data = patient.Cluster.ratio, aes(x = Type2, y = ratio, fill = Type1)) + 
    theme_bw()+
    geom_bar(stat= 'identity', position = 'fill',width = 0.85)+ #堆叠图，position = fill 表示堆叠图
    labs(x = '',y = 'Cell Ratio',fill = NULL)+ #定义坐标轴以及图例标题
    scale_fill_manual(values = color.len) +#自定义颜色，可通过`library(RColorBrewer);display.brewer.all()`来展示所有备选项
    scale_y_continuous(labels = c("0%","25%","50%", "75%", "100%")) +  ## 百分比坐标轴（需加载scales包）
    scale_x_discrete(labels = xlabs) +  
    theme(axis.text.x = element_text(angle = angle, hjust=0.3, size = 12), #x轴标签偏转45°，并下降0.5
          axis.text.y = element_text(size = 10),
          panel.grid = element_blank(),text = element_text(size = 12),
          legend.position = 'right',
          legend.key.height = unit(0.6,'cm'))  # 长宽比，以y / x表示
  print(p)
}
