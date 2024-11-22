Cell.integrated <- readRDS("F:/R/SLE lncRNA GSE181500/SystemicLupusErythematosus/RDS Files/Intergrated/4. Seurat.integrated.annotated.rds")

table(Cell.integrated@active.ident)

cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e","#4aef7b",
                "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", "#99db27", "#e07233", "#ff523f",
                "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6",
                "#d66551","#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,
                "#22547f", "#db5e92","#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,
                "#7b34c1" ,"#0cf29a","#d80fc1","#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5",
                "#925bea", "#63ff4f")

subset0 <- subcluster_view(seurat_obj = Cell.integrated, cluster_id = "B cells", resolution = 0.6, dims = 1:10)
table(subset0@active.ident)
DimPlot(Cell.integrated, reduction = "umap", label = T)
DimPlot(subset0, reduction = "umap", label = T)

Cell.integrated <- merge_subclusters(
  seurat_obj = Cell.integrated,
  subcluster_obj = subset0,
  new_labels = c("B0", "B0", "B0", "B0", "B1", "B1","B1"),
  new_id_column = "refined_clusters"
)


table(Cell.integrated@active.ident)

DimPlot(Cell.integrated, reduction = "umap", label = T)

saveRDS(Cell.integrated, "F:/R/SLE lncRNA GSE181500/SystemicLupusErythematosus/RDS Files/Intergrated/4. Seurat.integrated.TEST annotated.rds")

#########################################################################################################################

Cell.integrated.Annotated<- readRDS("F:/R/SLE lncRNA GSE181500/SystemicLupusErythematosus/RDS Files/Intergrated/4. Seurat.integrated.TEST annotated.rds")

conserved.gene <- scCEGs(Cell.integrated, grouping.var = "condition")

View(conserved.gene)

############################################################################################
####################################   Draw stacked bar plot  ####################################
####Calculate cell numbers for each cluster, based on clusters and groups/conditions
## extract meta data
md <- Cell.integrated.Annotated@meta.data %>% as.data.table
#将数据框或其他类型的数据对象转换为 data.table 对象。data.table 是R中一个高效的数据存储和操作框架。

# count the number of cells per unique combinations of "Sample" and "seurat_clusters"
Cell.integrated.Annotated_cluster.cell.numbers_1 <- md[, .N, by = c("condition", "refined_clusters")]
#理解不同条件下不同聚类中细胞的分布情况非常有用
view(Cell.integrated.Annotated_cluster.cell.numbers_1)

write.csv(Cell.integrated.Annotated_cluster.cell.numbers_1,
          "F:/R/Osteoarthritis/4. CSV files/3.1 Cell.integrated.Annotated_cluster.cell.numbers_1.csv", quote = F)


# with additional casting after the counting
Cell.integrated.Annotated_cluster.cell.numbers_2 <- md[, .N, by = c("condition", "refined_clusters")] %>%
  dcast(., condition ~ refined_clusters, value.var = "N")
#使用 dcast() 函数将数据从长格式转换为宽格式

view(Cell.integrated.Annotated_cluster.cell.numbers_2)

write.csv(Cell.integrated.Annotated_cluster.cell.numbers_2,
          "F:/R/Osteoarthritis/4. CSV files/3.1 Cell.integrated.Annotated_cluster.cell.numbers_2.csv", quote = F)
#########################################################################################################

#####stacked barplot
#Healthypare .csv file, containing 3 columns: Group (condition/Condition), sub-group(Cluster names/PodocyteSubset), Number (Cell numbers)

data <- Cell.integrated.Annotated_cluster.cell.numbers_1
view(data)

#draw a stacked barplot, vertical
P_stackedbarplot <- ggplot(data, aes(fill = refined_clusters, y = N , x = condition)) +
  geom_bar(position="fill", stat="identity")
#创建一个堆叠条形图（stacked bar plot），用于展示不同条件下各个聚类的细胞计数或比例。
#data：传递给ggplot()函数的数据集，包含至少三列：seurat_clusters（聚类标识）、N（计数或比例数值）、condition（条件标识）。
#geom_bar()：这是ggplot2中的一个几何对象（geom），用于绘制条形图。
#position = "fill"：指定条形的堆叠方式，这里是按填充颜色堆叠，即每个 condition 下的 seurat_clusters 将堆叠在一起。
#stat = "identity"：指定条形图的统计方法，"identity" 使用数据中直接给定的y值，而不是计算它们（如求和或平均）。
P_stackedbarplot

P_stackedbarplot + theme(panel.background = element_rect(fill = NA)) +
  theme(panel.grid.major = element_line(colour = NA))



#draw a stacked barplot, lateral
P_stackedbarplot_2 <- ggplot(data, aes(fill = refined_clusters, y = condition, x= N)) +
  geom_bar(position="fill", stat="identity")
#y = N：指定y轴的数值基于 N 列，这代表每个条形的高度或计数。
#x = condition：指定x轴的类别基于 condition 列，不同的条件将形成不同的条形。
P_stackedbarplot_2 + theme(panel.grid.major = element_line(colour = "lemonchiffon1"))

P_stackedbarplot_2 + theme(panel.background = element_rect(fill = NA)) +
  theme(panel.grid.major = element_line(colour = NA))+labs(title = "Frequencies",
                                                           x = "Frequency", y = NULL, fill = "Seurat_clusters",
                                                           size = 15)

############################################################################################
############################################################################################


