# ------------- Figure S3 --------------
#----Figure S3A----
#infercnv，Identify CNV levels at the gene level in epithelial cells
#Extract epithelial cells corresponding to module patients
rm(list = ls())
setwd("./t100553/wss")
sce <- readRDS("./sce_Annotated_umap.rds")
corrsce <- readRDS("./subsce.corr.rds")
epi <- subset(sce,celltype=="Epithelial")
epi$CELLTYPE <- "Epithelial"
count.epi <- as.data.frame(epi@assays$RNA$counts)
#Extract control cells from the module
DC <- subset(corrsce,CELLTYPE %in% c("DC"))
count.DC <- as.data.frame(DC@assays$RNA$counts)
NK <- subset(corrsce,CELLTYPE %in% c("NK"))
count.NK <- as.data.frame(NK@assays$RNA$counts)
#Merge expression matrix
dat <- cbind(count.epi,count.DC,count.NK)
colnames(dat) <- unique(colnames(dat))
groupinfo=data.frame(v1=colnames(dat),
                     v2=c(rep('epi',ncol(epi)),
                          rep('DC',ncol(DC)),
                          rep('NK',ncol(NK))))
library(AnnoProbe)
library(stringr)
library(infercnv)
library(gtools)

geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(start)),c(1,4:6)]
geneInfor=geneInfor[with(geneInfor, mixedorder(chr,decreasing = F)),]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)


dir.create("./Infercnv.DC.NK")
setwd("./Infercnv.DC.NK")
expFile='expFile.txt'
write.table(dat,file = expFile,sep = '\t',quote = F)
groupFiles='groupFiles.txt'
head(groupinfo)

write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)
rm(DC,corrsce,count.DC,count.epi,count.NK,sce,NK,epi)
infercnv_obj_pre = CreateInfercnvObject(raw_counts_matrix=expFile,
                                        annotations_file=groupFiles,
                                        delim="\t",
                                        gene_order_file= geneFile,
                                        min_max_counts_per_cell = c(100, +Inf),
                                        ref_group_names=c("DC","NK"))  ## This depends on the information in your group
infercnv_obj2 = infercnv::run(infercnv_obj_pre,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir= "./",  # dir is auto-created for storing outputs
                              denoise=T, #denoising
                              cluster_by_groups=F,   
                              analysis_mode = "subclusters", 
                              leiden_resolution = 0.0001,
                              hclust_method="ward.D2",
                              num_threads = 15,
                              HMM = F,
                              output_format="pdf",
                              plot_steps=F)

#----Figure S3B----
#Extract clustering results
setwd("./wss/Infercnv.DC.NK")
rm(list=ls())
infercnv_obj <- readRDS("./run.final.infercnv_obj")
#Identify malignant epithelial cells
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")

expr <- infercnv_obj@expr.data
DC_loc <- infercnv_obj@reference_grouped_cell_indices$DC
NK_loc <- infercnv_obj@reference_grouped_cell_indices$NK
epi_loc <- infercnv_obj@observation_grouped_cell_indices$epi

anno.df=data.frame(
  CB=c(colnames(expr)[DC_loc],colnames(expr)[NK_loc],colnames(expr)[epi_loc]),
  class=c(rep("DC",length(DC_loc)),rep("NK",length(NK_loc)),rep("epi",length(epi_loc))))


head(anno.df)

gn <- rownames(expr)
geneFile <- read.table("./geneFile.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]

#Cluster, 2 categories, extract results
set.seed(10000000)
kmeans.result <- kmeans(t(expr),5)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") 
kmeans_df_s=arrange(kmeans_df,kmeans_class) 
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) 
head(kmeans_df_s)

#Define annotations and color schemes for heat maps
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:7] #number of categories
names(color_v)=as.character(1:7)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("DC"="#dd738c","NK" = "#9cc53d","epi"="#6796c0"),kmeans_class=color_v))

#Below is the drawing
pdf("module2.epi.heatmap.pdf",width = 15,height = 10)

ht = Heatmap(t(expr)[rownames(kmeans_df_s),], #Maintain consistency between the CB order of drawing data and annotation CB order
             col = colorRamp2(c(0.4,1,1.6), c("#377EB8","#F0F0F0","#E41A1C")), 
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), 
             column_gap = unit(2, "mm"),
             
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
             
             top_annotation = top_anno,left_annotation = left_anno, #Add comments
             row_title = NULL,column_title = NULL,use_raster=T)
draw(ht, heatmap_legend_side = "right")
dev.off()

#The corresponding CB for each category is saved in the kmeans_df_2 data box
write.table(kmeans_df_s, file = "kmeans_df_s.txt", quote = FALSE, sep = '\t', row.names = T, col.names = T)

#Draw a violin plot to further examine the CNV levels of each cluster
expr1=expr-1
expr2=expr1 ^ 2
CNV_score=as.data.frame(colMeans(expr2))
colnames(CNV_score)="CNV_score"
CNV_score$name <- rownames(CNV_score)
cnv <- CNV_score[match(kmeans_df$CB,rownames(CNV_score)),]
cnv$group <- kmeans_df$kmeans_class
cnv$group <- factor(cnv$group,levels = c("1","2","3","4","5"))
colnames(cnv) <- c("score","name","group")
library(scales)

color_v=RColorBrewer::brewer.pal(8, "Set3")[1:8]
pdf("kmeans.pdf",width = 10,height = 8)

cnv <- cnv[rownames(cnv) != "P44-GSE202371_P51_C32f", ]
cnv$score <- rescale(cnv$score,to=c(0, 1))
cnv%>%ggplot(aes(group,score))+geom_violin(aes(fill=group),color="NA")+
  scale_fill_manual(values = color_v)+theme_bw()+scale_y_continuous(limits = c(0, 1))
dev.off()

#----Figure S3C----
#The visualization of consistency clustering indicators is generated by the code of GE identification process in Figure 3A

#----Figure S3D----
#The consistency clustering results are generated by the code of the GE recognition process in Figure 3A

#----Figure S3E----
#Scoring each driver subtype using the GE gene set
rm(list=ls())
corsce <- readRDS("./t100553/wss/subsce.corr.rds")
setwd("./t100553/wss/copykat/cancer state/")
GElist <- read.csv("cancer_GE.csv")

GElist <- split(GElist$gene,GElist$GE)
length(GElist) # get number of GEs
#score
corsce <- AddModuleScore_UCell(corsce,
                               features = GElist,
                               assay = "RNA")

#Draw a box line diagram
library(tidyverse)
library(rstatix)    
library(ggpubr)   
library(ggplot2)

# Ensure that cellmodule is of factor type and specify the order
corsce@meta.data$mutation <- factor(corsce@meta.data$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
#GE1

df <- na.omit(corsce@meta.data[, c("GE1_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")

# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE1_UCell, 0.25),
    Q3 = quantile(GE1_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE1_UCell < (Q1 - 1.5 * IQR) | GE1_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = mutation, y = GE1_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE1 Activity by Cell Module", x = NULL, y = "GE1_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size

#GE2

df <- na.omit(corsce@meta.data[, c("GE2_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")

# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE2_UCell, 0.25),
    Q3 = quantile(GE2_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE2_UCell < (Q1 - 1.5 * IQR) | GE2_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

ggplot(df, aes(x = mutation, y = GE2_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE2 Activity by Cell Module", x = NULL, y = "GE2_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size


#GE3
 
df <- na.omit(corsce@meta.data[, c("GE3_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")

# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE3_UCell, 0.25),
    Q3 = quantile(GE3_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE3_UCell < (Q1 - 1.5 * IQR) | GE3_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = mutation, y = GE3_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE3 Activity by Cell Module", x = NULL, y = "GE3_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size

#GE4
 
df <- na.omit(corsce@meta.data[, c("GE4_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")

# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE4_UCell, 0.25),
    Q3 = quantile(GE4_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE4_UCell < (Q1 - 1.5 * IQR) | GE4_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = mutation, y = GE4_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE4 Activity by Cell Module", x = NULL, y = "GE4_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size

#GE5
 
df <- na.omit(corsce@meta.data[, c("GE5_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")
 
# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE5_UCell, 0.25),
    Q3 = quantile(GE5_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE5_UCell < (Q1 - 1.5 * IQR) | GE5_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = mutation, y = GE5_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE5 Activity by Cell Module", x = NULL, y = "GE5_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size

#GE6
 
df <- na.omit(corsce@meta.data[, c("GE6_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")

# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE6_UCell, 0.25),
    Q3 = quantile(GE6_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE6_UCell < (Q1 - 1.5 * IQR) | GE6_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = mutation, y = GE6_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE6 Activity by Cell Module", x = NULL, y = "GE6_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size


 
#GE7
 
df <- na.omit(corsce@meta.data[, c("GE7_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")

# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE7_UCell, 0.25),
    Q3 = quantile(GE7_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE7_UCell < (Q1 - 1.5 * IQR) | GE7_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = mutation, y = GE7_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE7 Activity by Cell Module", x = NULL, y = "GE7_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size

 
#GE8
 
df <- na.omit(corsce@meta.data[, c("GE8_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")

# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE8_UCell, 0.25),
    Q3 = quantile(GE8_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE8_UCell < (Q1 - 1.5 * IQR) | GE8_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = mutation, y = GE8_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE8 Activity by Cell Module", x = NULL, y = "GE8_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size


 
#GE9
 
df <- na.omit(corsce@meta.data[, c("GE9_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")

# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE9_UCell, 0.25),
    Q3 = quantile(GE9_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE9_UCell < (Q1 - 1.5 * IQR) | GE9_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = mutation, y = GE9_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE9 Activity by Cell Module", x = NULL, y = "GE9_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size




 
#GE10
 
df <- na.omit(corsce@meta.data[, c("GE10_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")

# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE10_UCell, 0.25),
    Q3 = quantile(GE10_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE10_UCell < (Q1 - 1.5 * IQR) | GE10_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

ggplot(df, aes(x = mutation, y = GE10_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE10 Activity by Cell Module", x = NULL, y = "GE10_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size

 
#GE11
 
df <- na.omit(corsce@meta.data[, c("GE11_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")
 
# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE11_UCell, 0.25),
    Q3 = quantile(GE11_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE11_UCell < (Q1 - 1.5 * IQR) | GE11_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = mutation, y = GE11_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE11 Activity by Cell Module", x = NULL, y = "GE11_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size


#GE12
  
df <- na.omit(corsce@meta.data[, c("GE12_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")
# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE12_UCell, 0.25),
    Q3 = quantile(GE12_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE12_UCell < (Q1 - 1.5 * IQR) | GE12_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = mutation, y = GE12_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE12 Activity by Cell Module", x = NULL, y = "GE12_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size

#GE13
  
df <- na.omit(corsce@meta.data[, c("GE13_UCell", "mutation")])
df$mutation <- factor(df$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))
 
mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")
# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(mutation) %>%
  mutate(
    Q1 = quantile(GE13_UCell, 0.25),
    Q3 = quantile(GE13_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE13_UCell < (Q1 - 1.5 * IQR) | GE13_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = mutation, y = GE13_UCell, color = mutation)) +
  geom_boxplot(
    fill = NA,  # No filling color
    outlier.shape = NA,  # Remove default outliers
    width = 0.6, size = 1.2, show.legend = FALSE
  ) +
  geom_point(
    data = df_outliers_sampled,
    position = position_jitter(width = 0.2),
    size = 1.2, alpha = 0.6, show.legend = FALSE
  ) +
  scale_color_manual(values = mycolors) +
  
  labs(title = "GE13 Activity by Cell Module", x = NULL, y = "GE13_UCell Score") +
  theme_classic(base_size = 14) +
  theme(
    
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    axis.line = element_line(color = "black", size = 0.6)
  )+ stat_compare_means(method = "anova", #statistical method
                        aes(label = "p.format"), size = 5) #size






