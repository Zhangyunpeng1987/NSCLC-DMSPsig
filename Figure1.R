############################################################
## Project: Coordinated multicellular immune programs and drug targets revealed by single-cell analysis in driver-mutated NSCLC
##
## Purpose:
## Reproducible R workflow for single-cell analysis, TIME module identification, prognostic model construction, immune profiling, and drug sensitivity prediction.
##
## R version: 4.4.2
############################################################


# ------------- Path settings --------------
project_dir <- "./NSCLC-DMSPsig"
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
figure1_dir <- file.path(results_dir, "Figure1")
dir.create(figure1_dir, showWarnings = FALSE, recursive = TRUE)

# ------------- Figure 1 --------------
library(dplyr)
library(ggplot2)
library(Seurat)
#----Figure 1A---- 
sce <- readRDS(file.path(data_dir, "sce_annotated_umap.rds"))
# Display the number of patients for each driver mutation
# 1. Extract meta.data
df <- sce@meta.data %>%
  dplyr::select(patient, mutation)

# 2. Remove duplicates; each patient is counted only once within each mutation type
df_unique <- df %>%
  distinct(patient, mutation)

# 3. Count the number of patients for each mutation type
df_count <- df_unique %>%
  group_by(mutation) %>%
  summarise(n_patient = n()) %>%
  ungroup()

# 4. Optional: sort by the number of patients for better readability
df_count$mutation <- factor(
  df_count$mutation,
  levels = df_count$mutation[order(df_count$n_patient, decreasing = TRUE)]
)

# 5. Draw the bar plot
ggplot(df_count, aes(x = mutation, y = n_patient)) +
  geom_bar(stat = "identity", width = 0.7, fill = "lightblue") +
  theme_classic(base_size = 14) +
  labs(
    x = "Mutation type",
    y = "Number of patients"
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::label_number(accuracy = 1)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the number of cells for each driver mutation
df_cell_count <- df %>%
  group_by(mutation) %>%
  summarise(n_cells = n()) %>%
  ungroup()

# Draw the bar plot
ggplot(df_cell_count, aes(x = mutation, y = n_cells)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#6CB689") +
  theme_classic(base_size = 14) +
  labs(
    x = "Mutation type",
    y = "Number of cells"
  ) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::label_number(accuracy = 1)
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#----Figure 1B----
# Draw the circular UMAP plot using plot1cell
library(plot1cell)
library(Seurat)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(dplyr)
rm(list=ls())
setwd(figure1_dir)
sce <- readRDS(file.path(data_dir, "sce_annotated_umap.rds"))
circ_data <- prepare_circlize_data(sce, scale = 0.8 )
circ_data$celltype <- factor(circ_data$celltype,levels = c("T/NK","B","Myeloid","Epithelial","Fibroblast","Endothelial","Astrocyte"))
set.seed(1234)
#cluster_colors<-c("#aad09e","#91c7c2","#ebc1d8",'#d4c2db',"#f6e4e0","#99bcdd","#6bb3c0")#,"#4198b9","#B2DF8A",)
#cluster_colors <- c("#E6C5CE","#9589BA","#C63C55","#9FBEAB","#65729A","#F2D7A4","#F2B5A1")
#cluster_colors <- c("#EB4B59","#5FD490","#FAAC56","#9381BA","#EAA6C4","#3CBBB4","#239CC1")
cluster_colors <- c("#8fb4dc","#ffdd8e","#70cdbe","#ac99d2","#7ac3df","#eb7e60","#f5aa61")
cluster_colors <- c("#8fb4dc","#ffdd8e","#f5aa61","#ac99d2","#7ac3df","#eb7e60","#70cdbe")
#cluster_colors <- c("#ffdd8e","#f5aa61","#8fb4dc","#ac99d2","#7ac3df","#70cdbe","#eb7e60")
group_colors <- c("#8DD3C7",  "#BEBADA", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FB8072",  "#1F78B4", "#33A02C")
rep_colors <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
                '#c6adb0','#a5a9b0','#EAEAAE','#A8A8A8','#92699e',
                '#D19275','#df5734','#C8A8DA','#FFE5C0','#d4c2db',
                '#537eb7','#83ab8e','#ece399','#80B1D3','#EBC79E',
                '#b95055','#d5bb72','#FFB284','#e0cfda','#d8a0c0',
                '#e6b884','#FFCCBF','#fae3ae','#64a776','#cbdaa9',
                '#e9bebc','#a9c2cb','#D9D9F3','#63a3b8','#c8c7e1',
                '#d25774','#c49abc','#e6e2a3',"#9999FF",'#efd2c9',
                '#E9C2A6','#c4daec','#ebb1a4','#61bada','#b7deea',
                '#e29eaf','#4490c4',"#33CCCC","#66CCFF","#FF9966","#FF99CC")

#rand_color(length(names(table(sce$orig.ident))))
# Modify function parameters and increase the font size

plot_circlize(circ_data,do.label = T, pt.size = 0.1, col.use = cluster_colors ,bg.color = 'white', kde2d.n = 200, repel = T, label.cex = 0.8)
add_track(circ_data, group = "mutation", colors = group_colors, track_num = 2)
add_track(circ_data, group = "orig.ident",colors = rep_colors, track_num = 3)


#----Figure 1C----
# Draw UMAP plots and circular proportion plots for each driver subtype
table(sce$mutation)
EGFR <- subset(sce,mutation=="EGFR")
EGFR.BM <- subset(sce,mutation=="EGFR-BM") 
EGFR.co.mutation <- subset(sce,mutation=="EGFR-co-mutation")
KRAS <- subset(sce,mutation=="KRAS")
KRAS.co.mutation <- subset(sce,mutation=="KRAS-co-mutation")
ROS1 <- subset(sce,mutation=="ROS1")
TP53 <- subset(sce,mutation=="TP53")
ALK <- subset(sce,mutation=="ALK")
HER2 <- subset(sce,mutation=="HER2")
MET.BM <- subset(sce,mutation=="MET-BM")

#EGFR
EGFR <- ScaleData(EGFR)
EGFR <- RunPCA(EGFR, features = VariableFeatures(object = EGFR))
ElbowPlot(EGFR,ndims = 50)
# Calculate PCA
subcluster <- EGFR
calculate_pcs <- function(scobj) {
  # Extract standard deviations and calculate percentage contributions
  pct <- scobj[["pca"]]@stdev / sum(scobj[["pca"]]@stdev) * 100
  # Calculate cumulative contribution rates
  cumu <- cumsum(pct)
  # Identify the first PC where the cumulative contribution exceeds 90% and the individual PC contribution is below 5%
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Identify the first PC after consecutive PCs with a contribution change greater than 0.1%
  co2 <- if (length(pct) > 1) {
    sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  } else {
    # If only one PC is present, this criterion is not used by default
    NA
  }
  # Select the number of PCs; use the minimum of co1 and co2, or co1 if co2 is unavailable
  pcs <- ifelse(is.na(co2), co1, min(co1, co2, na.rm = TRUE))
  # Create a data frame containing the required information
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
  # Return the result list
  list(
    pcs = pcs,
    co1 = co1,
    co2 = co2,
    plot_df = plot_df
  )
}
pca <- calculate_pcs(subcluster)$pcs
pca

library(harmony)
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)
cluster_colors <- c("#8fb4dc","#ffdd8e","#f5aa61","#ac99d2","#7ac3df","#eb7e60","#70cdbe")
cluster_colors <- c("#eb7e60","#ffdd8e","#70cdbe","#ac99d2","#7ac3df","#f5aa61")
DimPlot(subcluster,label = T,reduction = "umap",group.by = "celltype",raster = F,cols = cluster_colors)
EGFR <- subcluster

#EGFR.BM
EGFR.BM <- ScaleData(EGFR.BM)
EGFR.BM <- RunPCA(EGFR.BM, features = VariableFeatures(object = EGFR.BM))
subcluster <- EGFR.BM
pca <- calculate_pcs(subcluster)$pcs
pca

library(harmony)
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)
table(subcluster$celltype)
cluster_colors <- c("#8fb4dc","#ffdd8e","#f5aa61","#ac99d2","#7ac3df","#eb7e60","#70cdbe")
cluster_colors <- c("#eb7e60","#ffdd8e","#70cdbe","#ac99d2","#7ac3df","#f5aa61","#8fb4dc")
DimPlot(subcluster,label = T,reduction = "umap",group.by = "celltype",raster = F,cols = cluster_colors)
EGFR.BM <- subcluster

#EGFR.co.mutation
EGFR.co.mutation <- ScaleData(EGFR.co.mutation)
EGFR.co.mutation <- RunPCA(EGFR.co.mutation, features = VariableFeatures(object = EGFR.co.mutation))
subcluster <- EGFR.co.mutation
pca <- calculate_pcs(subcluster)$pcs
pca
library(harmony)
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)

DimPlot(subcluster,label = T,reduction = "umap",group.by = "celltype",raster = F,cols = cluster_colors)
EGFR.co.mutation <- subcluster

#KRAS
KRAS <- ScaleData(KRAS)
KRAS <- RunPCA(KRAS, features = VariableFeatures(object = KRAS))
subcluster <- KRAS
pca <- calculate_pcs(subcluster)$pcs
pca
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)

DimPlot(subcluster,label = T,reduction = "umap",group.by = "celltype",raster = F,cols = cluster_colors)
KRAS <- subcluster

#KRAS-co-mutation
KRAS-co-mutation <- ScaleData(KRAS-co-mutation)
KRAS-co-mutation <- RunPCA(KRAS-co-mutation, features = VariableFeatures(object = KRAS-co-mutation))
subcluster <- KRAS-co-mutation
pca <- calculate_pcs(subcluster)$pcs
pca
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)

DimPlot(subcluster,label = T,reduction = "umap",group.by = "celltype",raster = F,cols = cluster_colors)
KRAS-co-mutation <- subcluster


#ROS1
ROS1 <- ScaleData(ROS1)
ROS1 <- RunPCA(ROS1, features = VariableFeatures(object = ROS1))
subcluster <- ROS1
pca <- calculate_pcs(subcluster)$pcs
pca
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)

DimPlot(subcluster,label = T,reduction = "umap",group.by = "celltype",raster = F,cols = cluster_colors)
ROS1 <- subcluster



#TP53
TP53 <- ScaleData(TP53)
TP53 <- RunPCA(TP53, features = VariableFeatures(object = TP53))
subcluster <- TP53
pca <- calculate_pcs(subcluster)$pcs
pca
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)

DimPlot(subcluster,label = T,reduction = "umap",group.by = "celltype",raster = F,cols = cluster_colors)
TP53 <- subcluster


#ALK
ALK <- ScaleData(ALK)
ALK <- RunPCA(ALK, features = VariableFeatures(object = ALK))
subcluster <- ALK
pca <- calculate_pcs(subcluster)$pcs
pca
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)

DimPlot(subcluster,label = T,reduction = "umap",group.by = "celltype",raster = F,cols = cluster_colors)
ALK <- subcluster

#HER2
HER2 <- ScaleData(HER2)
HER2 <- RunPCA(HER2, features = VariableFeatures(object = HER2))
subcluster <- HER2
pca <- calculate_pcs(subcluster)$pcs
pca
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)

DimPlot(subcluster,label = T,reduction = "umap",group.by = "celltype",raster = F,cols = cluster_colors)
HER2 <- subcluster


#MET-BM
MET-BM <- ScaleData(MET-BM)
MET-BM <- RunPCA(MET-BM, features = VariableFeatures(object = MET-BM))
subcluster <- MET-BM
pca <- calculate_pcs(subcluster)$pcs
pca
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)

DimPlot(subcluster,label = T,reduction = "umap",group.by = "celltype",raster = F,cols = cluster_colors)
MET-BM <- subcluster

## Draw a pie chart for mutation types
# Method 1
sce <- readRDS(file.path(data_dir, "sce_annotated_umap.rds"))
library(ggsci)
library(ggplot2)
library(ggforce)
data <- data.frame(table(sce@meta.data$mutation,sce@meta.data$orig.ident))
# Draw the number of cells for patients with each driver mutation
mutation <- c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53","MET-BM","HER2")
cellnumber <- c(40965,36041,5431,65791,13540,5990,15575,7954,3355,1918)
data <- data.frame(mutation,cellnumber)

mycol2 = c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")
data$mutation <- factor(data$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53","MET-BM","HER2"))
p1 <- 
  ggplot() +
  geom_arc_bar(data=data,stat = "pie",aes(x0=0,y0=0,r0=1,r=2,amount=cellnumber,fill=mutation))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+# Remove unnecessary ggplot background elements and axes  
  xlab("")+ylab('')+  
  scale_fill_manual(values = mycol2)  
p1

# Draw the cell proportion for each driver mutation
mutation <- c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53","MET-BM","HER2")
patientnumber <- c(11,10,2,12,3,1,3,1,1,1)
data2 <- data.frame(mutation,patientnumber)
data2$mutation <- factor(data$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53","MET-BM","HER2"))
p2 <- 
  ggplot() +
  geom_arc_bar(data=data2,stat = "pie",aes(x0=0,y0=0,r0=1,r=2,amount=patientnumber,fill=mutation))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+# Remove unnecessary ggplot background elements and axes  
  xlab("")+ylab('')+  
  scale_fill_manual(values = mycol2)
p2

#----Figure 1D----
# Draw cell distribution density plots
table(sce$mutation)
EGFR <- subset(sce,mutation=="EGFR")
EGFR.BM <- subset(sce,mutation=="EGFR-BM") 
EGFR.co.mutation <- subset(sce,mutation=="EGFR-co-mutation")
KRAS <- subset(sce,mutation=="KRAS")
KRAS.co.mutation <- subset(sce,mutation=="KRAS-co-mutation")
ROS1 <- subset(sce,mutation=="ROS1")
TP53 <- subset(sce,mutation=="TP53")
ALK <- subset(sce,mutation=="ALK")
HER2 <- subset(sce,mutation=="HER2")
MET.BM <- subset(sce,mutation=="MET-BM")

# Draw plots separately for EGFR, EGFR.BM, EGFR.co.mutation, KRAS, KRAS.co.mutation, ROS1, TP53, ALK, HER2, and MET.BM
umap_sce <- EGFR@reductions[["umap"]]@cell.embeddings
umap_sce <- as.data.frame(umap_sce)
xlims=c(min(sce@reductions[["umap"]]@cell.embeddings[,1])-1,max(sce@reductions[["umap"]]@cell.embeddings[,1])+1)
ylims=c(min(sce@reductions[["umap"]]@cell.embeddings[,2])-1,max(sce@reductions[["umap"]]@cell.embeddings[,2])+1)
library(viridis)
library(ggplot2)
RLI_plot=ggplot(umap_sce, aes(x=umap_1, y=umap_2) ) + xlim(xlims) + ylim(ylims)+
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)  +xlab('umap_1')+ylab('umap_2')+
  #scale_fill_distiller(palette=2, direction=0.1,expand = c(0, 0))+     # +ggtitle(iterm)
  #scale_fill_manual(values = mycol)+
  scale_fill_viridis(option='B',alpha = 1,direction=1) +  #ggtitle(iterm)+
  geom_point(aes(x=umap_1, y=umap_2), col='#FCFDBFFF',size=0.00001,alpha=0.5)+
  #scale_x_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    #legend.position='none',
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    #legend.position='none',
    plot.margin = margin(0,0,0,0,"cm")) 
RLI_plot


#----Figure 1E----
# Identify marker genes for each cell type and draw marker gene plots
Idents(sce) <- "celltype"
sce$seurat_clusters <- sce@active.ident

# Identify differentially expressed genes for each cluster and annotate cells
DefaultAssay(sce) <- "RNA"
all.markers  <- FindAllMarkers(sce, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
significant.markers  <- all.markers [all.markers $p_val_adj < 0.05, ]
celltype.deg_top5 <- significant.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(significant.markers, file = file.path(figure1_dir, "celltype.deg_5.significant.markers_p_0.05.csv"))# Save

# Display marker genes of each cell group using a dot plot
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(scales)
library(reshape2)

celltype.deg_top3 <- significant.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DotPlot(sce,
        features = split(celltype.deg_top3$gene, celltype.deg_top3$cluster),
        cols = c("white","#8385DA" )
) +
  RotatedAxis() + # From Seurat
  theme(
    panel.border = element_rect(color = "black"),
    panel.spacing = unit(1, "mm"),
    axis.title = element_blank(),
    axis.text.y = element_blank()
  )
#"#beb8dc","#77CBC6","#D47A9A"




#----Figure 1F----
# Stacked bar plot of cell proportions
library(ggpubr)
sce$celltype <- factor(sce$celltype,levels = c("Myeloid","B","T/NK","Epithelial","Fibroblast","Endothelial","Astrocyte"))
cluster_colors <- c("#8ECAC0","#B57EB2","#B1CB6C","#DC2177","#E4CD76","#E97974","#4067A4")
P1 <- 
  sce@meta.data %>%
  ggplot(aes(x = mutation, fill = celltype)) +
  theme_pubr(base_size = 5) +
  theme(plot.title = element_blank(),
        text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(3,"pt"),
        legend.position = "right") +
  #facet_grid(~mutation, scales = "free_x", space ="free_x") +
  #geom_bar( alpha = 0.9,colour="white")+
  geom_bar( alpha = 0.9,position = "fill")+
  scale_fill_manual(values = cluster_colors)+
  xlab("patient") + ylab("Fraction") +coord_flip()
P1

#----Figure 1G----
# Determine immune types of patients with different driver mutations
# Define major immune types based on T/B and Myeloid/NK cell proportions (innate/adaptive immunity)
library(dplyr)
TIME@meta.data <- TIME@meta.data %>% mutate(immunetypecell = case_when(
  CELLTYPE %in% c("CD84","CD8T","Treg","B") ~ "Adaptive",
  CELLTYPE %in% c("Macrophage","Mast","Neutrophil","DC") ~ "Innate",
  CELLTYPE %in% c("Endothelial","Fibroblast","Astrocyte","NK") ~"Stromal"
  
))
table(TIME$immunetypecell)
saveRDS(TIME,"./TIME.rds")
write.csv(table(TIME$immunetypecell,TIME$orig.ident),file = "./immune.type.csv")
immune.type <- read.csv("./immune.type.csv",row.names = 1)
immune.type <- t(immune.type)
ii <- immune.type[,-3]
immune.type <- prop.table(ii,1)
immune.type <- as.data.frame(immune.type)
immune.type$orig.ident <- rownames(immune.type)
immune.type <- immune.type %>% mutate(mutation = case_when(
  orig.ident %in% c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11") ~ "EGFR",
  orig.ident %in% c("P12","P13","P14","P15","P16","P17","P18","P19","P20","P21") ~ "EGFR-BM",
  orig.ident %in% c("P22","P23") ~ "EGFR-co-mutation",
  orig.ident %in% c("P24","P25","P26","P27","P28","P29","P30","P31","P32","P33","P34","P35") ~ "KRAS",
  orig.ident %in% c("P36","P37","P38") ~ "KRAS-co-mutation",
  orig.ident %in% c("P39") ~ "ALK",
  orig.ident %in% c("P40","P41","P42") ~ "ROS1",
  orig.ident %in% c("P43") ~ "TP53",
  orig.ident %in% c("P44") ~ "MET-BM",
  orig.ident %in% c("P45") ~ "HER2"))
immune.type$orig.ident <- NULL


immune.type1 <- immune.type
immune.type$Innate <- NULL
immune.type$patient <- rownames(immune.type)
immune.type$type <- "Adaptive"

immune.type1$Adaptive <- NULL
immune.type1$patient <- rownames(immune.type1)
immune.type1$type <- "Innate"
colnames(immune.type) <- c("percent","mutation","patient","type")
colnames(immune.type1) <- c("percent","mutation","patient","type")
DATA <- rbind(immune.type,immune.type1)


library(ggplot2)
library(ggsci)
library(ggthemes)
library(ggpubr)
library(ggprism)

# Plot
write.csv(DATA,file = "./immune.type.csv")
DATA <- read.csv(file="./immune.type.csv",row.names=1)
P <- ggplot(DATA, aes(x=factor(type), y=percent),color=mutation) +
  # Draw the boxplot
  geom_boxplot(aes(color=mutation),
               alpha=0.1,lwd=0.8)+ # Set transparency
  # Draw jittered points
  geom_jitter(aes(color = factor(mutation)), size = 1.5) +  # Map colors using aes()
  # Set colors
  scale_color_manual(values = pal_npg('nrc')(10))+
  scale_fill_manual(values = pal_npg('nrc')(10))+
  # Set theme
  theme_bw()+stat_boxplot(geom="errorbar",width=0.15,aes(color=mutation))+
  # Remove grid lines
  theme(panel.grid = element_blank())+facet_wrap(~mutation,scale="free")

P    
#t.test, wilcox.test, wilcox.test, wilcox.test, anova, anova, kruskal.test, kruskal.test
P+stat_compare_means(method = "t.test")+labs(title="immune type",
                                             x=NULL,y="Fraction") 

#----Figure 1H----
# Draw functional enrichment results for each cell type
# GSVA at the cell-subpopulation level
library(gplots)
library(ggplot2) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA) # BiocManager::install('GSVA')
library(GSEABase)
library(msigdbr)
library(pheatmap)
rm(list = ls())
sce <- readRDS(file.path(data_dir, "sce_annotated_umap.rds"))
TIME.sce <- readRDS(file.path(figure1_dir, "subsce.corr.rds"))
# Use hallmark gene sets from msigdbr
all_gene_sets <- msigdbr(species ="Homo sapiens",category = "H")
# Convert to the object required for GSVA analysis
gs <- split(all_gene_sets$gene_symbol,all_gene_sets$gs_name)
gs <- lapply(gs, unique)
gsc <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, gs, names(gs)))
gsc
geneset <- gsc
geneset # This is a GeneSetCollection object
library(Seurat)
# For all major cell types
all.celltype <- AggregateExpression(sce , group.by = "celltype", assays = "RNA") 
all.celltype <- all.celltype[[1]]
dim(all.celltype)
head(all.celltype) # The output is an integer matrix

X1 <- all.celltype 
## Convert X1 to as.matrix(X)
# Code for the older version of GSVA
es.max <- gsva(as.matrix(X1), geneset, mx.diff=FALSE, verbose=FALSE,  parallel.sz=8)
rownames(es.max) <- gsub("HALLMARK_","",rownames(es.max))
pheatmap(es.max[1:50, ],show_rownames = T,fontsize = 9,border = F,color = colorRampPalette(c("#63b7bd","#57bad7","#b8e8ee","white","#fdcec3","#ee654c","#d45241"))(50)) # Draw the plot

#----Figure 1I----
# Cell subclustering
rm(list=ls())
sce <- readRDS(file.path(data_dir, "sce_annotated_umap.rds"))
# Divide cells into lymphocytes, macrophages, endothelial cells, and stromal cells
# Lymphocytes (T/NK/B/Cycling)
#T/NK
T.NK.cell <- subset(sce,celltype == "T/NK")
#T.NK.cell <- subset(T.NK.cell,features=setdiff(rownames(T.NK.cell),delete.gene))
T.NK.cell <- ScaleData(T.NK.cell)
T.NK.cell <- RunPCA(T.NK.cell, features = VariableFeatures(object = T.NK.cell))
ElbowPlot(T.NK.cell,ndims = 50)

# Calculate PCA
subcluster <- T.NK.cell
calculate_pcs <- function(scobj) {
  # Extract standard deviations and calculate percentage contributions
  pct <- scobj[["pca"]]@stdev / sum(scobj[["pca"]]@stdev) * 100
  # Calculate cumulative contribution rates
  cumu <- cumsum(pct)
  # Identify the first PC where the cumulative contribution exceeds 90% and the individual PC contribution is below 5%
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Identify the first PC after consecutive PCs with a contribution change greater than 0.1%
  co2 <- if (length(pct) > 1) {
    sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1), decreasing = TRUE)[1] + 1
  } else {
    # If only one PC is present, this criterion is not used by default
    NA
  }
  # Select the number of PCs; use the minimum of co1 and co2, or co1 if co2 is unavailable
  pcs <- ifelse(is.na(co2), co1, min(co1, co2, na.rm = TRUE))
  # Create a data frame containing the required information
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
  # Return the result list
  list(
    pcs = pcs,
    co1 = co1,
    co2 = co2,
    plot_df = plot_df
  )
}
pca <- calculate_pcs(subcluster)$pcs
pca

library(harmony)
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)
library(clustree)
clustree(subcluster)
Idents(T.NK.) <- "RNA_snn_res.0.5"
subcluster$seurat_clusters <- subcluster@active.ident

DimPlot(subcluster,label = T,reduction = "umap",group.by = "RNA_snn_res.0.5",raster = F)

# Identify differentially expressed genes for each T/NK subcluster and annotate cells
DefaultAssay(subcluster) <- "RNA"
subcluster.markers  <- FindAllMarkers(subcluster, 
                                      only.pos = TRUE, 
                                      min.pct = 0.25, 
                                      logfc.threshold = 0.25)
subcluster.markers1  <- subcluster.markers [subcluster.markers $p_val_adj < 0.05, ]
subcluster.markers1_TOP30 <- subcluster.markers1 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
setwd(figure1_dir)
subcluster.markers1_TOP30 <- read.csv(file = file.path(figure1_dir, "T.NK_top30_markers_1.csv"),header = T,row.names = 1)
FeaturePlot(subcluster, features = c('CD3D', 'CD3E', 'CD8A', 'CD4','CD2'),raster=FALSE)
VlnPlot(subcluster,features = subcluster.markers1_TOP30$gene[361:380],raster = F,pt.size = 0)
# Examine the expression of known markers across clusters
# naive
naive <- c('CCR7','LEF1','SELL')
NKT <- AddModuleScore(T.NK.cell,features = list(naive),name = 'naive.scores')
p1 <- VlnPlot(NKT,features = 'naive.scores1',pt.size = 0,raster = F)+ labs(title = expression(atop('Naive T cells','(CCR7,LEF1,SELL)')))
p1
# NK
NK.genes <- c('NCAM1','NCR1','TYROBP','KLRD1','KLRF1')
NKT <- AddModuleScore(NKT,features = list(NK.genes),name = 'NK.scores')
p2 <- VlnPlot(NKT,features = 'NK.scores1',pt.size = 0,raster = F)+ labs(title = expression(atop('NK','(NCAM1,NCR1,TYROBP,KLRD1,KLRF1)')))
p2
# CD4 Treg
Treg.genes <- c('CD4','FOXP3','IL2RA','TNFRSF4')
NKT <- AddModuleScore(NKT,features = list(Treg.genes),name = 'Treg.scores')
p3 <- VlnPlot(NKT,features = 'Treg.scores1',pt.size = 0,raster = F)+ labs(title = expression(atop('Tregs','(CD4,FOXP3,IL2RA,TNFRSF4)')))
p3

# CD4 Tfh
cd4.Tfh.genes <- c('CD4','IL21','ICOS','CXCL13')
NKT <- AddModuleScore(NKT,features = list(cd4.Tfh.genes),name = 'CD4Tfh.scores')
p4 <- VlnPlot(NKT,features = 'CD4Tfh.scores1',pt.size = 0,raster = F)+ labs(title = expression(atop('CD4+ Tfh','(CD4,IL21,ICOS,CXCL13)')))
p4

# CD4 Tcm
cd4.Tcm.genes <- c('CD4','IL7R','ANXA1')
NKT <- AddModuleScore(NKT,features = list(cd4.Tcm.genes),name = 'CD4Tcm.scores')
p5 <- VlnPlot(NKT,features = 'CD4Tcm.scores1',pt.size = 0,raster = F)+ labs(title = expression(atop('CD4+ Tcm','(CD4,IL7R,ANXA1)')))
p5

# CD8 Tem
cd8.Tem.genes <- c('CD8A','CD8B','GZMK','CCL4')
NKT <- AddModuleScore(NKT,features = list(cd8.Tem.genes),name = 'CD8Tem.scores')
p6 <- VlnPlot(NKT,features = 'CD8Tem.scores1',pt.size = 0,raster = F)+ labs(title = expression(atop('CD8+ Tem','(CD8A,CD8B,GZMK,CCL4)')))
p6

# CD8 Trm
cd8.Trm.genes <- c('CD8A','CD8B','ZNF683','ITGA1')
NKT <- AddModuleScore(NKT,features = list(cd8.Trm.genes),name = 'CD8Trm.scores')
p7 <- VlnPlot(NKT,features = 'CD8Trm.scores1',pt.size = 0,raster = F)+ labs(title = expression(atop('CD8+ Trm','(CD8A,CD8B,ZNF683,ITGA1)')))
p7

# CD8 Teff
cd8.Teff.genes <- c('CD8A','CD8B','GZMH','CX3CR1','FGFBP2')
NKT <- AddModuleScore(NKT,features = list(cd8.Teff.genes),name = 'CD8Teff.scores')
p8 <- VlnPlot(NKT,features = 'CD8Teff.scores1',pt.size = 0,raster = F)+ labs(title = expression(atop('CD8+ Teff','(CD8A,CD8B,GZMH,CX3CR1,FGFBP2)')))
p8

# CD8 Tex
cd8.Tex.genes <- c('CD8A','CD8B','PDCD1','CTLA4','LAG3')
NKT <- AddModuleScore(NKT,features = list(cd8.Tex.genes),name = 'CD8Tex.scores')
p9 <- VlnPlot(NKT,features = 'CD8Tex.scores1',pt.size = 0,raster = F)+ labs(title = expression(atop('CD8+ Tex','(CD8A,CD8B,PDCD1,CTLA4,LAG3)')))
p9

# Annotate T/NK cell subclusters
new.cluster.ids <- c("0"="CD4T_01_IL7R", 
                     "1"="CD8T_01_CCL5",
                     "2"="CD4T_02_GNB2L1",
                     "3"="CD8T_02_GZMB",
                     "4"="CD4T_03_SCGB3A1",
                     "5"="Treg_FOXP3",
                     "6"="CD4T_04_SNORD3A",
                     "7"="CD4T_05_SFTPB",
                     "8"="NK_01_FGFBP2",
                     "9"="CD8T_03_MKI67",
                     "10"="CD4T_06_ACSL5",
                     "11"="CD8T_04_HLA-B",
                     "12"="CD4T_06_ACSL5"
)
T.NK.cell <- RenameIdents(subcluster, new.cluster.ids)                        
T.NK.cell$subcluster <- T.NK.cell@active.ident
T.NK.cell$subcluster <- factor(T.NK.cell$subcluster ,levels = c("CD4T_01_IL7R","CD4T_02_GNB2L1","CD4T_03_SCGB3A1","CD4T_04_SNORD3A","CD4T_05_SFTPB","CD4T_06_ACSL5","CD8T_01_CCL5","CD8T_02_GZMB","CD8T_03_MKI67","CD8T_04_HLA-B","Treg_FOXP3","NK_01_FGFBP2","Unidentified"))
#T.NK.cell <- readRDS("./T.NK.cell.rds")
#tcolor <- topo.colors(23)
tcolor <- c("#7eb6ed","#ffD65E","#d27ca0","#b4ffde",'#fae3ae',"#c9e2f6",'#be95db','#e16373','#c8c7e1','#ece366','#6585bc',"#b5ed98",'#5CB0C3','#80B1D3',"green","yellow","blue")
tcolor <- c("#7eb6ed","#d27ca0",'#fae3ae','#be95db','#e16373',"#ffD65E","#b4ffde",'#ece366',"#b5ed98","#c9e2f6",'#c8c7e1','#5CB0C3',"#99cb9d")
#tcolor <- C("#77c0a4","#c0afcd","#eeefb8","#cfe1cb","#88bec3","#ff7484","#ec8a55","#4486ab","#e8d79d","#bd9abd","#d78b95","#f5b576","#99cb9d")
DimPlot(T.NK.cell,group.by = "subcluster",raster = F,cols = tcolor,label = F)
DimPlot(T.NK.cell,group.by = "RNA_snn_res.0.5",raster = F,cols = tcolor,label = T)
saveRDS(T.NK.cell,file = "./T.NK.cell.rds")

#B
B.cell <- subset(sce,celltype == "B")

B.cell <- ScaleData(B.cell)
B.cell <- RunPCA(B.cell, features = VariableFeatures(object = B.cell))
ElbowPlot(B.cell,ndims = 50)

# Calculate PCA
subcluster <- B.cell
pca <- calculate_pcs(subcluster)$pcs
pca

library(harmony)
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)
library(clustree)
clustree(subcluster)
Idents(subcluster) <- "RNA_snn_res.0.4"
subcluster$seurat_clusters <- subcluster@active.ident

DimPlot(subcluster,label = T,reduction = "umap",group.by = "RNA_snn_res.0.4",raster = F)

# Identify differentially expressed genes for each B-cell subcluster and annotate cells
DefaultAssay(subcluster) <- "RNA"
subcluster.markers  <- FindAllMarkers(subcluster, 
                                      only.pos = TRUE, 
                                      min.pct = 0.25, 
                                      logfc.threshold = 0.25)
subcluster.markers1  <- subcluster.markers [subcluster.markers $p_val_adj < 0.05, ]
subcluster.markers1_TOP30 <- subcluster.markers1 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
setwd(figure1_dir)
write.csv(subcluster.markers1_TOP30, file = file.path(figure1_dir, "B_top30_markers_1.csv"))# Save
FeaturePlot(subcluster,features = c("CD19","CD79A","IGHG1","IGHA1","MZB1","MS4A1","CD20","CD79B","SDC1"),raster = F)
VlnPlot(subcluster,features = subcluster.markers1_TOP30$gene[421:440],raster = F,pt.size = 0)

# Annotate B-cell subclusters
new.cluster.ids <- c("0"="B_01_HLA-DRA", 
                     "1"="B_02_HLA.DQB1",
                     "2"="Plasma_05_IGHG_IGHGP",
                     "3"="B_03_CD3E",
                     "4"="Plasma_01_IGHA_MZB1",
                     "5"="Plasma_02_IGHA_IGLC3",
                     "6"="Plasma_06_IGHG_SFTPB",
                     "7"="B_04_PTMAP2",
                     "8"="B_05_FYN",
                     "9"="Plasma_03_IGHA_GUSBP11",
                     "10"="Plasma_07_IGHG_GRIFIN",
                     "11"="Plasma_08_IGHG_GUSBP11",
                     "12"="Plasma_04_IGHA_VPS13D",
                     "13"="Bn_06_GZMB",
                     "14"="Plasma_09_IGHG_RRM2"
)

B.cell <- RenameIdents(subcluster, new.cluster.ids)                        
B.cell$subcluster <- B.cell@active.ident
B.cell$subcluster <- factor(B.cell$subcluster ,levels = c("B_01_HLA-DRA","B_02_HLA.DQB1","B_03_CD3E","B_04_PTMAP2","B_05_FYN","Bn_06_GZMB","Plasma_01_IGHA_MZB1","Plasma_02_IGHA_IGLC3","Plasma_03_IGHA_GUSBP11","Plasma_04_IGHA_VPS13D","Plasma_05_IGHG_IGHGP",
                                                          "Plasma_06_IGHG_SFTPB","Plasma_07_IGHG_GRIFIN","Plasma_08_IGHG_GUSBP11","Plasma_09_IGHG_RRM2"))
#B.cell <- readRDS("./B.cell.rds")
bcolor <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#d6c2db','#c6adb0','#EAEAAE','#83ab8e','#cbdaa9','#92699e','#C8A8DA','#FFE5C0','#ece399','#D19275','#df6934','#537eb7')
#B.cell$subcluster <- factor(B.cell$subcluster,levels = c("B_01_RPS29","B_02_LINC01641","B_03_TNFRSF13B","B_04_Naive_TCL1A","B_05_CD27","B_06_S100A6","Plasma_01_IGHA1_IGHA2","Plasma_02_IGHA1_DERL3","Plasma_03_IGHA1_IL32","Plasma_04_IGHA1_FCER1G","Plasma_05_IGHA1_MT1E","Plasma_06_IGHG1_IGHG4","Plasma_07_IGHG1_SFTPA2","Plasma_08_IGHG1_IGHG2","Plasma_09_IGHG1_SYNJ1","Plasma_10_IGHG1_FER1L4"))
DimPlot(B.cell,group.by = "subcluster",raster = F,label = F,cols = bcolor)
DimPlot(B.cell,group.by = "RNA_snn_res.0.4",raster = F,label = T,cols = bcolor)
#save(B.cell,file = "B.cell.RData")
saveRDS(B.cell,"./B.cell.rds")

#Myeloid
Myeloid.cell <- subset(sce,celltype == "Myeloid")
Myeloid.cell <- ScaleData(Myeloid.cell)
Myeloid.cell <- RunPCA(Myeloid.cell, features = VariableFeatures(object = Myeloid.cell))
ElbowPlot(Myeloid.cell,ndims = 50)

# Calculate PCA
subcluster <- Myeloid.cell
pca <- calculate_pcs(subcluster)$pcs
pca

library(harmony)
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)
library(clustree)
clustree(subcluster)+coord_flip()+scale_color_manual(values = Mcolor)
Idents(subcluster) <- "RNA_snn_res.0.6"
subcluster$seurat_clusters <- subcluster@active.ident

DimPlot(subcluster,label = T,reduction = "umap",group.by = "RNA_snn_res.0.6",raster = F)

# Identify differentially expressed genes for each cell subcluster and annotate cells
DefaultAssay(subcluster) <- "RNA"
subcluster.markers  <- FindAllMarkers(subcluster, 
                                      only.pos = TRUE, 
                                      min.pct = 0.25, 
                                      logfc.threshold = 0.25)
subcluster.markers1  <- subcluster.markers [subcluster.markers $p_val_adj < 0.05, ]
#subcluster.markers1 <- read.csv(file.path(figure1_dir, "Myeloid_markers_1.csv"),header = T,row.names = 1)
subcluster.markers1_TOP30 <- subcluster.markers1 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
setwd(figure1_dir)
write.csv(subcluster.markers1_TOP30, file = file.path(figure1_dir, "Myeloid_top30_markers_1.csv"))# Save
VlnPlot(subcluster,features = c("KIT","HLA-DQA1","CD68","LILRA4","FCER1A","CD163","FOXP3","CSF3R","FCGR3B"),raster = F,pt.size = 0)
# Annotate myeloid cell subclusters
new.cluster.ids <- c("0"="Macrophage_01_LGMN", 
                     "1"="Macrophage_02_FCN1",
                     "2"="Neurophil_01_GOS2",
                     "3"="cDC_01_CLEC10A",
                     "4"="Macrophage_03_MARCO",
                     "5"="Macrophage_04_SPP1",
                     "6"="Macrophage_05_CSTB",
                     "7"="Mast_01_GNB2L1",
                     "8"="Macrophage_06_C3",
                     "9"="Macrophage_07_IGKC",
                     "10"="Macrophage_08_SCGB3A2",
                     "11"="Macrophage_09_MKI67",
                     "12"="cDC_02_CD207",
                     "13"="pDC_03_GZMB",
                     "14"="Mast_02_CHIT1",
                     "15"="Mast_03_PBK")


Myeloid.cell <- RenameIdents(subcluster, new.cluster.ids)                        
Myeloid.cell$subcluster <- Myeloid.cell@active.ident
Myeloid.cell$subcluster <- factor(Myeloid.cell$subcluster ,levels = c("Macrophage_01_LGMN","Macrophage_02_FCN1","Macrophage_03_MARCO","Macrophage_04_SPP1","Macrophage_05_CSTB","Macrophage_06_C3","Macrophage_07_IGKC","Macrophage_08_SCGB3A2","Macrophage_09_MKI67","Mast_01_GNB2L1","Mast_02_CHIT1","Mast_03_PBK","Neurophil_01_GOS2","cDC_01_CLEC10A","cDC_02_CD207","pDC_03_GZMB"))

Mcolor <- c("#93b8db","#a4cba8","#7eb4c6","#f4806e","#f5e886","#d6d5b7","#7898e1","#85bf67","#af87fe","#80c6d2","#f8cb7f","#a8d8b8",
            '#aec7e8','#5ac1b3','#57c3f3','#F58D93','#f7b6d2','#F7C394','#ade87c','#7382BC','#F0E442','#FD7014','#EAB67D','#58a4c9',
            '#DF75AE','#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175','#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3','#476D87', '#E95C59')
#Myeloid.cell <- readRDS("./Myeloid.cell.rds")
DimPlot(Myeloid.cell,group.by = "subcluster",raster = F,label = F,cols = Mcolor)
DimPlot(Myeloid.cell,group.by = "RNA_snn_res.0.6",raster = F,label = T,cols = Mcolor)
saveRDS(Myeloid.cell,file = file.path(figure1_dir, "Myeloid.cell.rds"))


#Endothelial
Endothelial <- subset(sce,celltype=="Endothelial")
Endothelial <- ScaleData(Endothelial)
Endothelial <- RunPCA(Endothelial, features = VariableFeatures(object = Endothelial))
ElbowPlot(Endothelial,ndims = 50)

# Calculate PCA
subcluster <- Endothelial
pca <- calculate_pcs(subcluster)$pcs
pca

library(harmony)
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)
library(clustree)
clustree(subcluster)
Idents(subcluster) <- "RNA_snn_res.0.2"
subcluster$seurat_clusters <- subcluster@active.ident

DimPlot(subcluster,label = T,reduction = "umap",group.by = "RNA_snn_res.0.2",raster = F)

# Identify differentially expressed genes for each cell subcluster and annotate cells
DefaultAssay(subcluster) <- "RNA"
subcluster.markers  <- FindAllMarkers(subcluster, 
                                      only.pos = TRUE, 
                                      min.pct = 0.25, 
                                      logfc.threshold = 0.25)
subcluster.markers1  <- subcluster.markers [subcluster.markers $p_val_adj < 0.05, ]
#subcluster.markers1 <- read.csv(file.path(figure1_dir, "Myeloid_markers_1.csv"),header = T,row.names = 1)
subcluster.markers1_TOP30 <- subcluster.markers1 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
setwd(figure1_dir)
write.csv(subcluster.markers1_TOP30, file = file.path(figure1_dir, "Endothelial_markers_1.csv"))# Save
VlnPlot(subcluster,features = subcluster.markers1_TOP30$gene[241:260],pt.size = 0,raster = F)
# Annotate cell subclusters
new.cluster.ids <- c("0"= "Endothelial_01_FCN3",
                     "1"= "Endothelial_02_PTMAP2",
                     "2"= "Endothelial_03_VCAM1",
                     "3"= "Endothelial_04_DKK2",
                     "4"= "Endothelial_05_IL1RL1",
                     "5"= "Endothelial_06_CD68",
                     "6"= "Endothelial_07_NOTCH3",
                     "7"= "Endothelial_08_CCL21",
                     "8"= "Endothelial_09_S100A14",
                     "9"= "Endothelial_02_PTMAP2")

Endothelial <- RenameIdents(subcluster, new.cluster.ids)                        
Endothelial$subcluster <- Endothelial@active.ident
Endothelial$subcluster <- factor(Endothelial$subcluster ,levels = c("Endothelial_01_FCN3","Endothelial_02_PTMAP2","Endothelial_03_VCAM1","Endothelial_04_DKK2","Endothelial_05_IL1RL1","Endothelial_06_CD68","Endothelial_07_NOTCH3","Endothelial_08_CCL21","Endothelial_09_S100A14"))

Ecolor <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3','#476D87', '#E95C59','#F7C394','#ade87c')

DimPlot(Endothelial,group.by = "subcluster",raster = F,label = F,cols = Ecolor)
DimPlot(Endothelial,group.by = "RNA_snn_res.0.2",raster = F,label = T,cols = Mcolor)
saveRDS(Endothelial,file = "./Endothelial.rds")

#Fibroblast
Fibroblast <- subset(sce,celltype=="Fibroblast")
Fibroblast <- ScaleData(Fibroblast)
Fibroblast <- RunPCA(Fibroblast, features = VariableFeatures(object = Fibroblast))
ElbowPlot(Fibroblast,ndims = 50)
# Calculate PCA
subcluster <- Fibroblast
pca <- calculate_pcs(subcluster)$pcs
pca

library(harmony)
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)
library(clustree)
clustree(subcluster)
Idents(subcluster) <- "RNA_snn_res.0.1"
subcluster$seurat_clusters <- subcluster@active.ident

DimPlot(subcluster,label = T,reduction = "umap",group.by = "RNA_snn_res.0.1",raster = F)

# Identify differentially expressed genes for each cell subcluster and annotate cells
DefaultAssay(subcluster) <- "RNA"
subcluster.markers  <- FindAllMarkers(subcluster, 
                                      only.pos = TRUE, 
                                      min.pct = 0.25, 
                                      logfc.threshold = 0.25)
subcluster.markers1  <- subcluster.markers [subcluster.markers $p_val_adj < 0.05, ]
#subcluster.markers1 <- read.csv(file.path(figure1_dir, "Myeloid_markers_1.csv"),header = T,row.names = 1)
subcluster.markers1_TOP30 <- subcluster.markers1 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
setwd(figure1_dir)
write.csv(subcluster.markers1_TOP30, file = file.path(figure1_dir, "Fibroblast_markers_1.csv"))# Save
VlnPlot(subcluster,features = subcluster.markers1_TOP30$gene[1:20],pt.size = 0,raster = F)
# Annotate cell subclusters
# Score alveolar, myofibroblast, and adventitial fibroblast marker gene sets
Alveolar <- list(c("BMP5","MYH10","ANGPT1","PLEKHH2","FGFR4","TCF21","MACF1","A2M","LIMCH1","RGCC"))
Myo <- list(c("EVA1A","LEF1","GJB2","KIAA1217","ADAM12","LOXL2","VCAM1","INHBA","CTHRC1","POSTN"))
Adventitial <- list(c("TNNT3","GPNMB","EFEMP1","SFRP1","FBLN1","SCARA5","OGN","CCDC80","PI16","MFAP5"))
Inscore <- AddModuleScore(subcluster,
                          features = Alveolar,
                          ctrl = 100,
                          name = "gene")


colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[24] <- 'Alveolar'
Inscore <- AddModuleScore(Inscore,
                          features = Myo,
                          ctrl = 100,
                          name = "gene")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[25] <- 'Myo'
Inscore <- AddModuleScore(Inscore,
                          features = Adventitial,
                          ctrl = 100,
                          name = "gene")

colnames(Inscore@meta.data)[26] <- 'Adventitial'

# A Seurat object named Inscore is constructed and visualized
VlnPlot(Inscore,features = 'Alveolar', 
        pt.size = 0, adjust = 2,group.by = "RNA_snn_res.0.2",cols = Fcolor,raster=FALSE)
VlnPlot(Inscore,features = 'Myo', 
        pt.size = 0, adjust = 2,group.by = "RNA_snn_res.0.2",cols = Fcolor,raster=FALSE)
VlnPlot(Inscore,features = 'Adventitial', 
        pt.size = 0, adjust = 2,group.by = "RNA_snn_res.0.2",cols = Fcolor,raster=FALSE)

new.cluster.ids <- c("0"= "Fibroblast_04_Alveolar_DPT",
                     "1"= "Fibroblast_01_Myo_HIGD1B",
                     "2"= "Fibroblast_02_Myo_COL10A1",
                     "3"= "Fibroblast_05_Alveolar_MYH11",
                     "4"= "Fibroblast_06_Alveolar_MS4A7",
                     "5"= "Fibroblast_03_Myo_MKI67",
                     "6"= "Unidentified")

Fibroblast <- RenameIdents(subcluster, new.cluster.ids)                        
Fibroblast$subcluster <- Fibroblast@active.ident
Fibroblast$subcluster <- factor(Fibroblast$subcluster ,levels = c("Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67","Fibroblast_04_Alveolar_DPT","Fibroblast_05_Alveolar_MYH11","Fibroblast_06_Alveolar_MS4A7","Unidentified"))

Fcolor <- c('#E5D2DD', '#F1BB72', '#F3B1A0', '#D6E7A3','lightblue','#F7C394','#ade87c', '#E95C59', '#57C3F3','#53A85F')

DimPlot(Fibroblast,group.by = "subcluster",raster = F,label = F,cols = Fcolor)
DimPlot(Fibroblast,group.by = "RNA_snn_res.0.1",raster = F,label = T,cols = Fcolor)
saveRDS(Fibroblast,file = "./Fibroblast.rds")


#Astrocyte
Astrocyte <- subset(sce,celltype=="Astrocyte")
Astrocyte <- ScaleData(Astrocyte)
Astrocyte <- RunPCA(Astrocyte, features = VariableFeatures(object = Astrocyte))
ElbowPlot(Astrocyte,ndims = 50)
# Calculate PCA
subcluster <- Astrocyte
pca <- calculate_pcs(subcluster)$pcs
pca
library(harmony)
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)
library(clustree)
clustree(subcluster)
Idents(subcluster) <- "RNA_snn_res.0.1"
subcluster$seurat_clusters <- subcluster@active.ident

DimPlot(subcluster,label = T,reduction = "umap",group.by = "RNA_snn_res.0.2",raster = F)

# Identify differentially expressed genes for each cell subcluster and annotate cells
DefaultAssay(subcluster) <- "RNA"
subcluster.markers  <- FindAllMarkers(subcluster, 
                                      only.pos = TRUE, 
                                      min.pct = 0.25, 
                                      logfc.threshold = 0.25)
subcluster.markers1  <- subcluster.markers [subcluster.markers $p_val_adj < 0.05, ]
#subcluster.markers1 <- read.csv(file.path(figure1_dir, "Myeloid_markers_1.csv"),header = T,row.names = 1)
subcluster.markers1_TOP30 <- subcluster.markers1 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
setwd(figure1_dir)
write.csv(subcluster.markers1_TOP30, file = file.path(figure1_dir, "Astrocyte_markers_1.csv"))# Save
VlnPlot(subcluster,features = subcluster.markers1_TOP30$gene[1:230],pt.size = 0,raster = F)
# Annotate cell subclusters
new.cluster.ids <- c("0"= "Astrocyte_01_CAV1",
                     "1"= "Astrocyte_02_MALAT1",
                     "2"= "Astrocyte_03_ROM1",
                     "3"= "Astrocyte_04_SCGB3A1",
                     "4"= "Astrocyte_05_HLA-DRB5",
                     "5"= "Astrocyte_06_CCL5")

Astrocyte <- RenameIdents(subcluster, new.cluster.ids)                        
Astrocyte$subcluster <- Astrocyte@active.ident


Acolor <- c('#f1b38a', '#f0db69', '#a2c246', '#f5cee0', '#b6d2b7','#cbd2e5')

DimPlot(Astrocyte,group.by = "subcluster",raster = F,label = F,cols = Acolor)
DimPlot(Astrocyte,group.by = "RNA_snn_res.0.1",raster = F,label = T,cols = Acolor)
saveRDS(Astrocyte,file = "./Astrocyte.rds")


# Stromal cells
Endothelial <- readRDS("./Endothelial.rds")
Fibroblast <- readRDS("./Fibroblast.rds")
Astrocyte <- readRDS("./Astrocyte.rds")
stromal <- merge(x=Endothelial,y=c(Fibroblast,Astrocyte))
stromal <- JoinLayers(stromal)
stromal <- ScaleData(stromal)
stromal <- RunPCA(stromal, features = VariableFeatures(object = stromal))
ElbowPlot(stromal,ndims = 50)
# Calculate PCA
subcluster <- stromal
pca <- calculate_pcs(subcluster)$pcs
pca
library(harmony)
subcluster <- RunHarmony(subcluster,reduction = "pca",group.by.vars = "dataset",reduction.save = "harmony")
subcluster <- subcluster %>% FindNeighbors(reduction = "harmony", dims = 1:pca) %>%
  FindClusters(resolution = seq(from = 0.1, 
                                to = 1.0, 
                                by = 0.1))%>%RunUMAP(reduction = "harmony", dims = 1:pca)
library(clustree)
clustree(subcluster)
Idents(subcluster) <- "RNA_snn_res.0.6"
Scolor <- c('#f1b38a', '#f0db69', '#a2c246', '#f5cee0', '#b6d2b7','#cbd2e5',
            "#F59B7B","#F0B0AD","#ABD7EC","#57AF37","#FCC41E",'#D9E9A9','#07C3F3',"#8D73BA",'#ade87c',"#33ABC1","#8FB4DC","#FFDD8E","#70CDBE","#AC99D2","#ED8828","#7AC3DF","#D9DEE7","#C74546","lightpink")

DimPlot(subcluster,group.by = "subcluster",raster = F,label = F,cols = Scolor) 
saveRDS(subcluster,file = "./stromal.rds")

# Merge T/NK, B, Myeloid, and Stromal cells
rm(list = ls())
setwd(figure1_dir)
T.NK <- readRDS("./T.NK.cell.rds")
B <- readRDS("./B.cell.rds")
Myeloid <- readRDS("./Myeloid.cell.rds")
Endothelial <- readRDS("./Endothelial.rds")
Fibroblast <- readRDS("./Fibroblast.rds")
Astrocyte <- readRDS("./Astrocyte.rds")
TIME <- merge(x = T.NK,y = c(B,Myeloid,Endothelial,Fibroblast,Astrocyte))
TIME <- JoinLayers(TIME)
saveRDS(TIME,"./TIME.rds")

# Further classify major cell types based on subclustering results, such as separating T.NK into T and NK, and Myeloid into Macrophage, Neutrophil, DC, and Mast cells
library(dplyr)
setwd(figure1_dir)
TIME <- readRDS("./TIME.rds")
unidentified <- c("Unidentified")
TIME <- TIME[,!TIME$subcluster %in% unidentified]
TIME@meta.data <- TIME@meta.data %>% mutate(CELLTYPE = case_when(
  subcluster %in% c("CD4T_01_IL7R","CD4T_02_GNB2L1","CD4T_03_SCGB3A1","CD4T_04_SNORD3A","CD4T_05_SFTPB","CD4T_06_ACSL5") ~ "CD4T",
  subcluster %in% c("CD8T_01_CCL5","CD8T_02_GZMB","CD8T_03_MKI67","CD8T_04_HLA-B") ~ "CD8T",
  subcluster %in% c("Treg_FOXP3") ~ "Treg",
  subcluster %in% c("NK_01_FGFBP2") ~ "NK",
  subcluster %in% c("B_01_HLA-DRA","B_02_HLA.DQB1","B_03_CD3E","B_04_PTMAP2","B_05_FYN","Bn_06_GZMB","Plasma_01_IGHA_MZB1","Plasma_02_IGHA_IGLC3","Plasma_03_IGHA_GUSBP11","Plasma_04_IGHA_VPS13D","Plasma_05_IGHG_IGHGP",
                    "Plasma_06_IGHG_SFTPB","Plasma_07_IGHG_GRIFIN","Plasma_08_IGHG_GUSBP11","Plasma_09_IGHG_RRM2") ~ "B",
  subcluster %in% c("Macrophage_01_LGMN","Macrophage_02_FCN1","Macrophage_03_MARCO","Macrophage_04_SPP1","Macrophage_05_CSTB","Macrophage_06_C3","Macrophage_07_IGKC","Macrophage_08_SCGB3A2","Macrophage_09_MKI67") ~ "Macrophage",
  subcluster %in% c("Neurophil_01_GOS2") ~ "Neurophil",
  subcluster %in% c("cDC_01_CLEC10A","cDC_02_CD207") ~ "DC",
  subcluster %in% c("Mast_01_GNB2L1","Mast_02_CHIT1","Mast_03_PBK") ~ "Mast",
  subcluster %in% c("Endothelial_01_FCN3","Endothelial_02_PTMAP2","Endothelial_03_VCAM1","Endothelial_04_DKK2","Endothelial_05_IL1RL1","Endothelial_06_CD68","Endothelial_07_NOTCH3","Endothelial_08_CCL21","Endothelial_09_S100A14") ~ "Endothelial",
  subcluster %in% c("Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67","Fibroblast_04_Alveolar_DPT","Fibroblast_05_Alveolar_MYH11","Fibroblast_06_Alveolar_MS4A7") ~ "Fibroblast",
  subcluster %in% c("Astrocyte_01_CAV1","Astrocyte_02_MALAT1","Astrocyte_03_ROM1","Astrocyte_04_SCGB3A1","Astrocyte_05_HLA-DRB5","Astrocyte_06_CCL5") ~ "Astrocyte"
))
table(TIME$CELLTYPE)

# Calculate differentially expressed genes for each subcluster
Idents(TIME) <- "subcluster"
DefaultAssay(TIME) <- "RNA"
TIME.markers <- FindAllMarkers(TIME, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25)
TIME.markers1 <- TIME.markers[TIME.markers$p_val_adj < 0.05, ]
# Top 5 genes for each subcluster
TIME.markers1_TOP5 <- TIME.markers1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

