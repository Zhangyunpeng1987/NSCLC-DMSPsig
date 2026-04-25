# ------------- Figure S1 --------------
# -------- Path settings --------
project_dir <- "./NSCLC-DMSPsig"
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
figure1_dir <- file.path(results_dir, "Figure1")
figure2_dir <- file.path(results_dir, "Figure2")
dir.create(figure1_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure2_dir, showWarnings = FALSE, recursive = TRUE)

#----Figure S1A----
# Load the annotated Seurat object for marker visualization
sce <- readRDS(file.path(data_dir, "sce_annotated_umap.rds"))
# Display marker genes on the UMAP
library(ggplot2)
library(pheatmap)
gene <- c("TRAC","NKG7","CD79A","CSF3R","CPA3","CD163","CLEC10A","VWF","COL1A2","EPCAM","PLP1")
color <- c('lightgrey','#a6261f')# Set colors  
FeaturePlot(sce, features = 'APOC1',cols = color, pt.size = 1)+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))# Add border 
#----Figure S1B----
# Display dataset, mutation, and patient information on UMAP
rm(list = setdiff(ls(), c("project_dir", "data_dir", "results_dir", "figure1_dir", "figure2_dir")))
library(Seurat)
sce <- readRDS(file.path(data_dir, "sce_annotated_umap.rds"))
generate_colors <- function(n) {
  colorspace::rainbow_hcl(n)
}

# Generate colors for datasets
set.seed(123)  # Ensure reproducibility
dataset <- unique(sce$dataset)
#dataset_colors <- setNames(generate_colors(length(dataset)), dataset)
library(RColorBrewer)
dataset_colors <- RColorBrewer::brewer.pal(7,"Set3")
dataset_colors <- c("#8DD3C7","#FFFFB3","#FDB462","#B3DE69","#BEBADA","#FB8072","#80B1D3")
sce$patient <- factor(sce$patient,levels = c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10",
                                             "P11","P12","P13","P14","P15","P16","P17","P18",
                                             "P19","P20","P21","P22","P23","P24","P25","P26",
                                             "P27","P28","P29","P30","P31","P32","P33","P34",
                                             "P35","P36","P37","P38","P39","P40","P41",
                                             "P42", "P43", "P44", "P45"))
sce$mutation <- factor(sce$mutation,
                       levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation",
                                  "ALK","ROS1","TP53","MET-BM","HER2"))
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
sce$mutation <- factor(sce$mutation, levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation",
                                                "ALK","ROS1","TP53","MET-BM","HER2"))

mutation_colors = c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")

DimPlot(sce,group.by = "dataset",label = F,raster = F,cols = dataset_colors)
DimPlot(sce,group.by = "mutation",label = F,raster = F,cols = mutation_colors)
DimPlot(sce,group.by = "patient",label = F,raster = F,cols = rep_colors)
#----Figure S1C----
library(ggpubr)
sce$celltype <- factor(sce$celltype,levels = c("Myeloid","B","T/NK","Epithelial","Fibroblast","Endothelial","Astrocyte"))
cluster_colors <- c("#70cdbe","#ffdd8e","#eb7e60","#ac99d2","#7ac3df","#f5aa61","#8fb4dc")
#cluster_colors <- c("#ac99d2","#7ac3df","#f5aa61","#ffdd8e","#eb7e60","#70cdbe","#8fb4dc")
P1 <- 
  sce@meta.data %>%
  ggplot(aes(x = orig.ident, fill = celltype)) +
  theme_pubr(base_size = 5) +
  theme(plot.title = element_blank(),
        text = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(3,"pt"),
        legend.position = "right") +
  facet_grid(~mutation, scales = "free_x", space ="free_x") +
  geom_bar( alpha = 0.9,position = "fill")+
  #geom_bar( alpha = 0.9,colour="white")+
  scale_fill_manual(values = cluster_colors)+
  xlab("patient") + ylab("Fraction") 
P1
saveRDS(sce, file = file.path(data_dir, "sce_annotated_umap.rds"))
#----Figure S1D----
# Draw tissue preference of cell populations across TIME subtypes
rm(list = setdiff(ls(), c("project_dir", "data_dir", "results_dir", "figure1_dir", "figure2_dir")))
TIME <- readRDS(file.path(figure1_dir, "subsce.corr.rds"))
# Phenotype enrichment of TIME clusters
TIME$mutation <- factor(TIME$mutation, levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","HER2","ALK","MET-BM","ROS1","TP53"))
res.chisq <- chisq.test(table(TIME$subcluster, TIME$mutation))
R.oe <- (res.chisq$observed) / (res.chisq$expected)
write.csv(R.oe, file = file.path(figure2_dir, "oe.csv"))
R.oe <- read.csv(file.path(figure2_dir, "oe.csv"), row.names = 1)
pdf(file = file.path(figure1_dir, "FigureS1D_TIME_cluster_tissue_preference.pdf"), width = 5.1, height = 15)
pheatmap(R.oe,
         cluster_rows = T,
         color = c("#FDE5CD", "#FBBE8B", "#F38F4B","#E5520C"),
         breaks = c(0, 1, 1.5, 3, max(R.oe)),
         cluster_cols = F,
         angle_col = 45,
         fontsize = 8,
         border_color = "white",
         display_numbers = matrix(ifelse(R.oe > 3, "+++", ifelse(R.oe > 1.5, "++", ifelse(R.oe > 1, "+", "+/-"))), nrow(R.oe)),
         number_color = "black"
)
dev.off()

