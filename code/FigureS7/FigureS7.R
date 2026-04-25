# ------------- Figure S7 --------------
#----Figure S7A----
library(Seurat)
library(ggplot2)
library(patchwork)

rm(list=ls());gc()
#1.ST validate prognostic marker----
visiumPath <- "./ST/P10-B1/"
visiumPath <- "./ST/P10-T1/"
visiumPath <- "./ST/P10-T2/"
visiumPath <- "./ST/P11-B1/"
visiumPath <- "./ST/P11-T1/"
visiumPath <- "./ST/P15-B1/"
visiumPath <- "./ST/P15-T1/"
visiumPath <- "./ST/P16-B1/"
visiumPath <- "./ST/P16-T1/"
visiumPath <- "./ST/P17-B1/"
visiumPath <- "./ST/P17-T1/"
visiumPath <- "./ST/P19-B1/"
visiumPath <- "./ST/P19-T1/"
visiumPath <- "./ST/P24-B1/"
visiumPath <- "./ST/P24-T1/"
visiumPath <- "./ST/P25-B1/"
visiumPath <- "./ST/P25-T1/"

visiumPath <- "./ST/D1/"
visiumPath <- "./ST/D2/"

seu <- Load10X_Spatial(
  data.dir = visiumPath,
  filename = "filtered_feature_bc_matrix.h5"
)

# Unify color
mycols <- c("#3E4A89","#26828E","#6DCD59","#B4DE2C","#FDE725")

# Encapsulated function
plot_gene <- function(seu, gene, title="P19-T2"){
  SpatialFeaturePlot(
    seu,
    features = gene,
    image.alpha = 0,
    slot = "counts",
    pt.size.factor = 2
  ) +
    scale_fill_gradientn(colors = mycols) +
    coord_fixed() +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "black", color = NA)
    )
}

# gene list
genes <- c("USF2","IGF1","SP1","ITGA6","SNAI1","ITGB4",
           "RELA","PPARG","RARRES2","PAX6","ITGB1",
           "NR1H4","CEBPB","RETN")

plots <- lapply(genes, function(g) plot_gene(seu, g))
plots[[1]]# View the first one
wrap_plots(plots, ncol = 4)

SpatialFeaturePlot(
  seu,
  features = c("USF2", "IGF1")
)

#Ligand-receptor co-localization----
library(SpaGene)
library(dplyr)

count <- seu@assays$Spatial$counts# count matrix
location <- GetTissueCoordinates(seu)# Coordinates
location <- location %>%
  rename(imagecol = x, imagerow = y)# Change colnames
rownames(location) <- location$cell# rownames=barcode

location <- location[, c("imagecol", "imagerow")]

new_df <- data.frame(
  ligand_gene_symbol = c("USF2", "SP1", "SNAI1", "RELA", "PPARG", "CEBPB","CEBPB"),
  receptor_gene_symbol = c("IGF1", "ITGA6", "ITGB4", "IGF1", "RARRES2", "RETN", "IGF1")
)
obj_lr <- SpaGene_LR(count, location, LRpair = new_df)
plot_lr_pair <- function(ligand, receptor, file_prefix = "P19-T2") {
  
  p <- plotLR(
    count,
    location,
    LRpair = c(ligand, receptor),
    alpha.min = 0.2,
    pt.size = 1
  ) +
    ggtitle(paste0(ligand, " - ", receptor)) +
    theme(
      plot.title = element_text(hjust = 0.5, color = "white"),
      panel.background = element_rect(fill = "black"),
      plot.background = element_rect(fill = "black"),
      panel.grid = element_blank(),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white")
    )
  
  return(p)
}
gc()
library(patchwork)

p_list <- list(
  plot_lr_pair("USF2", "IGF1"),
  plot_lr_pair("SP1", "ITGA6"),
  plot_lr_pair("SNAI1", "ITGB4"),
  plot_lr_pair("RELA", "IGF1"),
  plot_lr_pair("PPARG", "RARRES2"),
  plot_lr_pair("PAX6", "ITGB1"),
  plot_lr_pair("NR1H4", "RARRES2"),
  plot_lr_pair("CEBPB", "RETN"),
  plot_lr_pair("CEBPB", "IGF1")
)
setwd("./ST/")

gc()
pdf("./P24-T1.pdf", width = 17.5, height = 8)
wrap_plots(p_list, ncol = 3)
dev.off()
