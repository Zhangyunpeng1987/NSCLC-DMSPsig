# ------------- Figure 3 --------------
#----Figure 3A---- 
#Analysis of Cancer Cell Status
#CopyKat Identify malignant cells
rm(list = ls())
setwd("./t100553/wss")
sce <- readRDS("./sce_Annotated_umap.rds")
dir.create("./copykat/")
setwd("./copykat")

#corrsce <- readRDS("./subsce.corr.rds")
epi <- subset(sce,celltype=="Epithelial")
epi$CELLTYPE <- "Epithelial"
exp.rawdata <- as.matrix(epi@assays$RNA$counts)
library(copykat)
copykat.test <- copykat(rawmat = exp.rawdata,
                        id.type = "S",  # Gene ID type（gene Symbol: "S"，）
                        ngene.chr = 5,  # At least 5 genes per chromosome are used to calculate copy number
                        win.size = 25, # Each window contains at least 25 genes
                        KS.cut = 0.1, # Segmentation parameters（0.05-0.15）
                        sam.name = "test",
                        distance = "euclidean",  # Cluster distance type（euclidean / spearman / pearson）
                        norm.cell.names = NULL,  # Known normal cell names（default NULL）
                        output.seg = "FALSE",  # if output IGV visualization files
                        plot.genes = "TRUE",  # gene name displayed in the heatmap
                        genome = "hg20",  # Reference genome version
                        n.cores = 4)  # Number of parallel computing cores



save(copykat.test,file = "./copykat.test.RData")
test_copykat_prediction <- read.table("test_copykat_prediction.txt")
table(test_copykat_prediction$V2)

#Extract tumor cells from it
tumor <- subset(test_copykat_prediction,V2=="aneuploid")
tumor <- tumor$V1
malig.epi <- epi[,colnames(epi) %in% tumor]
saveRDS(malig.epi,file = "malig.epi.rds")

#Identify the state of cancer cells pmid3861409
#The method proposed using PMID38614094
library(Seurat)
library(BiocManager)
library(GEOquery) 
library(plyr)
library(dplyr) 
library(Matrix)
library(Seurat)
library(ggplot2)
library(cowplot) 
library(multtest)
library(msigdbr)
library(fgsea)
library(loomR)
library(clustree)
library(tibble)
library(SeuratData)
library(matrixStats)
library(sparseMatrixStats)
library(DESeq2)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(viridis)
library(gridExtra)
library(ggplotify)
library(multtest)
library(metap)
library(writexl)
library(Rcpp)
library(RcppZiggurat)
library(Rfast)
library(ggh4x)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationHub)
library(cola)
library(msigdbr)
library(UCell)
library(RColorBrewer)
rm(list = ls())
setwd("./t100553/wss/copykat")
cancer.epi <- readRDS("./malig.epi.rds")

# Functions
# DEG_remove_mito 
DEG_Remove_mito <- function (df){
  df_rm_mito <- df[!grepl("^MT-|^MT.",df$gene),]
  return(df_rm_mito)
}

# DEG_remove_heat 
DEG_Remove_heat <- function (df){
  df_rm_hsp <- df[!grepl("^HSP",rownames(df)),]
  return(df_rm_hsp)
}

# plotSimilarityMatrix 
# create a similarity matrix plot from a dataframe
# adapted from klic package (adjusted heatmap parameters)
plotSimilarityMatrix = function(X, y = NULL, 
                                min.val = 0, 
                                max.val = 1,
                                clusLabels = NULL, 
                                colX = NULL, colY = NULL, 
                                clr = FALSE, clc = FALSE, 
                                annotation_col = NULL, 
                                annotation_row = NULL, 
                                annotation_colors = NULL, 
                                myLegend = NULL, 
                                fileName = "posteriorSimilarityMatrix", 
                                savePNG = FALSE, 
                                semiSupervised = FALSE, 
                                showObsNames = FALSE) {
  
  if (!is.null(y)) {
    # Check if the rownames correspond to the ones in the similarity matrix
    check <- sum(1 - rownames(X) %in% row.names(y))
    if (check == 1)
      stop("X and y must have the same row names.")
  }
  
  if (!is.null(clusLabels)) {
    if (!is.integer(clusLabels))
      stop("Cluster labels must be integers.")
    
    n_clusters <- length(table(clusLabels))
    riordina <- NULL
    for (i in 1:n_clusters) {
      riordina <- c(riordina, which(clusLabels == i))
    }
    
    X <- X[riordina, riordina]
    y <- y[riordina, ]
    y <- as.data.frame(y)
  }
  
  if (savePNG)
    grDevices::png(paste(fileName, ".png", sep = ""))
  
  if (!is.null(y)) {
    ht <- ComplexHeatmap::pheatmap(X, legend = TRUE,  
                                   color = rev(brewer.pal(11, "RdBu")), 
                                   breaks = seq(min.val, max.val, length.out = 11),
                                   cluster_rows = clr, 
                                   cluster_cols = clc, 
                                   #annotation_col = y,
                                   show_rownames = showObsNames, 
                                   show_colnames = showObsNames, 
                                   drop_levels = TRUE, 
                                   treeheight_row = -1,
                                   treeheight_col = -1,
                                   #annotation_row = annotation_row, 
                                   #annotation_col = annotation_col, 
                                   annotation_colors = annotation_colors)
  } else {
    ht <- ComplexHeatmap::pheatmap(X, legend = TRUE,
                                   color = rev(brewer.pal(11, "RdBu")),  
                                   breaks = seq(min.val, max.val, length.out = 11),
                                   cluster_rows = clr, 
                                   cluster_cols = clc,
                                   show_rownames = showObsNames, 
                                   show_colnames = showObsNames,
                                   treeheight_row = -1,
                                   treeheight_col = -1,
                                   drop_levels = TRUE, 
                                   #annotation_row = annotation_row, 
                                   #annotation_col = annotation_col, 
                                   annotation_colors = annotation_colors)
  }
  
  if (savePNG)
    grDevices::dev.off()
  
  return(ht)
}

 
# calculate similarity matrix 
simil <- function(df, drop, file, method) {
  if (length(drop) > 0) {
    jc <- df[-drop, , drop = TRUE]
  }
  else {
    jc <- df
  }
  
  jc <- as.matrix(jc)
  
  if (method == "jaccard") {
    jc <- prabclus::jaccard(jc)
    jc <- 1 - jc
  }
  
  # calculate correlation matrix if method == "corr"
  else if (method == "corr") {
    jc <- Rfast::cora(jc)
  }
  
  saveRDS(jc, file = file) 
  
  v <- c(quantile(as.vector(jc), na.rm = TRUE), mean(as.vector(jc), na.rm = TRUE))
  names(v) <- c('0%','25%', '50%', '75%', '100%', 'mean')
  print(v)
  return(v)
}

# plot interactive similarity heatmap 

simil_plot <- function(a, min.val, max.val, annot) {
  jc <- readRDS(a) #read in similarity matrix
  
  # add annotation
  if (length(annot) > 0) {
    row_annot <- annot[rownames(jc), , drop = FALSE] 
    col_annot <- annot[colnames(jc), , drop = FALSE] 
    colors <- mako(n_distinct(annot)) 
    names(colors) <- base::unique(annot)[[1]]
    colors <- list(colors, colors)
    names(colors) <- c(as.name(names(annot)), as.name(names(annot)))
  } 
  else {
    row_annot <- NULL
    col_annot <- NULL
    colors <- NULL
  }
  
  hm <- plotSimilarityMatrix(jc, clr = TRUE, clc = TRUE, 
                             min.val = min.val, max.val = max.val,
                             annotation_row = row_annot, 
                             annotation_col = col_annot, 
                             annotation_colors = colors,
                             showObsNames = T) # plot full matrix
  return(hm)
}


simil_GE <- function(df, file, method) {
  jc <- df[GElist, , drop = FALSE] # subset to sc50 genes

  jc <- as.matrix(jc)
  
  if (method == "jaccard") {
    jc <- prabclus::jaccard(jc)
    jc <- 1 - jc
  }
  
  else if (method == "corr") {
    jc <- Rfast::cora(jc)
  }
  
  saveRDS(jc, file = file) 
  
  # return quantiles and mean similarity value
  v <- c(quantile(as.vector(jc), na.rm = TRUE), mean(as.vector(jc), na.rm = TRUE))
  names(v) <- c('0%','25%', '50%', '75%', '100%', 'mean')
  print(v)
  return(v)
}


library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
DefaultAssay(cancer.epi) <- "RNA"
cancer.epi <- NormalizeData(cancer.epi)
cancer.epi <- FindVariableFeatures(cancer.epi, nfeatures = 2000)
cancer.epi <- ScaleData(cancer.epi, features = rownames(cancer.epi))
meta <- FetchData(cancer.epi, vars = c("mutation", "patient"))
drop_patients <- meta %>% count(patient) %>% filter(n < 50)
meta <- meta %>% filter(!patient %in% drop_patients$patient)
meta <- meta %>% group_by(mutation, patient) %>% summarise(Nb = n(), .groups = "drop") %>% mutate(C = sum(Nb), percent = Nb / C * 100)
patients <- unique(meta$patient)
calc_pcs <- function(obj){
  pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- if (length(pct) > 1) sort(which(diff(pct) > 0.1), decreasing = TRUE)[1] + 1 else NA
  list(pcs = ifelse(is.na(co2), co1, min(co1, co2)))
}
pca <- calc_pcs(cancer.epi)$pcs
cancer.epi <- FindNeighbors(cancer.epi, dims = 1:pca)
cancer.epi <- FindClusters(cancer.epi, resolution = seq(0.01, 2, length.out = 15))
cancer.epi <- RunUMAP(cancer.epi, dims = 1:pca)
saveRDS(cancer.epi, "cancerepi_umap.rds")
j <- readRDS("cancerepi_umap.rds")
all_DGE <- do.call(rbind, lapply(colnames(j@meta.data)[c(9,12:19,22:28)], function(i){
  Idents(j) <- j@meta.data[, i]
  mk <- FindAllMarkers(j, only.pos = TRUE, min.cells.group = 50, min.diff.pct = 0.25, logfc.threshold = 0.25)
  mk$cluster_res <- paste0(mk$cluster, "_", strsplit(i, "res.")[[1]][2])
  mk
}))
all <- all_DGE %>% filter(p_val_adj < 0.05)
write.csv(all, "cancer_epi_DGE_unsupervised.csv")
deg <- read.csv("cancer_epi_DGE_unsupervised.csv", row.names = 1)
deg <- deg %>% filter(p_val_adj < 0.05, avg_log2FC > 0, pct.1 > 0.2) %>% group_by(cluster_res) %>% top_n(-200, p_val_adj) %>% top_n(200, avg_log2FC)
sig <- data.frame(gene = unique(deg$gene), row.names = unique(deg$gene))
for (i in unique(deg$cluster_res)) {
  tmp <- data.frame(gene = deg$gene[deg$cluster_res == i], v = 1)
  colnames(tmp) <- c("gene", i)
  rownames(tmp) <- tmp$gene
  sig <- left_join(sig, tmp, by = "gene")
}
sig[is.na(sig)] <- 0
sig <- sig[, -1]
counts <- deg %>% count(cluster_res)
drop <- which(colSums(sig) < 20)
sig <- sig[, -drop]
saveRDS(sig, "jaccard.rds")
jacc <- as.matrix(readRDS("jaccard.rds"))
drop <- c()
for (i in 1:ncol(jacc)){
  for (j in (i+1):ncol(jacc)){
    if (!is.na(jacc[i,j]) && jacc[i,j] > 0.95) drop <- c(drop, j)
  }
}
jacc <- jacc[-unique(drop), -unique(drop)]
saveRDS(jacc, "jaccard_noredundant.rds")
set.seed(123)
rh <- consensus_partition(jacc, top_value_method = "ATC", partition_method = "skmeans", p_sampling = 0.8, max_k = 13)
k <- 13
write.csv(get_classes(rh, k = k), "GE_classes.csv")
GE <- read.csv("GE_classes.csv")[,1:2]
colnames(GE) <- c("cluster_res", "GE")
deg_GE <- left_join(deg, GE, by = "cluster_res") %>% filter(!is.na(GE))
GE_gene <- deg_GE %>% group_by(GE, gene) %>% summarise(Nb = n(), .groups = "drop") %>% mutate(C = sum(Nb), percent = Nb / C * 100)
GE_gene <- GE_gene %>% group_by(GE) %>% top_n(350, Nb) %>% top_n(200, avg_log2FC)
write.csv(GE_gene, "cancer_GE.csv")
GElist <- split(GE_gene$gene, GE_gene$GE)
cancer.epi <- AddModuleScore_UCell(cancer.epi, features = GElist, assay = "RNA")
saveRDS(cancer.epi, "cancerepi_withGE.rds")
#GE heatmap
set.seed(123)
setwd("./t100553/wss/copykat/cancer state/")
# Unsupervised DGE generation
cancer.epi <- readRDS("./cancerepi_withGEs score.rds")
expdata <- t(cancer.epi@meta.data[,c(4,5,11,29:41)])
collapse_expdata <- as.data.frame(rownames(expdata))
for (i in patient) {
  subset <- expdata[,which(expdata[3,] == i)]
  subset <- subset[,sample.int(dim(subset)[2],min(dim(subset)[2],20000))]
  collapse_expdata <- cbind(collapse_expdata, subset)
}
labels <- rownames(collapse_expdata)[4:18]
collapse_zscore <- collapse_expdata[-c(1:3),]
collapse_zscore <- collapse_zscore[,-1]
collapse_zscore <- as.matrix(sapply(collapse_zscore, as.numeric))
collapse_zscore <- t(apply(collapse_zscore, 1, function(x) (x-mean(x))/sd(x)))
sort <- rbind(apply(collapse_zscore, 2, function(x) which.max(x)),
              apply(collapse_zscore, 2, function(x) max(x)), 
              collapse_zscore)
sort <- sort[,order(sort[2,],decreasing = T)] # sort by max z-score
sort <- sort[,order(sort[1,],decreasing = F)] #sort by GE
collapse_zscore <- sort[-c(1,2),]
collapse_expdata <- collapse_expdata[,colnames(collapse_zscore)]
library(pheatmap)
mutation_anno <- t(as.matrix(collapse_expdata[1,]))
colnames(mutation_anno) <- c("mutation")
dataset_anno <- t(as.matrix(collapse_expdata[2,]))
colnames(dataset_anno) <- c("dataset")
colors <- pal_npg("nrc")(7) 
collapse_zscore <- t(collapse_zscore)
anno_df <- data.frame(mutation = mutation_anno,
                      dataset = dataset_anno)
annotation_colors <- list(
  dataset = c(
    "GSE171145" = "#E64B35FF",
    "GSE202371" = "#4DBBD5FF",
    "GSE131907" = "#00A087FF",
    "CO.0121060.V1" = "#3C5488FF",
    "PMID35027529" = "#F39B7FFF",
    "GSE148071" = "#8491B4FF" ,
    "GSE136246" = "#91D1C2FF"
  ),
  mutation = c(
    "EGFR" = "#c6b7d4",
    "EGFR-BM" = "#d44e26",
    "EGFR-co-mutation" = "#e3a264",
    "KRAS" = "#6fc2d0",
    "KRAS-co-mutation" = "#6f9abf",
    "ALK" = "#a5c49b",
    "ROS1" = "#7266ac",
    "TP53" = "#FF9966",
    "MET-BM" = "#d84986",
    "HER2" = "#2d588e"
  )
)
pheatmap(t(collapse_zscore),cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = F,annotation_col  = anno_df,annotation_colors = annotation_colors,color = colorRampPalette(c("#1389CC","white","#FAB476","#F6885B","#d66c61","#d46258","#bd2b2e","#bc272c","#B8222E"))(1000))#,
pheatmap(t(collapse_zscore),cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = F,annotation_col  = anno_df,annotation_colors = annotation_colors,color = colorRampPalette(c("#A5CCE0","white","#F6C7A2","#FAB476","#F39869","#F6885B","#d66c61","#C84E4C","#D04B2B","#bd2b2e","#bc272c","#B8222E"))(1000))#,


#Define similarity function
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))  # Intersection size
  union <- length(union(set1, set2))            # Merge size
  return(intersection / union)                  # Similarity calculation
}

library(msigdbr)
msigdbr_species() #List some species

#Select gene set
human_KEGG = msigdbr(species = "Homo sapiens",category = "H" ) %>% dplyr::select(gs_name,gene_symbol)

H.geneset = human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)

# Initialize similarity matrix
num_sets_A <- length(GElist)
num_sets_B <- length(H.geneset)
similarity_matrix <- matrix(0, nrow = num_sets_A, ncol = num_sets_B)

# Calculate the Jaccard similarity between all gene set pairs
for (i in 1:num_sets_A) {
  for (j in 1:num_sets_B) {
    similarity_matrix[i, j] <- jaccard_similarity(GElist[[i]], H.geneset[[j]])
  }
}

# Set row and column names
rownames(similarity_matrix) <- names(GElist)
colnames(similarity_matrix) <- names(H.geneset)
rownames(similarity_matrix) <- paste0("GE",rownames(similarity_matrix))

# Output similarity matrix
print(similarity_matrix)
similarity_matrix <- data.frame(t(similarity_matrix))
rownames(similarity_matrix) <- gsub("HALLMARK_","",rownames(similarity_matrix))
pheatmap(similarity_matrix,color = colorRampPalette(c("white","#fdedea","#fbdedb","#f0b2aa","#dd7e74","#d66c61","#d46258","#bd2b2e","#bc272c"))(1000))

write_xlsx(similarity_matrix,
           path = "GE_H.geneset_jaccard.xlsx", 
           col_names = TRUE, 
           format_headers = TRUE)

library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

converted_gene_sets <- lapply(GElist, function(x) bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID)

go_results <- lapply(converted_gene_sets, function(x)
  enrichGO(x, OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod="BH"))

kegg_results <- lapply(converted_gene_sets, function(x)
  enrichKEGG(x, organism="hsa", pvalueCutoff=0.05, qvalueCutoff=0.05, pAdjustMethod="BH"))

save(go_results, kegg_results, file="go_kegg.RData")

extract_top <- function(res, n=10, type="GO"){
  df <- as.data.frame(res)
  if(type=="KEGG"){
    df <- df %>% arrange(qvalue) %>% head(n)
    df$term <- df$subcategory
  } else {
    df <- df %>% arrange(qvalue) %>% head(n)
    df$term <- df$Description
  }
  df
}

plot_func <- function(df, title, ylab){
  ggplot(df, aes(x=Count, y=reorder(term, Count), fill=-log10(qvalue))) +
    geom_bar(stat="identity", width=0.8) +
    scale_fill_distiller(palette="Blues", direction=1) +
    labs(x="Number of Gene", y=ylab, title=title) +
    theme_bw()
}

plot_list <- lapply(seq_along(go_results), function(i){
  go_df <- extract_top(go_results[[i]], type="GO")
  kegg_df <- extract_top(kegg_results[[i]], type="KEGG")
  list(
    GO = plot_func(go_df, paste0("go.GE", i), "GO"),
    KEGG = plot_func(kegg_df, paste0("kegg.GE", i), "KEGG")
  )
})

go_all <- bind_rows(lapply(seq_along(go_results), function(i){
  df <- extract_top(go_results[[i]], type="GO")
  df$GE <- paste0("GE", i)
  df
}))

write.csv(go_all, "./yang/wang/GO.csv", row.names=FALSE)

#Enrich various GE functions using the ClusterProfiler package
library(ggplot2)
library(clusterProfiler)
data(gcSample)
xx.go <- compareCluster(converted_gene_sets,fun="enrichGO", OrgDb="org.Hs.eg.db")
dotplot(xx.go, showCategory=5, includeAll=FALSE) 

xx.kegg <- compareCluster(
  geneClusters = converted_gene_sets,
  fun = "enrichKEGG",
  organism = "hsa"  
)
dotplot(xx.kegg, showCategory=5, includeAll=FALSE) 
save(xx.go,files = "./go.RData")
#----3、Calculate Jaccard correlation between GE module and NMF module of TIME cells
GElist <- readxl::read_xlsx("cancer_GE.xlsx")
GElist$GE <- paste0("GE",GElist$GE)
GElist <- split(GElist$gene,GElist$GE)
TIME.NMF <- read.csv("./TIME_NMF_markers2.csv")
TIME.NMFlist <- split(TIME.NMF$gene,TIME.NMF$cluster)

#Define similarity function
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))  
  union <- length(union(set1, set2))            
  return(intersection / union)               
}

# Initialize similarity matrix
num_sets_A <- length(GElist)
num_sets_B <- length(TIME.NMFlist)
similarity_matrix <- matrix(0, nrow = num_sets_A, ncol = num_sets_B)

# Calculate the Jaccard similarity between all gene set pairs
for (i in 1:num_sets_A) {
  for (j in 1:num_sets_B) {
    similarity_matrix[i, j] <- jaccard_similarity(GElist[[i]], TIME.NMFlist[[j]])
  }
}

rownames(similarity_matrix) <- names(GElist)
colnames(similarity_matrix) <- names(TIME.NMFlist)

print(similarity_matrix)
similarity_matrix <- data.frame(t(similarity_matrix))
pheatmap(similarity_matrix,color = colorRampPalette(c("#6c94c7","white","#fdedea","#fbdedb","#f0b2aa","#dd7e74","#d66c61","#d46258","#bd2b2e","#bc272c"))(1000))

write_xlsx(similarity_matrix,
           path = "GE_TIME.NMFlist_jaccard.xlsx", 
           col_names = TRUE, 
           format_headers = TRUE)

#Determine the cellular status to which each cancer epithelial cell belongs
setwd("./t100553/wss/copykat/cancer state/")
cancer.epi <- readRDS("./cancerepi_withGEs score.rds")
library(tidyverse)
cancer.epi@meta.data <- cancer.epi@meta.data %>% mutate(mutation = case_when(
  patient %in% c("P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11") ~ "EGFR",
  patient %in% c("P12","P13","P14","P15","P16","P17","P18","P19","P20","P21") ~ "EGFR-BM",
  patient %in% c("P22","P23") ~ "EGFR-co-mutation",
  patient %in% c("P24","P25","P26","P27","P28","P29","P30","P31","P32","P33","P34","P35") ~ "KRAS",
  patient %in% c("P36","P37","P38") ~ "KRAS-co-mutation",
  patient %in% c("P39") ~ "ALK",
  patient %in% c("P40","P41","P42") ~ "ROS1",
  patient %in% c("P43") ~ "TP53",
  patient %in% c("P44") ~ "MET-BM",
  patient %in% c("P45") ~ "HER2"))
df <- cancer.epi@meta.data[,c(29:41)]
colnames(df) <- gsub("_UCell","",colnames(df))
colnames(df) <- paste0("GE",colnames(df))
ge_cols <- grep("^GE", colnames(df), value = TRUE)  # Extract GE column names
df$type <- apply(df[, ge_cols], 1, function(row) names(which.max(row)))
cancer.epi$GEtype <- df$type
saveRDS(cancer.epi, file = "./wss/copykat/cancer state/cancer.epi_GE.rds")

#----Figure 3B----
# plot UMAP by mutation
cancer.epi <- readRDS("./wss/copykat/cancer state/cancer.epi_GE.rds")
cancer.epi$GEtype <- factor(cancer.epi$GEtype,levels=c("GE1","GE2","GE3","GE4","GE5","GE6","GE7","GE8","GE9","GE10","GE11","GE12","GE13"))
library(ggsci)
colors_10 <- pal_npg("nrc")(10) 
colors_13 <- colorRampPalette(colors_10)(13)
colors_13 <- c("#0073C2FF", "#438F49" , "#26ADAE", "#0F8C87" ,"#D25515","#DF86BC","#9887C8" ,"#699DC8","#91D1C2" ,"#C93430" ,"#DB873F" ,"#C9B8D1","#EFC000FF")
p1 <- DimPlot(cancer.epi, reduction = "umap", 
              label = F, 
              repel = TRUE, 
              raster = FALSE, 
              group.by = "GEtype",
              cols = colors_13) + SeuratAxes() 
p1
ggsave("./NSCLC driver genes/result diagram/figure3/copykat/cancer state/cancerepi_UMAP_GE.pdf", plot = p1, width = 6, height = 5.5)
cancer.epi$mutation <- factor(cancer.epi$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53 ","MET-BM","HER2"))

mycolors <- c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#d84986","#7266ac","#FF9966" )
p2 <- DimPlot(cancer.epi, reduction = "umap", 
              label = F, 
              repel = TRUE, 
              raster = FALSE, 
              group.by = "mutation",
              cols = mycolors) + SeuratAxes() 

p2
ggsave("E:/NSCLC driver genes/result diagram/figure3/copykat/cancer state/cancerepi_UMAP_mutation.pdf", plot = p2, width = 7, height = 5.5)

#----Figure 3C----
#Using the GE gene set to score each TIME and explore the association between oncogenes and TIME
corsce <- readRDS("./t100553/wss/subsce.corr.rds")
GElist <- read.csv("cancer_GE.order.csv")
GElist <- split(GElist$gene,GElist$GE)
length(GElist) 

library(UCell)
corsce <- AddModuleScore_UCell(corsce,
                               features = GElist,
                               assay = "RNA")
corsce$GE <- apply(corsce@meta.data, 1, function(row) names(corsce@meta.data)[which.max(row)])

#Draw a box line diagram
library(tidyverse)
library(rstatix)   
library(ggpubr)     
library(ggplot2)
df <- na.omit(corsce@meta.data[, c("GE1_UCell", "cellmodule")])

mycolors <- c("#70cdbe", "#80C1D7", "#968ab7", "#d586b3", "#eb7e60")
names(mycolors) <- unique(df$cellmodule)   

kruskal_p <- kruskal.test(GE1_UCell ~ cellmodule, data = df)$p.value
kruskal_label <- paste0("Kruskal-Wallis P = ", signif(kruskal_p, 3))


df_outliers_sampled <- df %>%
  group_by(cellmodule) %>%
  mutate(
    Q1 = quantile(GE1_UCell, 0.25),
    Q3 = quantile(GE1_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE1_UCell < (Q1 - 1.5 * IQR) | GE1_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()


ggplot(df, aes(x = cellmodule, y = GE1_UCell, color = cellmodule)) +
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
  ) + stat_compare_means(method = "anova", #statistical method
                         aes(label = "p.format"), size = 5) #size


#GE2
ggplot(corsce@meta.data, aes(x = cellmodule, y = GE2_UCell, fill = cellmodule)) +
  geom_boxplot(show.legend = FALSE,colour = "white", outlier.color = "black") +
  stat_compare_means(
    method = "wilcox.test",  # Non parametric testing (more robust)
    comparisons = comparisons,  # Automatically traverse all combinations
    label = "p.signif",      # display *** / ** / *
    hide.ns = TRUE,          # hide ns
    tip_length = 0.01,
    step.increase = 0.1      # Increase comment line spacing
  ) +
  scale_fill_manual(values = mycolors) +  # five colors
  labs(title = "GE2 Activity by module", x = NULL, y = "Score") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    legend.position = "right"
  )


#GE3
ggplot(corsce@meta.data, aes(x = cellmodule, y = GE3_UCell, fill = cellmodule)) +
  geom_boxplot(show.legend = FALSE,colour = "white", outlier.color = "black") +
  stat_compare_means(
    method = "wilcox.test",  
    comparisons = comparisons,  
    label = "p.signif",      
    hide.ns = TRUE,          
    tip_length = 0.01,
    step.increase = 0.1      
  ) +
  scale_fill_manual(values = mycolors) +   
  labs(title = "GE3 Activity by module", x = NULL, y = "Score") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    legend.position = "right"
  )


#GE4
ggplot(corsce@meta.data, aes(x = cellmodule, y = GE4_UCell, fill = cellmodule)) +
  geom_boxplot(show.legend = FALSE,colour = "white", outlier.color = "black") +
  stat_compare_means(
    method = "wilcox.test",  
    comparisons = comparisons,  
    label = "p.signif",      
    hide.ns = TRUE,          
    tip_length = 0.01,
    step.increase = 0.1      
  ) +
  scale_fill_manual(values = mycolors) +  
  labs(title = "GE4 Activity by module", x = NULL, y = "Score") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    legend.position = "right"
  )


#GE5
 
df <- na.omit(corsce@meta.data[, c("GE5_UCell", "cellmodule")])

 
mycolors <- c("#70cdbe", "#80C1D7", "#968ab7", "#d586b3", "#eb7e60")
names(mycolors) <- unique(df$cellmodule)   

 
kruskal_p <- kruskal.test(GE5_UCell ~ cellmodule, data = df)$p.value
kruskal_label <- paste0("Kruskal-Wallis P = ", signif(kruskal_p, 3))

# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(cellmodule) %>%
  mutate(
    Q1 = quantile(GE5_UCell, 0.25),
    Q3 = quantile(GE5_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE5_UCell < (Q1 - 1.5 * IQR) | GE5_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = cellmodule, y = GE5_UCell, color = cellmodule)) +
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
ggplot(corsce@meta.data, aes(x = cellmodule, y = GE6_UCell, fill = cellmodule)) +
  geom_boxplot(show.legend = FALSE,colour = "white", outlier.color = "black") +
  stat_compare_means(
    method = "wilcox.test",  
    comparisons = comparisons,  
    label = "p.signif",      
    hide.ns = TRUE,          
    tip_length = 0.01,
    step.increase = 0.1      
  ) +
  scale_fill_manual(values = mycolors) + 
  labs(title = "GE6 Activity by module", x = NULL, y = "Score") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    legend.position = "right"
  )


#GE7
ggplot(corsce@meta.data, aes(x = cellmodule, y = GE7_UCell, fill = cellmodule)) +
  geom_boxplot(show.legend = FALSE,colour = "white", outlier.color = "black") +
  stat_compare_means(
    method = "wilcox.test",  
    comparisons = comparisons,  
    label = "p.signif",      
    hide.ns = TRUE,          
    tip_length = 0.01,
    step.increase = 0.1      
  ) +
  scale_fill_manual(values = mycolors) +  
  labs(title = "GE7 Activity by module", x = NULL, y = "Score") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    legend.position = "right"
  )


#GE8
 
df <- na.omit(corsce@meta.data[, c("GE8_UCell", "cellmodule")])
mycolors <- c("#70cdbe", "#80C1D7", "#968ab7", "#d586b3", "#eb7e60")
names(mycolors) <- unique(df$cellmodule)   

kruskal_p <- kruskal.test(GE8_UCell ~ cellmodule, data = df)$p.value
kruskal_label <- paste0("Kruskal-Wallis P = ", signif(kruskal_p, 3))
# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(cellmodule) %>%
  mutate(
    Q1 = quantile(GE8_UCell, 0.25),
    Q3 = quantile(GE8_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE8_UCell < (Q1 - 1.5 * IQR) | GE8_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = cellmodule, y = GE8_UCell, color = cellmodule)) +
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
 
df <- na.omit(corsce@meta.data[, c("GE9_UCell", "cellmodule")])

mycolors <- c("#70cdbe", "#80C1D7", "#968ab7", "#d586b3", "#eb7e60")
names(mycolors) <- unique(df$cellmodule)  # 保证名字匹配

kruskal_p <- kruskal.test(GE9_UCell ~ cellmodule, data = df)$p.value
kruskal_label <- paste0("Kruskal-Wallis P = ", signif(kruskal_p, 3))
# Calculate outliers and draw a maximum of 100 for each group
df_outliers_sampled <- df %>%
  group_by(cellmodule) %>%
  mutate(
    Q1 = quantile(GE9_UCell, 0.25),
    Q3 = quantile(GE9_UCell, 0.75),
    IQR = Q3 - Q1,
    is_outlier = GE9_UCell < (Q1 - 1.5 * IQR) | GE9_UCell > (Q3 + 1.5 * IQR)
  ) %>%
  filter(is_outlier) %>%
  slice_head(n = 100) %>%
  ungroup()

 
ggplot(df, aes(x = cellmodule, y = GE9_UCell, color = cellmodule)) +
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
ggplot(corsce@meta.data, aes(x = cellmodule, y = GE10_UCell, fill = cellmodule)) +
  geom_boxplot(show.legend = FALSE,colour = "white", outlier.color = "black") +
  stat_compare_means(
    method = "wilcox.test",   
    comparisons = comparisons,   
    label = "p.signif",       
    hide.ns = TRUE,           
    tip_length = 0.01,
    step.increase = 0.1       
  ) +
  scale_fill_manual(values = mycolors) +   
  labs(title = "GE10 Activity by module", x = NULL, y = "Score") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    legend.position = "right"
  )


#GE11
ggplot(corsce@meta.data, aes(x = cellmodule, y = GE11_UCell, fill = cellmodule)) +
  geom_boxplot(show.legend = FALSE,colour = "white", outlier.color = "black") +
  stat_compare_means(
    method = "wilcox.test",   
    comparisons = comparisons,   
    label = "p.signif",       
    hide.ns = TRUE,           
    tip_length = 0.01,
    step.increase = 0.1       
  ) +
  scale_fill_manual(values = mycolors) +   
  labs(title = "GE11 Activity by module", x = NULL, y = "Score") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    legend.position = "right"
  )


#GE12
ggplot(corsce@meta.data, aes(x = cellmodule, y = GE12_UCell, fill = cellmodule)) +
  geom_boxplot(show.legend = FALSE,colour = "white", outlier.color = "black") +
  stat_compare_means(
    method = "wilcox.test",   
    comparisons = comparisons,   
    label = "p.signif",       
    hide.ns = TRUE,           
    tip_length = 0.01,
    step.increase = 0.1       
  ) +
  scale_fill_manual(values = mycolors) +   
  labs(title = "GE12 Activity by module", x = NULL, y = "Score") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    legend.position = "right"
  )


#GE13
ggplot(corsce@meta.data, aes(x = cellmodule, y = GE13_UCell, fill = cellmodule)) +
  geom_boxplot(show.legend = FALSE,colour = "white", outlier.color = "black") +
  stat_compare_means(
    method = "wilcox.test",   
    comparisons = comparisons,   
    label = "p.signif",       
    hide.ns = TRUE,           
    tip_length = 0.01,
    step.increase = 0.1       
  ) +
  scale_fill_manual(values = mycolors) +   
  labs(title = "GE13 Activity by module", x = NULL, y = "Score") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_text(size = 11),
    legend.position = "right"
  )
dir.create("./GE_TIME")

#----Figure 3D----
#Exploring the CNVs at the gene level of various cancer cell states based on the results of Infercnv
library(infercnv)
library(dplyr)
library(tidyverse)
library(circlize)
library(stringr)
setwd("./t100553/wss/Infercnv.DC.NK/")
cancer.epi <- readRDS("./wss/copykat/cancer state/cancer.epi_GE.rds")
cancer.epi$GEtype <- factor(cancer.epi$GEtype,levels=c("GE1","GE2","GE3","GE4","GE5","GE6","GE7","GE8","GE9","GE10","GE11","GE12","GE13"))
infercnv_obj <- readRDS("./t100553/wss/Infercnv.DC.NK/run.final.infercnv_obj")
group <- read.table("./t100553/wss/Infercnv.DC.NK/infercnv.observation_groupings.txt", header=T)
cancer.epi$cellid <- rownames(cancer.epi@meta.data)
cancerepi <- cancer.epi@meta.data[,c(43,42)]
table(cancerepi$GEtype)
#GE1
GE1 <- subset(cancerepi,GEtype =="GE1")
expr <- infercnv_obj@expr.data
expr.cancer <- expr[,colnames(expr) %in% GE1$cellid]

expr2=expr.cancer-1

exp_cnv1 <- data.frame(rowMeans(expr2))
exp_cnv1$name <- rownames(exp_cnv1)
gene <- infercnv_obj@gene_order
gene$name <- rownames(gene)
gene1 <- gene[match(exp_cnv1$name,gene$name),]
gene1$name <- NULL
gene1$cnv <- exp_cnv1$rowMeans.expr2.
gene1$name <- rownames(gene1)

gene_up  <- gene1 %>%  top_n(n = 50, wt = cnv)
gene_down  <- gene1 %>%  top_n(n = -50, wt = cnv)

gene2 <- rbind(gene_up,gene_down)

gene2 <- rbind(gene2,cnv)
write.csv(gene2,file = "./t100553/wss/copykat/CNV gene/GE1.csv")
cnv <- subset(gene1,rownames(gene1) %in% c("EGFR","MET","KRAS","TP53","BRAF","ALK","ROS1","HER2"))
gene_up  <- gene1 %>%  top_n(n = 5, wt = cnv)
gene_down  <- gene1 %>%  top_n(n = -5, wt = cnv)

gene2 <- rbind(gene_up,gene_down)

gene2 <- rbind(gene2,cnv)

library(circlize)
circos.clear() 
circos.par("start.degree" = 90)
circos.par("gap.degree" = rep(c(2, 2), 12), ADD = TRUE)
circos.initializeWithIdeogram(species = "hg38", plotType = c("axis", "labels"))

circos.initializeWithIdeogram(species = "hg38")
circos.genomicTrack(gene2, 
                    track.height = 0.15,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 1, col = 1)})

circos.genomicTrack(gene2, 
                    track.height = 0.08,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
                                         col = ifelse(value[[1]] > 0, "red", "green"),...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })



circos.genomicLabels(gene2, labels.column = 5, side = "inside",cex = 2.1)


#GE5
GE5 <- subset(cancerepi,GEtype =="GE5")
expr <- infercnv_obj@expr.data
expr.cancer <- expr[,colnames(expr) %in% GE5$cellid]


expr2=expr.cancer-1
exp_cnv1 <- data.frame(rowMeans(expr2))
exp_cnv1$name <- rownames(exp_cnv1)
gene <- infercnv_obj@gene_order
gene$name <- rownames(gene)
gene1 <- gene[match(exp_cnv1$name,gene$name),]
gene1$name <- NULL
gene1$cnv <- exp_cnv1$rowMeans.expr2.
gene1$name <- rownames(gene1)
cnv <- subset(gene1,rownames(gene1) %in% c("EGFR","MET","KRAS","TP53","BRAF","ALK","ROS1","HER2"))
gene_up  <- gene1 %>%  top_n(n = 50, wt = cnv)
gene_down  <- gene1 %>%  top_n(n = -50, wt = cnv)

gene2 <- rbind(gene_up,gene_down)

gene2 <- rbind(gene2,cnv)
write.csv(gene2,file = "./t100553/wss/copykat/CNV gene/GE5.csv")

gene_up  <- gene1 %>%  top_n(n = 5, wt = cnv)
gene_down  <- gene1 %>%  top_n(n = -5, wt = cnv)

gene2 <- rbind(gene_up,gene_down)
gene2 <- rbind(gene2,cnv)
 
 
library(circlize)
circos.clear() 
circos.par("start.degree" = 90)
circos.par("gap.degree" = rep(c(2, 2), 12), ADD = TRUE)
circos.initializeWithIdeogram(species = "hg38", plotType = c("axis", "labels"))
 
circos.initializeWithIdeogram(species = "hg38")
circos.genomicTrack(gene2, 
                    track.height = 0.15,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 1, col = 1)})

circos.genomicTrack(gene2, 
                    track.height = 0.08,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
                                         col = ifelse(value[[1]] > 0, "red", "green"),...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })



circos.genomicLabels(gene2, labels.column = 5, side = "inside",cex = 2.1)


#GE8
GE8 <- subset(cancerepi,GEtype =="GE8")
expr <- infercnv_obj@expr.data
expr.cancer <- expr[,colnames(expr) %in% GE8$cellid]


expr2=expr.cancer-1
 
exp_cnv1 <- data.frame(rowMeans(expr2))
exp_cnv1$name <- rownames(exp_cnv1)
gene <- infercnv_obj@gene_order
gene$name <- rownames(gene)
gene1 <- gene[match(exp_cnv1$name,gene$name),]
gene1$name <- NULL
gene1$cnv <- exp_cnv1$rowMeans.expr2.
gene1$name <- rownames(gene1)
 
cnv <- subset(gene1,rownames(gene1) %in% c("EGFR","MET","KRAS","TP53","BRAF","ALK","ROS1","HER2"))
gene_up  <- gene1 %>%  top_n(n = 50, wt = cnv)
gene_down  <- gene1 %>%  top_n(n = -50, wt = cnv)

gene2 <- rbind(gene_up,gene_down)

gene2 <- rbind(gene2,cnv)
write.csv(gene2,file = "./t100553/wss/copykat/CNV gene/GE8.csv")

gene_up  <- gene1 %>%  top_n(n = 5, wt = cnv)
gene_down  <- gene1 %>%  top_n(n = -5, wt = cnv)

gene2 <- rbind(gene_up,gene_down)
gene2 <- rbind(gene2,cnv)
 
 
library(circlize)
circos.clear()  
circos.par("start.degree" = 90)
circos.par("gap.degree" = rep(c(2, 2), 12), ADD = TRUE)
circos.initializeWithIdeogram(species = "hg38", plotType = c("axis", "labels"))
 
circos.initializeWithIdeogram(species = "hg38")
circos.genomicTrack(gene2, 
                    track.height = 0.15,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 1, col = 1)})

circos.genomicTrack(gene2, 
                    track.height = 0.08,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
                                         col = ifelse(value[[1]] > 0, "red", "green"),...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })



circos.genomicLabels(gene2, labels.column = 5, side = "inside",cex = 2.1)



#GE9
GE9 <- subset(cancerepi,GEtype =="GE9")
expr <- infercnv_obj@expr.data
expr.cancer <- expr[,colnames(expr) %in% GE9$cellid]
expr2=expr.cancer-1

exp_cnv1 <- data.frame(rowMeans(expr2))
exp_cnv1$name <- rownames(exp_cnv1)
gene <- infercnv_obj@gene_order
gene$name <- rownames(gene)
gene1 <- gene[match(exp_cnv1$name,gene$name),]
gene1$name <- NULL
gene1$cnv <- exp_cnv1$rowMeans.expr2.
gene1$name <- rownames(gene1)
 
cnv <- subset(gene1,rownames(gene1) %in% c("EGFR","MET","KRAS","TP53","BRAF","ALK","ROS1","HER2"))
gene_up  <- gene1 %>%  top_n(n = 50, wt = cnv)
gene_down  <- gene1 %>%  top_n(n = -50, wt = cnv)

gene2 <- rbind(gene_up,gene_down)

gene2 <- rbind(gene2,cnv)
write.csv(gene2,file = "./t100553/wss/copykat/CNV gene/GE9.csv")


gene_up  <- gene1 %>%  top_n(n = 5, wt = cnv)
gene_down  <- gene1 %>%  top_n(n = -5, wt = cnv)

gene2 <- rbind(gene_up,gene_down)
gene2 <- rbind(gene2,cnv)
 
 
library(circlize)
circos.clear()  
circos.par("start.degree" = 90)
circos.par("gap.degree" = rep(c(2, 2), 12), ADD = TRUE)
circos.initializeWithIdeogram(species = "hg38", plotType = c("axis", "labels"))
 
circos.initializeWithIdeogram(species = "hg38")
circos.genomicTrack(gene2, 
                    track.height = 0.15,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, pch = 16, cex = 1, col = 1)})

circos.genomicTrack(gene2, 
                    track.height = 0.08,
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
                                         col = ifelse(value[[1]] > 0, "red", "green"),...)
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                    })



circos.genomicLabels(gene2, labels.column = 5, side = "inside",cex = 2.1)


# Save all GE's top CNVs
all_ge_data <- list()

# Loop processing each GE group (GE1-GE13)
for (ge_num in 1:13) {
  ge_group <- paste0("GE", ge_num)
  
  GE_subset <- subset(cancerepi, GEtype == ge_group)
  
  expr <- infercnv_obj@expr.data
  expr.cancer <- expr[, colnames(expr) %in% GE_subset$cellid, drop = FALSE]
  
  expr_centered <- expr.cancer - 1
  gene_means <- rowMeans(expr_centered)
  
  gene_info <- infercnv_obj@gene_order
  cnv_data <- data.frame(
    Gene = rownames(gene_info),  
    gene_info,
    cnv_score = gene_means[match(rownames(gene_info), names(gene_means))]
  )
  rownames(cnv_data) <- NULL  
  
  top_gain <- cnv_data %>% 
    arrange(desc(cnv_score)) %>% 
    head(50) %>% 
    mutate(direction = "gain")
  
  top_loss <- cnv_data %>% 
    arrange(cnv_score) %>% 
    head(50) %>% 
    mutate(direction = "loss")
  
  ge_result <- rbind(top_gain, top_loss) %>% 
    mutate(GE_group = ge_group)  
  
  all_ge_data[[ge_group]] <- ge_result
}

final_combined_data <- do.call(rbind, all_ge_data) %>%
  select(Gene, everything())

write.csv(final_combined_data, 
          file = "GE1-GE13_combined_cnv_results.csv",
          row.names = FALSE)  

glimpse(final_combined_data)

#----Figure 3E----
#Survival analysis
#Perform survival analysis on 5 cell modules and identify modules with poor prognosis
library(survival)
library(ggpubr)
library(survminer)
library(tidyverse)
library(Seurat)
setwd("./t100553/wss")
dir.create("survival analysis")
setwd("./t100553/wss/survival analysis")
#Read in the processed clinical data and RNA seq data of LUAD and LUSC
clinical_LUAD <- read.csv(file = "./t100553/wss/survival analysis/Bulk data for LUAD and LUSC/clinical_LUAD.csv",row.names = 1)
clinical_LUSC <- read.csv(file = "./t100553/wss/survival analysis/Bulk data for LUAD and LUSC/clinical_LUSC.csv",row.names = 1)
exprset_LUAD <- read.csv(file = "./t100553/wss/survival analysis/Bulk data for LUAD and LUSC/exprset_LUAD.csv",row.names = 1)
exprset_LUSC <- read.csv(file = "./t100553/wss/survival analysis/Bulk data for LUAD and LUSC/exprset_LUSC.csv",row.names = 1)
clinical_LUAD$event=ifelse(clinical_LUAD$vital_status=='Alive',0,1)
clinical_LUSC$event=ifelse(clinical_LUSC$vital_status=='Alive',0,1)

#Merge data from LUAD and LUSC
clinical <- rbind(clinical_LUAD,clinical_LUSC)
exprset <- rbind(exprset_LUAD,exprset_LUSC)

clinical$time=round(clinical$os/30,2)
clinical <- clinical[clinical$time<36,]
clinical$time <- as.numeric(clinical$time)
exprset <- exprset[match(rownames(clinical),substr(rownames(exprset),1,12)),]

set.seed(111)
ind <- sample(nrow(exprset), nrow(exprset) * 0.5)
train.tpm <- exprset[ind, ]
test.tpm <- exprset[-ind, ]
train.tpm <- data.frame(t(train.tpm))
test.tpm <- data.frame(t(test.tpm))
train.meta <- clinical[ind, ]
test.meta <- clinical[-ind, ]
expr <- exprset
clin <- clinical
#Calculate the set of marker genes for each immune microenvironment subtype

#Read in subse.corr data
#Survival analysis of the first data
setwd("./t100553/wss/survival analysis")
exprset <- train.tpm
clinical <- train.meta

subsce.corr.markers1.up.csv <- read.csv(file = "./subsce.corr.markers1.up.csv",row.names = 1)
module1.marker <- subset(subsce.corr.markers1.up.csv,cluster=="module1")
module2.marker <- subset(subsce.corr.markers1.up.csv,cluster=="module2")
module3.marker <- subset(subsce.corr.markers1.up.csv,cluster=="module3")
module4.marker <- subset(subsce.corr.markers1.up.csv,cluster=="module4")
module5.marker <- subset(subsce.corr.markers1.up.csv,cluster=="module5")

module1.marker <- module1.marker %>% top_n(n = 30, wt = avg_log2FC)
module1.marker <- as.list(module1.marker$gene)
module2.marker <- module2.marker %>% top_n(n = 30, wt = avg_log2FC)
module2.marker <- as.list(module2.marker$gene)
module3.marker <- module3.marker %>% top_n(n = 30, wt = avg_log2FC)
module3.marker <- as.list(module3.marker$gene)
module4.marker <- module4.marker %>% top_n(n = 30, wt = avg_log2FC)
module4.marker <- as.list(module4.marker$gene)
module5.marker <- module5.marker %>% top_n(n = 30, wt = avg_log2FC)
module5.marker <- as.list(module5.marker$gene)

#Score the gene set, calculate the z-score of each patient for each immune microenvironment in the bulk data, and determine the dominant cell module
library(scales)
exp.1 <-  data.frame((colMeans(exprset[rownames(exprset) %in% module1.marker,])))
exp.2 <-  data.frame((colMeans(exprset[rownames(exprset) %in% module2.marker,])))
exp.3 <-  data.frame((colMeans(exprset[rownames(exprset) %in% module3.marker,])))
exp.4 <-  data.frame((colMeans(exprset[rownames(exprset) %in% module4.marker,])))
exp.5 <-  data.frame((colMeans(exprset[rownames(exprset) %in% module5.marker,])))
exp <- as.data.frame(cbind(exp.1,exp.2,exp.3,exp.4,exp.5))
colnames(exp) <- c("module1","module2","module3","module4","module5")
write.csv(exp,file = "./Survival analysis module rating_exp.csv")
#module1
#Calculate the optimal threshold
clinical$module1_score <- exp$module1
best_threshold_surv <- surv_cutpoint(clinical,
                                     time = "time",  # Survival time column name
                                     event = "event",      # Survival Event Listing
                                     variables = "module1_score",  # Need to find variable column names for threshold
                                     progressbar = TRUE)  # Display progress bar
module1_score <- ifelse((clinical$module1_score>best_threshold_surv$cutpoint$cutpoint),'high','low')
clinical$module1_score <- module1_score


#module2
#Calculate the optimal threshold
clinical$module2_score <- exp$module2
best_threshold_surv <- surv_cutpoint(clinical,
                                     time = "time",  # Survival time column name
                                     event = "event",      # Survival Event Listing
                                     variables = "module2_score",  # Need to find variable column names for threshold
                                     progressbar = TRUE)  # Display progress bar
module2_score <- ifelse((clinical$module2_score>best_threshold_surv$cutpoint$cutpoint),'high','low')
clinical$module2_score <- module2_score

#module3
clinical$module3_score <- exp$module3
best_threshold_surv <- surv_cutpoint(clinical,
                                     time = "time",  # Survival time column name
                                     event = "event",      # Survival Event Listing
                                     variables = "module3_score",  # Need to find variable column names for threshold
                                     progressbar = TRUE)  # Display progress bar
module3_score <- ifelse((clinical$module3_score>best_threshold_surv$cutpoint$cutpoint),'high','low')
clinical$module3_score <- module3_score

#module4
clinical$module4_score <- exp$module4
best_threshold_surv <- surv_cutpoint(clinical,
                                     time = "time",  # Survival time column name
                                     event = "event",      # Survival Event Listing
                                     variables = "module4_score",  # Need to find variable column names for threshold
                                     progressbar = TRUE)  # Display progress bar
module4_score <- ifelse((clinical$module4_score>best_threshold_surv$cutpoint$cutpoint),'high','low')
clinical$module4_score <- module4_score

#module5
clinical$module5_score <- exp$module5
best_threshold_surv <- surv_cutpoint(clinical,
                                     time = "time",  # Survival time column name
                                     event = "event",      # Survival Event Listing
                                     variables = "module5_score",  # Need to find variable column names for threshold
                                     progressbar = TRUE)  # Display progress bar
module5_score <- ifelse((clinical$module5_score>best_threshold_surv$cutpoint$cutpoint),'high','low')
clinical$module5_score <- module5_score

exprSet <- expr
meta = clinical
sfit1=survfit(Surv(time, event) ~ module1_score, data=meta)


ggsurvplot(sfit1,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change the color of the risk table based on stratification
           linetype = "strata", # Change line type based on layering
           surv.median.line = "hv", # Simultaneously display vertical and horizontal reference lines
           xlab ="Time in months", 
           ggtheme = theme_bw(), # Change the theme of ggplot2
           palette = c("#ff1f39","#4089c0"))#c9beff

sfit1=survfit(Surv(time, event) ~ module2_score, data=meta)


ggsurvplot(sfit1,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change the color of the risk table based on stratification
           linetype = "strata", # Change line type based on layering
           surv.median.line = "hv", # Simultaneously display vertical and horizontal reference lines
           xlab ="Time in months", 
           ggtheme = theme_bw(), # Change the theme of ggplot2
           palette = c("#ff1f39","#4089c0"))#c9beff


sfit1=survfit(Surv(time, event) ~ module3_score, data=meta)



ggsurvplot(sfit1,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change the color of the risk table based on stratification
           linetype = "strata", # Change line type based on layering
           surv.median.line = "hv", # Simultaneously display vertical and horizontal reference lines
           xlab ="Time in months", 
           ggtheme = theme_bw(), # Change the theme of ggplot2
           palette = c("#ff1f39","#4089c0"))#c9beff

sfit1=survfit(Surv(time, event) ~ module4_score, data=meta)


ggsurvplot(sfit1,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change the color of the risk table based on stratification
           linetype = "strata", # Change line type based on layering
           surv.median.line = "hv", # Simultaneously display vertical and horizontal reference lines
           xlab ="Time in months", 
           ggtheme = theme_bw(), # Change the theme of ggplot2
           palette = c("#ff1f39","#4089c0"))#c9beff

sfit1=survfit(Surv(time, event) ~ module5_score, data=meta)


ggsurvplot(sfit1,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change the color of the risk table based on stratification
           linetype = "strata", # Change line type based on layering
           surv.median.line = "hv", # Simultaneously display vertical and horizontal reference lines
           xlab ="Time in months", 
           ggtheme = theme_bw(), # Change the theme of ggplot2
           palette = c("#ff1f39","#4089c0"))#c9beff



