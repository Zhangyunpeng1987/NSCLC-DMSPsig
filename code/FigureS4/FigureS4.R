# ------------- Figure S4 --------------
#----Figure S4A----
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
#CM2
celltype_col <- c('#f1b38a', '#f0db69', '#a2c246', '#f5cee0', '#b6d2b7','#cbd2e5','#ad98c3','#be95db','#53A85F','#E5D2DD', '#F1BB72', '#F3B1A0',"#d6d5b7",'#D6E7A3',"#0073C2FF", "#438F49" ,"#699DC8","#91D1C2" ,"#DB873F" )
netVisual_circle(module2.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use = celltype_col)
netVisual_circle(module2.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use = celltype_col)

#----Figure S4B----
#CM5
celltype_col <- c('#E5D2DD','#23A23B','#F3B1A0','#57C3F3','#E95C59','#F7C394',"#0073C2FF","#C9B8D1","#93b8db","#7eb4c6","#f5e886",'#F1BB72')
netVisual_circle(module5.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use = celltype_col)
netVisual_circle(module5.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use = celltype_col)

#----Figure S4C----
# Identify the largest contribution signal 
library(RColorBrewer)
display.brewer.all()
module2_all <- readRDS("./wss/cellchat/module2/module2.rds")
celltype_col <- c('#f1b38a', '#f0db69', '#a2c246', '#f5cee0', '#b6d2b7','#cbd2e5','#ad98c3',"#b4ffde",'#53A85F','#E5D2DD', '#F1BB72', '#F3B1A0',"#0073C2FF", "#438F49" ,"#699DC8","#91D1C2" ,"#DB873F","#d6d5b7",'#D6E7A3')
names(celltype_col) <- c("Astrocyte_01_CAV1","Astrocyte_02_MALAT1","Astrocyte_03_ROM1","Astrocyte_04_SCGB3A1","Astrocyte_05_HLA-DRB5",
                         "Astrocyte_06_CCL5","B_04_PTMAP2","CD4T_04_SNORD3A","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67",
                         "GE1","GE2","GE8","GE9","GE11","Macrophage_06_C3","Neurophil_01_GOS2")
# Get celltype with high communication intensity
cell <- c("Astrocyte_05_HLA-DRB5","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67","GE1","GE2","GE8","GE9","GE11","Macrophage_06_C3")
module2_all <- subset(module2_all,subcluster %in% cell)
data.input <- module2_all@assays$RNA$data
meta = module2_all@meta.data
cell.use = rownames(meta)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "subcluster")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "subcluster")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

# Set database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling",key = "annotation")
cellchat@DB <- CellChatDB.use

# Subset the expression data to save fast computate
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 2 * 1024^3)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Expression data onto the PPI
cellchat <- projectData(cellchat, PPI.human)

cellchat@meta$subcluster <- factor(cellchat@meta$subcluster,levels = cell)
cellchat@idents <- cellchat@meta$subcluster
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Communication at the signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
celltype_col <- c('#b6d2b7','#53A85F','#E5D2DD', '#F1BB72', '#F3B1A0',"#0073C2FF", "#438F49" ,"#699DC8","#91D1C2" ,"#DB873F","#d6d5b7")
names(celltype_col) <- cell
# Calculate network centrality score
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
netAnalysis_signalingRole_network(cellchat, signaling = "CXCL", width = 8, height = 2.5, font.size = 6)
library(viridis)
# Function to change color
netAnalysis_signalingRole_heatmap2<- 
  function (object, signaling = NULL, pattern = c("outgoing", "incoming", 
                                                  "all"), slot.name = "netP", color.use = NULL, color.heatmap = "BuGn", 
            title = NULL, width = 10, height = 8, font.size = 8, font.size.title = 10, 
            cluster.rows = FALSE, cluster.cols = FALSE,color.heatmap.use = c('gray','red')){
    pattern <- match.arg(pattern)
    if (length(slot(object, slot.name)$centr) == 0) {
      stop("Please run `netAnalysis_computeCentrality` to compute the network centrality scores! ")
    }
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (i in 1:length(centr)) {
      outgoing[, i] <- centr[[i]]$outdeg
      incoming[, i] <- centr[[i]]$indeg
    }
    if (pattern == "outgoing") {
      mat <- t(outgoing)
      legend.name <- "Outgoing"
    }
    else if (pattern == "incoming") {
      mat <- t(incoming)
      legend.name <- "Incoming"
    }
    else if (pattern == "all") {
      mat <- t(outgoing + incoming)
      legend.name <- "Overall"
    }
    if (is.null(title)) {
      title <- paste0(legend.name, " signaling patterns")
    }
    else {
      title <- paste0(paste0(legend.name, " signaling patterns"), 
                      " - ", title)
    }
    if (!is.null(signaling)) {
      mat1 <- mat[rownames(mat) %in% signaling, , drop = FALSE]
      mat <- matrix(0, nrow = length(signaling), ncol = ncol(mat))
      idx <- match(rownames(mat1), signaling)
      mat[idx[!is.na(idx)], ] <- mat1
      dimnames(mat) <- list(signaling, colnames(mat1))
    }
    mat.ori <- mat
    mat <- sweep(mat, 1L, apply(mat, 1, max), "/", check.margin = FALSE)
    mat[mat == 0] <- NA
    if (is.null(color.use)) {
      color.use <- scPalette(length(colnames(mat)))
    }
    color.heatmap.use = color.heatmap.use
    df <- data.frame(group = colnames(mat))
    rownames(df) <- colnames(mat)
    names(color.use) <- colnames(mat)
    col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), 
                                        which = "column", show_legend = FALSE, show_annotation_name = FALSE, 
                                        simple_anno_size = grid::unit(0.2, "cm"))
    ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(mat.ori), 
                                                    border = FALSE, gp = gpar(fill = color.use, col = color.use)), 
                            show_annotation_name = FALSE)
    pSum <- rowSums(mat.ori)
    pSum.original <- pSum
    pSum <- -1/log(pSum)
    pSum[is.na(pSum)] <- 0
    idx1 <- which(is.infinite(pSum) | pSum < 0)
    if (length(idx1) > 0) {
      values.assign <- seq(0, 8, 
                           length.out = length(idx1))
      position <- sort(pSum.original[idx1], index.return = TRUE)$ix
      pSum[idx1] <- values.assign[match(1:length(idx1), position)]
    }
    ha1 = rowAnnotation(Strength = anno_barplot(pSum, border = FALSE), 
                        show_annotation_name = FALSE)
    if (min(mat, na.rm = T) == max(mat, na.rm = T)) {
      legend.break <- max(mat, na.rm = T)
    }
    else {
      legend.break <- c(round(min(mat, na.rm = T), digits = 1), 
                        round(max(mat, na.rm = T), digits = 1))
    }
    ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", 
                  name = "Relative strength", bottom_annotation = col_annotation, 
                  top_annotation = ha2, right_annotation = ha1, cluster_rows = cluster.rows, 
                  cluster_columns = cluster.rows, row_names_side = "left", 
                  row_names_rot = 0, row_names_gp = gpar(fontsize = font.size), 
                  column_names_gp = gpar(fontsize = font.size), width = unit(width, 
                                                                             "cm"), height = unit(height, "cm"), column_title = title, 
                  column_title_gp = gpar(fontsize = font.size.title), column_names_rot = 90, 
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, 
                                                              fontface = "plain"), title_position = "leftcenter-rot", 
                                              border = NA, at = legend.break, legend_height = unit(20, 
                                                                                                   "mm"), labels_gp = gpar(fontsize = 8), grid_width = unit(2, 
                                                                                                                                                            "mm")))
    return(ht1)
  }

library(ComplexHeatmap)

ht1 <- netAnalysis_signalingRole_heatmap2(cellchat, pattern = "outgoing",color.heatmap.use  = c("white","#75A9BF"),,color.use = celltype_col) 
ht2 <- netAnalysis_signalingRole_heatmap2(cellchat, pattern = "incoming",color.heatmap.use  = c("white","#75A9BF"),,color.use = celltype_col) 
setwd("./t100553/wss/cellchat/module2")
pdf(sprintf("in.out.pathway_key.pdf"), width = 10, height = 6)
ht1 + ht2
dev.off()

netVisual_aggregate(M2,signaling = )

save(cellchat,file = "./t100553/wss/cellchat/module2/module2.cellchat.key.RData")

#----Figure S4D----
# Identify the signals in module 5
setwd("./t100553/wss/cellchat/module5")
library(RColorBrewer)
display.brewer.all()
celltype_col <- c('#E5D2DD', '#F8BB88', '#F3B1A0','#57C3F3', '#E95C59','#F7C394',"#0073C2FF","#C9B8D1","#93b8db","#7eb4c6","#f5e886", '#F1BB72')
names(celltype_col) <- c("Endothelial_01_FCN3","Endothelial_03_VCAM1","Endothelial_04_DKK2","Endothelial_06_CD68","Endothelial_08_CCL21","Fibroblast_06_Alveolar_MS4A7",
                         "GE1","GE12","Macrophage_01_LGMN","Macrophage_03_MARCO","Macrophage_05_CSTB","Mast_02_CHIT1")
cell <- c("Endothelial_03_VCAM1","Endothelial_04_DKK2","Endothelial_06_CD68","Endothelial_08_CCL21","Fibroblast_06_Alveolar_MS4A7",
          "GE1","Macrophage_01_LGMN","Macrophage_03_MARCO","Mast_02_CHIT1")
celltype_col <- c( '#F8BB88', '#F3B1A0','#57C3F3', '#E95C59','#F7C394',"#0073C2FF","#93b8db","#7eb4c6",'#F1BB72')
names(celltype_col) <- cell
module5_all <- subset(module5_all,subcluster %in% cell)
data.input <- module5_all@assays$RNA$data
meta = module5_all@meta.data
cell.use = rownames(meta)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "subcluster")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "subcluster") 
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents))


CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling",key = "annotation")
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 2 * 1024^3)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- projectData(cellchat, PPI.human)

cellchat@meta$subcluster <- factor(cellchat@meta$subcluster,levels = cell)
cellchat@idents <- cellchat@meta$subcluster

ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# Calculate network centrality score
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = "CXCL", width = 8, height = 2.5, font.size = 6)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",color.heatmap = 'Blues',color.use = celltype_col) 
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",color.heatmap = 'Blues',color.use = celltype_col)
setwd("./t100553/wss/cellchat/module5")
pdf(sprintf("in.out.pathway_key.pdf"), width = 10, height = 6)
ht1 + ht2
dev.off()
save(cellchat,file = "module5.cellchat.key.RData")

#----Figure S4E----
library(scMetabolism)
library(Seurat)
library(AUCell)
library(GSEABase)
library(GSVA)
library(pheatmap)
library(dplyr)
sc.metabolism.Seurat <- function(obj, method = "AUCell", imputation = FALSE, ncores = 2,
                                 metabolism.type = "KEGG") {
  
  countexp <- as.data.frame(as.matrix(obj@assays$RNA$counts))
  
  gmtFile <- system.file(
    "data",
    ifelse(metabolism.type == "KEGG",
           "KEGG_metabolism_nc.gmt",
           "REACTOME_metabolism.gmt"),
    package = "scMetabolism"
  )
  
  countexp2 <- countexp
  
  if (method == "AUCell") {
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), nCores = ncores, plotStats = FALSE)
    geneSets <- getGmt(gmtFile)
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
    signature_exp <- as.data.frame(getAUC(cells_AUC))
  }
  
  obj@assays$METABOLISM$score <- signature_exp
  return(obj)
}

# read data
module2 <- readRDS("./t100553/wss/cellchat/module2/module2.rds")
module5 <- readRDS("./t100553/wss/cellchat/module5/module5.rds")

# add information
module2$type <- ifelse(startsWith(module2$subcluster, "GE"), "module2.epi", "module2.TIME")
module5$type <- ifelse(startsWith(module5$subcluster, "GE"), "module5.epi", "module5.TIME")

# merge
ss <- merge(module2, module5)
ss <- JoinLayers(ss)

# Score of metabolism
ss <- sc.metabolism.Seurat(ss, method = "AUCell", ncores = 2, metabolism.type = "KEGG")

df <- as.data.frame(t(ss@assays$METABOLISM$score))
df$type <- ss$type

avg_df <- aggregate(df[, -ncol(df)], list(df$type), mean)
rownames(avg_df) <- avg_df$Group.1
avg_df <- t(avg_df[, -1])

# plot
pdf("./all.metabolism.TIME_epi.pdf", width = 6, height = 10)
pheatmap(avg_df,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA,
         color = colorRampPalette(c("#40A429", "white", "#825695"))(200),
         fontsize = 6)
dev.off()
