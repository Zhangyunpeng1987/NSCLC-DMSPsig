# ------------- Figure 4 --------------
#----Figure 4A----
#Create CellChat objects for module2 and module5
library(CellChat)
library(patchwork)
library(tidyverse)
options(stringsAsFactors = FALSE)
corsce <- readRDS("./wss/subsce.corr.rds")
corsce@meta.data <- corsce@meta.data[,c(4,5,25,24,27)]
cancer.epi <- readRDS(file = "./wss/copykat/cancer state/cancer.epi_GE.rds")
cancer.epi@meta.data <- cancer.epi@meta.data %>% mutate(cellmodule = case_when(
  patient %in% c("P22","P43","P9","P1","P4","P25","P10","P27") ~ "module1",
  patient %in% c("P18","P15","P38","P14","P44","P12","P13","P3") ~ "module2",
  patient %in% c("P16","P42","P2","P23","P8","P41","P26","P20","P17","P45","P21","P36","P7","P11","P19","P39","P34","P5") ~ "module3",
  patient %in% c("P31","P28","P30","P29","P35","P32") ~ "module4",
  patient %in% c("P24","P40","P37","P33","P6") ~ "module5"))
cancer.epi$subcluster <- cancer.epi$GEtype
cancer.epi@meta.data <- cancer.epi@meta.data[,c(4,5,21,44,43)]

#Merge TIME cells and Cancer.epi
sce <- merge(corsce,cancer.epi)
sce <- JoinLayers(sce)

module2 <- subset(sce,cellmodule=="module2")
module5 <- subset(sce,cellmodule=="module5")
saveRDS(module2,file = "./module2.rds")
saveRDS(module5,file = "./module5.rds")
#Create separate CellChat objects
library(Seurat)
#module2
setwd("./t100553/wss/cellchat/module2/")
module2_all <- subset(sce,cellmodule == "module2")
table(module2_all$subcluster)
#Remove cell clusters with low cell counts
module2_all <- subset(module2_all,subcluster %in% c("Astrocyte_01_CAV1","Astrocyte_02_MALAT1","Astrocyte_03_ROM1","Astrocyte_04_SCGB3A1","Astrocyte_05_HLA-DRB5",
                                                    "Astrocyte_06_CCL5","B_04_PTMAP2","CD4T_04_SNORD3A","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67",
                                                    "GE1","GE2","GE8","GE9","GE11","Macrophage_06_C3","Neurophil_01_GOS2"))
data.input <- module2_all@assays$RNA$data
meta = module2_all@meta.data
cell.use = rownames(meta)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "subcluster")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "subcluster") 
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

#Set ligand-receptor interaction database
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
#Show database structure
dplyr::glimpse(CellChatDB$interaction)
#Perform cell communication analysis using a subset of CellChatDB
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling",key = "annotation")
#Set the database used in the object
cellchat@DB <- CellChatDB.use

#Subset signaling gene expression data
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 2 * 1024^3)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

ptm = Sys.time()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

#Project gene expression data onto PPI
cellchat <- projectData(cellchat, PPI.human)

#Infer cell-cell communication network
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
#Infer cell-cell communication at the signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#Calculate aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

#Visualize aggregated cell-cell communication network

cellchat@meta$subcluster <- factor(cellchat@meta$subcluster,levels = c("Astrocyte_01_CAV1","Astrocyte_02_MALAT1","Astrocyte_03_ROM1","Astrocyte_04_SCGB3A1","Astrocyte_05_HLA-DRB5",
                                                                       "Astrocyte_06_CCL5","B_04_PTMAP2","CD4T_04_SNORD3A","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67",
                                                                       "Macrophage_06_C3","Neurophil_01_GOS2","GE1","GE2","GE8","GE9","GE11"))
celltype_col <- c('#f1b38a', '#f0db69', '#a2c246', '#f5cee0', '#b6d2b7','#cbd2e5','#ad98c3','#be95db','#53A85F','#E5D2DD', '#F1BB72', '#F3B1A0',"#d6d5b7",'#D6E7A3',"#0073C2FF", "#438F49" ,"#699DC8","#91D1C2" ,"#DB873F" )

names(celltype_col) <- c("Astrocyte_01_CAV1","Astrocyte_02_MALAT1","Astrocyte_03_ROM1","Astrocyte_04_SCGB3A1","Astrocyte_05_HLA-DRB5",
                         "Astrocyte_06_CCL5","B_04_PTMAP2","CD4T_04_SNORD3A","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67",
                         "Macrophage_06_C3","Neurophil_01_GOS2","GE1","GE2","GE8","GE9","GE11")
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use = celltype_col)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use = celltype_col)
#ggsave("module2_cellchat1.pdf",width = 20,height = 20,units = "cm")

#Check signals sent from each cell group
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), label.edge= F, title.name = rownames(mat)[i],color.use = celltype_col)
}
#ggsave("module2_cellchat2.pdf",width = 10,height = 8,units = "cm")

#Calculate network centrality scores 
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
#Visualize calculated centrality scores using a heatmap 
netAnalysis_signalingRole_network(cellchat, signaling = "CXCL", width = 8, height = 2.5, font.size = 6)
gg1 <- netAnalysis_signalingRole_scatter(cellchat,pattern = "outgoing",color.use = celltype_col)
gg2 <- netAnalysis_signalingRole_scatter(cellchat,pattern = "incoming",color.use = celltype_col)
gg1
gg1 + gg2
save(cellchat,file = "./monocle2_cellchat.RData")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "all",color.heatmap = 'Purples',color.use = celltype_col)
save(cellchat,file = "./module2_cellchat_subcluster.RData")

#Expression of all genes for a specific signal across cell groups
mycol1 <- c("#33CC66","#33CCFF","#FF9966","#FF99CC","#FF9999","#6699FF","#66CC00")
mycol2 <- c("#33CCFF","#FF9966","#FF99CC","#FF9999","#6699FF","#66CC00")
mycol3 <- c("#99FF99","#33CCFF","#FF9966","#FF99CC","#6699FF","#9999FF","#00CCCC","#CCFF33","#66CC00")

plotGeneExpression(cellchat, signaling = "PSAP")#,color.use = mycol1)
mycolors2 <- rainbow(20)
p2 = plotGeneExpression(cellchat, signaling = "SPP1", type = "dot",color.use = mycolors2)
p2
#Identify signaling groups based on functional similarity
setwd("./t100553/wss/cellchat/module2")
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#Classification learning of the signaling networks for a single dataset
#Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
save(cellchat,file = "./module2_cellchat.RData")

#module5
dir.create("./t100553/wss/cellchat/module5/")
setwd("./t100553/wss/cellchat/module5/")
saveRDS(module5,file = "module5.rds")
module5_all <- readRDS(file = "./module5.rds")
table(module5_all$subcluster)
#Remove cell clusters with low cell counts
module5_all <- subset(module5_all,subcluster %in% c("Endothelial_01_FCN3","Endothelial_03_VCAM1","Endothelial_04_DKK2","Endothelial_06_CD68","Endothelial_08_CCL21","Fibroblast_06_Alveolar_MS4A7","GE1","GE12",
                                                    "Macrophage_01_LGMN","Macrophage_03_MARCO","Macrophage_05_CSTB","Mast_02_CHIT1"))
data.input <- module5_all@assays$RNA$data
meta = module5_all@meta.data

cell.use = rownames(meta)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "subcluster")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "subcluster") 
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) 

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
#Show database structure
dplyr::glimpse(CellChatDB$interaction)
#Perform cell communication analysis using a subset of CellChatDB
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling",key = "annotation")
#Set the database used in the object
cellchat@DB <- CellChatDB.use

# Subset signaling gene expression data
cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 2 * 1024^3)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

ptm = Sys.time()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

#Project gene expression data onto PPI
cellchat <- projectData(cellchat, PPI.human)

#Infer cell-cell communication network
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
#Infer cell-cell communication at the signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

#Visualize the aggregated cell-cell communication network
cellchat@meta$subcluster <- factor(cellchat@meta$subcluster,levels = c("Endothelial_01_FCN3","Endothelial_03_VCAM1","Endothelial_04_DKK2","Endothelial_06_CD68","Endothelial_08_CCL21","Fibroblast_06_Alveolar_MS4A7","GE1","GE12",
                                                                       "Macrophage_01_LGMN","Macrophage_03_MARCO","Macrophage_05_CSTB","Mast_02_CHIT1"))
celltype_col <- c('#E5D2DD','#23A23B','#F3B1A0','#57C3F3','#E95C59','#F7C394',"#0073C2FF","#C9B8D1","#93b8db","#7eb4c6","#f5e886",'#F1BB72')
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions",color.use = celltype_col)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use = celltype_col)
#ggsave("module2_cellchat1.pdf",width = 20,height = 20,units = "cm")

#Check signaling inputs from each cell group
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), label.edge= F, title.name = rownames(mat)[i],color.use = celltype_col)
}
#ggsave("module2_cellchat2.pdf",width = 10,height = 8,units = "cm")

#Show the number or strength of outgoing and incoming signaling for each cell group
library(ggplot2)
df.net <- subsetCommunication(cellchat)
source_counts <- table(df.net$source)
target_counts <- table(df.net$target)

data <- data.frame(
  CellType = c(names(source_counts), names(target_counts)),
  Count = c(source_counts, target_counts),
  Category = factor(rep(c("Source", "Target"), each = length(source_counts)))
)

ggplot(data, aes(x = CellType, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Cell Type Counts by Source and Target", x = "Cell Type", y = "Count")+
  scale_color_manual(values = c("blue","red"))+
  theme_bw()

#Calculate network centrality scores 
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
#Visualize the calculated centrality scores using a heatmap
netAnalysis_signalingRole_network(cellchat, signaling = "CXCL", width = 8, height = 2.5, font.size = 6)

gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat)
gg1
gg1 + gg2

netAnalysis_signalingRole_heatmap(cellchat, pattern = "all",color.heatmap = 'Purples',color.use = celltype_col)
save(cellchat,file = "./module5_cellchat.RData")

#Plot the total number of interactions for module2 and module5
load("./t100553/wss/cellchat/module2/module2_cellchat.RData")
module2.cellchat <- cellchat
load("./t100553/wss/cellchat/module5/module5_cellchat.RData")
module5.cellchat <- cellchat
n_interactions1 <- length(module2.cellchat@LR$LRsig$interaction_name)
n_interactions2 <- length(module5.cellchat@LR$LRsig$interaction_name)

#Compare the strengths
weights1 <- module2.cellchat@net$weight
total_strength1 <- sum(weights1, na.rm = TRUE)
weights2 <- module5.cellchat@net$weight
total_strength2 <- sum(weights2, na.rm = TRUE)

library(ggplot2)

df <- data.frame(
  CellChat = c("module2", "module5"),
  Interactions = c(n_interactions1, n_interactions2),
  Strength = c(total_strength1, total_strength2)
)

#Comparison of the number of interactions
ggplot(df, aes(x = CellChat, y = Interactions, fill = CellChat)) +
  geom_bar(stat = "identity") +
  labs(title = "Communication Number", y = "Count") +
  theme_bw() +
  scale_fill_manual(values = c("#80C1D7", "#eb7e60"))

#Comparison of communication strengths
ggplot(df, aes(x = CellChat, y = Strength, fill = CellChat)) +
  geom_bar(stat = "identity") +
  labs(title = "Communication Strength", y = "Strength") +
  theme_bw()+
  scale_fill_manual(values = c("#80C1D7", "#eb7e60"))

#----Figure 4B----
#Plot the communication strength between all cells in module2 and module5
pheatmap::pheatmap(module2.cellchat@net$count, border_color = "grey", 
                     cluster_cols = F, fontsize = 10, cluster_rows = F,
                     display_numbers = F,angle_col = 45,number_color="black",number_format = "%.0f",color = colorRampPalette(c("white","#B8DBE7", "#1E6FB1"))(100))

pheatmap::pheatmap(module5.cellchat@net$count, border_color = "grey", 
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = F,number_color="black",number_format = "%.0f",angle_col = 45,colorRampPalette(c("white","#F6DAD3", "#eb7e60"))(100))

#----Figure 4C----
#CM2
celltype_col <- c('#f1b38a', '#f0db69', '#a2c246', '#f5cee0', '#b6d2b7','#cbd2e5','#ad98c3','#be95db','#53A85F','#E5D2DD', '#F1BB72', '#F3B1A0',"#d6d5b7",'#D6E7A3',"#0073C2FF", "#438F49" ,"#699DC8","#91D1C2" ,"#DB873F" )
ht1 <- netAnalysis_signalingRole_scatter(module2.cellchat,color.use = celltype_col)
#CM5
celltype_col <- c('#E5D2DD','#23A23B','#F3B1A0','#57C3F3','#E95C59','#F7C394',"#0073C2FF","#C9B8D1","#93b8db","#7eb4c6","#f5e886",'#F1BB72')
ht2 <- netAnalysis_signalingRole_scatter(module5.cellchat,color.use = celltype_col)

#----Figure 4D----
#CM2
celltype_col <- c('#f1b38a', '#f0db69', '#a2c246', '#f5cee0', '#b6d2b7','#cbd2e5','#ad98c3','#be95db','#53A85F','#E5D2DD', '#F1BB72', '#F3B1A0',"#d6d5b7",'#D6E7A3',"#0073C2FF", "#438F49" ,"#699DC8","#91D1C2" ,"#DB873F" )
netVisual_aggregate(module2.cellchat, signaling = "SPP1",  vertex.receiver = c(1:5), layout = "hierarchy",color.use = celltype_col)
netVisual_aggregate(module2.cellchat, signaling = "MK",  vertex.receiver = c(1:5), layout = "hierarchy",color.use = celltype_col)
#CM5
celltype_col <- c('#E5D2DD','#23A23B','#F3B1A0','#57C3F3','#E95C59','#F7C394',"#0073C2FF","#C9B8D1","#93b8db","#7eb4c6","#f5e886",'#F1BB72')
netVisual_aggregate(module5.cellchat, signaling = "SPP1",  vertex.receiver = c(1:5), layout = "hierarchy",,color.use = celltype_col)
netVisual_aggregate(module5.cellchat, signaling = "MK",  vertex.receiver = c(1:5), layout = "hierarchy",color.use = celltype_col)

#----Figure 4E----
#Plot the ligand-receptor chord diagram for module2
table(module2.cellchat@meta$subcluster)
celltype_col <- c('#f1b38a', '#f0db69', '#a2c246', '#f5cee0', '#b6d2b7','#cbd2e5','#ad98c3','#be95db','#53A85F','#E5D2DD', '#F1BB72', '#F3B1A0',"#d6d5b7",'#D6E7A3',"#0073C2FF", "#438F49" ,"#699DC8","#91D1C2" ,"#DB873F" )
names(celltype_col) <- c("Astrocyte_01_CAV1","Astrocyte_02_MALAT1","Astrocyte_03_ROM1","Astrocyte_04_SCGB3A1","Astrocyte_05_HLA-DRB5",
                         "Astrocyte_06_CCL5","B_04_PTMAP2","CD4T_04_SNORD3A","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67",
                         "Macrophage_06_C3","Neurophil_01_GOS2","GE1","GE2","GE8","GE9","GE11")

#①MK GALECTIN VEGF 
col <- celltype_col[c(8:10,12,15,16,17,19)]
netVisual_chord_gene(module2.cellchat,signaling = c("MK","GALECTIN","VEGF"), sources.use = c(15,16,19) ,targets.use = c(8,9,10,12,17),lab.cex = 0.5,legend.pos.y = 30,color.use = col)
#②SPP1 PTN PSAP FGF 
col <- celltype_col[c(3:5,12:13,15,16,19)]
netVisual_chord_gene(module2.cellchat,signaling = c("SPP1","PTN","PSAP","FGF"), sources.use = c(3:5,13) ,targets.use = c(4,12,15,16,19),lab.cex = 0.5,legend.pos.y = 30,color.use = col)
#VISFATIN
col <- celltype_col[c(1,3,9,13:16,18,19)]
netVisual_chord_gene(module2.cellchat,signaling = c("GRN","VISFATIN","EGF"), sources.use = c(13,14,18) ,targets.use = c(1,3,9,15,16,19),lab.cex = 0.5,legend.pos.y = 30,color.use = col)

setwd("./t100553/wss/cellchat/module2/")
load("./wss/cellchat/module2/module2.cellchat.key.RData")
celltype_col <- c('#b6d2b7','#53A85F','#E5D2DD', '#F1BB72', '#F3B1A0',"#0073C2FF", "#438F49" ,"#699DC8","#91D1C2" ,"#DB873F","#d6d5b7")
cell <- c("Astrocyte_05_HLA-DRB5","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67","GE1","GE2","GE8","GE9","GE11","Macrophage_06_C3")
names(celltype_col) <- cell


#astro5 -> endo2 + fibro2/3
col <- celltype_col[c(1,2,4,5)]
netVisual_chord_gene(cellchat,signaling = c("PTN","SPP1","ANGPTL","ANGPT","PDGF","TGFb","BMP"), sources.use = c(1,2,4,5) ,targets.use = c(1,2,4,5),lab.cex = 0.5,legend.pos.y = 30,color.use = col)
#astro5 + endo2 + fibro1/2 +GE8 -> endo2
col <- celltype_col[c(1,2,3,4,8)]
netVisual_chord_gene(cellchat,signaling = c("VEGF","VISFATIN","ANGPTL","GDF","SEMA3","ANGPT","CALCR","BMP","APELIN"), sources.use = c(1,2,3,4,8) ,targets.use = 2,lab.cex = 0.5,legend.pos.y = 30,color.use = col)
#endo2 + fibro2 + GE11 + macro6 -> macro6
col <- celltype_col[c(2,4,10,11)]
netVisual_chord_gene(cellchat,signaling = c("GALECTIN","MIF","ANNEXIN","PROS","GAS","CXCL","CSF","TNF","TGFb","CCL","CX3C"),sources.use = c(2,4,10,11) ,targets.use = 11,lab.cex = 0.5,legend.pos.y = 30,color.use = col)


#Plot the ligand-receptor chord diagram for module5
table(module5.cellchat@meta$subcluster)
celltype_col <- c('#E5D2DD', '#F8BB88', '#F3B1A0','#57C3F3', '#E95C59','#F7C394',"#0073C2FF","#C9B8D1","#93b8db","#7eb4c6","#f5e886", '#F1BB72')
names(celltype_col) <- c("Endothelial_01_FCN3","Endothelial_03_VCAM1","Endothelial_04_DKK2","Endothelial_06_CD68","Endothelial_08_CCL21","Fibroblast_06_Alveolar_MS4A7",
                         "GE1","GE12","Macrophage_01_LGMN","Macrophage_03_MARCO","Macrophage_05_CSTB","Mast_02_CHIT1")

#①CXCL CCL GDF VEGF 
col <- celltype_col[c(1:3,5:10)]
netVisual_chord_gene(module5.cellchat,signaling = c("CXCL","CCL","GDF","VEGF"), sources.use = c(2,6:9,10) ,targets.use = c(1:3,5),lab.cex = 0.5,legend.pos.y = 30,color.use = col)

#②SPP1 GALECTIN MIF RESISTIN  
col <- celltype_col[c(2,4,6,7,9:12)]
netVisual_chord_gene(module5.cellchat,signaling = c("SPP1","GALECTIN","MIF","RESISTIN"), sources.use = c(2,4,6,7,9:11) ,targets.use = c(6,9,10,12),lab.cex = 0.5,legend.pos.y = 30,color.use = col)

#MK VISFATIN
col <- celltype_col[c(1:4,6,7,9,12)]
netVisual_chord_gene(module5.cellchat,signaling = c("MK","VISFATIN"), sources.use = c(2,4,6,7,9,12) ,targets.use = c(1,2,3,4,7),lab.cex = 0.5,legend.pos.y = 30,color.use = col)

table(cellchat@meta$subcluster)
celltype_col <- c('#F8BB88', '#F3B1A0','#57C3F3', '#E95C59','#F7C394',"#0073C2FF","#93b8db","#7eb4c6", '#F1BB72')
cell <- c("Endothelial_03_VCAM1","Endothelial_04_DKK2","Endothelial_06_CD68 ","Endothelial_08_CCL21","Fibroblast_06_Alveolar_MS4A7","GE1","Macrophage_01_LGMN","Macrophage_03_MARCO","Mast_02_CHIT1")
names(celltype_col) <- cell
#macro1/3 mast2 ge1 endo8 -> endo3
col <- celltype_col[c(1,4,6,7:9)]
netVisual_chord_gene(cellchat,signaling = c("CXCL","CCL","VISFATIN","VEGF","GDF","CALCR","ANGPT","OSM","TNF","IGF"), sources.use = c(4,6,7:9) ,targets.use = c(1),lab.cex = 0.5,legend.pos.y = 30,color.use = col)
#endo3/4/6/8 fibro6 ge1 macro1/3 -> macro1/3 mast2
col <- celltype_col[c(1:9)]
netVisual_chord_gene(cellchat,signaling = c("SPP1","MIF","GALECTIN","RESISTIN","COMPLEMENT","GAS","PROS"), sources.use = c(1:8) ,targets.use = c(7:9),lab.cex = 0.5,legend.pos.y = 30,color.use = col)
#GE1 ->  GE1 macro1/3 mast endo3/4/8 fibro6 
col <- celltype_col[c(1,2,4:9)]
netVisual_chord_gene(cellchat,signaling = c("MK","MIF","VEGF","UGRP1","GDF","EGF","KIT","PDGF"), sources.use = c(6) ,targets.use = c(1,2,4:9),lab.cex = 0.5,legend.pos.y = 30,color.use = col)


#----Figure 4F----
#HR analysis
load("./t100553/wss/cellchat/module5/JS_M_HR.RData")
M2 <- subset(M2_HR_result,p<0.05)
M2 <- subset(M2,cancer%in%c("LUAD","LUSC"))
M5 <- subset(M5_HR_result,p<0.05)
M5 <- subset(M5,cancer%in%c("LUAD","LUSC"))

library(GSVA)
library(TCGAplot)
library(tidyverse)
library(survival)
library(meta)
library(doParallel)  

registerDoParallel(cores = 20)
tpm <- get_all_tpm()
meta <- get_all_meta()
cancers <- unique(meta$Cancer)

dt <- read_csv("Figure5E_net_nmf.csv")
a <- dt %>% select(interaction_name_2, interaction_name) %>% separate_rows(interaction_name, sep = "_") %>% distinct()
genelist <- split(a$interaction_name, a$interaction_name_2)

final_pan <- foreach(geneset_name = names(genelist), 
                     .combine = bind_rows,
                     .packages = c("GSVA", "tidyverse", "survival", "meta", "TCGAplot"),
                     .export = c("tpm", "meta", "cancers")) %dopar% {
                       current_geneset <- genelist[[geneset_name]]
                       current_genelist <- list(current_geneset)
                       names(current_genelist) <- geneset_name
                       
                       cox_results <- list()
                       
                       for (cancer in cancers) {
                         exprSet <- subset(tpm, Group == "Tumor" & Cancer == cancer) %>% 
                           tibble::add_column(ID = stringr::str_sub(rownames(.), 1, 12), .before = "Cancer") %>% 
                           dplyr::filter(!duplicated(ID)) %>% 
                           tibble::remove_rownames() %>% 
                           tibble::column_to_rownames("ID") %>% 
                           dplyr::filter(rownames(.) %in% rownames(subset(meta, Cancer == cancer)))
                         
                         exprSet <- exprSet[, -(1:2)] %>% as.matrix() %>% t()
                         
                         gsvapar <- gsvaParam(exprData = exprSet, 
                                              geneSets = current_genelist, 
                                              kcdf = "Gaussian")
                         exprSet_gsva <- gsva(gsvapar)
                         
                         cl <- meta[colnames(exprSet_gsva), ]
                         cl$symbol <- exprSet_gsva[geneset_name, ]
                         
                         if (sum(!is.na(cl$time) & !is.na(cl$event)) > 0) {
                           m <- tryCatch(
                             coxph(Surv(time, event) ~ symbol + age, data = cl),
                             error = function(e) NULL
                           )
                           
                           if (!is.null(m)) {
                             beta <- coef(m)
                             se <- sqrt(diag(vcov(m)))
                             tmp <- round(cbind(
                               HR = exp(beta),
                               se = se,
                               lower = exp(beta - 1.96 * se),
                               upper = exp(beta + 1.96 * se),
                               p = 1 - pchisq((beta/se)^2, 1)
                             ), 3)
                             cox_results[[cancer]] <- tmp["symbol", ]
                           }
                         }
                       }
                       
                       if (length(cox_results) > 0) {
                         a <- do.call(rbind, cox_results) %>% 
                           as.data.frame() %>% 
                           rownames_to_column("cancer") %>% 
                           select(cancer, HR, se, lower, upper, p)
                        
                         meta_res <- metagen(log(a$HR), a$se, sm = "HR")
                         pan_row <- data.frame(
                           cancer = "Pan cancer",
                           HR = exp(meta_res$TE.random),
                           lower = exp(meta_res$lower.random),
                           upper = exp(meta_res$upper.random),
                           p = meta_res$pval.random,
                           stringsAsFactors = FALSE
                         )
                         
                         a <- a %>% select(cancer, HR, lower, upper, p)
                         
                         pan <- bind_rows(a, pan_row) %>% 
                           mutate(
                             celltype = geneset_name,
                             significance = case_when(
                               p <= 0.0001 ~ "****",
                               p <= 0.001 ~ "***",
                               p <= 0.01 ~ "**",
                               p <= 0.05 ~ "*",
                               TRUE ~ ""
                             ),
                             color = case_when(
                               HR > 1 ~ "Worse survival",
                               HR < 1 ~ "Better survival",
                               TRUE ~ "Neutral"
                             )
                           )
                         pan
                       } else {
                         NULL
                       }
                     }

stopImplicitCluster()

write.csv(final_pan, file = "TCGAplot_60_LRs_HR.csv")

final_result1 <- final_pan[final_pan$cancer == "Pan cancer", ]
final_result1 <- final_result1 %>% filter(p < 0.05)
table(final_result1$color) 
# Better survival  Worse survival 
# 9               12 
a <- dt[dt$interaction_name_2 %in% unique(final_result1$celltype), ]
write.csv(a, file = "./data1/yang/test/Results/Figure6/TCGA/pancancer_HightriskLRsgenes.csv")

#Heatmap of module-specific ligand-receptor pairs for HR-related driver genes
setwd("./t100553/wss/cellchat")
load("./genelist.RData")
load("./JS_M_HR.RData")
module2_all <- readRDS("./wss/cellchat/module2/module2.rds")
load("./wss/cellchat/module2/module2.cellchat.key.RData")

#Extract enriched ligand-receptor pairs for specified signaling pathways
enriched_LR <- extractEnrichedLR(cellchat, signaling = pathways)
pathway_gene_list <- list()

#Remove cell clusters with low cell counts
module2_all <- subset(module2_all,subcluster %in% c("Astrocyte_01_CAV1","Astrocyte_02_MALAT1","Astrocyte_03_ROM1","Astrocyte_04_SCGB3A1","Astrocyte_05_HLA-DRB5",
                                                    "Astrocyte_06_CCL5","B_04_PTMAP2","CD4T_04_SNORD3A","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67",
                                                    "GE1","GE2","GE8","GE9","GE11","Macrophage_06_C3","Neurophil_01_GOS2"))

module2_all$subcluster <- factor(module2_all$subcluster,levels = c("Astrocyte_01_CAV1","Astrocyte_02_MALAT1","Astrocyte_03_ROM1","Astrocyte_04_SCGB3A1","Astrocyte_05_HLA-DRB5",
                                                                   "Astrocyte_06_CCL5","B_04_PTMAP2","CD4T_04_SNORD3A","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67",
                                                                   "Macrophage_06_C3","Neurophil_01_GOS2","GE1","GE2","GE8","GE9","GE11"))
M2 <- subset(M2_HR_result,cancer == c("LUAD","LUSC"))
M2 <- subset(M2, p <0.05)
M2$LR <- gsub("\\(|\\)","",M2$celltype)
M2$LR <- gsub(" ", "",M2$LR)
M2$LR <- gsub("-","_",M2$LR)
M2$LR <- gsub("\\+","_",M2$LR)

library(stringr)
M2.name <- str_split(M2$LR, "_", simplify = FALSE)

names(M2.name) <- M2$LR
M2.name <- unique(M2.name)
names(M2.name) <- unique(M2$LR)

library(UCell)
library(Seurat)
module2_score <- AddModuleScore_UCell(module2_all,
                                      features = M2.name,
                                      assay = "RNA")


LR.score <- module2_score@meta.data[,c(1,4,6:15)]
colnames(LR.score) <- gsub("_UCell","",colnames(LR.score))

score <- data.frame(module2_score@meta.data[,c(6:15)])
library(pheatmap)
mutation_anno <- as.matrix(LR.score[,1])
colnames(mutation_anno) <- c("mutation")
subcluster_anno <- as.matrix(LR.score[,2])
colnames(subcluster_anno) <- c("subcluster")


anno_df <- data.frame(subcluster = subcluster_anno,mutation = mutation_anno)
rownames(anno_df) <- rownames(score)
annotation_colors <- list(
  subcluster = c(
    "Astrocyte_01_CAV1" = "#f1b38a",
    "Astrocyte_02_MALAT1" = "#f0db69",
    "Astrocyte_03_ROM1" = "#a2c246",
    "Astrocyte_04_SCGB3A1" = "#f5cee0",
    "Astrocyte_05_HLA-DRB5" = "#b6d2b7",
    "Astrocyte_06_CCL5" = "#cbd2e5",
    "B_04_PTMAP2" = "#ad98c3",
    "CD4T_04_SNORD3A" = "#be95db",
    "Endothelial_02_PTMAP2" = "#53A85F",
    "Fibroblast_01_Myo_HIGD1B" = "#E5D2DD",
    "Fibroblast_02_Myo_COL10A1" = "#F1BB72",
    "Fibroblast_03_Myo_MKI67" = "#F3B1A0",
    "Macrophage_06_C3" = "#d6d5b7",
    "Neurophil_01_GOS2" = "#D6E7A3",
    "GE1" = "#0073C2FF",
    "GE2" = "#438F49",
    "GE8" = "#699DC8",
    "GE9" = "#91D1C2",
    "GE11" = "#DB873F"
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
colnames(score) <-  gsub("_UCell","",colnames(score))
setwd("./t100553/wss/cellchat/module2")
pdf("./t100553/wss/cellchat/module2/module2.HR.LR.pdf",width = 14, height = 6)
pheatmap(t(score),cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = F,annotation_col  = anno_df,annotation_colors = annotation_colors,color = colorRampPalette(c("white", "#EEE2F2","#D5B7E0","#C297D3","#B077C5","#A362BC","#4B0082"))(1000))#
pheatmap(t(score),cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = F,annotation_col  = anno_df,annotation_colors = annotation_colors,color = colorRampPalette(c("white","#F1D8D1","#E0BAB3","#D7978A","#C6775E","#C27664"))(1000))#

dev.off()

#module5 
module5_all <- readRDS("./wss/cellchat/module5/module5.rds")
load("./wss/cellchat/module5/module5.cellchat.key.RData")

#Remove cell clusters with insufficient cell counts
module5_all <- subset(module5_all,subcluster %in% c("Endothelial_01_FCN3","Endothelial_03_VCAM1","Endothelial_04_DKK2","Endothelial_06_CD68","Endothelial_08_CCL21","Fibroblast_06_Alveolar_MS4A7","GE1","GE12",
                                                    "Macrophage_01_LGMN","Macrophage_03_MARCO","Macrophage_05_CSTB","Mast_02_CHIT1"))
module5_all$subcluster <- factor(module5_all$subcluster,levels = c("Endothelial_01_FCN3","Endothelial_03_VCAM1","Endothelial_04_DKK2","Endothelial_06_CD68","Endothelial_08_CCL21","Fibroblast_06_Alveolar_MS4A7","GE1","GE12",
                                                                   "Macrophage_01_LGMN","Macrophage_03_MARCO","Macrophage_05_CSTB","Mast_02_CHIT1"))
celltype_col <- c('#E5D2DD','#23A23B','#F3B1A0','#57C3F3','#E95C59','#F7C394',"#0073C2FF","#C9B8D1","#93b8db","#7eb4c6","#f5e886",'#F1BB72')
M5 <- subset(M5_HR_result,cancer == c("LUAD","LUSC"))
M5 <- subset(M5, p <0.05)
M5$LR <- gsub("\\(|\\)","",M5$celltype)
M5$LR <- gsub(" ", "",M5$LR)
M5$LR <- gsub("-","_",M5$LR)
M5$LR <- gsub("\\+","_",M5$LR)

library(stringr)
M5.name <- str_split(M5$LR, "_", simplify = FALSE)

names(M5.name) <- M5$LR
M5.name <- unique(M5.name)
names(M5.name) <- unique(M5$LR)

library(UCell)
library(Seurat)
module5_score <- AddModuleScore_UCell(module5_all,
                                      features = M5.name,
                                      assay = "RNA")


LR.score <- module5_score@meta.data[,c(1,4,6:11)]
colnames(LR.score) <- gsub("_UCell","",colnames(LR.score))

score <- data.frame(module5_score@meta.data[,c(6:11)])
library(pheatmap)
mutation_anno <- as.matrix(LR.score[,1])
colnames(mutation_anno) <- c("mutation")
subcluster_anno <- as.matrix(LR.score[,2])
colnames(subcluster_anno) <- c("subcluster")

anno_df <- data.frame(subcluster = subcluster_anno,mutation = mutation_anno)
rownames(anno_df) <- rownames(score)
annotation_colors <- list(
  subcluster = c(
    "Endothelial_01_FCN3" = "#E5D2DD",
    "Endothelial_03_VCAM1" = "#23A23B",
    "Endothelial_04_DKK2" = "#F3B1A0",
    "Endothelial_06_CD68" = "#57C3F3",
    "Endothelial_08_CCL21" = "#E95C59",
    "Fibroblast_06_Alveolar_MS4A7" = "#F7C394",
    "GE1" = "#0073C2FF",
    "GE12" = "#C9B8D1",
    "Macrophage_01_LGMN" = "#93b8db",
    "Macrophage_03_MARCO" = "#7eb4c6",
    "Macrophage_05_CSTB" = "#f5e886",
    "Mast_02_CHIT1" = "#F1BB72"
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
colnames(score) <-  gsub("_UCell","",colnames(score))
setwd("./t100553/wss/cellchat/module5")
pdf("./t100553/wss/cellchat/module5/module5.LR.HR.pdf",width = 14, height = 5)
pheatmap(t(score),cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = F,annotation_col  = anno_df,annotation_colors = annotation_colors,color = colorRampPalette(c("white", "#EEE2F2","#D5B7E0","#C297D3","#B077C5","#A362BC","#4B0082"))(1000))#
pheatmap(t(score),cluster_rows = F, cluster_cols = F,show_rownames = T, show_colnames = F,annotation_col  = anno_df,annotation_colors = annotation_colors,color = colorRampPalette(c("white","#F1D8D1","#E0BAB3","#D7978A","#C6775E","#C27664"))(1000))#
dev.off()

#----Figure 4G----
#The cellular mechanism diagram was drawn using Adobe Illustrator 2021
#----Figure 4H----
library(scMetabolism)
library(Seurat)
library(dplyr)
library(pheatmap)

sc.metabolism.Seurat <- function(obj, method="AUCell", imputation=F, ncores=2, metabolism.type="KEGG"){
  countexp <- as.data.frame(as.matrix(obj@assays$RNA$counts))
  gmtFile <- system.file("data","KEGG_metabolism_nc.gmt",package="scMetabolism")
  countexp2 <- countexp
  if(imputation){
    library(alra)
    countexp2 <- alra(as.matrix(countexp))[[3]]
  }
  if(method=="AUCell"){
    library(AUCell);library(GSEABase)
    cells_rankings <- AUCell_buildRankings(as.matrix(countexp2), nCores=ncores, plotStats=F)
    cells_AUC <- AUCell_calcAUC(getGmt(gmtFile), cells_rankings)
    signature_exp <- data.frame(getAUC(cells_AUC))
  }
  obj@assays$METABOLISM$score <- signature_exp
  obj
}
module2 <- readRDS("./module2/module2.rds")
module5 <- readRDS("./module5/module5.rds")
module2 <- subset(module2, subcluster %in% c("Astrocyte_01_CAV1","Astrocyte_02_MALAT1","Astrocyte_03_ROM1","Astrocyte_04_SCGB3A1","Astrocyte_05_HLA-DRB5","Astrocyte_06_CCL5","B_04_PTMAP2","CD4T_04_SNORD3A","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67","GE1","GE2","GE8","GE9","GE11","Macrophage_06_C3","Neurophil_01_GOS2"))
module5 <- subset(module5, subcluster %in% c("Endothelial_01_FCN3","Endothelial_03_VCAM1","Endothelial_04_DKK2","Endothelial_06_CD68","Endothelial_08_CCL21","Fibroblast_06_Alveolar_MS4A7","GE1","GE12","Macrophage_01_LGMN","Macrophage_03_MARCO","Macrophage_05_CSTB","Mast_02_CHIT1"))
module2$type <- ifelse(grepl("^GE",module2$subcluster),"module2.epi","module2.TIME")
module5$type <- ifelse(grepl("^GE",module5$subcluster),"module5.epi","module5.TIME")
ss <- merge(module2, module5)
ss <- JoinLayers(ss)
ss <- sc.metabolism.Seurat(ss, method="AUCell", imputation=F, ncores=2)
df <- as.data.frame(t(ss@assays$METABOLISM$score))
df$type <- ss$type
avg_df <- aggregate(.~type, data=df, mean)
rownames(avg_df) <- avg_df$type
avg_df <- as.data.frame(t(avg_df[,-1]))
pdf("all.metabolism.pdf", width=6, height=13)
pheatmap(avg_df,
         scale="row",
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         border_color=NA,
         color=colorRampPalette(c("#40A429","white","#825695"))(200),
         fontsize=6)
dev.off()
