# ------------- Figure 2 --------------
#----Figure 2A----
#co-occurrence analyses
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
setwd("./NSCLC/Figure/figure2/")
rm(list = ls())
TIME <- readRDS(file = "E:/NSCLC/Figure/figure1/TIME.rds")
sce <- TIME
#remove cell populations with low cell counts, low frequencies, or unidentified cell populations
length(table(sce$subcluster))
freq.low <- c("Unidentified")
freq.high <- sce$subcluster[!sce$subcluster %in% freq.low]
sce <- sce[,!sce$subcluster %in% freq.low]
#Calculate the frequency of cell subsets; the row labels are the samples, and the column labels are the cell subsets
write.csv(prop.table(table(sce$orig.ident,sce$subcluster),1),file = "./cell frequency cor.csv")
data.freq <- data.frame(read.csv(file="./cell frequency cor.csv"),row.names = 1)
data.freq <- data.freq[,!colnames(data.freq) %in% freq.low]#Remove cell clusters with low cell counts or frequencies
rowSums(data.freq)#Determine the proportion of the cell population in the sample
#Calculating Correlation
library(corrplot)
library(psych)
library(pheatmap)
#①corr.test(pearson)
data_test <- corr.test(data.freq, method= "pearson")#spearman/kendall
data_cor <- data_test[["r"]]
write.csv(data_cor,file = "./correlation of TIME cells.csv")
#pheatmap
library(circlize)
col_fun <- colorRamp2(
  c( -0.5,0, 1), 
  c("#638fa9","white","#b10928")
  #c("#01665e","white", "#8c510a")
)
pheatmap(data_cor,cluster_rows = T,cluster_cols = T,fontsize = 7,main = "correlation of TIME cells",border_color = NA,color = col_fun)


#Perform HR analysis on the various cell subsets in TIME
library(GSVA)
library(TCGAplot)
library(tidyverse)
library(survival)
library(meta)
library(doParallel)
#Set the number of parallel processing cores
registerDoParallel(cores = 15)
subcluster.marker <- read.csv(file="./yang/wang/TIME.markers1_TOP10.csv",row.names = 1)
tpm <- get_all_tpm()
meta <- get_all_meta()
cancers <- c("LUAD","LUSC")
genelist <- split(subcluster.marker$gene, subcluster.marker$cluster)
genelist <- genelist[-c(16,60)]
# Use a `foreach` loop to parallelize the outer loop
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
                         
                         # GSVA
                         gsvapar <- gsvaParam(exprData = exprSet, 
                                              geneSets = current_genelist, 
                                              kcdf = "Gaussian")
                         exprSet_gsva <- gsva(gsvapar)
                         
                         # Survival Analysis
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
                         
                         # Pan-cancer meta analysis
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
subcluster.hr <- subset(final_pan,cancer %in%c("LUAD","LUSC"))
subcluster.hr <- subcluster.hr %>% mutate(cellmodule = case_when(
  celltype %in% c("Bn_06_GZMB","pDC_03_GZMB","Plasma_09_IGHG_RRM2","Plasma_05_IGHG_IGHGP","Plasma_01_IGHA_MZB1","Plasma_02_IGHA_IGLC3","CD8T_02_GZMB",
                  "CD8T_03_MKI67","CD8T_01_CCL5","Treg_FOXP3","B_02_HLA.DQB1","CD4T_02_GNB2L1","CD4T_06_ACSL5","CD8T_04_HLA-B") ~ "module1",
  celltype %in% c("Astrocyte_03_ROM1","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67","B_04_PTMAP2","CD4T_04_SNORD3A","Astrocyte_05_HLA-DRB5","Neurophil_01_GOS2","Astrocyte_04_SCGB3A1","Macrophage_06_C3","Astrocyte_01_CAV1","Astrocyte_02_MALAT1","Astrocyte_06_CCL5") ~ "module2",
  celltype %in% c("NK_01_FGFBP2","CD4T_03_SCGB3A1","Endothelial_05_IL1RL1","Mast_01_GNB2L1","Mast_03_PBK","cDC_01_CLEC10A","cDC_02_CD207","Macrophage_02_FCN1","Macrophage_04_SPP1",
                  "Macrophage_09_MKI67","CD4T_05_SFTPB","Plasma_06_IGHG_SFTPB","Endothelial_09_S100A14","Macrophage_08_SCGB3A2","Fibroblast_04_Alveolar_DPT","Macrophage_07_IGKC") ~ "module3",
  celltype %in% c("CD4T_01_IL7R","B_03_CD3E","B_05_FYN","Fibroblast_05_Alveolar_MYH11","Endothelial_07_NOTCH3","B_01_HLA-DRA","Plasma_07_IGHG_GRIFIN","Plasma_08_IGHG_GUSBP11","Plasma_03_IGHA_GUSBP11","Plasma_04_IGHA_VPS13D") ~ "module4",
  celltype %in% c("Endothelial_03_VCAM1","Endothelial_08_CCL21","Endothelial_01_FCN3","Endothelial_04_DKK2","Fibroblast_06_Alveolar_MS4A7","Endothelial_06_CD68","Mast_02_CHIT1","Macrophage_03_MARCO","Macrophage_01_LGMN","Macrophage_05_CSTB") ~ "module5",
))

subcluster.hr$cellmodule <-  factor(subcluster.hr$cellmodule,levels = c("module1","module2","module3",'module4','module5'))
subcluster.hr <- subcluster.hr[order(subcluster.hr$cellmodule),]

write.csv(subcluster.hr,file = "/data/yang/wang//subcluster.hr.csv")

#----Figure 2B----
#Analysis of Gene Differential Expression Across Cell Modules
#Remove cell clusters that are less relevant and lie outside the module
rm(list = ls())
TIME <- readRDS(file = "./NSCLC/Figure/figure1/TIME.rds")
subsce.corr <- TIME
#corlow <- c("")
#subsce.corr <- subsce.corr[,!subsce.corr$subcluster %in% freq.low]
subsce.corr@meta.data <- subsce.corr@meta.data %>% mutate(cellmodule = case_when(
  subcluster %in% c("Bn_06_GZMB","pDC_03_GZMB","Plasma_09_IGHG_RRM2","Plasma_05_IGHG_IGHGP","Plasma_01_IGHA_MZB1","Plasma_02_IGHA_IGLC3","CD8T_02_GZMB",
                    "CD8T_03_MKI67","CD8T_01_CCL5","Treg_FOXP3","B_02_HLA.DQB1","CD4T_02_GNB2L1","CD4T_06_ACSL5","CD8T_04_HLA-B") ~ "module1",
  subcluster %in% c("Astrocyte_03_ROM1","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67","B_04_PTMAP2","CD4T_04_SNORD3A","Astrocyte_05_HLA-DRB5","Neurophil_01_GOS2","Astrocyte_04_SCGB3A1","Macrophage_06_C3","Astrocyte_01_CAV1","Astrocyte_02_MALAT1","Astrocyte_06_CCL5") ~ "module2",
  subcluster %in% c("NK_01_FGFBP2","CD4T_03_SCGB3A1","Endothelial_05_IL1RL1","Mast_01_GNB2L1","Mast_03_PBK","cDC_01_CLEC10A","cDC_02_CD207","Macrophage_02_FCN1","Macrophage_04_SPP1",
                    "Macrophage_09_MKI67","CD4T_05_SFTPB","Plasma_06_IGHG_SFTPB","Endothelial_09_S100A14","Macrophage_08_SCGB3A2","Fibroblast_04_Alveolar_DPT","Macrophage_07_IGKC") ~ "module3",
  subcluster %in% c("CD4T_01_IL7R","B_03_CD3E","B_05_FYN","Fibroblast_05_Alveolar_MYH11","Endothelial_07_NOTCH3","B_01_HLA-DRA","Plasma_07_IGHG_GRIFIN","Plasma_08_IGHG_GUSBP11","Plasma_03_IGHA_GUSBP11","Plasma_04_IGHA_VPS13D") ~ "module4",
  subcluster %in% c("Endothelial_03_VCAM1","Endothelial_08_CCL21","Endothelial_01_FCN3","Endothelial_04_DKK2","Fibroblast_06_Alveolar_MS4A7","Endothelial_06_CD68","Mast_02_CHIT1","Macrophage_03_MARCO","Macrophage_01_LGMN","Macrophage_05_CSTB") ~ "module5",
))
table(subsce.corr$cellmodule)
getwd()
saveRDS(subsce.corr,file = "./NSCLC/Figure/figure2/subsce.corr.rds")
subsce.corr <- readRDS("./subsce.corr.rds")
#Identify the differentially upregulated genes in each cell module
DefaultAssay(subsce.corr) <- "RNA"
Idents(subsce.corr) <- "cellmodule"
subsce.corr.markers.up  <- FindAllMarkers(subsce.corr, 
                                          min.pct = 0.25,
                                          only.pos = T,
                                          logfc.threshold = 0.25)
subsce.corr.markers1.up  <- subsce.corr.markers.up [subsce.corr.markers.up $p_val_adj < 0.05, ]
subsce.corr.markers1.up_TOP5 <- subsce.corr.markers1.up %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(subsce.corr.markers1.up,file = "./subsce.corr.markers1.up.csv")
subsce.corr.markers1.up <- read.csv("./subsce.corr.markers1.up.csv",row.names = 1,header = T)
subsce.corr.markers1.up$cluster <- factor(subsce.corr.markers1.up$cluster,levels = c("module1","module2","module3",'module4','module5'))
subsce.corr.markers1.top30 <- subsce.corr.markers1.up %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
subsce.corr.markers1.top30 <- subsce.corr.markers1.top30[order(subsce.corr.markers1.top30$cluster), ]
write.csv(subsce.corr.markers1.top30,file = "./NSCLC/Figure/figure2/subsce.corr.markers1.top30.csv")
#install.packages('devtools')
#devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
#DEG
DefaultAssay(subsce.corr) <- "RNA"
subsce.corr.markers  <- FindAllMarkers(subsce.corr, 
                                       min.pct = 0.25,
                                       only.pos = F,
                                       logfc.threshold = 0.25)
subsce.corr.markers1  <- subsce.corr.markers [subsce.corr.markers $p_val_adj < 0.05, ]
subsce.corr.markers1_TOP5 <- subsce.corr.markers1 %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.csv(subsce.corr.markers1,file = "./subsce.corr.markers1.csv")
#subsce.corr.markers1 <- read.csv("./subsce.corr.markers1.csv",row.names = 1)
subsce.corr.markers1$cluster <- factor(subsce.corr.markers1$cluster,levels = c("module1","module2","module3",'module4','module5'))
#save(subsce.corr,file = "F:/subsce.corr.RData")
mycolors <- c("#f9b3ad","#80C1D7","#70cdbe","#f1dba6","#CEB6E0")
mycolors <- c("#B1C1E3","#92D3EB","#C9E5D5","#E2EDBC","#FAE7F1")
mycolors <- c("#9B7EBD","#70cdbe","#46aec8","#6fc08d","#457baa")
mycolors <- c("#d72e2d","#ea9740","#316339","#7dbcd1","#313663")
mycolors <- c("#a5c469","#d586b3","#8ea1c9","#dd9f85","#88bdaa")
mycolors <- c("#70cdbe","#80C1D7","#968ab7","#d586b3","#eb7e60")
jjVolcano(diffData = subsce.corr.markers1,
          log2FC.cutoff = 0.25, 
          size  = 2.5, 
          fontface = 'italic', 
          aesCol = c('#8fb4dc','#eb7e60'), 
          tile.col = mycolors, 
          #col.type = "adjustP", 
          topGeneN = 5, 
          fontface = 'italic',
          polar = T)+ ylim(-7,12)

color =  colorRampPalette(c("white","#e6573f", "#d23229"))(100)
jjVolcano(diffData = subsce.corr.markers1,
          tile.col = mycolors,
          size  = 3.5,
          pSize = 0.9,
          topGeneN = 5,
          base_size = 20,
          fontface = 'italic',
          aesCol = c('#8fb4dc','#d23229'), 
          polar = T) + ylim(-12,10)#+xlim(-10,10)

#----Figure 2C----
#Functional enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(enrichplot)
group <- data.frame(
  gene = subsce.corr.markers.up$gene,
  module = subsce.corr.markers.up$cluster
)
#SYMBOL → ENTREZID
gene_id <- bitr(group$gene,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

mydf <- merge(gene_id, group, by.x = "SYMBOL", by.y = "gene")
#KEGG
kk <- compareCluster(
  ENTREZID ~ module,
  data = mydf,
  fun = "enrichKEGG",
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
dotplot(kk, showCategory = 10) +
  scale_shape_manual(values = 18)

#----Figure 2D----
#Gene set scoring
subsce.corr <- readRDS("./subsce.corr.rds")
DefaultAssay(subsce.corr) <- "RNA"
#Antigen-presenting cell co-inhibitory genes
APC_co_inhibition <- list(c("C10orf54","CD274","LGALS9","PDCD1LG2","PVRL3"))
#Co-stimulatory genes of antigen-presenting cells
APC_co_stimulation <- list(c("CD40","CD58","CD70","ICOSLG","SLAMF1","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF8","TNFSF9"))
#CC chemokine receptor genes
CCR <- list(c("CCL16","TPO","TGFBR2","CXCL2","CCL14","TGFBR3","IL11RA","CCL11","IL4I1","IL33","CXCL12","CXCL10","BMPER","BMP8A","CXCL11","IL21R","IL17B","TNFRSF9","ILF2","CX3CR1","CCR8","TNFSF12","CSF3","TNFSF4","BMP3","CX3CL1","BMP5","CXCR2","TNFRSF10D","BMP2","CXCL14",
              "CCL28","CXCL3","BMP6","CCL21","CXCL9","CCL23","IL6","TNFRSF18","IL17RD","IL17D","IL27","CCL7","IL1R1","CXCR4","CXCR2P1","TGFB1I1","IFNGR1","IL9R","IL1RAPL1","IL11","CSF1","IL20RA","IL25","TNFRSF4","IL18","ILF3","CCL20","TNFRSF12A","IL6ST","CXCL13","IL12B","TNFRSF8",
              "IL6R","BMPR2","IFNE","IL1RAPL2","IL3RA","BMP4","CCL24","TNFSF13B","CCR4","IL2RA","IL32","TNFRSF10C","IL22RA1","BMPR1A","CXCR5","CXCR3","IFNA8","IL17REL","IFNB1","IFNAR1","TNFRSF1B","CCL17","IFNL1","IL16","IL1RL1","ILK","CCL25","ILDR2","CXCR1","IL36RN","IL34","TGFB1","IFNG",
              "IL19","ILKAP","BMP2K","CCR10","ILDR1","EPO","CCR7","IL17C","IL23A","CCR5","IL7","EPOR","CCL13","IL2RG","IL31RA","TNFAIP6","IFNL2","BMP1","IL12RB1","TNFAIP8","IL4R","TNFRSF6B","TNFAIP8L1","TNFRSF10B","IFNL3","CCL5","CXCL6","CXCL1","CCR3","TNFSF11","CSF1R","IL21","IL1RAP","IL12RB2",
              "CCL1","IL17RA","CCR1","IL1RN","TNFRSF11B","TNFRSF14","IL13","IL2RB","BMP8B","CCL2","IL24","IL18RAP","TGFBI","TNFSF10","TNFRSF11A","CXCL5","IL5RA","TNFSF9","IL1RL2","TNFRSF13C","IL36G","IL15RA","TNFRSF21","CXCL8","IL22RA2","TNFAIP8L2","IL18R1","IFNLR1","CXCR6","CCL3L3","TNFRSF1A","IL17RE","IFNGR2",
              "IL17RC","TNFAIP8L3","ILVBL","TGFBRAP1","CCL4L1","CSF2RA","CCRN4L","CCL26","TNFAIP1","CCRL2","IFNA10","TNFRSF17","IFNA13","IL20","IL18BP",'CCL3L1',"TNFSF12-TNFSF13","IL5","IL23R","IL26","TNF","TGFA","CSF2","IL1F10","CXCL17","TNFSF13","IFNA4","IL37","IL12A","IL7R","IFNA1","IL1A","IL4","IL2",
              "CCL22","CSF3R","IL10","IFNK","TGFB2","IL1R2","IL1B","IL17F","IL27RA","IL15","TNFSF8","IL36B","XCL1","CXCL16","TNFRSF19","IL3","CCL3","IFNA2","BMPR1B","IFNA21","TNFSF18","CCL8","IL17RB","TNFRSF25","IL22","IL10RB","IFNAR2","CCL18","IFNA16","CSF2RB","IL36A","TNFAIP3","IL13RA2","IL13RA1","CCR9","TNFRSF10A",
              "IFNA7","IFNW1","XCL2","TNFSF14","CCR2","BMP15","BMP10","CCL15-CCL14","TGFBR1","IFNA5","BMP7","IFNA14","IL20RB","IL10RA","IFNA17","CCR6","TGFB3","CCL15","CCL4","CCL27","TNFRSF13B","TNFAIP2","IL31","IL17A","TNFSF15","CCL19","IFNA6","IL9"))
#Checkpoint genes
Check_point <- list(c("IDO1","LAG3","CTLA4","TNFRSF9","ICOS","CD80","PDCD1LG2","TIGIT","CD70","TNFSF9","ICOSLG","KIR3DL1","CD86","PDCD1","LAIR1","TNFRSF8","TNFSF15","TNFRSF14","IDO2","CD276","CD40","TNFRSF4","TNFSF14","HHLA2","CD244",
                      "CD274","HAVCR2","CD27","BTLA","LGALS9","TMIGD2","CD28","CD48","TNFRSF25","CD40LG","ADORA2A","VTCN1","CD160","CD44","TNFSF18","TNFRSF18","BTNL2","C10orf54","CD200R1","TNFSF4","CD200","NRP1"))
#Genes encoding cell lysis activity
Cytolytic_activity <- list(c("PRF1","GZMA"))

#Human leukocyte antigen genes
HLA <- list(c("HLA-E","HLA-DPB2","HLA-C","HLA-J","HLA-DQB1","HLA-DQB2","HLA-DQA2","HLA-DQA1","HLA-A","HLA-DMA","HLA-DOB","HLA-DRB1","HLA-H","HLA-B","HLA-DRB5","HLA-DOA","HLA-DPB1","HLA-DRA","HLA-DRB6","HLA-L","HLA-F","HLA-G","HLA-DMB","HLA-DPA1"))

#pro-inflammatory genes
Inflammation_promoting <- list(c("CCL5","CD19","CD8B","CXCL10","CXCL13","CXCL9","GNLY","GZMB","IFNG","IL12A","IL12B","IRF1","PRF1","STAT1","TBX21"))

#Major Histocompatibility Complex Class I Gene
MHC_class_I <- list(c("B2M","HLA-A","TAP1"))

#anti-inflammatory genes
Parainflammation <- list(c("CXCL10","PLAT","CCND1","LGMN","PLAUR","AIM2","MMP7","ICAM1","MX2","CXCL9","ANXA1","TLR2","PLA2G2D","ITGA2","MX1","HMOX1",
                           "CD276","TIRAP","IL33","PTGES","TNFRSF12A","SCARB1","CD14","BLNK","IFIT3",'RETNLB',"IFIT2","ISG15","OAS2","REL","OAS3","CD44","PPARG","BST2","OAS1","NOX1","PLA2G2A","IFIT1",'IFITM3',"IL1RN"))

#T-cell co-inhibitory genes
T_cell_co_inhibition <- list(c("CD160","CD244","CD274","CTLA4","HAVCR2","LAG3","LAIR1","TIGIT"))

#T-cell costimulation
T_cell_co_stimulation <- list(c("CD2","CD226","CD27","CD28","CD40LG","ICOS","SLAMF1","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF14")) 

#Type I interferon
Type_I_IFN_Reponse <- list(c("DDX4","IFIT1","IFIT2","IFIT3","IRF7","ISG20","MX1","MX2","RSAD2","TNFSF10"))

#Type II interferon
Type_II_IFN_Reponse <- list(c("GPR146","SELP","AHR"))

#T-cell status genes
T_naiveness <- list(c("CCR6","CCR7","TCF7","SELL","LEF1","IL7R"))
T_cytotoxicity <- list(c("GZMA","GZMB","GZMH","GZMK","PRF1","CXCL13","GNLY","IFNG","NKG7"))
T_exhaustion <- list(c("PDCD1","CTLA4","HAVCR2","LAG3","TIGIT","LAYN"))
T_proliferation <- list(c("MKI67","TOP2A","CD3D","CD3E"))

#gene list
library(readxl)
genelist <- "./GeneList1.xlsX"
genes <- read_excel(genelist)
#Antigen Presentation and Processing Gene Set
Antigen_Processing_and_Presentation <- subset(genes,Category=="Antigen_Processing_and_Presentation") 
Antigen_Processing_and_Presentation <- list(Antigen_Processing_and_Presentation$Symbol)

#TCR signaling
TCRsignalingPathway <- subset(genes,Category=="TCRsignalingPathway") 
TCRsignalingPathway <- list(TCRsignalingPathway$Symbol)

#TNF (Tumor Necrosis Factor) gene set
TNF_Family_Members <- subset(genes,Category==c("TNF_Family_Members"))
TNF_Family_Members_Receptors <- subset(genes,Category==c("TNF_Family_Members_Receptors"))
TNF_Family <- rbind(TNF_Family_Members,TNF_Family_Members_Receptors)
TNF_Family <- list(TNF_Family)

#Natural killer cell cytotoxicity
NaturalKiller_Cell_Cytotoxicity <- subset(genes,Category=="NaturalKiller_Cell_Cytotoxicity")
NaturalKiller_Cell_Cytotoxicity <- list(NaturalKiller_Cell_Cytotoxicity)

#Cytokines
Cytokines <- subset(genes,Category=="Cytokines")
Cytokines <- list(Cytokines)

#BCRSignalingPathway
BCRSignalingPathway <- subset(genes,Category=="BCRSignalingPathway")
BCRSignalingPathway <- list(BCRSignalingPathway)

#Interleukins
Interleukins <- subset(genes,Category=="Interleukins")
Interleukins <- list(Interleukins)

#TGFb_Family_Member
TGFb_Family_Member <- subset(genes,Category=="TGFb_Family_Member")
TGFb_Family_Member <- list(TGFb_Family_Member)

#TIMELASER_marker
TIMELASER_marker <- read.csv(file = "./TIMELASER_marker.csv")
TIME_IA <- subset(TIMELASER_marker, TIMELASER.subtype=="TIME-IA")
IA.gene <- list(TIME_IA$Gene)
TIME_ISM <- subset(TIMELASER_marker, TIMELASER.subtype=="TIME-ISM")
ISM.gene <- list(TIME_ISM$Gene)
TIME_ISS <- subset(TIMELASER_marker, TIMELASER.subtype=="TIME-ISS")
ISS.gene <- list(TIME_ISS$Gene)
TIME_IE <- subset(TIMELASER_marker, TIMELASER.subtype=="TIME-IE")
IE.gene <- list(TIME_IE$Gene)
TIME_IR <- subset(TIMELASER_marker, TIMELASER.subtype=="TIME-IR")
IR.gene <- list(TIME_IR$Gene)

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

DefaultAssay(subsce.corr) <- "RNA"
gene_sets <- list(
  IA = IA.gene,
  ISM = ISM.gene,
  ISS = ISS.gene,
  IE = IE.gene,
  IR = IR.gene,
  APC_co_inhibition = APC_co_inhibition,
  APC_co_stimulation = APC_co_stimulation,
  CCR = CCR,
  Check_point = Check_point,
  Cytolytic_activity = Cytolytic_activity,
  HLA = HLA,
  Inflammation_promoting = Inflammation_promoting,
  MHC_class_I = MHC_class_I,
  Parainflammation = Parainflammation,
  T_cell_co_inhibition = T_cell_co_inhibition,
  T_cell_co_stimulation = T_cell_co_stimulation,
  Type_I_IFN_Reponse = Type_I_IFN_Reponse,
  Type_II_IFN_Reponse = Type_II_IFN_Reponse,
  T_naiveness = T_naiveness,
  T_cytotoxicity = T_cytotoxicity,
  T_exhaustion = T_exhaustion,
  T_proliferation = T_proliferation,
  Antigen_Processing_and_Presentation = Antigen_Processing_and_Presentation,
  TCRsignalingPathway = TCRsignalingPathway,
  TNF_Family = TNF_Family,
  NaturalKiller_Cell_Cytotoxicity = NaturalKiller_Cell_Cytotoxicity,
  Cytokines = Cytokines,
  BCRSignalingPathway = BCRSignalingPathway,
  Interleukins = Interleukins,
  TGFb_Family_Member = TGFb_Family_Member
)
#AddModuleScore
for (nm in names(gene_sets)) {
  subsce.corr <- AddModuleScore(
    subsce.corr,
    features = gene_sets[[nm]],
    ctrl = 100,
    name = "score"
  )
}
# score1~scoreN
score_cols <- grep("^score", colnames(subsce.corr@meta.data), value = TRUE)
colnames(subsce.corr@meta.data)[match(score_cols, colnames(subsce.corr@meta.data))] <- names(gene_sets)
plot_list <- list()
mycolors <- c("#70cdbe","#80C1D7","#968ab7","#d586b3","#eb7e60")
for (nm in names(gene_sets)) {
  
  df <- subsce.corr@meta.data %>%
    mutate(cellmodule = factor(cellmodule))
  
  # Kruskal-Wallis
  kw <- kruskal.test(df[[nm]] ~ df$cellmodule)
  
  p_label <- paste0(
    "Kruskal-Wallis: χ²(", kw$parameter, ") = ",
    format(kw$statistic, digits = 3),
    ifelse(kw$p.value < 0.001,
           ", p < 0.001",
           paste0(", p = ", format(kw$p.value, digits = 3)))
  )
  p <- ggplot(df, aes(x = cellmodule, y = .data[[nm]], fill = cellmodule)) +
    geom_violin(width = 0.8, alpha = 0.85, scale = "width") +
    geom_boxplot(width = 0.15, fill = "white", outlier.size = 1) +
    scale_fill_manual(values = mycolors) +
    theme_bw() +
    annotate(
      "text",
      x = 3,
      y = Inf,
      label = p_label,
      vjust = 1.5,
      size = 4,
      fontface = "italic"
    ) +
    labs(title = nm)
  
  plot_list[[nm]] <- p
}

#----Figure 2E---- 
#Determine the TIME subtype for each patient
write.csv(prop.table(table(subsce.corr$subcluster,subsce.corr$orig.ident),1),file = "./patient.type.csv")
ALL <- data.frame(read.csv(file="./NSCLC/Figure/figure2/patient.type.csv"))#,row.names = 1))
ALL <- ALL %>% mutate(cellmodule = case_when(
  X %in% c("Bn_06_GZMB","pDC_03_GZMB","Plasma_09_IGHG_RRM2","Plasma_05_IGHG_IGHGP","Plasma_01_IGHA_MZB1","Plasma_02_IGHA_IGLC3","CD8T_02_GZMB",
           "CD8T_03_MKI67","CD8T_01_CCL5","Treg_FOXP3","B_02_HLA.DQB1","CD4T_02_GNB2L1","CD4T_06_ACSL5","CD8T_04_HLA.B") ~ "module1",
  X %in% c("Astrocyte_03_ROM1","Endothelial_02_PTMAP2","Fibroblast_01_Myo_HIGD1B","Fibroblast_02_Myo_COL10A1","Fibroblast_03_Myo_MKI67","B_04_PTMAP2","CD4T_04_SNORD3A",
           "Astrocyte_05_HLA.DRB5","Neurophil_01_GOS2","Astrocyte_04_SCGB3A1","Macrophage_06_C3","Astrocyte_01_CAV1","Astrocyte_02_MALAT1","Astrocyte_06_CCL5") ~ "module2",
  X %in% c("NK_01_FGFBP2","CD4T_03_SCGB3A1","Endothelial_05_IL1RL1","Mast_01_GNB2L1","Mast_03_PBK","cDC_01_CLEC10A","cDC_02_CD207","cDC_02_CD207","Macrophage_02_FCN1","Macrophage_04_SPP1",
           "Macrophage_09_MKI67","CD4T_05_SFTPB","Plasma_06_IGHG_SFTPB","Endothelial_09_S100A14","Macrophage_08_SCGB3A2","Fibroblast_04_Alveolar_DPT","Macrophage_07_IGKC") ~ "module3",
  X %in% c("CD4T_01_IL7R","B_03_CD3E","B_05_FYN","Fibroblast_05_Alveolar_MYH11","Endothelial_07_NOTCH3","B_01_HLA.DRA","Plasma_07_IGHG_GRIFIN","Plasma_08_IGHG_GUSBP11","Plasma_03_IGHA_GUSBP11","Plasma_04_IGHA_VPS13D") ~ "module4",
  X %in% c("Endothelial_03_VCAM1","Endothelial_08_CCL21","Endothelial_01_FCN3","Endothelial_04_DKK2","Fibroblast_06_Alveolar_MS4A7","Endothelial_06_CD68","Mast_02_CHIT1","Macrophage_03_MARCO","Macrophage_01_LGMN","Macrophage_05_CSTB") ~ "module5",
))

which.max ()
rownames(ALL) <- ALL$X
ALL$X <-NULL

library(stringr)
module1 <- subset(ALL,cellmodule=="module1")
module1$cellmodule <- NULL
module_freq <- data.frame(colSums(module1))
module2 <- subset(ALL,cellmodule=="module2")
module2$cellmodule <- NULL
module_freq[[2]] <- data.frame(colSums(module2))
module3 <- subset(ALL,cellmodule=="module3")
module3$cellmodule <- NULL
module_freq[[3]] <- data.frame(colSums(module3))
module4 <- subset(ALL,cellmodule=="module4")
module4$cellmodule <- NULL
module_freq[[4]] <- data.frame(colSums(module4))
module5 <- subset(ALL,cellmodule=="module5")
module5$cellmodule <- NULL
module_freq[[5]] <- data.frame(colSums(module5))
module_freq[[6]] <- data.frame(rowSums(module_freq))
colnames(module_freq) <- c("module1","module2","module3","module4","module5","sum_freq")
module_freq$module1 <- module_freq$module1/module_freq$sum_freq
module_freq$module2 <- module_freq$module2/module_freq$sum_freq
module_freq$module3 <- module_freq$module3/module_freq$sum_freq
module_freq$module4 <- module_freq$module4/module_freq$sum_freq
module_freq$module5 <- module_freq$module5/module_freq$sum_freq
module_freq[[6]] <- NULL
module_freq$patient <- rownames(module_freq)
row_anno <- module_freq$patient
row_anno <- data.frame(row_anno)
colnames(row_anno) <- "patient"
#row_anno$patient <- rownames(row_anno)
#row_anno$patient <- NULL
row_anno <- row_anno %>% mutate(mutation = case_when(
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
rownames(row_anno) <- row_anno$patient
row_anno$mutation<- factor(row_anno$mutation)
rownames(module_freq) <- module_freq$patient
row_anno$patient <- NULL
module_freq$patient <- NULL
#factor(row_anno$mutation)
row_anno <- data.frame(t(row_anno))
module_freq <- data.frame(t(module_freq))
#row_anno$mutation = list(c("EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR","EGFR-co-mutation","EGFR-co-mutation","EGFR-co-mutation","KRAS","KRAS","KRAS","KRAS","KRAS","KRAS","KRAS","KRAS","KRAS","KRAS","KRAS","KRAS","KRAS","KRAS-co-mutation","KRAS-co-mutation","KRAS-co-mutation","ALK","ROS1","ROS1","ROS1","TP53","MET","BRAF"))
module_freq
colSums(module_freq)
row_anno <- data.frame(t(row_anno))
row_anno$mutation <- factor(row_anno$mutation,levels = c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53","MET-BM","HER2"))
library("pheatmap")
ann_colors=list(class=c(L='#009933',R='#CC33CC',F='#FDDCA9'))
mycol2 = list(mutation = c("EGFR"="#c6b7d4","EGFR-BM"="#d44e26","EGFR-co-mutation"="#e3a264","KRAS"="#6fc2d0","KRAS-co-mutation"="#6f9abf","ALK"="#a5c49b","ROS1"="#7266ac","TP53"="#FF9966","MET-BM"="#d84986","HER2"="#2d588e"))
#names(mycol2) <- unique(row_anno$mutation)
module_freq <- data.frame(t(module_freq))
patientorder <- c("P22","P43","P9","P1","P4","P25","P10","P27",
                  "P18","P15","P38","P14","P44","P12","P13","P3",
                  "P16","P42","P2","P23","P8","P41","P26","P20","P17","P45","P21","P36","P7","P11","P19","P39","P34","P5",
                  "P31","P28","P30","P29","P35","P32",
                  "P24","P40","P37","P33","P6")
df_sorted <- module_freq[patientorder,]


df_sorted <- data.frame(t(df_sorted))
#color =  colorRampPalette(c("white","#e6573f", "#d23229"))(100)
color = colorRampPalette(c("white","#0172be"))(100)
pheatmap(df_sorted,
         color = color,
         annotation_col  = row_anno,
         annotation_colors = mycol2,
         display_numbers = F,
         border="white", 
         cluster_cols = F, 
         cluster_rows = F,
         treeheight_col = 40, 
         treeheight_row = 45,
         border_color = NA)            

#----Figure 2F---- 
#Plotting the Organizational Preference Map for the TIME Subtype
subsce.corr <- readRDS(file = "./NSCLC/Figure/figure2/subsce.corr.rds")
ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}
divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    } 
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}
metadata <- subsce.corr@meta.data

meta_filt <- metadata[metadata$mutation %in% c("ALK","EGFR","EGFR-BM","EGFR-co-mutation","HER2","KRAS","KRAS-co-mutation","MET-BM","ROS1","TP53"),]
meta_filt$mutation <- factor(as.vector(meta_filt$mutation),levels=c("ALK","EGFR","EGFR-BM","EGFR-co-mutation","HER2","KRAS","KRAS-co-mutation","MET-BM","ROS1","TP53"))
summary <- table(meta_filt[,c('cellmodule','mutation')])

# ro/e
roe <- as.data.frame(ROIE(summary))
library("pheatmap")
pheatmap(roe, display_numbers = TRUE,number_color = "black",cluster_row = FALSE,cluster_col = FALSE,
         color = colorRampPalette(c("#e7eb8e","#cd782d","firebrick3"))(50)) 

#Create a lollipop chart showing organizational preferences for the TIME subtype
roe0 <- as.data.frame(t(roe))
roe0$group <- rownames(roe0)

library(ggpubr)
mycol2 = c("#c6b7d4","#d44e26","#e3a264","#6fc2d0","#6f9abf","#a5c49b","#7266ac","#FF9966","#d84986","#2d588e")
#module1
ggdotchart(roe0, x = "group", y = "module1",
           color = "group",   # Color by groups
           sorting = "descending",  # Sort value in descending order
           add = "segments", # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 11, # Large dot size
           label = round(roe0$module1,2),  # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),  # Adjust label parameters
           ggtheme = theme_pubr())+ 
  geom_hline(yintercept = 1, linetype = 2, color = "black")+
  theme(legend.position='none')+ylab("Ro/e")+ggtitle("Ro/e-module1")+
  theme(axis.text.x = element_text(angle = 45,vjust = 1))+
  scale_color_manual(values = c("EGFR"="#c6b7d4","EGFR-BM"="#d44e26","EGFR-co-mutation"="#e3a264","KRAS"="#6fc2d0","KRAS-co-mutation"="#6f9abf","ALK"="#a5c49b","ROS1"="#7266ac","TP53"="#FF9966","MET-BM"="#d84986","HER2"="#2d588e"))

#module2
ggdotchart(roe0, x = "group", y = "module2",
           color = "group",   # Color by groups
           sorting = "descending",  # Sort value in descending order
           add = "segments", # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 11, # Large dot size
           label = round(roe0$module2,2),  # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),  # Adjust label parameters
           ggtheme = theme_pubr())+ 
  geom_hline(yintercept = 1, linetype = 2, color = "black")+
  theme(legend.position='none')+ylab("Ro/e")+ggtitle("Ro/e-module2")+
  theme(axis.text.x = element_text(angle = 45,vjust = 1))+
  scale_color_manual(values = c("EGFR"="#c6b7d4","EGFR-BM"="#d44e26","EGFR-co-mutation"="#e3a264","KRAS"="#6fc2d0","KRAS-co-mutation"="#6f9abf","ALK"="#a5c49b","ROS1"="#7266ac","TP53"="#FF9966","MET-BM"="#d84986","HER2"="#2d588e"))

#module3
ggdotchart(roe0, x = "group", y = "module3",
           color = "group",   # Color by groups
           sorting = "descending",  # Sort value in descending order
           add = "segments", # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 11, # Large dot size
           label = round(roe0$module3,2),  # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),  # Adjust label parameters
           ggtheme = theme_pubr())+ 
  geom_hline(yintercept = 1, linetype = 2, color = "black")+
  theme(legend.position='none')+ylab("Ro/e")+ggtitle("Ro/e-module3")+
  theme(axis.text.x = element_text(angle = 45,vjust = 1))+
  scale_color_manual(values = c("EGFR"="#c6b7d4","EGFR-BM"="#d44e26","EGFR-co-mutation"="#e3a264","KRAS"="#6fc2d0","KRAS-co-mutation"="#6f9abf","ALK"="#a5c49b","ROS1"="#7266ac","TP53"="#FF9966","MET-BM"="#d84986","HER2"="#2d588e"))

#module4
ggdotchart(roe0, x = "group", y = "module4",
           color = "group",   # Color by groups
           sorting = "descending",  # Sort value in descending order
           add = "segments", # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 11, # Large dot size
           label = round(roe0$module4,2),  # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),  # Adjust label parameters
           ggtheme = theme_pubr())+ 
  geom_hline(yintercept = 1, linetype = 2, color = "black")+
  theme(legend.position='none')+ylab("Ro/e")+ggtitle("Ro/e-module4")+
  theme(axis.text.x = element_text(angle = 45,vjust = 1))+
  scale_color_manual(values = c("EGFR"="#c6b7d4","EGFR-BM"="#d44e26","EGFR-co-mutation"="#e3a264","KRAS"="#6fc2d0","KRAS-co-mutation"="#6f9abf","ALK"="#a5c49b","ROS1"="#7266ac","TP53"="#FF9966","MET-BM"="#d84986","HER2"="#2d588e"))

#module5
ggdotchart(roe0, x = "group", y = "module5",
           color = "group",   # Color by groups
           sorting = "descending",  # Sort value in descending order
           add = "segments", # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 11, # Large dot size
           label = round(roe0$module5,2),  # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),  # Adjust label parameters
           ggtheme = theme_pubr())+ 
  geom_hline(yintercept = 1, linetype = 2, color = "black")+
  theme(legend.position='none')+ylab("Ro/e")+ggtitle("Ro/e-module5")+
  theme(axis.text.x = element_text(angle = 45,vjust = 1))+
  scale_color_manual(values = c("EGFR"="#c6b7d4","EGFR-BM"="#d44e26","EGFR-co-mutation"="#e3a264","KRAS"="#6fc2d0","KRAS-co-mutation"="#6f9abf","ALK"="#a5c49b","ROS1"="#7266ac","TP53"="#FF9966","MET-BM"="#d84986","HER2"="#2d588e"))


#Plotting a scatter plot of gene set scores
library(ggplot2)
library(ggsankey)
subsce.corr <- readRDS("./NSCLC/Figure/figure2/subsce.corr.rds")
DefaultAssay(subsce.corr) <- "RNA"
immune.inhibitor <- list(c("ADORA2A","BTLA","CD160","CD244","CD274","CD96","CSF1R","CTLA4","HAVCR2","IDO1","IL10","IL10RB","KDR","KIR2DL1","KIR2DL3","LAG3","LGALS9","PDCD1","PDCD1LG2","TGFB1","TGFBR1","TIGIT","VTCN1"))
Inscore <- AddModuleScore(subsce.corr,
                          features = immune.inhibitor,
                          ctrl = 100,
                          name = "AUCELL")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[28] <- 'immune.inhibitor'

immune.stimulator <- list(c("CD27","CD276","CD28","CD40","CD40LG","CD48","CD70","CD80","CD86","CXCL12","CXCR4","ENTPD1","HHLA2","ICOS","ICOSLG","IL2RA","IL6","IL6R","KLRC1","KLRK1","LTA","MICB","NT5E","PVR","RAET1E","TMIGD2","TNFRSF13B","TNFRSF13C","TNFRSF14","TNFRSF17","TNFRSF18","TNFRSF25","TNFRSF4","TNFRSF8","TNFRSF9","TNFSF13","TNFSF13B","TNFSF14","TNFSF15","TNFSF18","TNFSF4","TNFSF9","ULBP1"))
Inscore <- AddModuleScore(Inscore,
                          features = immune.stimulator,
                          ctrl = 100,
                          name = "AUCELL")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[29] <- 'immune.stimulator'

Chemokine.receptor <- list(c("CCR1","CCR2","CCR3","CCR4","CCR5","CCR6","CCR7","CCR8","CCR9","CCR10", "CXCR1","CXCR2","CXCR3","CXCR4","CXCR5","CXCR6","XCR1","CX3R1"))
Inscore <- AddModuleScore(Inscore,
                          features = Chemokine.receptor,
                          ctrl = 100,
                          name = "AUCELL")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[30] <- 'Chemokine.receptor'

Chemokine <- list(c("CCL1","CCL2","CCL3","CCL4","CCL5","CCL7","CCL8","CCL11","CCL13","CCL14","CCL15","CCL16","CCL17","CCL18","CCL19","CCL20","CCL21","CCL22","CCL23","CCL24","CCL25","CCL26","CCL28","CX3CL1","CXCL1","CXCL2","CXCL3","CXCL5","CXCL6","CXCL8","CXCL9","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14","CXCL16","CXCL17"))
Inscore <- AddModuleScore(Inscore,
                          features = Chemokine,
                          ctrl = 100,
                          name = "AUCELL")
colnames(Inscore@meta.data)
colnames(Inscore@meta.data)[31] <- 'Chemokine'

#----Figure 2G----
#The 14 malicious pathways in the decouple package are scored on CM1-CM5
library(Seurat)
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyverse)
corsce <- readRDS("./yang/wang/NSCLC/subsce.corr.rds")
net <- get_progeny(organism = 'human', top = 500)
mat <- as.matrix(corsce@assays$RNA@layers$data)
rownames(mat) <- rownames(corsce)
colnames(mat) <- colnames(corsce)
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                
                .mor='weight', minsize = 5)

acts
corsce[['pathwaysmlm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject()

DefaultAssay(corsce) <- 'pathwaysmlm'

corsce <- ScaleData(corsce)
corsce@assays$pathwaysmlm@data <- corsce@assays$pathwaysmlm@scale.data
#module
Idents(corsce) <- "cellmodule"
df <- t(as.matrix(corsce@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(cellmodule = Idents(corsce)) %>%  
  pivot_longer(cols = -cellmodule, names_to = "pathway", values_to = "score") %>%
  group_by(cellmodule, pathway) %>%
  summarise(mean_score = mean(score), .groups = "drop")

selected_pathways <- c("Androgen","EGFR","Estrogen","Hypoxia","JAK-STAT","MAPK","NFkB","p53","PI3K","TGFb","TNFa","Trail","VEGF","WNT")

df_filtered <- df %>%
  filter(pathway %in% selected_pathways)

top_acts_mat <- df_filtered %>%
  pivot_wider(names_from = pathway, values_from = mean_score) %>%
  column_to_rownames('cellmodule') %>%
  as.matrix()
library(pheatmap)
module.top_acts_mat <- top_acts_mat
pheatmap(module.top_acts_mat,
         scale = "row", 
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("#A5CCE0", "white", "#d66c61"))(100),
         main = "Pathway Activity across Cell Modules")


#mutation
Idents(corsce) <- "mutation"
df <- t(as.matrix(corsce@assays$pathwaysmlm@data)) %>%
  as.data.frame() %>%
  mutate(mutation = Idents(corsce)) %>%
  pivot_longer(cols = -mutation, names_to = "pathway", values_to = "score") %>%
  group_by(mutation, pathway) %>%
  summarise(mean_score = mean(score), .groups = "drop")

selected_pathways <- c("Androgen","EGFR","Estrogen","Hypoxia","JAK-STAT","MAPK","NFkB","p53","PI3K","TGFb","TNFa","Trail","VEGF","WNT") 

df_filtered <- df %>%
  filter(pathway %in% selected_pathways)

top_acts_mat <- df_filtered %>%
  pivot_wider(names_from = pathway, values_from = mean_score) %>%
  column_to_rownames('mutation') %>%
  as.matrix()
library(pheatmap)

pheatmap(top_acts_mat,
         scale = "row", 
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border="white",
         border_color = NA,
         color = colorRampPalette(c("#A5CCE0", "white", "#d66c61"))(100),
         main = "Pathway Activity across Cell Modules")

df.module <- data.frame(module.top_acts_mat)
df.mutation <- data.frame(top_acts_mat)
write.csv(df.module,file = "/data/yang/wang/NSCLC/df.module.csv")
write.csv(df.mutation,file = "/data/yang/wang/NSCLC/df.mutation.csv")

#----Figure 2H----
#Drawing Sankey
#Classify regulon_group expression values into high and low categories
df <- Inscore@meta.data[,c(4,27:31)]
df$score <- ifelse(df$immune.inhibitor >= median(df$immune.inhibitor, na.rm = TRUE), "High", "Low")
df$score <- ifelse(df$immune.stimulator >= median(df$immune.stimulator, na.rm = TRUE), "High", "Low")
df$score <- ifelse(df$Chemokine.receptor >= median(df$Chemokine.receptor, na.rm = TRUE), "High", "Low")
df$score <- ifelse(df$Chemokine >= median(df$Chemokine, na.rm = TRUE), "High", "Low")
write.csv(df,file = "./NSCLC/Figure/figure2/df.csv")
library(ggsankey)
library(tidyverse)

level <- unique(df$score)
mutation <- unique(df$mutation)
module <- unique(df$cellmodule)
level_colors <- c(
  "Low"    = "#ADD4A5",
  "High"      = "#C9B8DA"
)
mutation_colors <- c(
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

module_colors <- c(
  "module1" = "#70cdbe",
  "module2" = "#80C1D7",
  "module3" = "#968ab7",
  "module4" = "#d586b3",
  "module5" = "#eb7e60")
level <- names(level_colors)
mutation <- names(mutation_colors)

module <- names(module_colors)

color_palette <- c(module_colors,level_colors,mutation_colors)
df$cellmodule <- factor(df$cellmodule,levels = c("module1","module2","module3","module4","module5"))
df_sankey <- df %>%
  select(cellmodule,score, mutation) %>%
  make_long(cellmodule,score, mutation)

ggplot(df_sankey, aes(x = x,
                      next_x = next_x,
                      node = node,
                      next_node = next_node,
                      fill = node,
                      label = node)) +
  scale_fill_manual(values = color_palette) +
  geom_sankey(flow.alpha = 0.6, smooth = 6, width = 0.15) +
  geom_sankey_text(size = 3, color = "black") +
  theme_void() +
  theme(legend.position = "none")