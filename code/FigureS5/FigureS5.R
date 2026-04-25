# ------------- Figure S5 --------------
#----Figure S5A----
module2 <- readRDS("./yang/wang/figure5/module2.rds")
module5 <- readRDS("./yang/wang/figure5/module5.rds")
M2.epi <- subset(module2,CELLTYPE %in% c("Epithelial"))
M5.epi <- subset(module5,CELLTYPE %in% c("Epithelial"))
M2.5.epi <- merge(M2.epi,M5.epi)
M2.5.epi <- JoinLayers(M2.5.epi)
Idents(M2.5.epi) <- "cellmodule"
net <- get_collectri(organism='human', split_complexes=FALSE)
mat <- GetAssayData(M2.5.epi, slot = "data")
mat <- as.matrix(mat)
plan("multisession",workers = 5)
acts <- run_ulm(mat, net,minsize = 5,
                .source='source', .target='target',.mor='mor')
M2.5.epi[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', 
              names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
DefaultAssay(object = M2.5.epi) <- "tfsulm"
options(future.globals.maxSize = 50 * 1024^3)  
M2.5.epi <- ScaleData(M2.5.epi)
M2.5.epi@assays$tfsulm@data <- M2.5.epi@assays$tfsulm@scale.data

sce <- M2.5.epi
#subcluster
Idents(sce) <- sce@meta.data$cellmodule

n_tfs_per_group <- 20  
df <- t(as.matrix(sce@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(sce)) %>%
  pivot_longer(cols = -cluster, 
               names_to = "source",
               values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

#Identify the most active TFs in CMs
top_tfs_list <- df %>%
  group_by(cluster) %>%  
  arrange(desc(mean), .by_group = TRUE) %>% 
  slice_head(n = n_tfs_per_group) %>%  
  ungroup() %>%
  pull(source) %>% 
  unique()  


#Construct a heatmap matrix
top_acts_mat <- df %>%
  filter(source %in% top_tfs_list) %>%
  pivot_wider(id_cols = 'cluster', 
              names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()

palette_length = 100
my_color = colorRampPalette(c("#0F7B9F","#418AAA","#8CB4C9","#ABC7D6","white","#FEC8B8","#F4896C","#DF482C","#D72F13"))(palette_length)

my_breaks <- c(seq(-1, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 1, length.out=floor(palette_length/2)))

pdf(file = "module.epi.tf.pdf",height = 13,width = 6)
pheatmap(t(top_acts_mat), 
         border_color = NA,
         color=my_color, 
         breaks = my_breaks,
         angle_col = 45)
dev.off()
saveRDS(M2.5.epi,file = "./yang/wang/figure5/M2.5.epi.rds")

#module2.epi
tf_data <- t(as.matrix(sce@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(sce)) %>%
  pivot_longer(cols = -cluster, 
               names_to = "TF",
               values_to = "score") %>%
  group_by(cluster, TF) %>%
  summarise(mean_score = mean(score)) %>%
  filter(mean_score > 0) 

M2 <- subset(tf_data,cluster == "module2")
M2 <- M2[order(M2$mean_score,decreasing = T),]
M2$TF <- factor(M2$TF, levels = M2$TF, ordered = TRUE)
M2$top5 <- ifelse(rank(-M2$mean_score) <= 5, "Top5", "Other")
M2$top10 <- ifelse(rank(-M2$mean_score) <= 10, "Top10", "Other")
library(ggplot2)
library(ggrepel)

ggplot(data = M2, aes(x = reorder(TF, -mean_score), y = mean_score)) +
  geom_point(aes(color = top10), size = 3) +
  geom_point(data = subset(M2, top10 == "Top10"), color = "red", size = 3) +
  geom_text_repel(
    data = subset(M2, top10 == "Top10"),
    aes(label = TF),
    color = "red",
    size = 4,
    nudge_x = 0.3,
    segment.color = "red",
    segment.size = 0.2,
    max.overlaps = Inf, 
    box.padding = 0.5,   
    point.padding = 0.5  
  ) +
  scale_color_manual(values = c("Top10" = "red", "Other" = "gray60"), guide = "none") +
  labs(x = "Transcription Factor (TF)", y = "Module2.EPI Expression") +
  theme(
    text = element_text(family = "Helvetica", size = 12),
    axis.title = element_text(size = 12),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.line = element_line(linewidth = 0.2, color = "black"),
    axis.ticks = element_line(linewidth = 0.2, color = "black"),
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.1),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 20, 10, 10)  
  )

#module5.epi
M5 <- subset(tf_data,cluster == "module5")
M5 <- M5[order(M5$mean_score,decreasing = T),]
M5$TF <- factor(M5$TF, levels = M5$TF, ordered = TRUE)
M5$top5 <- ifelse(rank(-M5$mean_score) <= 5, "Top5", "Other")
M5$top10 <- ifelse(rank(-M5$mean_score) <= 10, "Top10", "Other")
library(ggplot2)
library(ggrepel) 

ggplot(data = M5, aes(x = reorder(TF, -mean_score), y = mean_score)) +
  geom_point(aes(color = top10), size = 3) +
  geom_point(data = subset(M5, top10 == "Top10"), color = "red", size = 3) +
  geom_text_repel(
    data = subset(M5, top10 == "Top10"),
    aes(label = TF),
    color = "red",
    size = 4,
    nudge_x = 0.3,
    segment.color = "red",
    segment.size = 0.2,
    max.overlaps = Inf,  
    box.padding = 0.5,  
    point.padding = 0.5 
  ) +
  scale_color_manual(values = c("Top10" = "red", "Other" = "gray60"), guide = "none") +
  labs(x = "Transcription Factor (TF)", y = "Module5.EPI Expression") +
  theme(
    text = element_text(family = "Helvetica", size = 12),
    axis.title = element_text(size = 12),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.line = element_line(linewidth = 0.2, color = "black"),
    axis.ticks = element_line(linewidth = 0.2, color = "black"),
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.1),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 20, 10, 10)  
  )

#----Figure S5B----
rm(list=ls())
net <- get_collectri(organism='human', split_complexes=FALSE)
M2.5.epi <- readRDS(file = "./yang/wang/figure5/M2.5.epi.rds")
sce <- M2.5.epi
#----subcluster
Idents(sce) <- sce@meta.data$cellmodule
df <- t(as.matrix(sce@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(sce)) %>%
  pivot_longer(cols = -cluster, 
               names_to = "source",
               values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))

df <- subset(df,mean >0)
table(df$cluster)

R2 <- subset(df, cluster== "module2")
R5 <- subset(df, cluster== "module5")
library(data.table)
net2 <- merge(net, R2["source"], by = "source")
net5 <- merge(net, R5["source"], by = "source")

#Read in the key CNV genes associated with cancer cell states from Result 3
GE1 <- read.csv(file = "./yang/wang/figure5/CNV gene/GE1.csv",row.names = 1)
GE5<- read.csv(file = "./yang/wang/figure5/CNV gene/GE5.csv",row.names = 1)
GE8 <- read.csv(file = "./yang/wang/figure5/CNV gene/GE8.csv",row.names = 1)
GE9 <- read.csv(file = "./yang/wang/figure5/CNV gene/GE9.csv",row.names = 1)
GE <- rbind(GE1,GE5)
GE <- rbind(GE,GE8)
GE <- rbind(GE,GE9)
GE$CNV <- ifelse(GE$cnv>0,"up","down")
sorted_GE <- GE[order(GE$cnv), ]
subset_GE <- rbind(head(sorted_GE, 10), tail(sorted_GE, 10))
mutation <- subset(GE,rownames(GE) %in% c("EGFR","MET","KRAS","TP53","BRAF","ALK","ROS1","HER2"))
ge <- rbind(subset_GE,mutation)
ge <- ge[,c(5,6)]
ge <- unique(ge)

library(stringr)
GE.name <- ge$name

#Extract the regulatory network of CNV-related genes.
net2.hr <- subset(net2,target %in% GE.name)
net5.hr <- subset(net5,target %in% GE.name)
write.csv(net2.hr,file = "/data/yang/wang/figure5/module2/net2.epi.csv")
write.csv(net5.hr,file = "/data/yang/wang/figure5/module5/net5.epi.csv")
#The network diagram was generated using Cytoscape software.

#----Figure S5C----
#ROC
library(timeROC)
library(ROCit)
#train
train$risk.score <- risk_score 
roc_data <- timeROC(T = train$time,
                    delta = train$event,
                    marker = train$risk.score,
                    cause=1, times = c(12, 24, 36))  # 1year、2year、3year

plot(roc_data, time=12, title=FALSE)
lines(roc_data$FP[,2], roc_data$TP[,2], col="blue")  
lines(roc_data$FP[,3], roc_data$TP[,3], col="green") 
legend("bottomright", 
       legend=paste(c("1-Year","2-Year","3-Year"), 
                    "AUC =", round(roc_data$AUC,2)),
       col=c("red","blue","green"), lty=1)

#test
library(timeROC)
roc_data <- timeROC(T = test$time,
                    delta = test$event,
                    marker = test$risk.score,
                    cause=1, times = c(12, 24, 36))   # 1year、2year、3year

plot(roc_data, time=12, title=FALSE)
lines(roc_data$FP[,2], roc_data$TP[,2], col="blue") 
lines(roc_data$FP[,3], roc_data$TP[,3], col="green") 
legend("bottomright", 
       legend=paste(c("1-Year","2-Year","3-Year"), 
                    "AUC =", round(roc_data$AUC,2)),
       col=c("red","blue","green"), lty=1)

#----Figure S5D----
library(survival)
library(ggrisk)
library(ROCit)
genes <- c("PAX6-ITGB1","CEBPB-IGF1","CEBPB-RETN","NR1H4-RARRES2","PPARG-RARRES2","RELA-IGF1","SNAI1-ITGB4","SP1-ITGA6","USF2-IGF1")
#train
fit_train <- coxph(Surv(time,event)~age+gender+stage,data=train)
ggrisk(fit_train,cutoff.value="median",cutoff.y=2,heatmap.genes=genes,color.A=c(low="#2E9FDF",high="#E7B800"),color.B=c(code.0="#2E9FDF",code.1="#E7B800"),color.C=c(low="#2E9FDF",median="white",high="#E7B800"))
ggrisk(fit_train,cutoff.value="roc",cutoff.y=2,heatmap.genes=genes,color.A=c(low="#2E9FDF",high="#E7B800"),color.B=c(code.0="#2E9FDF",code.1="#E7B800"),color.C=c(low="#2E9FDF",median="white",high="#E7B800"))
#test
fit_test <- coxph(Surv(time,event)~age+gender+stage,data=test)
test$risk.score <- risk_score
ggrisk(fit_test,cutoff.value="median",cutoff.y=2,heatmap.genes=genes,color.A=c(low="#2E9FDF",high="#E7B800"),color.B=c(code.0="#2E9FDF",code.1="#E7B800"),color.C=c(low="#2E9FDF",median="white",high="#E7B800"))
ggrisk(fit_test,cutoff.value="roc",cutoff.y=2,heatmap.genes=genes,color.A=c(low="#2E9FDF",high="#E7B800"),color.B=c(code.0="#2E9FDF",code.1="#E7B800"),color.C=c(low="#2E9FDF",median="white",high="#E7B800"))

#----Figure S5E----
#Sankey diagram plotting
#Classify regulon_group expression values into high and low groups.
df <- Mscore@meta.data[,c(1,3,4,5,8:17)]
df$rs_group <- ifelse(df$rs_score >= median(df$rs_score, na.rm = TRUE), "High", "Low")
write.csv(df,file = "/data/yang/wang/figure5/df.csv")
library(ggsankey)
library(tidyverse)
level <- unique(df_sankey$node[df_sankey$x == "rs_group"])

level_colors <- c(
  "Low"    = "#ADD4A5",
  "High"      = "#C9B8DA"
)

subcluster <- unique(df_sankey$node[df_sankey$x == "subcluster"])
subcluster_colors <- c(
  "Astrocyte_01_CAV1" = "#f1b38a",
  "Astrocyte_02_MALAT1" = "#f0db69",
  "Astrocyte_03_ROM1" = "#a2c246",
  "Astrocyte_04_SCGB3A1" = "#f5cee0",
  "Astrocyte_05_HLA-DRB5" = "#b6d2b7",
  "Astrocyte_06_CCL5" = "#cbd2e5",
  "B_04_PTMAP2" = "#ad98c3",
  "CD4T_04_SNORD3A" = "#b4ffde",
  "Endothelial_01_FCN3" = "#E5D2DD",
  "Endothelial_02_PTMAP2" = "#53A85F",
  "Endothelial_03_VCAM1" = "#F1BB72",
  "Endothelial_04_DKK2" = "#F3B1A0",
  "Endothelial_06_CD68" = "#57C3F3",
  "Endothelial_08_CCL21" = "#E95C59",
  "Fibroblast_01_Myo_HIGD1B" = "#E5D2DD",
  "Fibroblast_02_Myo_COL10A1" = "#F1BB72",
  "Fibroblast_03_Myo_MKI67" = "#F3B1A0",
  "Fibroblast_06_Alveolar_MS4A7" = "#F7C394",
  "Macrophage_01_LGMN" = "#93b8db",
  "Macrophage_03_MARCO" = "#7eb4c6",
  "Macrophage_05_CSTB" = "#f5e886",
  "Macrophage_06_C3" = "#d6d5b7",
  "Mast_02_CHIT1" = "#D6E7A3",
  "Neurophil_01_GOS2" = "#F1BB72"
)

subcluster <- names(subcluster_colors)

module <- unique(df_sankey$node[df_sankey$x == "cellmodule"])
cellmodule_color <- c("module2" = "#80C1D7","module5" ="#eb7e60")
module <- names(cellmodule_color)

color_palette <- c(level_colors, subcluster_colors, mutation_colors,cellmodule_color)
df_sankey <- df %>%
  select(rs_group, subcluster, mutation,cellmodule) %>%
  make_long(rs_group, subcluster, mutation,cellmodule)

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
