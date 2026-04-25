# ------------- Figure 5 --------------
#----Figure 5A----
library(Seurat)
library(dplyr)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)
library(decoupleR)
library(ggplot2)
library(patchwork)
library(OmnipathR)
# Infer transcription factor activity in CM2 and CM5 cells
module2 <- readRDS("./data/yang/wang/figure5/module2.rds")
module5 <- readRDS("./data/yang/wang/figure5/module5.rds")
M2 <- subset(module2,CELLTYPE %in% c("Astrocyte","B","CD4T","Endothelial","Fibroblast","Macrophage","Neurophil"))
M5 <- subset(module5,CELLTYPE %in% c("Endothelial","Fibroblast","Macrophage","Mast"))
M2.5 <- merge(M2,M5)
M2.5 <- JoinLayers(M2.5)

Idents(M2.5) <- "cellmodule"
net <- get_collectri(organism='human', split_complexes=FALSE)
mat <- GetAssayData(M2.5, slot = "data")
mat <- as.matrix(mat)
plan("multisession",workers = 5)
acts <- run_ulm(mat, net,minsize = 5,
                .source='source', .target='target',.mor='mor')


M2.5[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', 
              names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
DefaultAssay(object = M2.5) <- "tfsulm"
options(future.globals.maxSize = 50 * 1024^3)  
M2.5 <- ScaleData(M2.5)
M2.5@assays$tfsulm@data <- M2.5@assays$tfsulm@scale.data

sce <- M2.5
#----subcluster
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



# Identify top active transcription factors in each group
top_tfs_list <- df %>%
  group_by(cluster) %>% 
  arrange(desc(mean), .by_group = TRUE) %>% 
  slice_head(n = n_tfs_per_group) %>% 
  ungroup() %>%
  pull(source) %>% 
  unique() 

top_acts_mat <- df %>%
  filter(source %in% top_tfs_list) %>%
  pivot_wider(id_cols = 'cluster', 
              names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()
saveRDS(M2.5,file = "./data/yang/wang/figure5/M2.5.rds")

#module2
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
  
  labs(x = "Transcription Factor (TF)", y = "Module2 Expression") +
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

#module5
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
  
  labs(x = "Transcription Factor (TF)", y = "Module5 Expression") +
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

#----Figure 5B----
# Extract top TFs from module 2 and module 5 for network construction
# Select TFs significantly enriched in module 2 and module 5 by chi-square test
tf <- df$source
tf <- intersect(tf, rownames(M2.5))
bulk_matrix <- 
  Seurat::AggregateExpression(
    object = M2.5,
    assays = "RNA",               
    features = tf,                  
    group.by = "cellmodule",       
    slot = "counts",              
    # slot = "data",              
    return.seurat = FALSE,        
    aggregation.fun = sum          
    # aggregation.fun = mean     
  )$RNA
bulk_matrix <- as.data.frame(bulk_matrix)
mat <- bulk_matrix[, c("module2", "module5")]
results <- t(apply(mat, 1, function(x) {
  x <- as.numeric(x)
  total <- sum(x)
  
  
  if(total == 0) {
    return(c(chi_sq = NA, p_val = NA))
  }
  
  
  test <- chisq.test(x, p = p_expected)
  
  
  c(chi_sq = test$statistic, p_val = test$p.value)
}))

results_df <- data.frame(
  gene = rownames(mat),
  module2_count = mat[,"module2"],
  module5_count = mat[,"module5"],
  chi_sq = results[, 1],
  p_val = results[, 2],
  stringsAsFactors = FALSE
)

results_df$adj_pval <- p.adjust(results_df$p_val, method = "BH")

results_df <- results_df[order(-results_df$chi_sq), ]
results_df <- subset(results_df,p_val <0.05)
write.csv(results_df,file = "./data/yang/wang/figure5/results_tf.csv")

head(results_df)
df <- subset(df,source %in% results_df$gene)

R2 <- subset(df, cluster== "module2")
R5 <- subset(df, cluster== "module5")
library(data.table)
net2 <- merge(net, R2["source"], by = "source")
net5 <- merge(net, R5["source"], by = "source")


load("./data/yang/wang/figure5/JS_M_HR.RData")
M2 <- subset(M2_HR_result,cancer == c("LUAD","LUSC"))
M2 <- subset(M2, p <0.05)
M2$LR <- gsub("\\(|\\)","",M2$celltype)
M2$LR <- gsub(" ", "",M2$LR)
M2$LR <- gsub("-","_",M2$LR)
M2$LR <- gsub("\\+","_",M2$LR)
library(stringr)
M2.name <- str_split(M2$LR, "_", simplify = FALSE)
M2.name <- unlist(M2.name)
M2.name <- unique(M2.name)

M5 <- subset(M5_HR_result,cancer == c("LUAD","LUSC"))
M5 <- subset(M5, p <0.05)
M5$LR <- gsub("\\(|\\)","",M5$celltype)
M5$LR <- gsub(" ", "",M5$LR)
M5$LR <- gsub("-","_",M5$LR)
M5$LR <- gsub("\\+","_",M5$LR)
M5.name <- str_split(M5$LR, "_", simplify = FALSE)
M5.name <- unlist(M5.name)
M5.name <- unique(M5.name)

net2.hr <- subset(net2,target %in% M2.name)
net5.hr <- subset(net5,target %in% M5.name)

data <- net2.hr
colnames(data) <- c("TFs","Genes","mor")

edges <- data.frame(from = net2.hr$source,to = net2.hr$target)
show <- c(unique(net2.hr$source),unique(net2.hr$target))
vertices <- data.frame(name = unique(c(net2.hr$source,net2.hr$target)), type = 'Targets', show_name = NA, size = 1)

vertices$type[match(unique(net2.hr$source),vertices$name) ] <- 'TFs'
vertices$color <- vertices$type
vertices$color[match(show,vertices$name) ] <- 'black'
vertices$size[match(unique(net2.hr$source),vertices$name) ] <- 2
vertices$show_name[match(show,vertices$name) ] <- show

ggraph_data <- igraph::graph_from_data_frame(d = edges, 
                                             vertices = vertices,
                                             directed = T)


#----Figure 5C----
# HR analysis of top regulators in module 2 and module 5
library(GSVA)
library(TCGAplot)
library(tidyverse)
library(survival)
library(meta)
library(doParallel) 

registerDoParallel(cores = 15)

#----module2
net2.hr <- read.csv("./data/yang/wang/figure5/module2/net2.hr.csv",row.names = 1)
net5.hr <- read.csv("./data/yang/wang/figure5/module5/net5.hr.csv",row.names = 1)

tpm <- get_all_tpm()
meta <- get_all_meta()
cancers <- unique(meta$Cancer)
dt <- net2.hr
dt$name <- paste0(dt$source,sep = "-",dt$target)
a <- apply(dt, 1, function(row) c(row["source"], row["target"]))
genelist <- lapply(seq_len(ncol(a)), function(i) a[, i])
names(genelist) <- dt$name

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
save(genelist,file = "./data/yang/wang/figure5/module2/net2.hr.list.RData")


write.csv(final_pan, file = "./data/yang/wang/figure5/module2/TCGAplot_net2.tf_HR.csv")
final_result1 <- final_pan[final_pan$cancer != "Pan cancer", ]
final_result1 <- final_result1 %>% filter(p < 0.05)
lung <- final_result1 %>% filter(cancer == c("LUAD","LUSC"))
table(final_result1$color) 
# Better survival  Worse survival 
# 9               12 
write.csv(final_result1, file = "./data/yang/wang/figure5/module2/TCGAplot_net2.tf.sigure_HR.csv")

# LASSO regression on HR-related regulators to build a prognostic model
# Train–validation split (70%/30%)

setwd("./data/yang/wang/NSCLC/figure5/lassco/")
library(TCGAplot)
load("./data/yang/wang/figure5/genelist.RData")
tpm <- get_all_tpm()
meta <- get_all_meta()
lung.tpm <- subset(tpm,Cancer %in% c("LUAD","LUSC"))
lung.tpm <- subset(lung.tpm,Group =="Tumor")
lung.meta <- subset(meta,Cancer %in% c("LUAD","LUSC"))
lung.meta <- na.omit(lung.meta)
lung.tpm <- lung.tpm[match(rownames(lung.meta),substr(rownames(lung.tpm),1,12)),]
rownames(lung.tpm) <- substr(rownames(lung.tpm),1,12)

# Calculate regulon activity scores in the lung TPM dataset
library(GSVA)
tpm <- lung.tpm[,-c(1,2)]
tpm <- t(tpm)
param <- gsvaParam(exprData = tpm, 
                   geneSets = genelist, 
                   kcdf = "Gaussian")

gsva_result <- gsva(param) 
gsva_result <- t(gsva_result)
rownames(gsva_result) <- gsub("\\.","-",rownames(gsva_result))
tpm <- gsva_result[match(rownames(lung.meta),rownames(gsva_result)),]

set.seed(321)
ind <- sample(nrow(tpm), nrow(tpm) * 0.7)
train.tpm <- tpm[ind, ]
test.tpm <- tpm[-ind, ]
train.meta <- lung.meta[ind, ]
test.meta <- lung.meta[-ind, ]
save(train.tpm,train.meta,test.tpm,test.meta,file = "./data/yang/wang/figure5/lassco/data.RData")
save(train.tpm,file = "./data/yang/wang/figure5/lassco/data.train.tpm.RData")

#Perform Cox regression analysis on regulons in the training dataset.
library(survival)
td <- cbind(train.meta,train.tpm)
pFilter=0.05
outResult=data.frame()
sigGenes=c("event","time")
for(i in colnames(td[,7:ncol(td)])){
  tdcox <- coxph(Surv(time, event) ~ td[,i], data = td)
  tdcoxSummary = summary(tdcox)
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"]
  if(pvalue<pFilter){
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}


library(survival)
train.tpm <- train.tpm[,colnames(train.tpm) %in% outResult$id]
train <- cbind(train.meta,train.tpm)
surv_obj <- with(train, Surv(time, event))

library(glmnet)

x <- train.tpm

fit1 <- glmnet(x, surv_obj, family = "cox", alpha = 1)
plot(fit1, xvar = "lambda", label = TRUE)

fit <- cv.glmnet(x, surv_obj, family = "cox", alpha = 1)
coef(fit,s=fit$lambda.min)
plot(fit, xvar = "lambda", label = TRUE)

#----Figure 5D----
best_lambda <- fit$lambda.min

coef_lasso <- coef(fit1, s = best_lambda)
selected <- coef_lasso[coef_lasso[,1] != 0, , drop = FALSE]
save(selected,file = "./data/yang/wang/figure5/lassco/selected.coef.RData")
selected_features <- rownames(selected)
print(selected_features)
coeffident <- as.matrix(selected)
coeffident <- data.frame(coeffident)
coeffident$regulor <- rownames(coeffident)
write.csv(coeffident,file = "./data/yang/wang/figure5/lassco/coeffident.csv")

risk_score <- x[, selected_features] %*% selected[selected_features, 1]

risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")

surv_fit <- survfit(surv_obj ~ risk_group)

library(survminer)
ggsurvplot(surv_fit, data = data.frame(risk_group),
           pval = TRUE, risk.table = TRUE,conf.int = T,palette = 
             c("#E7B800", "#2E9FDF"),
           title = "LASSO Prognostic Model based on TF-target pairs-train")

# Validate model performance in the validation dataset
test.tpm <- test.tpm[,colnames(test.tpm) %in% outResult$id]
test <- cbind(test.meta,test.tpm)
surv_obj2 <- with(test, Surv(time, event))
x1 <- test.tpm

test$min<-predict(fit1,newx = x1,s=fit$lambda.min)



# Calculate the C-index in the validation dataset
model1<-coxph(Surv(time,event)~min,data=test)
summary(model1)

risk_score <- x1[, selected_features] %*% selected[selected_features, 1]
save(c(selected,selected_features),file = "./data/yang/wang/figure5/lassco/Risk score coefficients.RData")
risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")

surv_fit2 <- survfit(surv_obj2 ~ risk_group)

library(survminer)
ggsurvplot(surv_fit2, data = data.frame(risk_group),
           pval = TRUE, risk.table = TRUE,conf.int = T,palette = 
             c("#E7B800", "#2E9FDF"),
           title = "LASSO Prognostic Model based on TF-target pairs-test")

#----Figure 5E----
# Plot a forest plot of model genes in the training dataset
train.tpm1 <- train.tpm[,colnames(train.tpm) %in% selected_features]
td <- cbind(train.meta,train.tpm1)
pFilter=1
outResult=data.frame()
sigGenes=c("event","time")
for(i in colnames(td[,7:ncol(td)])){
  tdcox <- coxph(Surv(time, event) ~ td[,i], data = td)
  tdcoxSummary = summary(tdcox)
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"]
  if(pvalue<pFilter){
    sigGenes=c(sigGenes,i)
    outResult=rbind(outResult,
                    cbind(id=i,
                          HR=tdcoxSummary$conf.int[,"exp(coef)"],
                          L95CI=tdcoxSummary$conf.int[,"lower .95"],
                          H95CI=tdcoxSummary$conf.int[,"upper .95"],
                          pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

all_vars <- c("age", "gender", "stage") 

univ_results <- data.frame()
for (var in all_vars) {
  fmla <- reformulate(termlabels = var, response = "Surv(time, event)")
  fit  <- coxph(fmla, data = td)
  su   <- summary(fit)
  
  result <- data.frame(
    Variable = var,
    HR       = su$conf.int[, "exp(coef)"],
    Lower95  = su$conf.int[, "lower .95"],
    Upper95  = su$conf.int[, "upper .95"],
    pValue   = su$coefficients[, "Pr(>|z|)"]
  )
  univ_results <- rbind(univ_results, result)
}

colnames(univ_results) <- colnames(outResult2)
outResult <- rbind(outResult,univ_results)
outResult$id[12:14] <- c("stageII","stageIII","stageIV") 

library(forestplot)
gene <- outResult[c(10:14,1:9),]
write.csv(gene,file = "./data/yang/wang/figure5/lassco/train.forest.csv")
gene$HR_CI <- sprintf("%.4s (%.4s-%.4s)", gene$HR, gene$L95CI, gene$H95CI)

gene$pvalue_formatted <- as.numeric(gene$pvalue_formatted)
gene$pvalue_formatted <- round(gene$pvalue_formatted,)
gene$pvalue_formatted <- sprintf("%.9f", gene$pvalue)

gene$significance <- ifelse(gene$pvalue < 0.001, "***",
                            ifelse(gene$pvalue < 0.01, "**",
                                   ifelse(gene$pvalue < 0.05, "*", "")))


dat <- rbind(
  c("Variable", NA, NA, NA, "HR (95% CI)", "P value", "Significance"),
  cbind(
    as.character(gene$id),
    gene$HR,
    gene$L95CI,
    gene$H95CI,
    gene$HR_CI,
    gene$pvalue_formatted,
    gene$significance
  )
)

forestplot(
  labeltext = dat[, c(1, 5, 6, 7)],  
  mean = as.numeric(c(NA, gene$HR)), 
  lower = as.numeric(c(NA, gene$L95CI)), 
  upper = as.numeric(c(NA, gene$H95CI)), 
  zero = 1,                          
  boxsize = 0.3,                     
  graph.pos = 3,                     
  xticks = c(0.5, 1, 1.5, 2, 2.5),   
  is.summary = c(TRUE, rep(FALSE, nrow(gene))),  
  
  txt_gp = fpTxtGp(
    label = gpar(cex = 1.1),       
    ticks = gpar(cex = 1.0),       
    xlab = gpar(cex = 1.1)         
  ),
  
  hrzl_lines = list(
    "1" = gpar(lty = 1, lwd = 2),   
    "2" = gpar(lty = 1, lwd = 2),   
    "16" = gpar(lty = 1, lwd = 2)    
  ),
  
  col = fpColors(
    box = "#FF0000",           
    lines = "black",        
    zero = "darkred"           
  ),
  
  lwd.zero = 2,        
  lwd.ci = 2,            
  lty.ci = 1,            
  ci.vertices = TRUE,    
  ci.vertices.height = 0.1, 
  
  title = "Gene Interaction Hazard Ratios",
  xlab = "Hazard Ratio (95% CI)",
  
  grid = structure(c(0.5, 1, 1.5, 2), 
                   gp = gpar(lty = 2, col = "#CCCCCC")),
  
  clip = c(0.5, 2.5)  
)

# Perform Cox regression analysis on regulons in the validation dataset
test.tpm1 <- test.tpm[,colnames(test.tpm) %in% selected_features]
td <- cbind(test.meta,test.tpm1)
td$gender <- as.factor(td$gender)
td$stage <- as.factor(td$stage)
pFilter=1
outResult2=data.frame()
sigGenes=c("event","time")
for(i in colnames(td[,7:ncol(td)])){
  tdcox <- coxph(Surv(time, event) ~ td[,i], data = td)
  tdcoxSummary = summary(tdcox)
  pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"]
  if(pvalue<pFilter){
    sigGenes=c(sigGenes,i)
    outResult2=rbind(outResult2,
                     cbind(id=i,
                           HR=tdcoxSummary$conf.int[,"exp(coef)"],
                           L95CI=tdcoxSummary$conf.int[,"lower .95"],
                           H95CI=tdcoxSummary$conf.int[,"upper .95"],
                           pvalue=tdcoxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}


all_vars <- c("age", "gender", "stage") 

univ_results <- data.frame()
for (var in all_vars) {
  fmla <- reformulate(termlabels = var, response = "Surv(time, event)")
  fit  <- coxph(fmla, data = td)
  su   <- summary(fit)
  
  result <- data.frame(
    Variable = var,
    HR       = su$conf.int[, "exp(coef)"],
    Lower95  = su$conf.int[, "lower .95"],
    Upper95  = su$conf.int[, "upper .95"],
    pValue   = su$coefficients[, "Pr(>|z|)"]
  )
  univ_results <- rbind(univ_results, result)
}

colnames(univ_results) <- colnames(outResult2)
outResult2 <- rbind(outResult2,univ_results)
outResult2$id[12:14] <- c("stageII","stageIII","stageIV") 

# Plot a forest plot of model genes in the validation dataset
library(forestplot)
gene2 <- outResult2[c(10:14,1:9),]
write.csv(gene2,file ="./data/yang/wang/figure5/lassco/test.forest.csv")
gene2$HR_CI <- sprintf("%.4s (%.4s-%.4s)", gene2$HR, gene2$L95CI, gene2$H95CI)

gene2$pvalue_formatted <- sprintf("%.7s", gene2$pvalue)

gene2$significance <- ifelse(gene2$pvalue < 0.001, "***",
                             ifelse(gene2$pvalue < 0.01, "**",
                                    ifelse(gene2$pvalue < 0.05, "*", "")))


dat <- rbind(
  c("Variable", NA, NA, NA, "HR (95% CI)", "P value", "Significance"),
  cbind(
    as.character(gene2$id),
    gene2$HR,
    gene2$L95CI,
    gene2$H95CI,
    gene2$HR_CI,
    gene2$pvalue_formatted,
    gene2$significance
  )
)

forestplot(
  labeltext = dat[, c(1, 5, 6, 7)],  
  mean = as.numeric(c(NA, gene2$HR)), 
  lower = as.numeric(c(NA, gene2$L95CI)), 
  upper = as.numeric(c(NA, gene2$H95CI)), 
  zero = 1,                          
  boxsize = 0.3,                     
  graph.pos = 3,                     
  xticks = c(0.5, 1, 1.5, 2, 2.5),   
  is.summary = c(TRUE, rep(FALSE, nrow(gene2))),  
  
  txt_gp = fpTxtGp(
    label = gpar(cex = 1.1),       
    ticks = gpar(cex = 1.0),       
    xlab = gpar(cex = 1.1)         
  ),
  
  hrzl_lines = list(
    "1" = gpar(lty = 1, lwd = 2),    
    "2" = gpar(lty = 1, lwd = 2),    
    "16" = gpar(lty = 1, lwd = 2)    
  ),
  
  col = fpColors(
    box = "#FF0000",           
    lines = "black",         
    zero = "darkred"          
  ),
  
  lwd.zero = 2,          
  lwd.ci = 2,            
  lty.ci = 1,            
  ci.vertices = TRUE,    
  ci.vertices.height = 0.1, 
  
  title = "Gene Interaction Hazard Ratios",
  xlab = "Hazard Ratio (95% CI)",
  
  grid = structure(c(0.5, 1, 1.5, 2), 
                   gp = gpar(lty = 2, col = "#CCCCCC")),
  
  clip = c(0.5, 2.5)  
)

#----Figure 5F----
# Lollipop plot of coefficients for each variable
library(ggplot2)
library(dplyr)
library(forcats) 
coffident_sorted <- coeffident %>%
  arrange(X1) %>%            
  mutate(regulor = fct_inorder(regulor)) 

ggplot(coffident_sorted, aes(x = X1, y = regulor)) +
  geom_segment(
    aes(x = 0, xend = X1, y = regulor, yend = regulor),
    color = "gray70",
    size = 1
  ) +
  geom_point(
    aes(color = X1, size = abs(X1)), 
    show.legend = TRUE
  ) +
  scale_color_gradient2(
    low = "#6FA9B1", 
    mid = "#FACF7F",
    high = "#D03A4F",
    midpoint = median(coffident_sorted$X1),
    name = "Value"
  ) +
  scale_size_continuous(
    range = c(2, 8), 
    guide = FALSE   
  ) +
  labs(
    title = "Model coefficient",
    x = "coefficient",
    y = NULL,
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "italic")
  )

#----Figure 5G----
# Perform immune infiltration analysis on the nine regulons in the model
library(TCGAplot)
coef <- read.csv("./coeffident.csv",row.names = 1)
selected_features <- rownames(coef)
tpm <- get_all_tpm()
meta <- get_all_meta()
tpm <- subset(tpm,Cancer %in% c("LUAD","LUSC"))
tpm <- subset(tpm,Group %in% c("Tumor"))
immu_ratio <- get_immu_ratio()
rownames(immu_ratio) <- gsub("[- ]", ".", rownames(immu_ratio)) 
immu_ratio <- immu_ratio[match(rownames(tpm),rownames(immu_ratio)),]
rownames(immu_ratio) <- substr(rownames(immu_ratio), 1, 12)

library(psych)
genelist <- list(
  c("PAX6", "ITGB1"),
  c("CEBPB", "IGF1"),
  c("CEBPB", "RETN"),
  c("NR1H4", "RARRES2"),
  c("PPARG", "RARRES2"),
  c("RELA", "IGF1"),
  c("SNAI1", "ITGB4"),
  c("SP1", "ITGA6"),
  c("USF2", "IGF1")
)
names(genelist) <- selected_features
ss <- unique(c("PAX6", "ITGB1","CEBPB", "IGF1","CEBPB", "RETN","NR1H4", "RARRES2","PPARG", "RARRES2","RELA", "IGF1","SNAI1", "ITGB4","SP1", "ITGA6","USF2", "IGF1"))

library(GSVA)
tpm <- tpm[,-c(1,2)]
tpm <- t(tpm)
tpm <- tpm[rownames(tpm) %in% ss,]
tpm <- as.matrix(tpm[, colSums(is.na(tpm)) == 0])
param <- gsvaParam(exprData = tpm, 
                   geneSets = genelist, 
                   kcdf = "Gaussian")

gsva_result <- gsva(param) 

tpm <- data.frame(t(gsva_result))
colnames(tpm) <- gsub("\\.","-",colnames(tpm))
tpm <- as.matrix(tpm)
risk_score <- tpm[, selected_features] %*% selected[selected_features, 1]
regulor.tpm <- data.frame(tpm)
regulor.tpm$RS <- risk_score

dd <- corr.test(immu_ratio,regulor.tpm,method = "pearson")
cor <- dd$r
pval <- dd$p

library(pheatmap)
library(reshape2)
cor_long <- melt(cor, varnames = c("Immune", "Regulator"), value.name = "cor")
pval_long <- melt(pval, varnames = c("Immune", "Regulator"), value.name = "p.value")

data <- merge(cor_long, pval_long, by = c("Immune", "Regulator"))

data$pstar <- ifelse(data$p.value < 0.001, "***",
                     ifelse(data$p.value < 0.01, "**",
                            ifelse(data$p.value < 0.05, "*", "")))

library(ggplot2)
write.csv(data,"./data.csv")
ggplot(data, aes(x = Regulator, y = Immune)) +
  geom_tile(aes(fill = cor), color = "grey90", size = 0.4) +
  scale_fill_gradient2(
    low = "#A4B9DD", mid = "white", high = "#CDA3CB", midpoint = 0,
    name = "Correlation"
  ) +
  geom_text(aes(label = pstar), color = "black", size = 4.5) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(caption = "* p < 0.05, ** p < 0.01, *** p < 0.001")
color <- c("#3B97CD","#DCB125")

#----Figure 5H----
# Score each driver gene using the risk score
library(UCell)
library(GSVA)
library(cowplot)
library(ggpubr)
M2.5 <- readRDS("./data/yang/wang/figure5/M2.5.rds")
coeffident <- read.csv("./data/yang/wang/figure5/lassco/coeffident.csv",row.names = 1)
rs.gene <- coeffident$regulor
gene_pairs <- strsplit(rs.gene, "-")
names(gene_pairs) <- rs.gene
Mscore <- AddModuleScore_UCell(M2.5,
                               features = gene_pairs,
                               assay = "RNA")

colnames(Mscore@meta.data)[8:16] <- gsub("_UCell", "", colnames(Mscore@meta.data)[8:16])

needed_columns <- coeffident$regulor
rs_scores <- sapply(1:nrow(Mscore@meta.data), function(i) {
  cell_values <- Mscore@meta.data[i, needed_columns]
  sum(cell_values * coeffident$X1)
})

Mscore@meta.data$rs_score <- rs_scores
saveRDS(Mscore,file = "./data/yang/wang/figure5/Msore.tf.rds")
custom_colors <- c(
  "EGFR"             = "#c6b7d4",
  "EGFR-BM"          = "#d44e26",
  "EGFR-co-mutation" = "#e3a264",
  "KRAS"             = "#6fc2d0",
  "KRAS-co-mutation" = "#6f9abf",
  "ALK"              = "#a5c49b",
  "ROS1"             = "#7266ac",
  "TP53"             = "#FF9966",
  "MET-BM"           = "#d84986",
  "HER2"             = "#2d588e"
)


group_levels <- c("EGFR","EGFR-BM","EGFR-co-mutation","KRAS","KRAS-co-mutation","ALK","ROS1","TP53","MET-BM","HER2")
Mscore$mutation <- factor(Mscore$mutation, levels = group_levels)
group_count <- length(group_levels)
y_max <- max(Mscore$rs_score, na.rm = TRUE)

plot_data <- Mscore@meta.data
stat.test <- plot_data %>%
  wilcox_test(rs_score ~ mutation, comparisons = comparisons_list, p.adjust.method = "none")

stat.test.signif <- stat.test %>% filter(p < 0.05)

p <- ggplot(plot_data, aes(x = mutation, y = rs_score, fill = mutation)) +
  geom_violin(width = 0.8, alpha = 0.85, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 1.0, outlier.alpha = 0.5, outlier.color = "gray30") +
  scale_fill_manual(values = custom_colors) +
  theme_bw()

if (nrow(stat.test.signif) > 0) {
  for (i in 1:nrow(stat.test.signif)) {
    group1 <- stat.test.signif$group1[i]
    group2 <- stat.test.signif$group2[i]
    p_val <- stat.test.signif$p[i]
    signif_label <- ifelse(p_val < 0.001, " ** *",
                           ifelse(p_val < 0.01, " ** ",
                                  ifelse(p_val < 0.05, "*", "")))
    y_pos <- max(plot_data$rs_score) * (1 + 0.1 * i)
    p <- p + geom_signif(
      comparisons = list(c(group1, group2)),
      annotations = signif_label,
      y_position = y_pos,
      tip_length = 0.01,
      vjust = 0.5
    )
  }
}

p

kw_test <- kruskal.test(rs_score ~ mutation, data = plot_data)
p_label <- paste0(
  "Kruskal-Wallis: χ²(", kw_test$parameter, ") = ",
  format(kw_test$statistic, digits = 3),
  ifelse(kw_test$p.value < 0.001, ", p < 0.001", 
         paste0(", p = ", format(kw_test$p.value, digits = 3)))
)

p <- ggplot(plot_data, aes(x = mutation, y = rs_score, fill = mutation)) +
  geom_violin(width = 0.8, alpha = 0.85, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white", outlier.size = 1.0, outlier.alpha = 0.5, outlier.color = "gray30") +
  scale_fill_manual(values = custom_colors) +
  theme_bw()+annotate(
    "text",
    x = 5,      
    y = Inf,          
    label = p_label,
    hjust = 0.5, 
    vjust = 1.5,      
    size = 4.5,
    fontface = "italic"
  ) 

#----Figure 5I----
# Association between the nine regulons in the model and mutation-derived cell types
library(networkD3)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggalluvial)
df <- Mscore@meta.data[,c(1,3,4,8:17)]
write.csv(df,file = "./data/yang/wang/figure5/sangshen.csv")
df <- read.csv(".Figure/figure5/sangshen.csv",row.names = 1)

library(tidyverse)
library(ggsankey)
df_long <- df %>%
  pivot_longer(cols = 4:13, names_to = "Regulator", values_to = "expression")
df_grouped <- df_long %>%
  group_by(Regulator, CELLTYPE, mutation) %>%
  summarise(mean_expr = mean(expression), .groups = "drop")
df_sankey <- df_grouped %>%
  mutate(mean_expr = ifelse(mean_expr < 0, 0, mean_expr)) %>%
  make_long(Regulator, CELLTYPE, mutation, value = mean_expr)
mutation_colors <- c(
  "EGFR"             = "#c6b7d4",
  "EGFR-BM"          = "#d44e26",
  "EGFR-co-mutation" = "#e3a264",
  "KRAS"             = "#6fc2d0",
  "KRAS-co-mutation" = "#6f9abf",
  "ALK"              = "#a5c49b",
  "ROS1"             = "#7266ac",
  "TP53"             = "#FF9966",
  "MET-BM"           = "#d84986",
  "HER2"             = "#2d588e"
)

celltype_colors <- c(
  "Neurophil" = "#D6E7A3","Mast" = "#53A85F", "Macrophage"= "#f5e886", "Fibroblast"="#7ac3df","Endothelial"="#f5aa61","CD4T"="#b4ffde","B"="#ffdd8e","Astrocyte"="#8fb4dc")
library(RColorBrewer)
regulators <- unique(df_sankey$node[df_sankey$x == "Regulator"])
regulator_colors <- brewer.pal(10, "Set3")
regulator_colors <- c(
  "PAX6-ITGB1"    = "#8DD3C7",
  "CEBPB-IGF1"    = "#FFFFB3",
  "CEBPB-RETN"    = "#BEBADA",
  "NR1H4-RARRES2" = "#FB8072",
  "PPARG-RARRES2" = "#80B1D3",
  "RELA-IGF1"     = "#FDB462",
  "SNAI1-ITGB4"   = "#B3DE69",
  "SP1-ITGA6"     = "#FCCDE5",
  "USF2-IGF1"     = "#D9D9D9",
  "rs_score"      = "#BC80BD"
)
CELLTYPE <- unique(df_sankey$node[df_sankey$x == "CELLTYPE"])
regulators <- names(regulator_colors)
mutations <- names(mutation_colors)  
CELLTYPE <- names(celltype_colors)
all_colors <- c(regulator_colors, celltype_colors, mutation_colors)
ggplot(df_sankey, aes(x = x,
                      next_x = next_x,
                      node = node,
                      next_node = next_node,
                      value = value,
                      fill = node,
                      label = node)) +
  scale_fill_manual(values = all_colors) +  
  geom_sankey(flow.alpha = 0.5, smooth = 6, width = 0.15) +
  geom_sankey_text(size = 3, color = "black") +
  theme_void() +
  theme(legend.position = "none")
