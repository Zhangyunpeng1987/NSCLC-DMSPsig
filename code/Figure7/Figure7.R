# ------------- Figure 7 --------------
#----Figure 7A----
library(oncoPredict)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(DescTools)
library(testthat)
read_data <- function(path) {
  data <- readRDS(path)
  cat("Dimensions:", dim(data), "\n")
  print(data[1:5, 1:5])
  return(data)
}
# read data
GDSC2_Expr <- read_data("./DataFiles_1/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res <- read_data("./DataFiles_1/DataFiles/Training Data/GDSC2_Res.rds")
CTRP2_Expr <- read_data("./DataFiles_1/DataFiles/Training Data/CTRP2_Expr (TPM, not log transformed).rds")
CTRP2_Res <- read_data("./DataFiles_1/DataFiles/Training Data/CTRP2_Res.rds")

# Calculate phenotype
load("./TCGA_exp_meta.RData")
exp_meta <- exp_meta[exp_meta$Cancer %in% c("LUAD", "LUSC"), ]
sig.genes <- c("ITGA6", "IGF1", "ITGB1", "ITGB4", "PPARG", "RELA", "CEBPB", "SNAI1")
exp_meta <- exp_meta[, colnames(exp_meta) %in% sig.genes]

calcPhenotype(
  trainingExprData = CTRP2_Expr,
  trainingPtype = CTRP2_Res,
  testExprData = as.matrix(exp_meta),
  batchCorrect = 'eb',
  minNumSamples = 20,
  printOutput = TRUE,
  removeLowVaryingGenes = 0.2,
  removeLowVaringGenesFrom = "homogenizeData"
)
sen <- read.csv("./calcPhenotype_Output/DrugPredictions.csv", 
                row.names = 1, header = TRUE, check.names = FALSE)
core_drugs <- c("AGK-2", "lovastatin", "LY-2157299", "NSC30930", "saracatinib", "FSC231", "JW-55")
dat <- sen[, colnames(sen) %in% c("AGK-2", "lovastatin", "LY-2157299", "NSC30930", "saracatinib", "FSC231", "JW-55")]

DescTools::AllIdentical(
  rownames(exp_meta),
  rownames(dat)
) %>% testthat::expect_true()

### Merge
dd1 <- cbind(exp_meta, dat)

library(ggpubr)
ggplot(dd1, aes(x = IGF1, y = lovastatin)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = "top"
  ) +
  labs(y = "lovastatin sensitivity", x = "Expression of SNAI1") +
  theme_bw()
ggsave("./calcPhenotype_Output/lovastatin_IGF1.pdf", 
       width = 5, height = 5)