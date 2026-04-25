# ------------- Figure S6 --------------
#----Figure S6A----
library(survival)
library(survminer)
library(patchwork)
load("./figure5/lassco/data.RData")
coefficient <- read.csv("./figure5/lassco/coeffident.csv")

coef_vec <- setNames(coefficient$X1, coefficient$regulor)

calc_df <- function(tpm, meta){
  tpm <- as.matrix(tpm)
  selected_features <- intersect(colnames(tpm), names(coef_vec))
  risk_score <- tpm[, selected_features] %*% coef_vec[selected_features]
  
  df <- as.data.frame(tpm)
  df$RS <- as.numeric(risk_score)
  df$sample <- rownames(df)
  
  meta$sample <- rownames(meta)
  merge(df, meta, by = "sample")
}

plot_km <- function(data, title){
  data$RS_group <- ifelse(data$RS >= median(data$RS, na.rm = TRUE), "High", "Low")
  data$RS_group <- factor(data$RS_group, levels = c("Low","High"))
  
  fit <- survfit(Surv(time, event) ~ RS_group, data = data)
  
  ggsurvplot(
    fit,
    data = data,
    pval = TRUE,
    risk.table = TRUE,
    palette = c("#2E9FDF","#E7B800"),
    legend.title = "RS",
    legend.labs = c("Low","High"),
    title = title
  )
}

run_survival <- function(df, prefix){
  plots <- list(
    LUAD = plot_km(subset(df, Cancer == "LUAD"), paste0(prefix," LUAD")),
    LUSC = plot_km(subset(df, Cancer == "LUSC"), paste0(prefix," LUSC")),
    Stage_I = plot_km(subset(df, stage == "I"), paste0(prefix," Stage I")),
    Stage_II = plot_km(subset(df, stage == "II"), paste0(prefix," Stage II")),
    Stage_III = plot_km(subset(df, stage == "III"), paste0(prefix," Stage III")),
    Stage_IV = plot_km(subset(df, stage == "IV"), paste0(prefix," Stage IV"))
  )
  
  res <- arrange_ggsurvplots(
    plots, print = FALSE,
    ncol = 3, nrow = 2,
    risk.table.height = 0.25
  )
  
  pdf(paste0("./Survival/", prefix, "_survival.pdf"),
      width = 14, height = 10.5)
  grid::grid.draw(res)
  dev.off()
}

df_train <- calc_df(train.tpm, train.meta)
run_survival(df_train, "train")

#----Figure S6B----
df_test  <- calc_df(test.tpm, test.meta)
run_survival(df_test, "test")