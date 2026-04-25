# ------------- Figure S2 --------------
#----Figure S2A----
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
#Extract the top 200 upregulated genes for each module
modules <- paste0("module", 1:5)
subsce.corr.markers.up <- read.csv("./NSCLC/Figure/figure2/subsce.corr.markers.up.csv")
get_top_genes <- function(m) {
  subset(subsce.corr.markers.up, cluster == m) %>%
    group_by(cluster) %>%
    top_n(n = 200, wt = avg_log2FC) %>%
    pull(gene)
}

gene_list <- lapply(modules, get_top_genes)
names(gene_list) <- modules

#GO enrichment analysis
run_go <- function(genes, module_name) {
  
  ego <- enrichGO(
    gene = bitr(genes,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  # Top10
  go <- as.data.frame(ego) %>%
    filter(pvalue < 0.05) %>%
    head(10) %>%
    arrange(desc(qvalue))
  if (nrow(go) == 0) return(NULL)
  
  go$Description <- factor(go$Description, levels = go$Description)
  p <- ggplot(go, aes(x = Count, y = Description, fill = -log10(qvalue))) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_distiller(palette = "Blues", direction = 1) +
    labs(
      title = paste0("GO BP enrichment - ", module_name),
      x = "Gene count",
      y = NULL
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )
  return(p)
}
plot_list <- list()
for (m in modules) {
  plot_list[[m]] <- run_go(gene_list[[m]], m)
}

#----Figure S2B----
#The result figure was generated using the code for Figure 2D