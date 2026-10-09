source("R/bmd_data.R")
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
library(ggplot2)

ralox_file <- "data/raw_processed/GSE207739_BMDExpress_ready_Raloxifene.tsv.gz"
res_file <- "data/raw_processed/GSE207739_BMDExpress_ready_Resorcinol.tsv.gz"

make_pca <- function(file, chemical) {
  
  dat <- read.delim(
    gzfile(file),
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  dose_row <- dat[1, -1]
  expr <- dat[-1, ]
  
  gene_ids <- expr[[1]]
  expr <- expr[, -1]
  
  expr[] <- lapply(expr, as.numeric)
  rownames(expr) <- gene_ids
  
  sample_names <- colnames(expr)
  
  parts <- strsplit(sample_names, "_")
  
  group <- sapply(parts, function(x) x[length(x) - 1])
  replicate <- sapply(parts, function(x) x[length(x)])
  
  dose <- as.numeric(dose_row)
  
  metadata <- data.frame(
    sample = sample_names,
    chemical = chemical,
    group = group,
    dose = dose,
    replicate = replicate,
    stringsAsFactors = FALSE
  )
  
  log_expr <- log2(expr + 1)
  
  keep <- apply(log_expr, 1, var, na.rm = TRUE) > 0
  log_expr <- log_expr[keep, ]
  
  pca <- prcomp(t(log_expr), scale. = TRUE)
  
  var_exp <- round(
    100 * pca$sdev^2 / sum(pca$sdev^2),
    1
  )
  
  pca_df <- data.frame(
    sample = sample_names,
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    metadata[, -1]
  )
  
  list(
    data = pca_df,
    variance = var_exp
  )
}

ralox_pca <- make_pca(
  ralox_file,
  "Raloxifene"
)

res_pca <- make_pca(
  res_file,
  "Resorcinol"
)

ralox_pca$data$replicate <- factor(
  ralox_pca$data$replicate,
  levels = c("A", "B", "C", "D", "E")
)

res_pca$data$replicate <- factor(
  res_pca$data$replicate,
  levels = c("A", "B", "C", "D", "E")
)

ralox_pca$data$group <- factor(
  ralox_pca$data$group,
  levels = c("CC20", "10m4", "10m3", "10m2", "10m1", "EC20")
)

res_pca$data$group <- factor(
  res_pca$data$group,
  levels = c("CC20", "10m4", "10m3", "10m2", "10m1", "EC20")
)

ralox_plot <- ggplot(
  ralox_pca$data,
  aes(
    x = PC1,
    y = PC2,
    colour = group,
    shape = replicate
  )
) +
  geom_point(size = 4, alpha = 0.9) +
  scale_shape_manual(
    values = c(
      "A" = 16,
      "B" = 17,
      "C" = 15,
      "D" = 18,
      "E" = 8
    )
  ) +
  labs(
    title = "Raloxifene PCA",
    subtitle = "Shape indicates replicate",
    x = paste0("PC1 (", ralox_pca$variance[1], "%)"),
    y = paste0("PC2 (", ralox_pca$variance[2], "%)"),
    colour = "Dose group",
    shape = "Replicate"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

res_plot <- ggplot(
  res_pca$data,
  aes(
    x = PC1,
    y = PC2,
    colour = group,
    shape = replicate
  )
) +
  geom_point(size = 4, alpha = 0.9) +
  scale_shape_manual(
    values = c(
      "A" = 16,
      "B" = 17,
      "C" = 15,
      "D" = 18,
      "E" = 8
    )
  ) +
  labs(
    title = "Resorcinol PCA",
    subtitle = "Shape indicates replicate",
    x = paste0("PC1 (", res_pca$variance[1], "%)"),
    y = paste0("PC2 (", res_pca$variance[2], "%)"),
    colour = "Dose group",
    shape = "Replicate"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  "results/figures/raloxifene_PCA_replicates.png",
  ralox_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "results/figures/resorcinol_PCA_replicates.png",
  res_plot,
  width = 8,
  height = 6,
  dpi = 300
)

print(ralox_plot)
print(res_plot)