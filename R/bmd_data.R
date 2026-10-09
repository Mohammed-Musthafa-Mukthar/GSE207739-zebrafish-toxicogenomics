# GSE207739 zebrafish toxicogenomics
# first look at the BMDExpress data


# files

raloxifene_file <- file.path(
  "data",
  "raw_processed",
  "GSE207739_BMDExpress_ready_Raloxifene.tsv.gz"
)

resorcinol_file <- file.path(
  "data",
  "raw_processed",
  "GSE207739_BMDExpress_ready_Resorcinol.tsv.gz"
)


# both files 

file.exists(raloxifene_file)
file.exists(resorcinol_file)


# read data

raloxifene_raw <- read.delim(
  gzfile(raloxifene_file),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

resorcinol_raw <- read.delim(
  gzfile(resorcinol_file),
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# quick look

dim(raloxifene_raw)
dim(resorcinol_raw)

head(raloxifene_raw[, 1:6])
head(resorcinol_raw[, 1:6])


# seperate dose row and gene expression

prepare_bmd_data <- function(data, chemical) {
  
  sample_names <- colnames(data)[-1]
  
  dose_row <- data[data$SampleID == "Dose", -1, drop = FALSE]
  
  doses <- as.numeric(dose_row[1, ])
  
  # remove dose row
  gene_data <- data[data$SampleID != "Dose", ]
  
  gene_ids <- gene_data$SampleID
  
  expr <- as.matrix(gene_data[, -1])
  
  storage.mode(expr) <- "numeric"
  
  rownames(expr) <- gene_ids
  
  
  # get group and replicate from sample name
  sample_parts <- strsplit(sample_names, "_")
  
  group <- sapply(
    sample_parts,
    function(x) x[2]
  )
  
  replicate <- sapply(
    sample_parts,
    function(x) x[3]
  )
  
  
  metadata <- data.frame(
    sample_id = sample_names,
    chemical = chemical,
    group = group,
    dose_uM = doses,
    replicate = replicate,
    stringsAsFactors = FALSE
  )
  
  metadata$treatment <- ifelse(
    metadata$dose_uM == 0,
    "Carrier control",
    "Exposed"
  )
  
  
  return(
    list(
      expression = expr,
      metadata = metadata
    )
  )
}


# prepare both chemicals

raloxifene <- prepare_bmd_data(
  raloxifene_raw,
  "Raloxifene"
)

resorcinol <- prepare_bmd_data(
  resorcinol_raw,
  "Resorcinol"
)


# number of genes and samples

cat(
  "\nRaloxifene:",
  nrow(raloxifene$expression),
  "genes,",
  ncol(raloxifene$expression),
  "samples\n"
)

cat(
  "Resorcinol:",
  nrow(resorcinol$expression),
  "genes,",
  ncol(resorcinol$expression),
  "samples\n"
)


# sample info

print(raloxifene$metadata)
print(resorcinol$metadata)


# number of reps at each dose

table(raloxifene$metadata$dose_uM)

table(resorcinol$metadata$dose_uM)


# missing values?

sum(is.na(raloxifene$expression))

sum(is.na(resorcinol$expression))


# make output folders

dir.create(
  "results",
  showWarnings = FALSE
)

dir.create(
  file.path("results", "figures"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path("results", "tables"),
  recursive = TRUE,
  showWarnings = FALSE
)


# save metadata

sample_metadata <- rbind(
  raloxifene$metadata,
  resorcinol$metadata
)

write.csv(
  sample_metadata,
  file.path(
    "results",
    "tables",
    "sample_metadata.csv"
  ),
  row.names = FALSE
)


# PCA

run_pca <- function(dataset, chemical) {
  
  expr <- dataset$expression
  metadata <- dataset$metadata
  
  # log transform
  log_expr <- log2(expr + 1)
  
  # remove genes with no variation
  gene_var <- apply(
    log_expr,
    1,
    var,
    na.rm = TRUE
  )
  
  keep <- gene_var > 0 & !is.na(gene_var)
  
  pca_data <- log_expr[
    keep,
    ,
    drop = FALSE
  ]
  
  cat(
    "\n",
    chemical,
    ":",
    nrow(pca_data),
    "genes used for PCA\n"
  )
  
  
  pca <- prcomp(
    t(pca_data),
    center = TRUE,
    scale. = FALSE
  )
  
  
  variance <- (
    pca$sdev^2 /
      sum(pca$sdev^2)
  ) * 100
  
  
  pca_table <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    stringsAsFactors = FALSE
  )
  
  
  # match metadata with PCA samples
  pca_table <- merge(
    pca_table,
    metadata,
    by = "sample_id",
    sort = FALSE
  )
  
  pca_table <- pca_table[
    match(
      metadata$sample_id,
      pca_table$sample_id
    ),
  ]
  
  
  write.csv(
    pca_table,
    file.path(
      "results",
      "tables",
      paste0(
        tolower(chemical),
        "_PCA_coordinates.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # dose groups in increasing order
  group_order <- unique(
    metadata$group[
      order(metadata$dose_uM)
    ]
  )
  
  symbols <- as.integer(
    factor(
      pca_table$group,
      levels = group_order
    )
  )
  
  
  # save PCA plot
  
  png(
    file.path(
      "results",
      "figures",
      paste0(
        tolower(chemical),
        "_PCA.png"
      )
    ),
    width = 1400,
    height = 1100,
    res = 160
  )
  
  
  plot(
    pca_table$PC1,
    pca_table$PC2,
    pch = symbols,
    cex = 1.5,
    xlab = paste0(
      "PC1 (",
      round(variance[1], 1),
      "%)"
    ),
    ylab = paste0(
      "PC2 (",
      round(variance[2], 1),
      "%)"
    ),
    main = paste(
      chemical,
      "PCA"
    )
  )
  
  
  # letters are the biological reps
  text(
    pca_table$PC1,
    pca_table$PC2,
    labels = pca_table$replicate,
    pos = 3,
    cex = 0.7
  )
  
  
  legend(
    "topright",
    legend = group_order,
    pch = seq_along(group_order),
    title = "Dose group",
    cex = 0.8
  )
  
  
  dev.off()
  
  
  cat(
    chemical,
    "PC1:",
    round(variance[1], 2),
    "%\n"
  )
  
  cat(
    chemical,
    "PC2:",
    round(variance[2], 2),
    "%\n"
  )
  
  
  return(pca)
}


# run PCA for both

raloxifene_pca <- run_pca(
  raloxifene,
  "Raloxifene"
)

resorcinol_pca <- run_pca(
  resorcinol,
  "Resorcinol"
)


cat("\nDone\n")


# PCA coloured by replicate

plot_pca_rep <- function(pca, metadata, chemical) {
  
  variance <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  
  reps <- factor(metadata$replicate)
  
  png(
    file.path(
      "results",
      "figures",
      paste0(tolower(chemical), "_PCA_replicates.png")
    ),
    width = 1400,
    height = 1100,
    res = 160
  )
  
  plot(
    pca$x[, 1],
    pca$x[, 2],
    pch = as.integer(reps) + 15,
    cex = 1.5,
    xlab = paste0("PC1 (", round(variance[1], 1), "%)"),
    ylab = paste0("PC2 (", round(variance[2], 1), "%)"),
    main = paste(chemical, "PCA by replicate")
  )
  
  text(
    pca$x[, 1],
    pca$x[, 2],
    labels = metadata$group,
    pos = 3,
    cex = 0.65
  )
  
  legend(
    "topright",
    legend = levels(reps),
    pch = seq_along(levels(reps)) + 15,
    title = "Replicate"
  )
  
  dev.off()
}


plot_pca_rep(
  raloxifene_pca,
  raloxifene$metadata,
  "Raloxifene"
)

plot_pca_rep(
  resorcinol_pca,
  resorcinol$metadata,
  "Resorcinol"
)