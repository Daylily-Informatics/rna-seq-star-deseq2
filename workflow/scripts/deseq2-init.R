log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(stringr)
  library(cli)
  library(DESeq2)
})

parallel <- FALSE
if (snakemake@threads > 1) {
  suppressPackageStartupMessages({
    library(BiocParallel)
  })
  register(MulticoreParam(snakemake@threads))
  parallel <- TRUE
}

counts <- read.table(
  snakemake@input[["counts"]],
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)
counts <- counts[!grepl("^N_", rownames(counts)), , drop = FALSE]

if (ncol(counts) == 0) {
  cli_abort("Counts matrix is empty")
}

samples_path <- snakemake@config[["samples"]]
coldata <- read.table(
  samples_path,
  header = TRUE,
  sep = "\t",
  row.names = "sample_name",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

if (!all(colnames(counts) %in% rownames(coldata))) {
  missing <- setdiff(colnames(counts), rownames(coldata))
  cli_abort("Sample metadata missing for: {missing}")
}

coldata <- coldata[colnames(counts), , drop = FALSE]

vof <- snakemake@config[["diffexp"]][["variables_of_interest"]]
if (is.null(vof)) {
  vof <- list()
}

for (name in names(vof)) {
  if (!name %in% colnames(coldata)) {
    cli_abort("Column '{name}' required by variables_of_interest is missing in samples.tsv")
  }
  base_level <- vof[[name]][["base_level"]]
  coldata[[name]] <- factor(coldata[[name]])
  if (!is.null(base_level) && base_level %in% levels(coldata[[name]])) {
    coldata[[name]] <- relevel(coldata[[name]], base = base_level)
  }
}

batch_effects <- snakemake@config[["diffexp"]][["batch_effects"]]
if (is.null(batch_effects)) {
  batch_effects <- character()
} else if (is.character(batch_effects)) {
  batch_effects <- batch_effects[nzchar(batch_effects)]
} else {
  batch_effects <- unlist(batch_effects)
  batch_effects <- batch_effects[nzchar(batch_effects)]
}

for (effect in batch_effects) {
  if (!effect %in% colnames(coldata)) {
    cli_abort("Batch effect column '{effect}' is missing in samples.tsv")
  }
  coldata[[effect]] <- factor(coldata[[effect]])
}

if ("patient" %in% colnames(coldata) && !"patient" %in% batch_effects) {
  if (length(unique(coldata[["patient"]])) > 1) {
    coldata[["patient"]] <- factor(coldata[["patient"]])
    batch_effects <- union(batch_effects, "patient")
  }
}

design_formula <- snakemake@config[["diffexp"]][["model"]]
if (is.null(design_formula)) {
  design_formula <- ""
}
design_formula <- str_trim(design_formula)

if (design_formula == "") {
  terms <- character()
  if (length(batch_effects) > 0) {
    terms <- c(terms, batch_effects)
  }
  if (length(names(vof)) > 0) {
    terms <- c(terms, names(vof))
  }
  if (length(terms) == 0) {
    design_formula <- "~ 1"
  } else {
    design_formula <- paste("~", paste(unique(terms), collapse = " + "))
  }
}

cli_inform(c("Design formula" = design_formula))

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = coldata,
  design = as.formula(design_formula)
)

keep <- rowSums(counts(dds)) >= 10
if (!any(keep)) {
  cli_warn("No genes with at least 10 counts; keeping all rows")
} else {
  dds <- dds[keep, ]
}

dds <- estimateSizeFactors(dds)

if (ncol(counts(dds)) >= 2) {
  dds <- DESeq(dds, parallel = parallel)
}

saveRDS(dds, file = snakemake@output[["dds"]])

norm_counts <- counts(dds, normalized = TRUE)
norm_counts <- as.data.frame(norm_counts)
norm_counts <- cbind(gene = rownames(norm_counts), norm_counts)
write.table(
  norm_counts,
  file = snakemake@output[["normcounts"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  norm_counts,
  file = snakemake@output[["normalized"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

sink(type = "message")
sink()
close(log)
