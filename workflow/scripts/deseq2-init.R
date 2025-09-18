# deseq2-init.R — robust init with defensive design + logging

# ---- logging to Snakemake log ----
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")
.onexit <- function() {
  sink(type = "message"); sink(); close(log)
}
# ensure cleanup even on error
reg.finalizer(environment(), function(...) try(.onexit(), silent = TRUE), onexit = TRUE)

suppressPackageStartupMessages({
  library(stringr)
  library(cli)
  library(DESeq2)
})

# ---- parallel setup (POSIX safe; MulticoreParam not on Windows) ----
parallel <- FALSE
if (snakemake@threads > 1) {
  suppressPackageStartupMessages({ library(BiocParallel) })
  if (.Platform$OS.type != "windows") {
    register(MulticoreParam(workers = snakemake@threads))
  } else {
    register(SnowParam(workers = snakemake@threads, type = "SOCK"))
  }
  parallel <- TRUE
  cli_inform("BiocParallel registered with {snakemake@threads} workers (parallel={parallel}).")
}

# ---- read counts (STAR ReadsPerGene-like) ----
counts <- read.table(
  snakemake@input[["counts"]],
  header = TRUE, sep = "\t", row.names = 1,
  check.names = FALSE, quote = "", comment.char = ""
)
# drop STAR N_* rows
counts <- counts[!grepl("^N_", rownames(counts)), , drop = FALSE]

if (ncol(counts) == 0) cli_abort("Counts matrix is empty")

# ---- read sample metadata ----
samples_path <- snakemake@config[["samples"]]
coldata <- read.table(
  samples_path,
  header = TRUE, sep = "\t", row.names = "sample_name",
  check.names = FALSE, quote = "", comment.char = ""
)

# ensure all count columns have metadata
if (!all(colnames(counts) %in% rownames(coldata))) {
  missing <- setdiff(colnames(counts), rownames(coldata))
  cli_abort("Sample metadata missing for: {paste(missing, collapse=', ')}")
}
# order coldata to match counts
coldata <- coldata[colnames(counts), , drop = FALSE]

# ---- variables of interest (VOI) ----
vof <- snakemake@config[["diffexp"]][["variables_of_interest"]]
if (is.null(vof)) vof <- list()

for (name in names(vof)) {
  if (!name %in% colnames(coldata)) {
    cli_abort("Column '{name}' required by variables_of_interest is missing in samples.tsv")
  }
  coldata[[name]] <- factor(coldata[[name]])
  base_level <- vof[[name]][["base_level"]]
  if (!is.null(base_level) && base_level %in% levels(coldata[[name]])) {
    # FIX: use ref= (not base=)
    coldata[[name]] <- relevel(coldata[[name]], ref = base_level)
  }
}

# ---- batch effects (sanitize shapes) ----
batch_effects <- snakemake@config[["diffexp"]][["batch_effects"]]
if (is.null(batch_effects)) {
  batch_effects <- character()
} else if (is.character(batch_effects)) {
  batch_effects <- batch_effects[nzchar(batch_effects)]
} else {
  batch_effects <- unlist(batch_effects, use.names = FALSE)
  batch_effects <- batch_effects[nzchar(batch_effects)]
}

for (effect in batch_effects) {
  if (!effect %in% colnames(coldata)) {
    cli_abort("Batch effect column '{effect}' is missing in samples.tsv")
  }
  coldata[[effect]] <- factor(coldata[[effect]])
}

# opportunistic inclusion of 'patient' as batch if multi-patient
if ("patient" %in% colnames(coldata) && !"patient" %in% batch_effects) {
  if (length(unique(coldata[["patient"]])) > 1) {
    coldata[["patient"]] <- factor(coldata[["patient"]])
    batch_effects <- union(batch_effects, "patient")
    cli_inform("Added 'patient' to batch_effects (multiple patients detected).")
  }
}

# ---- design builder (defensive) ----
design_formula <- snakemake@config[["diffexp"]][["model"]]
if (is.null(design_formula)) design_formula <- ""
design_formula <- str_trim(design_formula)

if (design_formula == "") {
  terms <- character()
  if (length(batch_effects) > 0) terms <- c(terms, batch_effects)
  if (length(names(vof)) > 0)   terms <- c(terms, names(vof))
  terms <- unique(terms)

  # drop degenerate terms with <2 observed levels
  if (length(terms) > 0) {
    keep_terms <- vapply(terms, function(tn) {
      if (!tn %in% colnames(coldata)) return(FALSE)
      x <- coldata[[tn]]
      nlev <- if (is.factor(x)) nlevels(x) else length(unique(x))
      nlev >= 2
    }, logical(1))
    dropped <- terms[!keep_terms]
    if (length(dropped) > 0) {
      cli_warn("Dropping non-varying terms from design: {paste(dropped, collapse=', ')}")
    }
    terms <- terms[keep_terms]
  }

  design_formula <- if (length(terms) == 0) "~ 1" else paste("~", paste(terms, collapse = " + "))
} else {
  # optional soft warning on degenerate user-specified terms
  parsed_terms <- setdiff(all.vars(stats::terms(stats::as.formula(design_formula))), "1")
  if (length(parsed_terms) > 0) {
    deg <- parsed_terms[parsed_terms %in% colnames(coldata) &
      vapply(parsed_terms, function(tn) length(unique(coldata[[tn]])) < 2, logical(1))]
    if (length(deg) > 0) {
      cli_warn("Design includes non-varying terms: {paste(deg, collapse=', ')}; DESeq2 may drop them.")
    }
  }
}

cli_inform(c("Design formula" = design_formula))

# print observed levels for each term in the final design (helps debug)
terms_now <- setdiff(all.vars(terms(as.formula(design_formula))), "1")
if (length(terms_now) > 0) {
  for (tn in terms_now) {
    if (tn %in% colnames(coldata)) {
      lv <- if (is.factor(coldata[[tn]])) levels(coldata[[tn]]) else sort(unique(coldata[[tn]]))
      cli_inform("Term '{tn}' levels: {paste(lv, collapse=', ')}")
    }
  }
} else {
  cli_inform("No design terms (intercept-only). Will perform normalization only.")
}

# ---- build DESeq2 dataset ----
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = coldata,
  design = as.formula(design_formula)
)

# filter low-count genes (sum < 10 across all samples)
keep <- rowSums(counts(dds)) >= 10
if (!any(keep)) {
  cli_warn("No genes with total counts ≥ 10; keeping all rows (no filtering applied).")
} else {
  dds <- dds[keep, ]
  cli_inform("Kept {nrow(dds)} genes after filtering (sum counts ≥ 10).")
}

# size factors always OK
dds <- estimateSizeFactors(dds)

# ---- safe gate for DESeq() ----
run_deseq <- function(dds, design_formula) {
  if (ncol(counts(dds)) < 2) return(FALSE)
  if (identical(trimws(design_formula), "~ 1")) return(FALSE)
  terms_now <- setdiff(all.vars(stats::terms(stats::as.formula(design_formula))), "1")
  if (length(terms_now) == 0) return(FALSE)
  varying <- vapply(terms_now, function(tn) length(unique(colData(dds)[[tn]])) >= 2, logical(1))
  any(varying)
}

if (run_deseq(dds, design_formula)) {
  cli_inform("Running DESeq() …")
  dds <- DESeq(dds, parallel = parallel)
} else {
  cli_inform("Skipping DESeq(): design has no testable effects (normalization only).")
}

# ---- outputs ----
saveRDS(dds, file = snakemake@output[["dds"]])

norm_counts <- counts(dds, normalized = TRUE)
norm_counts <- as.data.frame(norm_counts, check.names = FALSE)
norm_counts <- cbind(gene = rownames(norm_counts), norm_counts)

write.table(
  norm_counts,
  file = snakemake@output[["normcounts"]],
  sep = "\t", quote = FALSE, row.names = FALSE
)
# legacy alias
write.table(
  norm_counts,
  file = snakemake@output[["normalized"]],
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ---- tidy logging teardown ----
.onexit()
