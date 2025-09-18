# deseq2-init.R — robust init with saturated-design / zero-variance guards

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")
.onexit <- function(){ sink(type="message"); sink(); close(log) }
reg.finalizer(environment(), function(...) try(.onexit(), silent=TRUE), onexit=TRUE)

suppressPackageStartupMessages({
  library(stringr)
  library(cli)
  library(DESeq2)
})

# ---------- parallel ----------
parallel <- FALSE
if (snakemake@threads > 1) {
  suppressPackageStartupMessages({ library(BiocParallel) })
  if (.Platform$OS.type != "windows") {
    register(MulticoreParam(workers = snakemake@threads))
  } else {
    register(SnowParam(workers = snakemake@threads, type="SOCK"))
  }
  parallel <- TRUE
  cli_inform("BiocParallel workers = {snakemake@threads}")
}

# ---------- counts ----------
counts <- read.table(
  snakemake@input[["counts"]],
  header=TRUE, sep="\t", row.names=1,
  check.names=FALSE, quote="", comment.char=""
)
counts <- counts[!grepl("^N_", rownames(counts)), , drop=FALSE]
if (ncol(counts) == 0) cli_abort("Counts matrix is empty")

# quick duplication/identity checks
dup_names <- duplicated(colnames(counts))
if (any(dup_names)) cli_warn("Duplicate sample names: {paste(unique(colnames(counts)[dup_names]), collapse=', ')}")
if (ncol(counts) >= 2) {
  ref <- counts[,1,drop=FALSE]; ident <- vapply(2:ncol(counts), function(j) isTRUE(all(ref[,1]==counts[,j])), TRUE)
  if (all(ident)) cli_warn("All columns identical to the first — DE impossible.")
}

# ---------- metadata ----------
samples_path <- snakemake@config[["samples"]]
coldata <- read.table(
  samples_path,
  header=TRUE, sep="\t", row.names="sample_name",
  check.names=FALSE, quote="", comment.char=""
)
if (!all(colnames(counts) %in% rownames(coldata))) {
  missing <- setdiff(colnames(counts), rownames(coldata))
  cli_abort("Sample metadata missing for: {paste(missing, collapse=', ')}")
}
coldata <- coldata[colnames(counts), , drop=FALSE]

# ---------- variables of interest ----------
vof <- snakemake@config[["diffexp"]][["variables_of_interest"]]
if (is.null(vof)) vof <- list()
for (name in names(vof)) {
  if (!name %in% colnames(coldata)) cli_abort("Missing VOI column '{name}' in samples.tsv")
  coldata[[name]] <- factor(coldata[[name]])
  base_level <- vof[[name]][["base_level"]]
  if (!is.null(base_level) && base_level %in% levels(coldata[[name]])) {
    coldata[[name]] <- relevel(coldata[[name]], ref=base_level)  # <- correct arg
  }
}

# ---------- batch effects ----------
batch_effects <- snakemake@config[["diffexp"]][["batch_effects"]]
if (is.null(batch_effects)) {
  batch_effects <- character()
} else if (is.character(batch_effects)) {
  batch_effects <- batch_effects[nzchar(batch_effects)]
} else {
  batch_effects <- unlist(batch_effects, use.names=FALSE)
  batch_effects <- batch_effects[nzchar(batch_effects)]
}
for (effect in batch_effects) {
  if (!effect %in% colnames(coldata)) cli_abort("Missing batch column '{effect}' in samples.tsv")
  coldata[[effect]] <- factor(coldata[[effect]])
}

# opportunistic patient as batch
if ("patient" %in% colnames(coldata) && !"patient" %in% batch_effects) {
  if (length(unique(coldata[["patient"]])) > 1) {
    coldata[["patient"]] <- factor(coldata[["patient"]])
    batch_effects <- union(batch_effects, "patient")
    cli_inform("Added 'patient' to batch_effects")
  }
}

# ---------- design (defensive) ----------
design_formula <- snakemake@config[["diffexp"]][["model"]]
if (is.null(design_formula)) design_formula <- ""
design_formula <- str_trim(design_formula)

if (design_formula == "") {
  terms <- unique(c(batch_effects, names(vof)))
  if (length(terms) > 0) {
    keep_terms <- vapply(terms, function(tn){
      if (!tn %in% colnames(coldata)) return(FALSE)
      x <- coldata[[tn]]
      nlev <- if (is.factor(x)) nlevels(x) else length(unique(x))
      nlev >= 2
    }, logical(1))
    if (any(!keep_terms)) cli_warn("Dropping non-varying terms: {paste(terms[!keep_terms], collapse=', ')}")
    terms <- terms[keep_terms]
  }
  design_formula <- if (length(terms)==0) "~ 1" else paste("~", paste(terms, collapse=" + "))
} else {
  parsed <- setdiff(all.vars(stats::terms(stats::as.formula(design_formula))), "1")
  if (length(parsed) > 0) {
    deg <- parsed[parsed %in% colnames(coldata) & vapply(parsed, function(tn) length(unique(coldata[[tn]]))<2, TRUE)]
    if (length(deg) > 0) cli_warn("Design has non-varying terms: {paste(deg, collapse=', ')}")
  }
}
cli_inform(c("Design formula" = design_formula))

# show levels
terms_now <- setdiff(all.vars(terms(as.formula(design_formula))), "1")
for (tn in terms_now) {
  if (tn %in% colnames(coldata)) {
    lv <- if (is.factor(coldata[[tn]])) levels(coldata[[tn]]) else sort(unique(coldata[[tn]]))
    cli_inform("Term '{tn}' levels: {paste(lv, collapse=', ')}")
  }
}

# ---------- build DESeqDataSet ----------
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = coldata,
  design    = as.formula(design_formula)
)

# filter low-count genes
keep <- rowSums(counts(dds)) >= 10
if (any(keep)) {
  dds <- dds[keep, ]
  cli_inform("Kept {nrow(dds)} genes after filtering (sum ≥ 10).")
} else {
  cli_warn("No genes meet sum ≥ 10; skipping filter.")
}

dds <- estimateSizeFactors(dds)

# ---------- preflight: variance + saturated design ----------
has_gene_variation <- function(dds) {
  if (ncol(counts(dds)) < 2) return(FALSE)
  sum(apply(counts(dds), 1, function(x) var(as.numeric(x)) > 0), na.rm=TRUE) > 0
}
nzv <- if (ncol(counts(dds)) >= 2) sum(apply(counts(dds),1,function(x) var(as.numeric(x))>0)) else 0
cli_inform("Genes with nonzero variance across samples: {nzv}")

is_saturated <- function(design_formula, coldata) {
  X <- try(model.matrix(as.formula(design_formula), data=coldata), silent=TRUE)
  if (inherits(X, "try-error")) return(TRUE)
  n <- nrow(X); p <- ncol(X)
  cli_inform("Model matrix: n={n}, p={p}, residual df={n - p}")
  n <= p
}

should_run_deseq <- function(dds, design_formula) {
  if (ncol(counts(dds)) < 2) return(FALSE)
  if (identical(trimws(design_formula), "~ 1")) return(FALSE)
  if (!has_gene_variation(dds)) return(FALSE)
  if (is_saturated(design_formula, as.data.frame(colData(dds)))) return(FALSE)
  TRUE
}

if (should_run_deseq(dds, design_formula)) {
  cli_inform("Running DESeq() …")
  dds <- DESeq(dds, parallel=parallel)
} else {
  cli_inform("Skipping DESeq(): no testable effects and/or saturated design and/or no gene variance (normalization only).")
}

# ---------- outputs ----------
saveRDS(dds, file = snakemake@output[["dds"]])

norm_counts <- counts(dds, normalized=TRUE)
norm_counts <- as.data.frame(norm_counts, check.names=FALSE)
norm_counts <- cbind(gene=rownames(norm_counts), norm_counts)

write.table(norm_counts, file=snakemake@output[["normcounts"]], sep="\t", quote=FALSE, row.names=FALSE)
write.table(norm_counts, file=snakemake@output[["normalized"]],  sep="\t", quote=FALSE, row.names=FALSE)

.onexit()
