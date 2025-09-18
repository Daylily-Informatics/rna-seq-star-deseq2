# scripts/deseq2.R — robust per-contrast runner

# ---- logging ----
log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")
.onexit <- function(){ sink(type="message"); sink(); close(log) }
on.exit(.onexit(), add = TRUE)

suppressPackageStartupMessages({
  library(cli)
  library(DESeq2)
})

# ---- parallel ----
parallel <- FALSE
if (!is.null(snakemake@threads) && snakemake@threads > 1) {
  suppressPackageStartupMessages({ library(BiocParallel) })
  if (.Platform$OS.type != "windows") {
    register(MulticoreParam(workers = snakemake@threads))
  } else {
    register(SnowParam(workers = snakemake@threads, type = "SOCK"))
  }
  parallel <- TRUE
  cli_inform("BiocParallel workers = {snakemake@threads}")
}

# ---- IO ----
dds_path   <- snakemake@input[[1]]
tsv_out    <- snakemake@output[["table"]]
maplot_out <- snakemake@output[["ma_plot"]]

dds <- readRDS(dds_path)
if (ncol(dds) < 2) cli_abort("Need at least 2 samples; got {ncol(dds)}")

# ---- build 'contrast' from config ----
cc <- snakemake@config[["diffexp"]][["contrasts"]][[ snakemake@wildcards[["contrast"]] ]]
if (is.null(cc)) cli_abort("contrast '{snakemake@wildcards[['contrast']]}' not found in config")

if (is.list(cc) && length(cc) == 2) {
  voi <- cc[["variable_of_interest"]]
  loi <- cc[["level_of_interest"]]
  if (is.null(voi) || is.null(loi)) {
    cli_abort("For list-form contrast, need 'variable_of_interest' and 'level_of_interest'")
  }
  if (!(voi %in% names(snakemake@config[["diffexp"]][["variables_of_interest"]]))) {
    cli_abort("VOI '{voi}' not present under diffexp.variables_of_interest in config.yaml")
  }
  base <- snakemake@config[["diffexp"]][["variables_of_interest"]][[voi]][["base_level"]]
  if (is.null(base)) cli_abort("VOI '{voi}' lacks 'base_level' in config")
  contrast <- c(voi, loi, base)
} else if (is.character(cc) && length(cc) == 1) {
  # advanced expression, e.g., 'list(c("condition","tumor","normal"), c("batch","B","A"))'
  contrast <- eval(parse(text = cc))
} else {
  cli_abort("Unsupported contrast config for '{snakemake@wildcards[['contrast']]}'")
}

# helpers
write_empty_outputs <- function(reason) {
  cli_warn("{reason} — emitting empty outputs.")
  header <- c("gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  out <- setNames(data.frame(matrix(nrow = 0, ncol = length(header))), header)
  write.table(out, file = tsv_out, sep = "\t", quote = FALSE, row.names = FALSE)
  svg(maplot_out, width = 7, height = 5); plot.new(); title(main = reason); dev.off()
}

is_saturated <- function(formula_obj, coldata) {
  X <- try(model.matrix(formula_obj, data = as.data.frame(coldata)), silent = TRUE)
  if (inherits(X, "try-error")) return(TRUE)
  n <- nrow(X); p <- ncol(X)
  cli_inform("Model matrix (subset): n={n}, p={p}, residual df={n - p}")
  n <= p
}

has_gene_variation <- function(dds) {
  if (ncol(counts(dds)) < 2) return(FALSE)
  sum(apply(counts(dds), 1, function(x) var(as.numeric(x)) > 0), na.rm = TRUE) > 0
}

# ---- subset to the two levels in this contrast (first element can be a factor name or a true results() list) ----
if (is.character(contrast) && length(contrast) == 3) {
  fact <- contrast[[1]]; num <- contrast[[2]]; den <- contrast[[3]]
  if (!fact %in% colnames(colData(dds))) {
    write_empty_outputs(sprintf("Factor '%s' not in colData", fact)); quit(save="no")
  }
  keep <- colData(dds)[[fact]] %in% c(num, den)
  dds2 <- dds[, keep, drop = FALSE]
  colData(dds2)[[fact]] <- droplevels(colData(dds2)[[fact]])
  present <- levels(colData(dds2)[[fact]])
  cli_inform("Contrast: {fact}: {num} vs {den}; kept {sum(keep)} samples; levels present: {paste(present, collapse=', ')}")
  if (!all(c(num, den) %in% present)) {
    write_empty_outputs("Missing one or both contrast levels in subset"); quit(save="no")
  }
  # ensure denominator is reference
  if (is.factor(colData(dds2)[[fact]]) && den %in% levels(colData(dds2)[[fact]])) {
    colData(dds2)[[fact]] <- relevel(colData(dds2)[[fact]], ref = den)
  }
  design_formula <- design(dds2)
} else {
  # complex contrast list: we can't easily infer which factor/levels; use all samples
  dds2 <- dds
  design_formula <- design(dds2)
  cli_inform("Advanced contrast; running on all samples (cannot pre-subset safely).")
}

# ---- run DE or short-circuit ----
if (identical(trimws(format(design_formula)), "~ 1")) {
  write_empty_outputs("Intercept-only design in subset"); quit(save="no")
}
if (!has_gene_variation(dds2)) {
  write_empty_outputs("No gene variance in subset"); quit(save="no")
}
if (is_saturated(design_formula, colData(dds2))) {
  write_empty_outputs("Saturated design (n <= p)"); quit(save="no")
}

cli_inform("Running DESeq() on subset …")
dds2 <- DESeq(dds2, parallel = parallel)

# ---- results + shrink ----
res <- try(results(dds2, contrast = contrast, parallel = parallel), silent = TRUE)
if (inherits(res, "try-error")) {
  write_empty_outputs(sprintf("results() failed: %s", as.character(res))); quit(save="no")
}

# Prefer apeglm if available, otherwise skip shrinkage gracefully
shrunken <- FALSE
if (requireNamespace("apeglm", quietly = TRUE)) {
  suppressPackageStartupMessages(library(apeglm))
  res <- try(lfcShrink(dds2, contrast = contrast, res = res, type = "apeglm"), silent = TRUE)
  if (!inherits(res, "try-error")) {
    shrunken <- TRUE
  } else {
    cli_warn("lfcShrink(apeglm) failed; continuing without shrinkage.")
    res <- results(dds2, contrast = contrast, parallel = parallel)
  }
} else {
  cli_warn("Package 'apeglm' not installed; continuing without shrinkage.")
}

res_df <- as.data.frame(res)
res_df <- cbind(gene = rownames(res_df), res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

# ---- outputs ----
write.table(res_df, file = tsv_out, sep = "\t", quote = FALSE, row.names = FALSE)

svg(maplot_out, width = 7, height = 5)
plotMA(res, ylim = c(-5, 5))
ttl <- if (is.character(contrast) && length(contrast) == 3) {
  sprintf("MA: %s (%s vs %s)%s",
          contrast[[1]], contrast[[2]], contrast[[3]],
          if (shrunken) " [apeglm]" else "")
} else {
  sprintf("MA: custom contrast%s", if (shrunken) " [apeglm]" else "")
}
title(main = ttl)
dev.off()
