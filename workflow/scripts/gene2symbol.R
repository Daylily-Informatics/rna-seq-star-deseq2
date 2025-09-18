log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(conflicted)
library(biomaRt)
library(tidyverse)
# useful error messages upon aborting
library("cli")

# this variable holds a mirror name until
# useEnsembl succeeds ("www" is last, because
# of very frequent "Internal Server Error"s)
mart <- "useast"
rounds <- 0
while (class(mart)[[1]] != "Mart") {
  mart <- tryCatch(
    {
      # done here, because error function does not
      # modify outer scope variables, I tried
      if (mart == "www") rounds <- rounds + 1
      # equivalent to useMart, but you can choose
      # the mirror instead of specifying a host
      biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = str_c(snakemake@params[["species"]], "_gene_ensembl"),
        mirror = mart
      )
    },
    error = function(e) {
      # change or make configurable if you want more or
      # less rounds of tries of all the mirrors
      if (rounds >= 3) {
        cli_abort(
          str_c(
            "Have tried all 4 available Ensembl biomaRt mirrors ",
            rounds,
            " times. You might have a connection problem, or no mirror is responsive.\n",
            "The last error message was:\n",
            message(e)
          )
        )
      }
      # hop to next mirror
      mart <- switch(
        mart,
        useast = "uswest",
        uswest = "asia",
        asia = "www",
        www = {
          # wait before starting another round through the mirrors,
          # hoping that intermittent problems disappear
          Sys.sleep(30)
          "useast"
        }
      )
    }
  )
}




# ---- load counts ----
df <- read.table(
  snakemake@input[["counts"]],
  sep = "\t", header = TRUE, check.names = FALSE,
  quote = "", comment.char = ""
)

# explicit override or heuristic detection
gene_col <- snakemake@params[["gene_col"]]
if (is.null(gene_col) || !nzchar(gene_col) || !(gene_col %in% colnames(df))) {
  candidates <- c("gene","Gene","gene_id","Geneid","GeneID","ENSEMBL","Ensembl","EnsemblID")
  gene_col <- intersect(candidates, colnames(df))[1]
  if (is.na(gene_col)) {
    # pick a char-like column that looks like ENS*
    char_cols <- names(df)[vapply(df, function(x) is.character(x) || is.factor(x), logical(1))]
    ens_like <- char_cols[vapply(df[char_cols], function(x) any(grepl("^ENS[A-Z]*\\d+", as.character(x))), logical(1))]
    gene_col <- if (length(ens_like)) ens_like[1] else NA
  }
}
if (is.na(gene_col)) {
  cli::cli_abort(c(
    "!"="Could not identify a gene column.",
    "i"="Columns present: {toString(colnames(df))}",
    "i"="Set params[['gene_col']] or include a column named 'gene'/'gene_id'."
  ))
}

# normalize to 'gene' column
df <- dplyr::rename(df, gene = !!gene_col)
df$gene <- as.character(df$gene)

# STAR ReadsPerGene summary lines? drop if present
df <- df[!df$gene %in% c("N_unmapped","N_multimapping","N_noFeature","N_ambiguous"), , drop=FALSE]

# strip ENS version suffixes like '.7'
df$gene <- sub("\\.\\d+$", "", df$gene)

# drop blank/NA
df <- dplyr::filter(df, !is.na(gene) & gene != "")
if (nrow(df) == 0) {
  cli::cli_abort("After parsing, no usable gene IDs found in column '{gene_col}'.")
}

uniq_ids <- unique(df$gene)

# ---- biomaRt ----
g2g <- biomaRt::getBM(
  attributes = c("ensembl_gene_id","external_gene_name"),
  filters    = "ensembl_gene_id",
  values     = uniq_ids,
  mart       = mart
)

# ---- join + prefer symbols ----
annotated <- df |>
  dplyr::left_join(g2g, by = c("gene" = "ensembl_gene_id")) |>
  dplyr::mutate(gene = dplyr::if_else(is.na(external_gene_name) | external_gene_name == "", gene, external_gene_name)) |>
  dplyr::select(-external_gene_name)

write.table(annotated, snakemake@output[["symbol"]], sep = "\t", row.names = FALSE, quote = FALSE)

sink(type = "message")
sink()
close(log)
