#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(conflicted)
  library(biomaRt)
  library(tidyverse)
  library(cli)
})

parse_args <- function(args) {
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument '%s'. Expected '--key value' pairs.", key))
    }
    if (i == length(args)) {
      stop(sprintf("No value provided for argument '%s'", key))
    }
    value <- args[[i + 1]]
    parsed[[substring(key, 3)]] <- value
    i <- i + 2
  }
  parsed
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("counts", "output", "species")
missing <- setdiff(required, names(args))
if (length(missing) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")))
}

counts_path <- args[["counts"]]
output_path <- args[["output"]]
species <- args[["species"]]

mart <- "useast"
rounds <- 0
while (class(mart)[[1]] != "Mart") {
  mart <- tryCatch(
    {
      if (mart == "www") rounds <- rounds + 1
      biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = str_c(species, "_gene_ensembl"),
        mirror = mart
      )
    },
    error = function(e) {
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
      mart <<- switch(
        mart,
        useast = "uswest",
        uswest = "asia",
        asia = "www",
        www = {
          Sys.sleep(30)
          "useast"
        }
      )
    }
  )
}

df <- read.table(counts_path, sep='\t', header=TRUE)

g2g <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = df$gene,
  mart = mart
)

annotated <- merge(df, g2g, by.x="gene", by.y="ensembl_gene_id")
annotated$gene <- ifelse(annotated$external_gene_name == '', annotated$gene, annotated$external_gene_name)
annotated$external_gene_name <- NULL
write.table(annotated, output_path, sep='\t', row.names=FALSE, quote=FALSE)
