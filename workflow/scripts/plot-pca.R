suppressPackageStartupMessages({
  library("DESeq2")
})

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

svg(snakemake@output[[1]])
if (ncol(dds) < 2) {
  plot.new()
  text(0.5, 0.5, "PCA requires at least two samples")
  dev.off()
  sink(type = "message")
  sink()
  close(log)
  quit(save = "no")
}

# obtain normalized counts
counts <- rlog(dds, blind = FALSE)
plotPCA(counts, intgroup = snakemake@wildcards[["variable"]])
dev.off()

sink(type = "message")
sink()
close(log)
