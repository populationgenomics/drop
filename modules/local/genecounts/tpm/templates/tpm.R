#!/usr/bin/env Rscript
# CPG: gene-level TPM from DROP's merged raw counts, in the long format seqr's TPM upload reads
# (sample_id, gene_id, TPM). See docs/CPG_METAMIST_SEQR_DESIGN.md (T3.3).

suppressPackageStartupMessages({
    library(data.table)
    library(SummarizedExperiment)
})

prefix <- ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix')

total_counts <- readRDS("${counts}")
counts <- assay(total_counts)

# union-exon length per gene: rowRanges holds each gene's exons from the GTF
gene_lengths <- sum(width(reduce(rowRanges(total_counts))))
stopifnot(identical(names(gene_lengths), rownames(counts)))

rpk <- counts / gene_lengths * 1e3
library_rpk <- colSums(rpk)
if (any(library_rpk == 0)) {
    stop(paste("No counts for:", paste(colnames(counts)[library_rpk == 0], collapse = ", ")))
}
tpm <- t(t(rpk) / library_rpk) * 1e6

# genes with no reads in any sample are 0 everywhere; seqr would only store them as zeros
tpm <- tpm[rowSums(counts) > 0, , drop = FALSE]

fwrite(
    data.table(
        sample_id = rep(colnames(tpm), each = nrow(tpm)),
        gene_id = rep(rownames(tpm), times = ncol(tpm)),
        TPM = signif(as.vector(tpm), 6)
    ),
    file = paste0(prefix, ".tsv.gz"),
    sep = "\t", quote = FALSE, compress = "gzip"
)

## VERSIONS FILE
writeLines(
    c(
        '"${task.process}":',
        paste('    r-base:', strsplit(version[['version.string']], ' ')[[1]][3]),
        paste('    r-data.table:', as.character(packageVersion('data.table'))),
        paste('    bioconductor-summarizedexperiment:', as.character(packageVersion('SummarizedExperiment')))
    ),
'versions.yml')