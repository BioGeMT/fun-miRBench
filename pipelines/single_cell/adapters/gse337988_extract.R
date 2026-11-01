# Cell metadata and authors' guide calls (perturbationIndicator) of one GSE337988 sub-library object, without SingleCellExperiment.
# Usage: Rscript gse337988_extract.R <object_se.rds> <out dir>
suppressPackageStartupMessages({library(S4Vectors); library(Matrix); library(SummarizedExperiment)})
args <- commandArgs(TRUE)
slots <- attributes(readRDS(args[1]))
out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write.table(as.data.frame(slots$colData), file.path(out, "cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
indicator <- attr(attr(attr(slots$int_colData@listData$altExps@listData$perturbationExperiment, "se"), "assays"), "data")@listData[[1]]
writeMM(as(as(indicator, "dMatrix"), "CsparseMatrix"), file.path(out, "guides.mtx"))
writeLines(rownames(indicator), file.path(out, "guides.tsv"))
writeLines(colnames(indicator), file.path(out, "guide_cells.tsv"))
writeLines(dimnames(attr(slots$assays, "data")@listData[[1]])[[1]], file.path(out, "genes.tsv"))
