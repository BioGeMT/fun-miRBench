# Raw counts and metadata of a Seurat object, read without Seurat (attributes only), for selected cells.
# Usage: Rscript seurat_counts.R <object.rds|object.RData.gz> <assay> <meta column> <file with one kept value per line, or -> <out dir> [object name in RData]
# SEURAT_LAYER=data takes the log-normalised layer instead (no raw layer in the object).
# Writes <out dir>/counts.mtx (genes x cells), features.tsv, barcodes.tsv, meta.tsv.
args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(Matrix))
path <- args[1]; assay <- args[2]; column <- args[3]; keep_file <- args[4]; out <- args[5]
if (grepl("\\.rds(\\.gz)?$", path, ignore.case = TRUE)) {
  object <- readRDS(path)
} else {
  env <- new.env(); load(path, envir = env)
  object <- get(if (length(args) >= 6) args[6] else ls(env)[1], envir = env)
}
slots <- attributes(object)
meta <- slots$meta.data
piece <- slots$assays[[assay]]
parts <- attributes(piece)
layer <- Sys.getenv("SEURAT_LAYER", "counts")
counts <- if (!is.null(parts$layers)) parts$layers[[layer]] else parts[[layer]]
if (is.null(counts)) stop("no counts layer")
if (is.null(rownames(counts))) {
  rownames(counts) <- attr(parts$features, "dimnames")[[1]]
  colnames(counts) <- attr(parts$cells, "dimnames")[[1]]
}
if (keep_file != "-") {
  wanted <- readLines(keep_file)
  chosen <- rownames(meta)[as.character(meta[[column]]) %in% wanted]
  chosen <- intersect(chosen, colnames(counts))
} else {
  chosen <- colnames(counts)
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
sub <- counts[, chosen, drop = FALSE]
if (layer == "counts") stopifnot(all(sub@x == round(sub@x)))
writeMM(sub, file.path(out, "counts.mtx"))
writeLines(rownames(sub), file.path(out, "features.tsv"))
writeLines(colnames(sub), file.path(out, "barcodes.tsv"))
write.table(cbind(cell = chosen, meta[chosen, , drop = FALSE]), file.path(out, "meta.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(sub), "genes", ncol(sub), "cells\n")
