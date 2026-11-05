# Size factors from the comparison's control samples only.
#
# DESeq2's median-of-ratios formula, with the reference profile taken as the geometric mean of the control samples
# instead of all samples: a sample's factor is exp(median(log(count) - log reference)) over the genes counted above 0 in
# every control and in the sample. No factor depends on the treated samples, so the control profiles do not either.
# DESeq2 also rescales its factors to a geometric mean of 1 over all samples; that step is left out for the same reason.

SIZE_FACTOR_METHOD <- "median of ratios to the control samples"

# counts: genes x samples; control_columns: the control samples among colnames(counts).
control_size_factors <- function(counts, control_columns) {
  log_reference <- rowMeans(log(counts[, control_columns, drop = FALSE]))
  out <- data.frame(sample_id = colnames(counts), size_factor = NA_real_, reference_genes_used = 0L)
  for (j in seq_len(ncol(counts))) {
    used <- is.finite(log_reference) & counts[, j] > 0
    out$reference_genes_used[j] <- sum(used)
    out$size_factor[j] <- exp(median(log(counts[used, j]) - log_reference[used]))
  }
  if (!all(is.finite(out$size_factor))) {
    stop(sprintf("No gene is counted above 0 in every control sample and in: %s",
                 paste(out$sample_id[!is.finite(out$size_factor)], collapse = ", ")), call. = FALSE)
  }
  out
}
