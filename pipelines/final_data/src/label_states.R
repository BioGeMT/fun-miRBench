# Final per-gene table for one comparison, from its saved DESeq2 analysis (complete_analysis.rds).
#
# Every gene of the benchmark gene list gets one row. Genes the fit did not include (zero counts in
# every sample) are kept, with no statistics and the reason.
#
#  1. Standard test against 0 and the one-sided tests against each margin d (OE: altHypothesis
#     "less", KO: "greater") are run on the fitted object restricted to the listed genes, so DESeq2's
#     independent filtering and the Benjamini-Hochberg FDR use the listed genes only. Cook's
#     filtering is DESeq2's default.
#  2. apeglm shrinkage uses the full fit (its prior is estimated from all fitted genes).
#  3. FSOS ("false sign or small") s-values for each d of the grid, over the listed genes: the per-gene
#     probability that the change is smaller than d or of the other sign, pnorm((d - |shrunken log2FC|) /
#     shrunken SD), which is what lfcShrink(lfcThreshold = d) gives, then apeglm::svalue.
#  4. States per d of the grid, on the shrunken oriented change o (-log2FC for OE, +log2FC for KO):
#       responsive                eligible, o > 0 and FSOS s-value < 0.05;
#       no_response_above_margin  eligible and o + 1.959964 * shrunken SD < d (whole 95% interval below d);
#       unresolved                everything else.
#     Eligible: finite log2FC, finite lfcSE > 0 and a finite p-value in the standard test, and a finite
#     shrunken log2FC with shrunken SD > 0.
#  5. Raw states at d = 0.5 and 1.0 (state_raw_*), kept as a check: responsive if the threshold-test FDR < 0.05,
#     no response if o + 1.959964 * lfcSE < d on the raw log2FC.
#
# Usage: Rscript label_states.R <complete_analysis.rds> <gene_list.tsv[.gz]> <OE|KO> <output_prefix>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4, file.exists(args[[1]]), file.exists(args[[2]]), args[[3]] %in% c("OE", "KO"))
options(warn = 1)
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})
analysis <- readRDS(args[[1]])
genes <- read.delim(args[[2]], stringsAsFactors = FALSE)
perturbation <- args[[3]]
prefix <- args[[4]]
MARGINS <- c(d05 = 0.5, d10 = 1.0)
D_GRID <- c(d00 = 0, d01 = 0.1, d025 = 0.25, d05 = 0.5, d075 = 0.75, d10 = 1.0, d15 = 1.5, d20 = 2.0, d30 = 3.0)
Z95 <- qnorm(0.975)

dds <- analysis$dds
coef_name <- "condition_treated_vs_control"
stopifnot(identical(levels(dds$condition), c("control", "treated")), coef_name %in% resultsNames(dds))
rownames(dds) <- sub("\\.[0-9]+$", "", rownames(dds))
stopifnot(!anyDuplicated(rownames(dds)), !anyDuplicated(genes$gene_id))

events <- data.frame(step = character(), kind = character(), message = character())
capture <- function(step, expr) {
  withCallingHandlers(expr,
    warning = function(w) {
      events[nrow(events) + 1, ] <<- c(step, "warning", conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      events[nrow(events) + 1, ] <<- c(step, "message", trimws(conditionMessage(m)))
      invokeRestart("muffleMessage")
    })
}

listed <- intersect(genes$gene_id, rownames(dds))
sub_dds <- dds[listed, ]
std <- results(sub_dds, name = coef_name, alpha = 0.05)
std_raw <- results(sub_dds, name = coef_name, cooksCutoff = FALSE, independentFiltering = FALSE)
shrunk_list <- capture("apeglm", lfcShrink(dds, coef = coef_name, type = "apeglm", svalue = TRUE,
                                           returnList = TRUE, quiet = FALSE))
shrunk <- shrunk_list$res[listed, ]
diag <- shrunk_list$fit$diag[listed, , drop = FALSE]

norm <- counts(sub_dds, normalized = TRUE)
raw <- counts(sub_dds)
ctrl <- sub_dds$condition == "control"
trt <- sub_dds$condition == "treated"

fit <- data.frame(
  gene_id = listed,
  baseMean = std$baseMean,
  logFC = std$log2FoldChange,
  lfcSE = std$lfcSE,
  stat = std$stat,
  PValue = std$pvalue,
  FDR = std$padj,
  control_mean_normalized_count = rowMeans(norm[, ctrl, drop = FALSE]),
  treated_mean_normalized_count = rowMeans(norm[, trt, drop = FALSE]),
  control_nonzero_samples = rowSums(raw[, ctrl, drop = FALSE] > 0),
  treated_nonzero_samples = rowSums(raw[, trt, drop = FALSE] > 0),
  logFC_shrunk = shrunk$log2FoldChange,
  lfcSE_shrunk = shrunk$lfcSE,
  svalue = shrunk$svalue,
  apeglm_conv = diag[, "conv"],
  stringsAsFactors = FALSE
)
fit$analysis_status <- "tested_with_adjusted_pvalue"
fit$analysis_status[!is.na(std$pvalue) & is.na(std$padj)] <- "independent_filtering"
fit$analysis_status[is.na(std$pvalue) & !is.na(std_raw$pvalue)] <- "cooks_outlier"
fit$analysis_status[is.na(std$pvalue) & is.na(std_raw$pvalue)] <- "pvalue_missing_other"

alt <- if (perturbation == "OE") "less" else "greater"
sign <- if (perturbation == "OE") -1 else 1
threshold_info <- list()
for (tag in names(MARGINS)) {
  d <- MARGINS[[tag]]
  thr <- capture(paste0("threshold_", tag),
                 results(sub_dds, name = coef_name, lfcThreshold = d, altHypothesis = alt, alpha = 0.05))
  stopifnot(identical(rownames(thr), listed))
  fit[[paste0("thr_pvalue_", tag)]] <- thr$pvalue
  fit[[paste0("thr_padj_", tag)]] <- thr$padj
  threshold_info[[tag]] <- list(lfcThreshold = d, altHypothesis = alt, filterThreshold = metadata(thr)$filterThreshold,
                                n_tested = sum(!is.na(thr$padj)))
}

out <- merge(genes[, c("gene_id", "gene_name")], fit, by = "gene_id", all.x = TRUE, sort = FALSE)
out <- out[match(genes$gene_id, out$gene_id), ]
out$analysis_status[is.na(out$analysis_status)] <- "not_fitted_all_zero"
out$oriented_logFC <- sign * out$logFC
out$oriented_ci95_low <- out$oriented_logFC - Z95 * out$lfcSE
out$oriented_ci95_high <- out$oriented_logFC + Z95 * out$lfcSE
out$oriented_logFC_shrunk <- sign * out$logFC_shrunk
out$oriented_ci95_low_shrunk <- out$oriented_logFC_shrunk - Z95 * out$lfcSE_shrunk
out$oriented_ci95_high_shrunk <- out$oriented_logFC_shrunk + Z95 * out$lfcSE_shrunk
out$eligible <- is.finite(out$logFC) & is.finite(out$lfcSE) & out$lfcSE > 0 & is.finite(out$PValue) &
  is.finite(out$logFC_shrunk) & is.finite(out$lfcSE_shrunk) & out$lfcSE_shrunk > 0
out$ineligible_reason <- ifelse(out$eligible, "",
  ifelse(out$analysis_status %in% c("not_fitted_all_zero", "cooks_outlier", "pvalue_missing_other"),
         out$analysis_status, "no_usable_estimate"))
for (tag in names(D_GRID)) {
  d <- D_GRID[[tag]]
  small <- pnorm((d - abs(out$logFC_shrunk[out$eligible])) / out$lfcSE_shrunk[out$eligible])
  sval <- rep(NA_real_, nrow(out))
  sval[out$eligible] <- apeglm::svalue(small)
  responsive <- out$eligible & out$oriented_logFC_shrunk > 0 & sval < 0.05
  no_response <- out$eligible & out$oriented_ci95_high_shrunk < d
  stopifnot(!any(responsive & no_response))
  state <- rep("unresolved", nrow(out))
  state[responsive] <- "responsive"
  state[no_response] <- "no_response_above_margin"
  out[[paste0("svalue_", tag)]] <- sval
  out[[paste0("state_", tag)]] <- state
}
for (tag in names(MARGINS)) {
  d <- MARGINS[[tag]]
  padj <- out[[paste0("thr_padj_", tag)]]
  state <- rep("unresolved", nrow(out))
  state[out$eligible & !is.na(padj) & padj < 0.05] <- "responsive"
  state[out$eligible & out$oriented_ci95_high < d] <- "no_response_above_margin"
  stopifnot(!any(out$eligible & !is.na(padj) & padj < 0.05 & out$oriented_ci95_high < d))
  out[[paste0("state_raw_", tag)]] <- state
}

gz <- gzfile(paste0(prefix, ".tsv.gz"), "w")
write.table(out, gz, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
close(gz)
write.table(events, file = paste0(prefix, "_events.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
capture.output(list(design = design(dds), coefficient = coef_name, perturbation = perturbation,
                    n_listed = nrow(genes), n_listed_fitted = length(listed),
                    standard_test = list(alpha = 0.05, filterThreshold = metadata(std)$filterThreshold),
                    threshold_tests = threshold_info, prior_info = priorInfo(shrunk_list$res),
                    nonconverged_map_rows = sum(!is.na(out$apeglm_conv) & out$apeglm_conv != 0),
                    session = sessionInfo()),
               file = paste0(prefix, "_metadata.txt"))
