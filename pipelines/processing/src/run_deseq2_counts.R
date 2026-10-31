args <- commandArgs(trailingOnly = TRUE)

read_arg <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    stop(sprintf("Missing required argument: %s", flag), call. = FALSE)
  }
  args[[idx + 1]]
}

read_table_auto <- function(path) {
  lower <- tolower(path)
  if (grepl("\\.csv(\\.gz)?$", lower)) {
    return(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

split_columns <- function(value) {
  columns <- trimws(unlist(strsplit(value, ",", fixed = TRUE)))
  columns[nzchar(columns)]
}

counts_path <- read_arg("--counts")
gene_id_column <- read_arg("--gene-id-column")
control_columns <- split_columns(read_arg("--control-columns"))
treated_columns <- split_columns(read_arg("--treated-columns"))
output_path <- read_arg("--output")

suppressPackageStartupMessages(library(DESeq2))

counts_df <- read_table_auto(counts_path)
all_columns <- c(control_columns, treated_columns)

if (!(gene_id_column %in% colnames(counts_df))) {
  stop(sprintf("Gene ID column not found: %s", gene_id_column), call. = FALSE)
}
missing_columns <- setdiff(all_columns, colnames(counts_df))
if (length(missing_columns) > 0) {
  stop(
    sprintf("Count matrix is missing columns: %s", paste(missing_columns, collapse = ", ")),
    call. = FALSE
  )
}

gene_ids <- trimws(as.character(counts_df[[gene_id_column]]))
valid_rows <- nzchar(gene_ids) & !is.na(gene_ids)
counts_df <- counts_df[valid_rows, , drop = FALSE]
gene_ids <- gene_ids[valid_rows]

counts_only <- counts_df[, all_columns, drop = FALSE]
for (column in all_columns) {
  counts_only[[column]] <- as.numeric(counts_only[[column]])
}
if (any(is.na(as.matrix(counts_only)))) {
  stop("Count matrix contains non-numeric values in sample columns.", call. = FALSE)
}
if (any(as.matrix(counts_only) < 0)) {
  stop("Count matrix contains negative values.", call. = FALSE)
}

count_matrix <- as.matrix(round(counts_only))
rownames(count_matrix) <- gene_ids
count_matrix <- rowsum(count_matrix, group = rownames(count_matrix), reorder = FALSE)

condition <- factor(
  c(rep("control", length(control_columns)), rep("treated", length(treated_columns))),
  levels = c("control", "treated")
)
col_data <- data.frame(row.names = all_columns, condition = condition)

design_formula <- ~ condition
if ("--sample-covariates" %in% args) {
  covariates <- read_table_auto(read_arg("--sample-covariates"))
  if (!all(c("sample_id", "batch") %in% colnames(covariates)) ||
      anyDuplicated(covariates$sample_id) || !setequal(covariates$sample_id, all_columns)) {
    stop("Sample covariates must provide one batch per biological sample.", call. = FALSE)
  }
  batches <- covariates$batch[match(all_columns, covariates$sample_id)]
  if (anyNA(batches) || any(!nzchar(trimws(batches)))) {
    stop("Missing sample batch.", call. = FALSE)
  }
  col_data$batch <- factor(batches)
  if (nlevels(col_data$batch) > 1) design_formula <- ~ batch + condition
}
design_matrix <- model.matrix(design_formula, col_data)
if (qr(design_matrix)$rank < ncol(design_matrix) || nrow(design_matrix) <= ncol(design_matrix)) {
  stop("DE design is confounded or has no residual replication.", call. = FALSE)
}

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = design_formula)
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treated", "control"))
raw_counts <- counts(dds, normalized = FALSE)
normalized_counts <- counts(dds, normalized = TRUE)
control_mean_normalized_count <- rowMeans(
  normalized_counts[, control_columns, drop = FALSE]
)

out <- data.frame(
  gene_id = sub("\\.[0-9]+$", "", rownames(res)),
  baseMean = res$baseMean,
  logFC = res$log2FoldChange,
  control_mean_normalized_count = control_mean_normalized_count[rownames(res)],
  lfcSE = res$lfcSE,
  stat = res$stat,
  PValue = res$pvalue,
  FDR = res$padj,
  raw_control_mean = rowMeans(raw_counts[, control_columns, drop = FALSE]),
  raw_treated_mean = rowMeans(raw_counts[, treated_columns, drop = FALSE]),
  normalized_control_mean = rowMeans(normalized_counts[, control_columns, drop = FALSE]),
  normalized_treated_mean = rowMeans(normalized_counts[, treated_columns, drop = FALSE]),
  stringsAsFactors = FALSE
)
out <- out[order(out$PValue, na.last = TRUE), ]

write.table(out, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE)

# Keep uncertainty and all input genes available for a future benchmark redesign.
# The canonical DE output and raw, unshrunk estimates above are unchanged.
output_dir <- dirname(output_path)
native <- data.frame(gene_id = rownames(res), as.data.frame(res), dispersion = dispersions(dds))
write.table(native, file = file.path(output_dir, "native_de_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
all_normalized <- sweep(count_matrix, 2, sizeFactors(dds), "/")
write.table(data.frame(gene_id = rownames(all_normalized), all_normalized, check.names = FALSE),
            file = file.path(output_dir, "normalized_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample_id = all_columns, condition = condition,
                       size_factor = sizeFactors(dds)),
            file = file.path(output_dir, "sample_size_factors.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
control_values <- all_normalized[, control_columns, drop = FALSE]
treated_values <- all_normalized[, treated_columns, drop = FALSE]
summary <- data.frame(
  gene_id = rownames(all_normalized),
  control_n = length(control_columns), treated_n = length(treated_columns),
  control_mean = rowMeans(control_values), control_sd = apply(control_values, 1, sd),
  control_sem = apply(control_values, 1, sd) / sqrt(length(control_columns)),
  treated_mean = rowMeans(treated_values), treated_sd = apply(treated_values, 1, sd),
  treated_sem = apply(treated_values, 1, sd) / sqrt(length(treated_columns)),
  control_raw_mean = rowMeans(count_matrix[, control_columns, drop = FALSE]),
  treated_raw_mean = rowMeans(count_matrix[, treated_columns, drop = FALSE]),
  all_zero = rowSums(count_matrix) == 0
)
write.table(summary, file = file.path(output_dir, "expression_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


write_tsv <- function(value, filename) {
  write.table(value, file.path(output_dir, filename), sep = "\t",
              quote = FALSE, row.names = FALSE, na = "NA")
}

# These diagnostic results explain NA values. They do not replace primary results,
# and do not reverse any automatic count replacement done during DESeq().
diagnostic_results <- results(dds, contrast = c("condition", "treated", "control"),
                              cooksCutoff = FALSE, independentFiltering = FALSE)
status <- rep("tested_with_adjusted_pvalue", nrow(res))
status[!is.na(res$pvalue) & is.na(res$padj)] <- "independent_filtering"
status[is.na(res$pvalue) & !is.na(diagnostic_results$pvalue)] <- "cooks_results_filter"
status[is.na(res$pvalue) & is.na(diagnostic_results$pvalue)] <- "pvalue_missing_other"
all_status <- rep("all_zero_not_fitted", nrow(count_matrix))
all_status[match(rownames(res), rownames(count_matrix))] <- status
all_native <- native[match(rownames(count_matrix), native$gene_id), , drop = FALSE]
all_native$gene_id <- rownames(count_matrix)
all_native$analysis_status <- all_status
write_tsv(all_native, "all_gene_statistics.tsv")
write_tsv(data.frame(gene_id = rownames(dds), as.data.frame(mcols(dds)),
                     check.names = FALSE), "model_gene_diagnostics.tsv")
write_tsv(data.frame(gene_id = rownames(diagnostic_results),
                     as.data.frame(diagnostic_results)),
          "diagnostic_results_without_results_filters.tsv")

summary$control_median <- apply(control_values, 1, median)
summary$treated_median <- apply(treated_values, 1, median)
summary$control_min <- apply(control_values, 1, min)
summary$control_max <- apply(control_values, 1, max)
summary$treated_min <- apply(treated_values, 1, min)
summary$treated_max <- apply(treated_values, 1, max)
summary$control_nonzero_samples <- rowSums(count_matrix[, control_columns, drop = FALSE] > 0)
summary$treated_nonzero_samples <- rowSums(count_matrix[, treated_columns, drop = FALSE] > 0)
summary$control_raw_sd <- apply(count_matrix[, control_columns, drop = FALSE], 1, sd)
summary$treated_raw_sd <- apply(count_matrix[, treated_columns, drop = FALSE], 1, sd)
summary$analysis_status <- all_status
write_tsv(summary, "expression_summary.tsv")
sample_summary <- data.frame(sample_id = all_columns, as.data.frame(colData(dds)),
                              raw_library_sum = colSums(count_matrix),
                              genes_nonzero = colSums(count_matrix > 0),
                              normalized_library_sum = colSums(all_normalized),
                              check.names = FALSE)
write_tsv(sample_summary, "sample_statistics.tsv")
saveRDS(list(dds = dds, results = res,
             diagnostic_results_without_results_filters = diagnostic_results,
             all_input_counts = count_matrix, sample_metadata = col_data,
             control_columns = control_columns, treated_columns = treated_columns,
             session_info = sessionInfo()),
        file.path(output_dir, "complete_analysis.rds"), compress = "gzip")
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
capture.output(list(design = design(dds), contrast = c("condition", "treated", "control"),
                    results_metadata = metadata(res), assays = assayNames(dds),
                    coefficient_names = resultsNames(dds),
                    note = "Descriptive control normalization uses both groups; SD/SEM are not lfcSE."),
               file = file.path(output_dir, "analysis_metadata.txt"))
