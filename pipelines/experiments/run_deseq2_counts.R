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

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = col_data, design = ~ condition)
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
