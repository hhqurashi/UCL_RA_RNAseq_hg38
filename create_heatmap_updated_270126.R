#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(pheatmap)
})

## -----------------------------
## USER EDITS (D280 + P1 exclusion)
## -----------------------------
PROJ_OUT_DIR <- "/mnt/scratch/hqurashi/13_Kelly_RNA_Seq/outputs/D280"

# DESeq2 outputs (from updated run_DESeq2.R)
DESEQ2_BASE <- file.path(PROJ_OUT_DIR, "DESeq2", "DESeq2_pairs_hg38_tximport_shrinkLFC_gene_names")

A6P_DE_SIG      <- file.path(DESEQ2_BASE, "A6_vs_P", "A6_vs_P_DE_sig_padj_lfc.csv")
A6P_DE_SIG_CHRX <- file.path(DESEQ2_BASE, "A6_vs_P", "A6_vs_P_DE_sig_padj_lfc_chrX.csv")

# Optional: used only for putting * on the H column vs P (if present)
# Leave "" to auto-detect.
PVSH_DE_RESULTS <- ""  # e.g. "/.../P_vs_H/P_vs_H_DE_results.csv"

# TPM file from nf-core/rnaseq outputs
TPM_TSV <- file.path(PROJ_OUT_DIR, "star_salmon", "salmon.merged.gene_tpm.tsv")

# Exclude outlier replicate from TPM mean calculations
EXCLUDE_SAMPLES <- c("P1")

OUT_BASE <- DESEQ2_BASE
OUT_DIR  <- file.path(OUT_BASE, "heatmaps_A6_vs_P_top100_meanTPM_excluding_P1")
if (!dir.exists(OUT_DIR)) {
  ok <- dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  if (!ok) stop("Failed to create output directory: ", OUT_DIR)
}

TOP_N      <- 100
ALPHA      <- 0.05
LFC_CUTOFF <- 1
STAR_REQUIRE_LFC <- TRUE  # TRUE = (padj < ALPHA & |LFC|>=LFC_CUTOFF), FALSE = padj-only

## Plot sizing / text
ROW_FONTSIZE  <- 7
COL_FONTSIZE  <- 12
STAR_FONTSIZE <- 10
CELLHEIGHT_PT <- 12
CELLWIDTH_PT  <- 60
SCALE         <- 1.4

## -----------------------------
## Helpers
## -----------------------------
warn_if_ensg <- function(x, what = "gene IDs") {
  if (any(grepl("^ENSG", x))) {
    message("NOTE: Found ENSG IDs in ", what, " (this is fine for the updated nf-core reference outputs).")
  }
}

padj_fix <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 1
  pmax(x, .Machine$double.xmin)
}

sig_call <- function(padj, lfc) {
  if (STAR_REQUIRE_LFC) {
    (padj < ALPHA) & (abs(lfc) >= LFC_CUTOFF)
  } else {
    (padj < ALPHA)
  }
}

read_auto <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(path, data.table = FALSE)))
  }
  df <- tryCatch(
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (!is.null(df) && ncol(df) > 1) return(df)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

infer_group_from_col <- function(colname) {
  if (grepl("^A6_", colname)) return("A6")
  if (grepl("^H",  colname))  return("H")
  if (grepl("^P",  colname))  return("P")
  NA_character_
}

row_zscore <- function(m) {
  mu  <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  (m - mu) / sdv
}

## -----------------------------
## Load DE sig list: A6 vs P (already padj+LFC filtered)
## -----------------------------
de_sig <- read_auto(A6P_DE_SIG)
need_cols <- c("gene_id", "log2FoldChange", "padj")
miss <- setdiff(need_cols, colnames(de_sig))
if (length(miss) > 0) stop("A6_vs_P_DE_sig_padj_lfc.csv missing columns: ", paste(miss, collapse = ", "))

de_sig$gene <- as.character(de_sig$gene_id)
warn_if_ensg(de_sig$gene, "A6_vs_P sig gene_id")
de_sig$lfc  <- suppressWarnings(as.numeric(de_sig$log2FoldChange))
de_sig$padj <- padj_fix(de_sig$padj)

## -----------------------------
## chrX gene list (from chrX sig file)
## -----------------------------
chrx <- read_auto(A6P_DE_SIG_CHRX)
if (!("gene_id" %in% colnames(chrx))) stop("chrX sig file missing gene_id column: ", A6P_DE_SIG_CHRX)

chrx_genes <- unique(as.character(chrx$gene_id))
chrx_genes <- chrx_genes[!is.na(chrx_genes) & nzchar(chrx_genes)]
warn_if_ensg(chrx_genes, "chrX gene_id list")

## -----------------------------
## Optional DE results for P vs H (for stars on H column vs P)
## -----------------------------
auto_candidates <- c(
  file.path(dirname(A6P_DE_SIG), "..", "P_vs_H", "P_vs_H_DE_results.csv"),
  file.path(dirname(A6P_DE_SIG), "..", "H_vs_P", "H_vs_P_DE_results.csv")
)
auto_candidates <- normalizePath(auto_candidates, winslash = "/", mustWork = FALSE)

if (!nzchar(PVSH_DE_RESULTS)) {
  cand <- auto_candidates[file.exists(auto_candidates)][1]
  if (!is.na(cand)) PVSH_DE_RESULTS <- cand else PVSH_DE_RESULTS <- ""
}

de_pvh <- NULL
if (nzchar(PVSH_DE_RESULTS) && file.exists(PVSH_DE_RESULTS)) {
  tmp <- read_auto(PVSH_DE_RESULTS)
  if (all(c("gene_id", "log2FoldChange", "padj") %in% colnames(tmp))) {
    tmp$gene <- as.character(tmp$gene_id)
    warn_if_ensg(tmp$gene, "P_vs_H/H_vs_P gene_id")
    tmp$lfc  <- suppressWarnings(as.numeric(tmp$log2FoldChange))
    tmp$padj <- padj_fix(tmp$padj)
    de_pvh <- tmp
  } else {
    warning("PVSH_DE_RESULTS found but lacks required columns; H-vs-P stars will be blank.")
  }
} else {
  message("No P_vs_H/H_vs_P DE results found; H column will have no stars.")
}

## -----------------------------
## Load TPM + set up sample columns (excluding P1)
## -----------------------------
tpm <- read_auto(TPM_TSV)
if (!all(c("gene_id", "gene_name") %in% colnames(tpm))) {
  stop("TPM file must have columns: gene_id, gene_name, then sample columns.")
}
tpm$gene_id   <- as.character(tpm$gene_id)
tpm$gene_name <- as.character(tpm$gene_name)
warn_if_ensg(tpm$gene_id, "TPM gene_id")

sample_cols_all <- setdiff(colnames(tpm), c("gene_id", "gene_name"))
sample_cols <- setdiff(sample_cols_all, EXCLUDE_SAMPLES)

dropped <- setdiff(sample_cols_all, sample_cols)
if (length(dropped) > 0) message("Excluding samples from TPM means: ", paste(dropped, collapse = ", "))

groups <- vapply(sample_cols, infer_group_from_col, character(1))
keep <- !is.na(groups)
sample_cols <- sample_cols[keep]
groups <- groups[keep]

if (length(sample_cols) == 0) stop("No TPM sample columns detected by prefixes A6_/H*/P* (after exclusions).")

group_levels <- c("H", "A6", "P")
for (g in group_levels) {
  if (!any(groups == g)) stop("No samples detected for group: ", g, " (after exclusions).")
}

mean_tpm_for_genes <- function(gene_vec) {
  df <- tpm[tpm$gene_id %in% gene_vec, c("gene_id", "gene_name", sample_cols), drop = FALSE]
  df <- df[match(gene_vec, df$gene_id), , drop = FALSE]

  mat <- as.matrix(df[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- df$gene_id

  group_means <- sapply(group_levels, function(g) {
    cols <- sample_cols[groups == g]
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  })

  logm <- log2(group_means + 1)
  zmat <- row_zscore(logm)

  list(mean_tpm = group_means, log2_mean = logm, rowz = zmat)
}

make_stars_vs_P <- function(gene_vec) {
  stars <- matrix("", nrow = length(gene_vec), ncol = length(group_levels),
                  dimnames = list(gene_vec, group_levels))

  # A6 vs P: genes are from sig list -> always "*"
  stars[, "A6"] <- "*"

  # H vs P: star if significant (if available)
  if (!is.null(de_pvh)) {
    idx <- match(gene_vec, de_pvh$gene)
    padj <- de_pvh$padj[idx]
    lfc  <- de_pvh$lfc[idx]
    stars[, "H"] <- ifelse(sig_call(padj, lfc), "*", "")
  } else {
    stars[, "H"] <- ""
  }

  stars[, "P"] <- ""
  stars
}

## -----------------------------
## Convert ENSG rownames to SYMBOL labels for display
## -----------------------------
map_rows_to_symbols <- function(mat, stars_mat = NULL) {
  ensg_ids <- rownames(mat)
  sym <- as.character(tpm$gene_name[match(ensg_ids, tpm$gene_id)])
  sym[is.na(sym) | !nzchar(sym)] <- ensg_ids[is.na(sym) | !nzchar(sym)]
  sym_unique <- make.unique(sym)

  mat2 <- mat
  rownames(mat2) <- sym_unique

  out <- list(mat = mat2, stars = stars_mat)

  if (!is.null(stars_mat)) {
    stars2 <- stars_mat
    rownames(stars2) <- sym_unique
    out$stars <- stars2
  }

  # helpful mapping table
  map_df <- data.frame(
    gene_id = ensg_ids,
    gene_name = sym,
    gene_name_unique = sym_unique,
    stringsAsFactors = FALSE
  )

  out$map_df <- map_df
  out
}

## -----------------------------
## Heatmap runner (writes CSVs + PNG + PDF)
## -----------------------------
run_heatmap <- function(gene_vec, tag) {
  gene_vec <- unique(gene_vec[!is.na(gene_vec) & nzchar(gene_vec)])
  if (length(gene_vec) == 0) stop(tag, ": gene list is empty.")

  ex <- mean_tpm_for_genes(gene_vec)
  zmat  <- ex$rowz
  stars <- make_stars_vs_P(rownames(zmat))

  # Display labels as gene symbols (from TPM gene_name)
  disp <- map_rows_to_symbols(zmat, stars)
  zmat_disp  <- disp$mat
  stars_disp <- disp$stars

  excl_tag <- if (length(dropped) > 0) paste0("_TPM_excl_", paste(dropped, collapse = "_")) else ""
  out_prefix <- file.path(OUT_DIR, paste0("A6_vs_P_top", TOP_N, "_", tag, excl_tag))

  # write tables
  write.csv(data.frame(gene_id = rownames(zmat), stringsAsFactors = FALSE),
            paste0(out_prefix, "_genes.csv"), row.names = FALSE)
  write.csv(disp$map_df, paste0(out_prefix, "_gene_id_to_name_map.csv"), row.names = FALSE)

  write.csv(ex$mean_tpm,  paste0(out_prefix, "_meanTPM.csv"))
  write.csv(ex$log2_mean, paste0(out_prefix, "_log2_meanTPM.csv"))
  write.csv(zmat,         paste0(out_prefix, "_rowZ_log2_meanTPM.csv"))            # ENSG rownames
  write.csv(zmat_disp,    paste0(out_prefix, "_rowZ_log2_meanTPM_SYMBOLS.csv"))    # SYMBOL rownames
  write.csv(stars,        paste0(out_prefix, "_stars_vs_P.csv"))                   # ENSG rownames
  write.csv(stars_disp,   paste0(out_prefix, "_stars_vs_P_SYMBOLS.csv"))           # SYMBOL rownames

  cols <- colorRampPalette(c("blue", "white", "red"))(100)

  n_genes <- nrow(zmat_disp)
  base_w_px <- 1100
  base_h_px <- 1600
  res_png   <- 220
  base_h_in <- base_h_px / res_png
  required_h_in <- (n_genes * CELLHEIGHT_PT / 72) + 4.5
  scale2 <- max(SCALE, required_h_in / base_h_in)

  png_w_px <- round(base_w_px * scale2)
  png_h_px <- round(base_h_px * scale2)

  main_title <- paste0(
    "Top ", TOP_N, " significant DEGs: A6 vs P (", tag, ")\n",
    "Ranked by |log2FC| from A6_vs_P_DE_sig_padj_lfc.csv\n",
    "Values: row z-score of log2(mean TPM + 1) across H/A6/P\n",
    "* = significantly different vs P\n",
    "TPM means exclude: ", ifelse(length(dropped) > 0, paste(dropped, collapse = ", "), "none")
  )

  # PNG
  png(paste0(out_prefix, ".png"), width = png_w_px, height = png_h_px, res = res_png)
  pheatmap(
    zmat_disp,
    color = cols,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    display_numbers = stars_disp,
    number_color = "black",
    fontsize_number = STAR_FONTSIZE,
    fontsize_row = ROW_FONTSIZE,
    fontsize_col = COL_FONTSIZE,
    cellheight = CELLHEIGHT_PT,
    cellwidth  = CELLWIDTH_PT,
    main = main_title
  )
  dev.off()

  # PDF
  pdf(paste0(out_prefix, ".pdf"), width = 7 * scale2, height = 10 * scale2)
  pheatmap(
    zmat_disp,
    color = cols,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    border_color = NA,
    display_numbers = stars_disp,
    number_color = "black",
    fontsize_number = STAR_FONTSIZE,
    fontsize_row = ROW_FONTSIZE,
    fontsize_col = COL_FONTSIZE,
    cellheight = CELLHEIGHT_PT,
    cellwidth  = 14,
    main = main_title
  )
  dev.off()

  message("Wrote heatmap: ", tag, " -> ", out_prefix, ".png/.pdf")
}

## -----------------------------
## Build TOP_N gene lists
## -----------------------------
# Includes chrX: top by abs(LFC), then padj
de_inclX <- de_sig[order(-abs(de_sig$lfc), de_sig$padj), , drop = FALSE]
top_genes_inclX <- head(unique(de_inclX$gene), TOP_N)
if (length(top_genes_inclX) < TOP_N) warning("Only found ", length(top_genes_inclX), " sig genes for includes-chrX heatmap.")

# Excludes chrX: drop chrX genes then rank
de_noX <- de_sig[!(de_sig$gene %in% chrx_genes), , drop = FALSE]
if (nrow(de_noX) == 0) stop("After excluding chrX, there are 0 significant A6 vs P genes left.")
de_noX <- de_noX[order(-abs(de_noX$lfc), de_noX$padj), , drop = FALSE]
top_genes_exclX <- head(unique(de_noX$gene), TOP_N)
if (length(top_genes_exclX) < TOP_N) warning("Only found ", length(top_genes_exclX), " sig genes for excludes-chrX heatmap.")

## -----------------------------
## Run BOTH heatmaps
## -----------------------------
run_heatmap(top_genes_inclX, "includes_chrX")
run_heatmap(top_genes_exclX, "excludes_chrX")

message("Done. Outputs in: ", OUT_DIR)