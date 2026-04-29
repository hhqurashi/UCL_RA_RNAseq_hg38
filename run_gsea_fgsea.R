#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

## -----------------------------
## USER EDITS (D280 + exclude P1)
## -----------------------------
PROJ_OUT_DIR <- "/mnt/scratch/hqurashi/13_Kelly_RNA_Seq/outputs/D280"

# Where the DESeq2 pair folders live (from your updated run_DESeq2.R)
DESEQ2_BASE <- file.path(PROJ_OUT_DIR, "DESeq2", "DESeq2_pairs_hg38_tximport_shrinkLFC_gene_names")

# TPM from nf-core output
TPM_TSV <- file.path(PROJ_OUT_DIR, "star_salmon", "salmon.merged.gene_tpm.tsv")

# Exclude outlier replicate from TPM means / z-scores (affects TOP100 + ALLSIG heatmapZ split)
EXCLUDE_SAMPLES <- c("P1")

# A6 vs P significant (padj+LFC) + chrX sig list
A6P_DE_SIG      <- file.path(DESEQ2_BASE, "A6_vs_P", "A6_vs_P_DE_sig_padj_lfc.csv")
A6P_DE_SIG_CHRX <- file.path(DESEQ2_BASE, "A6_vs_P", "A6_vs_P_DE_sig_padj_lfc_chrX.csv")

# Full DE results (all genes) + chrX-only results (for preranked filtering)
A6P_DE_ALL      <- file.path(DESEQ2_BASE, "A6_vs_P", "A6_vs_P_DE_results.csv")
A6P_DE_ALL_CHRX <- file.path(DESEQ2_BASE, "A6_vs_P", "A6_vs_P_DE_results_chrX.csv")

HVP_DE_ALL      <- file.path(DESEQ2_BASE, "H_vs_P", "H_vs_P_DE_results.csv")
HVP_DE_ALL_CHRX <- file.path(DESEQ2_BASE, "H_vs_P", "H_vs_P_DE_results_chrX.csv")

OUT_BASE <- DESEQ2_BASE

# ===== chrX EXCLUDED outputs (as before) =====
OUT_GSEA_EXCL_TOP100   <- file.path(OUT_BASE, "GSEA_top100_A6vsP_excl_chrX_split_by_HEATMAP_Z")
OUT_GSEA_EXCL_ALLSIG   <- file.path(OUT_BASE, "GSEA_allSig_A6vsP_excl_chrX_split_by_HEATMAP_Z")
OUT_GSEA_EXCL_PRERANK  <- file.path(OUT_BASE, "GSEA_preranked_allGenes_excl_chrX_A6vsP_and_HvsP")

# ===== chrX INCLUDED counterparts =====
OUT_GSEA_INCL_TOP100   <- file.path(OUT_BASE, "GSEA_top100_A6vsP_incl_chrX_split_by_HEATMAP_Z")
OUT_GSEA_INCL_ALLSIG   <- file.path(OUT_BASE, "GSEA_allSig_A6vsP_incl_chrX_split_by_HEATMAP_Z")
OUT_GSEA_INCL_PRERANK  <- file.path(OUT_BASE, "GSEA_preranked_allGenes_incl_chrX_A6vsP_and_HvsP")

TOP_N <- 100

## fgsea params (for SMALL universes; keep unchanged)
MIN_SIZE   <- 3
MAX_SIZE   <- 500
TOP_DOT_N  <- 15
ALPHA_GSEA <- 0.05

## --- SPLIT LOGIC CONTROLS (heatmap z split) ---
Z_EPS <- 0
CONCORDANT_REQUIRES_NONZERO <- FALSE

## fgsea params for FULL preranked runs
PRERANK_MIN_SIZE  <- 15
PRERANK_MAX_SIZE  <- 500
PRERANK_TOP_DOT_N <- 15

set.seed(1)

## -----------------------------
## Helpers
## -----------------------------
dir_create_or_stop <- function(d) {
  if (!dir.exists(d)) {
    ok <- dir.create(d, recursive = TRUE, showWarnings = FALSE)
    if (!ok) stop("Failed to create directory: ", d)
  }
}

read_auto <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  as.data.frame(data.table::fread(path, data.table = FALSE))
}

strip_ens_version <- function(x) sub("\\..*$", "", as.character(x))

padj_fix <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 1
  pmax(x, .Machine$double.xmin)
}

infer_group_from_col <- function(colname) {
  if (grepl("^A6_", colname)) return("A6")
  if (grepl("^H",  colname))  return("H")
  if (grepl("^P",  colname))  return("P")
  NA_character_
}

to_pathway_list <- function(msig_tbl) split(msig_tbl$gene_symbol, msig_tbl$gs_name)

get_msigdb_pathways <- function(species = "Homo sapiens") {
  list(
    HALLMARK = msigdbr(species = species, collection = "H"),
    REACTOME = msigdbr(species = species, collection = "C2", subcollection = "CP:REACTOME"),
    GO_BP    = msigdbr(species = species, collection = "C5", subcollection = "GO:BP"),
    GO_CC    = msigdbr(species = species, collection = "C5", subcollection = "GO:CC"),
    GO_MF    = msigdbr(species = species, collection = "C5", subcollection = "GO:MF")
  )
}

make_rank_stat <- function(log2fc, lfcse, pval) {
  log2fc <- suppressWarnings(as.numeric(log2fc))
  lfcse  <- suppressWarnings(as.numeric(lfcse))
  pval   <- suppressWarnings(as.numeric(pval))

  z <- log2fc / lfcse
  bad <- is.na(z) | is.infinite(z)

  if (any(bad)) {
    p <- pval
    p[is.na(p)] <- 1
    p <- pmax(p, .Machine$double.xmin)
    z_fallback <- sign(log2fc) * (-log10(p))
    z[bad] <- z_fallback[bad]
  }

  z[is.na(z) | is.infinite(z)] <- 0
  z
}

row_zscore <- function(m) {
  mu  <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  (m - mu) / sdv
}

plot_dot <- function(res, out_png, title, top_n = TOP_DOT_N, alpha = ALPHA_GSEA) {
  res <- res[is.finite(res$NES) & !is.na(res$padj), , drop = FALSE]
  if (nrow(res) == 0) return(invisible(NULL))

  res <- res[order(res$padj, res$pval), , drop = FALSE]
  res$padj <- pmax(res$padj, .Machine$double.xmin)
  res$is_sig <- res$padj < alpha
  res$fdr_log10 <- -log10(res$padj)

  top <- head(res, top_n)
  if (nrow(top) == 0) return(invisible(NULL))

  top$pathway_pretty <- gsub("_", " ", top$pathway)
  top$pathway_pretty <- factor(top$pathway_pretty, levels = rev(top$pathway_pretty))

  p <- ggplot(top, aes(x = NES, y = pathway_pretty, size = fdr_log10, colour = is_sig)) +
    geom_point() +
    scale_size_continuous(range = c(4, 12), name = "FDR (-log10)") +
    scale_colour_manual(
      values = c(`TRUE` = "black", `FALSE` = "grey60"),
      labels = c(`TRUE` = paste0("FDR < ", alpha), `FALSE` = "Not significant"),
      name = NULL
    ) +
    labs(title = title, x = "NES", y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 9),
      plot.margin = margin(12, 18, 12, 12)
    )

  ggsave(
    out_png, p,
    width = 14,
    height = max(6, 0.28 * nrow(top) + 3),
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )
}

collapse_duplicates_maxabs <- function(stats_named) {
  df <- data.frame(
    gene = names(stats_named),
    stat = as.numeric(stats_named),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$gene) & df$gene != "" & !is.na(df$stat), , drop = FALSE]
  df <- df[order(abs(df$stat), decreasing = TRUE), , drop = FALSE]
  df <- df[!duplicated(df$gene), , drop = FALSE]
  out <- df$stat
  names(out) <- df$gene
  sort(out, decreasing = TRUE)
}

## -----------------------------
## Load TPM (needed for heatmapZ split + ENSG->SYMBOL mapping)
## -----------------------------
tpm <- read_auto(TPM_TSV)
if (!all(c("gene_id", "gene_name") %in% colnames(tpm))) {
  stop("TPM_TSV must have gene_id and gene_name columns + samples.")
}

tpm$gene_id   <- strip_ens_version(tpm$gene_id)
tpm$gene_name <- as.character(tpm$gene_name)

# Mapping ENSG -> SYMBOL (fallback to ENSG if missing)
id2sym <- setNames(tpm$gene_name, tpm$gene_id)

map_ids_to_symbols <- function(ids) {
  ids <- strip_ens_version(ids)
  sym <- unname(id2sym[ids])
  sym <- as.character(sym)
  sym[is.na(sym) | !nzchar(sym)] <- ids[is.na(sym) | !nzchar(sym)]
  sym
}

## -----------------------------
## TPM sample columns (exclude P1 for means / z-scores)
## -----------------------------
sample_cols <- setdiff(colnames(tpm), c("gene_id", "gene_name"))
if (length(EXCLUDE_SAMPLES) > 0) {
  present <- intersect(sample_cols, EXCLUDE_SAMPLES)
  if (length(present) > 0) {
    message("Excluding TPM samples from means/z: ", paste(present, collapse = ", "))
    sample_cols <- setdiff(sample_cols, present)
  } else {
    message("Requested TPM exclusions not found in TPM columns: ", paste(EXCLUDE_SAMPLES, collapse = ", "))
  }
}

groups <- vapply(sample_cols, infer_group_from_col, character(1))
keep <- !is.na(groups)
sample_cols <- sample_cols[keep]
groups <- groups[keep]

if (length(sample_cols) == 0) stop("No TPM sample columns detected by prefixes A6_/H*/P* (after exclusions).")

group_levels <- c("H", "A6", "P")
for (g in group_levels) {
  if (!any(groups == g)) stop("No TPM samples detected for group: ", g)
}

compute_mean_logm_z <- function(gene_vec, label_for_msgs = "") {
  gene_vec <- strip_ens_version(gene_vec)
  gene_vec <- gene_vec[!is.na(gene_vec) & nzchar(gene_vec)]
  gene_vec <- gene_vec[!duplicated(gene_vec)]

  tpm_sub <- tpm[tpm$gene_id %in% gene_vec, c("gene_id", sample_cols), drop = FALSE]
  missing <- setdiff(gene_vec, tpm_sub$gene_id)
  if (length(missing) > 0) {
    warning(label_for_msgs, ": ", length(missing), " genes not found in TPM gene_id; dropping them from this run.")
  }

  genes_keep <- gene_vec[gene_vec %in% tpm_sub$gene_id]
  tpm_sub <- tpm_sub[match(genes_keep, tpm_sub$gene_id), , drop = FALSE]

  mat <- as.matrix(tpm_sub[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- tpm_sub$gene_id

  mean_by_group <- sapply(group_levels, function(g) {
    cols <- sample_cols[groups == g]
    rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  })

  logm <- log2(mean_by_group + 1)
  zmat <- row_zscore(logm)

  list(mean_by_group = mean_by_group, logm = logm, zmat = zmat)
}

split_by_heatmapZ <- function(zmat) {
  zH  <- zmat[, "H"]
  zA6 <- zmat[, "A6"]

  dirH  <- sign(zH)
  dirA6 <- sign(zA6)

  dirH[abs(zH) < Z_EPS]   <- 0
  dirA6[abs(zA6) < Z_EPS] <- 0

  if (CONCORDANT_REQUIRES_NONZERO) {
    concordant <- (dirH != 0) & (dirA6 != 0) & (dirH == dirA6)
  } else {
    concordant <- (dirH == dirA6)
  }

  list(
    concordant = concordant,
    zH = zH,
    zA6 = zA6,
    zP = zmat[, "P"],
    dirH = dirH,
    dirA6 = dirA6
  )
}

## -----------------------------
## fgsea runners
## -----------------------------
msig_sets <- get_msigdb_pathways("Homo sapiens")

run_one_small_universe <- function(base_dir, tag, stats_named, title_suffix = "") {
  out_tag <- file.path(base_dir, tag)
  dir_create_or_stop(out_tag)

  if (length(stats_named) < 5) {
    warning(basename(base_dir), " | ", tag, ": only ", length(stats_named),
            " ranked genes; fgsea may return few/no pathways.")
  }

  for (set_name in names(msig_sets)) {
    msig_tbl <- msig_sets[[set_name]]
    pathways <- to_pathway_list(msig_tbl)

    out_set <- file.path(out_tag, set_name)
    dir_create_or_stop(out_set)

    res <- fgseaMultilevel(
      pathways = pathways,
      stats    = stats_named,
      minSize  = MIN_SIZE,
      maxSize  = MAX_SIZE,
      eps      = 0
    )
    res <- as.data.frame(res)
    res <- res[order(res$padj, res$pval), , drop = FALSE]

    # flatten list-cols for CSV
    for (cn in names(res)) {
      if (is.list(res[[cn]])) {
        res[[cn]] <- vapply(res[[cn]], function(x) paste(x, collapse = ";"), character(1))
      }
    }

    write.csv(res, file.path(out_set, paste0(tag, "_", set_name, "_fgsea.csv")), row.names = FALSE)

    plot_dot(
      res,
      file.path(out_set, paste0(tag, "_", set_name, "_dotplot_top", TOP_DOT_N, ".png")),
      title = paste0(tag, " | ", set_name, " | fgsea (heatmapZ split; Z_EPS=", Z_EPS, ")", title_suffix),
      top_n = TOP_DOT_N,
      alpha = ALPHA_GSEA
    )
  }
}

run_fgsea_sets_prerank <- function(stats_named, base_dir, tag, min_size, max_size, top_dot_n, alpha, title_suffix = "") {
  out_tag <- file.path(base_dir, tag)
  dir_create_or_stop(out_tag)

  write.csv(
    data.frame(gene = names(stats_named), rank = as.numeric(stats_named), stringsAsFactors = FALSE),
    file.path(out_tag, paste0(tag, "_ranked_genes.csv")),
    row.names = FALSE
  )

  res_all <- list()

  for (set_name in names(msig_sets)) {
    msig_tbl <- msig_sets[[set_name]]
    pathways <- to_pathway_list(msig_tbl)

    out_set <- file.path(out_tag, set_name)
    dir_create_or_stop(out_set)

    res <- fgseaMultilevel(
      pathways = pathways,
      stats    = stats_named,
      minSize  = min_size,
      maxSize  = max_size,
      eps      = 0
    )
    res <- as.data.frame(res)
    res <- res[order(res$padj, res$pval), , drop = FALSE]

    for (cn in names(res)) {
      if (is.list(res[[cn]])) {
        res[[cn]] <- vapply(res[[cn]], function(x) paste(x, collapse = ";"), character(1))
      }
    }

    res$set <- set_name
    write.csv(res, file.path(out_set, paste0(tag, "_", set_name, "_fgsea.csv")), row.names = FALSE)

    plot_dot(
      res,
      file.path(out_set, paste0(tag, "_", set_name, "_dotplot_top", top_dot_n, ".png")),
      title = paste0(tag, " | ", set_name, " | fgsea (preranked all genes)", title_suffix),
      top_n = top_dot_n,
      alpha = alpha
    )

    res_all[[set_name]] <- res
  }

  do.call(rbind, res_all)
}

## -----------------------------
## Load DE sig + chrX list
## -----------------------------
de_sig <- read_auto(A6P_DE_SIG)
need <- c("gene_id", "log2FoldChange", "padj")
miss <- setdiff(need, colnames(de_sig))
if (length(miss) > 0) stop("Missing columns in A6P_DE_SIG: ", paste(miss, collapse = ", "))

de_sig$gene <- strip_ens_version(de_sig$gene_id)   # MUST match TPM gene_id
de_sig$padj <- padj_fix(de_sig$padj)
de_sig$lfc  <- suppressWarnings(as.numeric(de_sig$log2FoldChange))

# optional symbol column if present
if ("label" %in% colnames(de_sig)) {
  de_sig$symbol <- as.character(de_sig$label)
} else {
  de_sig$symbol <- map_ids_to_symbols(de_sig$gene)
}
de_sig$symbol[is.na(de_sig$symbol) | !nzchar(de_sig$symbol)] <- de_sig$gene[is.na(de_sig$symbol) | !nzchar(de_sig$symbol)]

chrx <- read_auto(A6P_DE_SIG_CHRX)
if (!("gene_id" %in% colnames(chrx))) stop("chrX file missing gene_id column: ", A6P_DE_SIG_CHRX)
chrx_genes <- strip_ens_version(chrx$gene_id)
chrx_genes <- unique(chrx_genes[!is.na(chrx_genes) & nzchar(chrx_genes)])

# chrX excluded sig list
de_noX <- de_sig[!(de_sig$gene %in% chrx_genes), , drop = FALSE]
if (nrow(de_noX) == 0) stop("After excluding chrX, 0 genes remain in A6 vs P sig list.")
de_noX <- de_noX[order(-abs(de_noX$lfc), de_noX$padj), , drop = FALSE]

# chrX included sig list
de_inclX <- de_sig[order(-abs(de_sig$lfc), de_sig$padj), , drop = FALSE]

## -----------------------------
## Make ranking stats for SMALL-universe fgsea
## IMPORTANT: stats must be named by SYMBOLS (msigdbr gene_symbol),
## even though gene sets for z-split come from ENSG IDs.
## -----------------------------
make_stats_from_sig_tbl <- function(genes, sig_tbl) {
  genes <- strip_ens_version(genes)
  sig_tbl$gene <- strip_ens_version(sig_tbl$gene)

  idx <- match(genes, sig_tbl$gene)
  de_sub <- sig_tbl[idx, , drop = FALSE]

  bad <- is.na(de_sub$gene) | !nzchar(de_sub$gene)
  if (any(bad)) {
    warning("Dropping ", sum(bad), " genes that could not be matched back to DE table (ID mismatch).")
    de_sub <- de_sub[!bad, , drop = FALSE]
  }

  p <- padj_fix(de_sub$padj)
  stat <- sign(de_sub$lfc) * (-log10(p))

  # name by SYMBOL so fgsea matches msigdbr pathways
  sym <- de_sub$symbol
  sym[is.na(sym) | !nzchar(sym)] <- map_ids_to_symbols(de_sub$gene)[is.na(sym) | !nzchar(sym)]

  names(stat) <- sym
  stat <- stat[!is.na(stat)]
  stat <- stat[stat != 0]

  collapse_duplicates_maxabs(stat)
}

## -----------------------------
## Setup output directories
## -----------------------------
dir_create_or_stop(OUT_GSEA_EXCL_TOP100)
dir_create_or_stop(OUT_GSEA_EXCL_ALLSIG)
dir_create_or_stop(OUT_GSEA_EXCL_PRERANK)

dir_create_or_stop(OUT_GSEA_INCL_TOP100)
dir_create_or_stop(OUT_GSEA_INCL_ALLSIG)
dir_create_or_stop(OUT_GSEA_INCL_PRERANK)

## =====================================================================
## A) TOP100 heatmapZ split — chrX EXCLUDED
## =====================================================================
top_genes_excl <- head(unique(de_noX$gene), TOP_N)

mz <- compute_mean_logm_z(top_genes_excl, label_for_msgs = "TOP100_excl_chrX")
zmat <- mz$zmat
mean_by_group <- mz$mean_by_group
logm <- mz$logm

sp <- split_by_heatmapZ(zmat)
concordant <- sp$concordant

genes_conc <- rownames(zmat)[which(concordant)]
genes_disc <- rownames(zmat)[which(!concordant)]

message("TOP100 (excl chrX) total: ", nrow(zmat))
message("  Concordant:  ", length(genes_conc))
message("  Discordant: ", length(genes_disc))

split_df <- data.frame(
  gene = rownames(zmat),
  meanTPM_H  = mean_by_group[, "H"],
  meanTPM_A6 = mean_by_group[, "A6"],
  meanTPM_P  = mean_by_group[, "P"],
  log2mean_H  = logm[, "H"],
  log2mean_A6 = logm[, "A6"],
  log2mean_P  = logm[, "P"],
  z_H  = sp$zH,
  z_A6 = sp$zA6,
  z_P  = sp$zP,
  dir_H  = sp$dirH,
  dir_A6 = sp$dirA6,
  group = ifelse(concordant, "CONCORDANT", "DISCORDANT"),
  stringsAsFactors = FALSE
)

write.csv(split_df, file.path(OUT_GSEA_EXCL_TOP100, "top100_split_debug_heatmapZ.csv"), row.names = FALSE)
write.csv(data.frame(gene = genes_conc), file.path(OUT_GSEA_EXCL_TOP100, "genes_concordant.csv"), row.names = FALSE)
write.csv(data.frame(gene = genes_disc), file.path(OUT_GSEA_EXCL_TOP100, "genes_discordant.csv"), row.names = FALSE)

stats_all_top100 <- make_stats_from_sig_tbl(rownames(zmat), de_noX)
stats_conc <- stats_all_top100[map_ids_to_symbols(genes_conc)]
stats_disc <- stats_all_top100[map_ids_to_symbols(genes_disc)]
stats_conc <- stats_conc[!is.na(stats_conc)]
stats_disc <- stats_disc[!is.na(stats_disc)]

write.csv(data.frame(gene = names(stats_conc), rank = as.numeric(stats_conc)),
          file.path(OUT_GSEA_EXCL_TOP100, "ranked_genes_concordant.csv"), row.names = FALSE)
write.csv(data.frame(gene = names(stats_disc), rank = as.numeric(stats_disc)),
          file.path(OUT_GSEA_EXCL_TOP100, "ranked_genes_discordant.csv"), row.names = FALSE)

run_one_small_universe(OUT_GSEA_EXCL_TOP100, "CONCORDANT", stats_conc, title_suffix = " | chrX excluded")
run_one_small_universe(OUT_GSEA_EXCL_TOP100, "DISCORDANT", stats_disc, title_suffix = " | chrX excluded")

## =====================================================================
## B) ALL-SIGNIFICANT heatmapZ split — chrX EXCLUDED
## =====================================================================
all_sig_excl <- unique(de_noX$gene)

mz_all <- compute_mean_logm_z(all_sig_excl, label_for_msgs = "ALLSIG_excl_chrX")
zmat_all <- mz_all$zmat
mean_by_group_all <- mz_all$mean_by_group
logm_all <- mz_all$logm

sp_all <- split_by_heatmapZ(zmat_all)
conc_all <- sp_all$concordant

genes_conc_all <- rownames(zmat_all)[which(conc_all)]
genes_disc_all <- rownames(zmat_all)[which(!conc_all)]

message("ALLSIG (excl chrX) total (after TPM dropouts): ", nrow(zmat_all))
message("  Concordant:  ", length(genes_conc_all))
message("  Discordant: ", length(genes_disc_all))

split_df_all <- data.frame(
  gene = rownames(zmat_all),
  meanTPM_H  = mean_by_group_all[, "H"],
  meanTPM_A6 = mean_by_group_all[, "A6"],
  meanTPM_P  = mean_by_group_all[, "P"],
  log2mean_H  = logm_all[, "H"],
  log2mean_A6 = logm_all[, "A6"],
  log2mean_P  = logm_all[, "P"],
  z_H  = sp_all$zH,
  z_A6 = sp_all$zA6,
  z_P  = sp_all$zP,
  dir_H  = sp_all$dirH,
  dir_A6 = sp_all$dirA6,
  group = ifelse(conc_all, "CONCORDANT", "DISCORDANT"),
  stringsAsFactors = FALSE
)

write.csv(split_df_all, file.path(OUT_GSEA_EXCL_ALLSIG, "allSig_split_debug_heatmapZ.csv"), row.names = FALSE)
write.csv(data.frame(gene = genes_conc_all), file.path(OUT_GSEA_EXCL_ALLSIG, "genes_concordant.csv"), row.names = FALSE)
write.csv(data.frame(gene = genes_disc_all), file.path(OUT_GSEA_EXCL_ALLSIG, "genes_discordant.csv"), row.names = FALSE)

stats_allsig <- make_stats_from_sig_tbl(rownames(zmat_all), de_noX)
stats_conc_all <- stats_allsig[map_ids_to_symbols(genes_conc_all)]
stats_disc_all <- stats_allsig[map_ids_to_symbols(genes_disc_all)]
stats_conc_all <- stats_conc_all[!is.na(stats_conc_all)]
stats_disc_all <- stats_disc_all[!is.na(stats_disc_all)]

write.csv(data.frame(gene = names(stats_conc_all), rank = as.numeric(stats_conc_all)),
          file.path(OUT_GSEA_EXCL_ALLSIG, "ranked_genes_concordant.csv"), row.names = FALSE)
write.csv(data.frame(gene = names(stats_disc_all), rank = as.numeric(stats_disc_all)),
          file.path(OUT_GSEA_EXCL_ALLSIG, "ranked_genes_discordant.csv"), row.names = FALSE)

run_one_small_universe(OUT_GSEA_EXCL_ALLSIG, "CONCORDANT", stats_conc_all, title_suffix = " | allSig | chrX excluded")
run_one_small_universe(OUT_GSEA_EXCL_ALLSIG, "DISCORDANT", stats_disc_all, title_suffix = " | allSig | chrX excluded")

## =====================================================================
## A2) TOP100 heatmapZ split — chrX INCLUDED
## =====================================================================
top_genes_incl <- head(unique(de_inclX$gene), TOP_N)

mz_i <- compute_mean_logm_z(top_genes_incl, label_for_msgs = "TOP100_incl_chrX")
zmat_i <- mz_i$zmat
mean_by_group_i <- mz_i$mean_by_group
logm_i <- mz_i$logm

sp_i <- split_by_heatmapZ(zmat_i)
conc_i <- sp_i$concordant

genes_conc_i <- rownames(zmat_i)[which(conc_i)]
genes_disc_i <- rownames(zmat_i)[which(!conc_i)]

message("TOP100 (incl chrX) total: ", nrow(zmat_i))
message("  Concordant:  ", length(genes_conc_i))
message("  Discordant: ", length(genes_disc_i))

split_df_i <- data.frame(
  gene = rownames(zmat_i),
  meanTPM_H  = mean_by_group_i[, "H"],
  meanTPM_A6 = mean_by_group_i[, "A6"],
  meanTPM_P  = mean_by_group_i[, "P"],
  log2mean_H  = logm_i[, "H"],
  log2mean_A6 = logm_i[, "A6"],
  log2mean_P  = logm_i[, "P"],
  z_H  = sp_i$zH,
  z_A6 = sp_i$zA6,
  z_P  = sp_i$zP,
  dir_H  = sp_i$dirH,
  dir_A6 = sp_i$dirA6,
  group = ifelse(conc_i, "CONCORDANT", "DISCORDANT"),
  stringsAsFactors = FALSE
)

write.csv(split_df_i, file.path(OUT_GSEA_INCL_TOP100, "top100_split_debug_heatmapZ.csv"), row.names = FALSE)
write.csv(data.frame(gene = genes_conc_i), file.path(OUT_GSEA_INCL_TOP100, "genes_concordant.csv"), row.names = FALSE)
write.csv(data.frame(gene = genes_disc_i), file.path(OUT_GSEA_INCL_TOP100, "genes_discordant.csv"), row.names = FALSE)

stats_top100_i <- make_stats_from_sig_tbl(rownames(zmat_i), de_inclX)
stats_conc_i <- stats_top100_i[map_ids_to_symbols(genes_conc_i)]
stats_disc_i <- stats_top100_i[map_ids_to_symbols(genes_disc_i)]
stats_conc_i <- stats_conc_i[!is.na(stats_conc_i)]
stats_disc_i <- stats_disc_i[!is.na(stats_disc_i)]

write.csv(data.frame(gene = names(stats_conc_i), rank = as.numeric(stats_conc_i)),
          file.path(OUT_GSEA_INCL_TOP100, "ranked_genes_concordant.csv"), row.names = FALSE)
write.csv(data.frame(gene = names(stats_disc_i), rank = as.numeric(stats_disc_i)),
          file.path(OUT_GSEA_INCL_TOP100, "ranked_genes_discordant.csv"), row.names = FALSE)

run_one_small_universe(OUT_GSEA_INCL_TOP100, "CONCORDANT", stats_conc_i, title_suffix = " | chrX included")
run_one_small_universe(OUT_GSEA_INCL_TOP100, "DISCORDANT", stats_disc_i, title_suffix = " | chrX included")

## =====================================================================
## B2) ALL-SIGNIFICANT heatmapZ split — chrX INCLUDED
## =====================================================================
all_sig_incl <- unique(de_inclX$gene)

mz_all_i <- compute_mean_logm_z(all_sig_incl, label_for_msgs = "ALLSIG_incl_chrX")
zmat_all_i <- mz_all_i$zmat
mean_by_group_all_i <- mz_all_i$mean_by_group
logm_all_i <- mz_all_i$logm

sp_all_i <- split_by_heatmapZ(zmat_all_i)
conc_all_i <- sp_all_i$concordant

genes_conc_all_i <- rownames(zmat_all_i)[which(conc_all_i)]
genes_disc_all_i <- rownames(zmat_all_i)[which(!conc_all_i)]

message("ALLSIG (incl chrX) total (after TPM dropouts): ", nrow(zmat_all_i))
message("  Concordant:  ", length(genes_conc_all_i))
message("  Discordant: ", length(genes_disc_all_i))

split_df_all_i <- data.frame(
  gene = rownames(zmat_all_i),
  meanTPM_H  = mean_by_group_all_i[, "H"],
  meanTPM_A6 = mean_by_group_all_i[, "A6"],
  meanTPM_P  = mean_by_group_all_i[, "P"],
  log2mean_H  = logm_all_i[, "H"],
  log2mean_A6 = logm_all_i[, "A6"],
  log2mean_P  = logm_all_i[, "P"],
  z_H  = sp_all_i$zH,
  z_A6 = sp_all_i$zA6,
  z_P  = sp_all_i$zP,
  dir_H  = sp_all_i$dirH,
  dir_A6 = sp_all_i$dirA6,
  group = ifelse(conc_all_i, "CONCORDANT", "DISCORDANT"),
  stringsAsFactors = FALSE
)

write.csv(split_df_all_i, file.path(OUT_GSEA_INCL_ALLSIG, "allSig_split_debug_heatmapZ.csv"), row.names = FALSE)
write.csv(data.frame(gene = genes_conc_all_i), file.path(OUT_GSEA_INCL_ALLSIG, "genes_concordant.csv"), row.names = FALSE)
write.csv(data.frame(gene = genes_disc_all_i), file.path(OUT_GSEA_INCL_ALLSIG, "genes_discordant.csv"), row.names = FALSE)

stats_allsig_i <- make_stats_from_sig_tbl(rownames(zmat_all_i), de_inclX)
stats_conc_all_i <- stats_allsig_i[map_ids_to_symbols(genes_conc_all_i)]
stats_disc_all_i <- stats_allsig_i[map_ids_to_symbols(genes_disc_all_i)]
stats_conc_all_i <- stats_conc_all_i[!is.na(stats_conc_all_i)]
stats_disc_all_i <- stats_disc_all_i[!is.na(stats_disc_all_i)]

write.csv(data.frame(gene = names(stats_conc_all_i), rank = as.numeric(stats_conc_all_i)),
          file.path(OUT_GSEA_INCL_ALLSIG, "ranked_genes_concordant.csv"), row.names = FALSE)
write.csv(data.frame(gene = names(stats_disc_all_i), rank = as.numeric(stats_disc_all_i)),
          file.path(OUT_GSEA_INCL_ALLSIG, "ranked_genes_discordant.csv"), row.names = FALSE)

run_one_small_universe(OUT_GSEA_INCL_ALLSIG, "CONCORDANT", stats_conc_all_i, title_suffix = " | allSig | chrX included")
run_one_small_universe(OUT_GSEA_INCL_ALLSIG, "DISCORDANT", stats_disc_all_i, title_suffix = " | allSig | chrX included")

## =====================================================================
## C) PRERANKED GSEA on ALL GENES
##   C1: chrX EXCLUDED + concordance + visuals
##   C2: chrX INCLUDED + concordance + visuals
## =====================================================================
get_chrX_gene_set <- function(chrX_df) {
  g <- character(0)
  if ("gene_id" %in% names(chrX_df)) g <- c(g, strip_ens_version(chrX_df$gene_id))
  if ("label"   %in% names(chrX_df)) g <- c(g, as.character(chrX_df$label))
  g <- unique(g[!is.na(g) & nzchar(g)])
  g
}

pick_gene_symbol <- function(df) {
  if ("label" %in% names(df)) {
    lab <- as.character(df$label)
    if (any(nzchar(lab))) return(lab)
  }
  # fallback: map ENSG -> symbol
  map_ids_to_symbols(df$gene_id)
}

make_preranked_vector_excl_chrX <- function(de_path, chrx_path, tag) {
  de <- read_auto(de_path)
  chrx <- read_auto(chrx_path)

  need <- c("log2FoldChange", "lfcSE", "pvalue")
  miss <- setdiff(need, names(de))
  if (length(miss) > 0) stop(tag, ": missing columns in DE file: ", paste(miss, collapse = ", "))

  gene_sym <- pick_gene_symbol(de)
  chrX_genes <- get_chrX_gene_set(chrx)

  keep <- !is.na(gene_sym) & nzchar(gene_sym) & !(gene_sym %in% chrX_genes)
  de <- de[keep, , drop = FALSE]
  gene_sym <- gene_sym[keep]

  stat <- make_rank_stat(de$log2FoldChange, de$lfcSE, de$pvalue)
  names(stat) <- gene_sym
  stat <- stat[!is.na(stat)]
  stat <- stat[stat != 0]
  stat <- collapse_duplicates_maxabs(stat)

  if (length(stat) < 1000) warning(tag, ": ranked list has only ", length(stat), " genes after filtering; check mapping.")
  stat
}

make_preranked_vector_incl_chrX <- function(de_path, tag) {
  de <- read_auto(de_path)

  need <- c("log2FoldChange", "lfcSE", "pvalue")
  miss <- setdiff(need, names(de))
  if (length(miss) > 0) stop(tag, ": missing columns in DE file: ", paste(miss, collapse = ", "))

  gene_sym <- pick_gene_symbol(de)
  keep <- !is.na(gene_sym) & nzchar(gene_sym)
  de <- de[keep, , drop = FALSE]
  gene_sym <- gene_sym[keep]

  stat <- make_rank_stat(de$log2FoldChange, de$lfcSE, de$pvalue)
  names(stat) <- gene_sym
  stat <- stat[!is.na(stat)]
  stat <- stat[stat != 0]
  stat <- collapse_duplicates_maxabs(stat)

  if (length(stat) < 1000) warning(tag, ": ranked list has only ", length(stat), " genes; check mapping.")
  stat
}

make_concordance_tables_and_visuals <- function(res_A6P, res_HVP, out_dir, title_prefix) {
  keep_cols <- c("set", "pathway", "NES", "pval", "padj", "size")

  res_A6P2 <- res_A6P[, intersect(keep_cols, names(res_A6P)), drop = FALSE]
  res_HVP2 <- res_HVP[, intersect(keep_cols, names(res_HVP)), drop = FALSE]

  colnames(res_A6P2) <- sub("^NES$",  "NES_A6vsP",  colnames(res_A6P2))
  colnames(res_A6P2) <- sub("^pval$", "pval_A6vsP", colnames(res_A6P2))
  colnames(res_A6P2) <- sub("^padj$", "padj_A6vsP", colnames(res_A6P2))
  colnames(res_A6P2) <- sub("^size$", "size_A6vsP", colnames(res_A6P2))

  colnames(res_HVP2) <- sub("^NES$",  "NES_HvsP",  colnames(res_HVP2))
  colnames(res_HVP2) <- sub("^pval$", "pval_HvsP", colnames(res_HVP2))
  colnames(res_HVP2) <- sub("^padj$", "padj_HvsP", colnames(res_HVP2))
  colnames(res_HVP2) <- sub("^size$", "size_HvsP", colnames(res_HVP2))

  merged <- merge(res_A6P2, res_HVP2, by = c("set", "pathway"), all = TRUE)

  merged$sign_A6 <- sign(merged$NES_A6vsP)
  merged$sign_H  <- sign(merged$NES_HvsP)

  merged$concordance <- ifelse(
    is.na(merged$NES_A6vsP) | is.na(merged$NES_HvsP),
    "MISSING",
    ifelse(merged$sign_A6 == merged$sign_H, "CONCORDANT", "DISCORDANT")
  )

  merged$sig_A6 <- !is.na(merged$padj_A6vsP) & (merged$padj_A6vsP < ALPHA_GSEA)
  merged$sig_H  <- !is.na(merged$padj_HvsP)  & (merged$padj_HvsP  < ALPHA_GSEA)

  merged$sig_class <- ifelse(
    merged$sig_A6 & merged$sig_H, "SIG_BOTH",
    ifelse(merged$sig_A6 | merged$sig_H, "SIG_ONE", "NOT_SIG")
  )

  merged <- merged[order(merged$set, merged$sig_class, pmin(merged$padj_A6vsP, merged$padj_HvsP, na.rm = TRUE)), ]

  write.csv(merged, file.path(out_dir, "pathway_concordance_all.csv"), row.names = FALSE)
  write.csv(merged[merged$sig_class != "NOT_SIG", , drop = FALSE],
            file.path(out_dir, "pathway_concordance_sigOnly.csv"), row.names = FALSE)

  vis_dir <- file.path(out_dir, "COMMON_CONCORDANT_visuals")
  dir_create_or_stop(vis_dir)

  cc <- merged[
    merged$sig_class == "SIG_BOTH" &
      merged$concordance == "CONCORDANT" &
      is.finite(merged$NES_A6vsP) & is.finite(merged$NES_HvsP),
    , drop = FALSE
  ]
  write.csv(cc, file.path(vis_dir, "common_concordant_SIG_BOTH_ALLSETS.csv"), row.names = FALSE)

  if (nrow(cc) == 0) {
    warning(title_prefix, ": No pathways found with SIG_BOTH + CONCORDANT. Nothing to plot.")
    return(invisible(merged))
  }

  TOP_COMMON <- 30
  set_order <- c("GO_BP", "GO_CC", "GO_MF", "HALLMARK", "REACTOME")
  sets_present <- intersect(set_order, unique(cc$set))
  sets_present <- c(sets_present, setdiff(unique(cc$set), sets_present))

  for (set_name in sets_present) {
    cc_set <- cc[cc$set == set_name, , drop = FALSE]
    if (nrow(cc_set) == 0) next

    set_dir <- file.path(vis_dir, set_name)
    dir_create_or_stop(set_dir)

    cc_set$min_padj <- pmin(cc_set$padj_A6vsP, cc_set$padj_HvsP, na.rm = TRUE)
    cc_set <- cc_set[order(cc_set$min_padj, cc_set$pathway), , drop = FALSE]

    write.csv(cc_set, file.path(set_dir, paste0("common_concordant_SIG_BOTH_", set_name, "_ALL.csv")), row.names = FALSE)

    n_take <- min(TOP_COMMON, nrow(cc_set))
    top_cc <- head(cc_set, n_take)
    write.csv(top_cc, file.path(set_dir, paste0("common_concordant_SIG_BOTH_", set_name, "_TOP", n_take, "_minFDR.csv")),
              row.names = FALSE)

    top_cc$pathway_pretty <- gsub("_", " ", top_cc$pathway)
    top_cc$pathway_pretty <- factor(top_cc$pathway_pretty, levels = rev(top_cc$pathway_pretty))

    p_sc <- ggplot(top_cc, aes(x = NES_A6vsP, y = NES_HvsP)) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_point(alpha = 0.8) +
      ggrepel::geom_text_repel(aes(label = pathway_pretty), size = 3, max.overlaps = Inf) +
      theme_bw(base_size = 12) +
      labs(
        title = paste0(title_prefix, " | ", set_name, " | Common concordant (SIG_BOTH) | Top ", n_take, " by min FDR"),
        x = "NES (A6 vs P)",
        y = "NES (H vs P)"
      )
    ggsave(file.path(set_dir, paste0("NES_scatter_common_concordant_", set_name, "_top", n_take, ".png")),
           p_sc, width = 10, height = 8, dpi = 300)

    p_db <- ggplot(top_cc, aes(y = pathway_pretty)) +
      geom_vline(xintercept = 0) +
      geom_segment(aes(x = NES_A6vsP, xend = NES_HvsP, yend = pathway_pretty), alpha = 0.6) +
      geom_point(aes(x = NES_A6vsP, shape = "A6 vs P"), size = 2) +
      geom_point(aes(x = NES_HvsP,  shape = "H vs P"),  size = 2) +
      theme_bw(base_size = 12) +
      labs(
        title = paste0(title_prefix, " | ", set_name, " | Common concordant (SIG_BOTH) | Top ", n_take, " by min FDR"),
        x = "NES",
        y = NULL,
        shape = NULL
      ) +
      theme(panel.grid = element_blank())
    ggsave(file.path(set_dir, paste0("NES_dumbbell_common_concordant_", set_name, "_top", n_take, ".png")),
           p_db, width = 12, height = max(6, 0.28 * nrow(top_cc) + 3), dpi = 300)

    pp <- as.character(top_cc$pathway_pretty)
    long_cc <- rbind(
      data.frame(pathway_pretty = pp, contrast = "A6 vs P", NES = top_cc$NES_A6vsP),
      data.frame(pathway_pretty = pp, contrast = "H vs P",  NES = top_cc$NES_HvsP)
    )
    long_cc$pathway_pretty <- factor(long_cc$pathway_pretty, levels = levels(top_cc$pathway_pretty))
    long_cc$contrast <- factor(long_cc$contrast, levels = c("A6 vs P", "H vs P"))

    p_hm <- ggplot(long_cc, aes(x = contrast, y = pathway_pretty, fill = NES)) +
      geom_tile() +
      theme_bw(base_size = 12) +
      labs(
        title = paste0(title_prefix, " | ", set_name, " | NES heatmap | Top ", n_take, " by min FDR"),
        x = NULL, y = NULL
      ) +
      theme(panel.grid = element_blank())
    ggsave(file.path(set_dir, paste0("NES_heatmap_common_concordant_", set_name, "_top", n_take, ".png")),
           p_hm, width = 7, height = max(6, 0.28 * nrow(top_cc) + 3), dpi = 300)
  }

  sc <- merged[is.finite(merged$NES_A6vsP) & is.finite(merged$NES_HvsP), , drop = FALSE]
  if (nrow(sc) > 0) {
    p <- ggplot(sc, aes(x = NES_A6vsP, y = NES_HvsP)) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_point(aes(shape = concordance), alpha = 0.6) +
      theme_bw(base_size = 12) +
      labs(
        title = paste0(title_prefix, " | Pathway NES concordance (all pathways)"),
        x = "NES (A6 vs P)",
        y = "NES (H vs P)"
      )
    ggsave(file.path(out_dir, "NES_scatter_A6vsP_vs_HvsP.png"),
           p, width = 8, height = 6, dpi = 300)
  }

  invisible(merged)
}

## ---- C1: preranked chrX excluded ----
message("\n--- Preranked GSEA on ALL genes (chrX excluded) ---")
stats_A6P_excl <- make_preranked_vector_excl_chrX(A6P_DE_ALL, A6P_DE_ALL_CHRX, tag = "A6_vs_P")
stats_HVP_excl <- make_preranked_vector_excl_chrX(HVP_DE_ALL, HVP_DE_ALL_CHRX, tag = "H_vs_P")

res_A6P_excl <- run_fgsea_sets_prerank(
  stats_named = stats_A6P_excl,
  base_dir    = OUT_GSEA_EXCL_PRERANK,
  tag         = "A6_vs_P",
  min_size    = PRERANK_MIN_SIZE,
  max_size    = PRERANK_MAX_SIZE,
  top_dot_n   = PRERANK_TOP_DOT_N,
  alpha       = ALPHA_GSEA,
  title_suffix = " | chrX excluded"
)

res_HVP_excl <- run_fgsea_sets_prerank(
  stats_named = stats_HVP_excl,
  base_dir    = OUT_GSEA_EXCL_PRERANK,
  tag         = "H_vs_P",
  min_size    = PRERANK_MIN_SIZE,
  max_size    = PRERANK_MAX_SIZE,
  top_dot_n   = PRERANK_TOP_DOT_N,
  alpha       = ALPHA_GSEA,
  title_suffix = " | chrX excluded"
)

make_concordance_tables_and_visuals(
  res_A6P = res_A6P_excl,
  res_HVP = res_HVP_excl,
  out_dir = OUT_GSEA_EXCL_PRERANK,
  title_prefix = "Preranked (chrX excluded)"
)

## ---- C2: preranked chrX included ----
message("\n--- Preranked GSEA on ALL genes (chrX included) ---")
stats_A6P_incl <- make_preranked_vector_incl_chrX(A6P_DE_ALL, tag = "A6_vs_P_incl_chrX")
stats_HVP_incl <- make_preranked_vector_incl_chrX(HVP_DE_ALL, tag = "H_vs_P_incl_chrX")

res_A6P_incl <- run_fgsea_sets_prerank(
  stats_named = stats_A6P_incl,
  base_dir    = OUT_GSEA_INCL_PRERANK,
  tag         = "A6_vs_P",
  min_size    = PRERANK_MIN_SIZE,
  max_size    = PRERANK_MAX_SIZE,
  top_dot_n   = PRERANK_TOP_DOT_N,
  alpha       = ALPHA_GSEA,
  title_suffix = " | chrX included"
)

res_HVP_incl <- run_fgsea_sets_prerank(
  stats_named = stats_HVP_incl,
  base_dir    = OUT_GSEA_INCL_PRERANK,
  tag         = "H_vs_P",
  min_size    = PRERANK_MIN_SIZE,
  max_size    = PRERANK_MAX_SIZE,
  top_dot_n   = PRERANK_TOP_DOT_N,
  alpha       = ALPHA_GSEA,
  title_suffix = " | chrX included"
)

make_concordance_tables_and_visuals(
  res_A6P = res_A6P_incl,
  res_HVP = res_HVP_incl,
  out_dir = OUT_GSEA_INCL_PRERANK,
  title_prefix = "Preranked (chrX included)"
)

message("\nDone.")
message("Top100 (excl chrX) outputs in:   ", OUT_GSEA_EXCL_TOP100)
message("AllSig (excl chrX) outputs in:   ", OUT_GSEA_EXCL_ALLSIG)
message("Preranked (excl chrX) outputs in:", OUT_GSEA_EXCL_PRERANK)
message("Top100 (incl chrX) outputs in:   ", OUT_GSEA_INCL_TOP100)
message("AllSig (incl chrX) outputs in:   ", OUT_GSEA_INCL_ALLSIG)
message("Preranked (incl chrX) outputs in:", OUT_GSEA_INCL_PRERANK)
message("TPM exclusions used for mean/z:  ", ifelse(length(EXCLUDE_SAMPLES) > 0, paste(EXCLUDE_SAMPLES, collapse = ", "), "none"))