#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(EnhancedVolcano)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

## -----------------------------
## USER EDITS (D280 + P1 outlier handling)
## -----------------------------
PROJ_OUT_DIR    <- "/mnt/scratch/hqurashi/13_Kelly_RNA_Seq/outputs/D280"
STAR_SALMON_DIR <- file.path(PROJ_OUT_DIR, "star_salmon")
TX2GENE_TSV     <- file.path(STAR_SALMON_DIR, "salmon.merged.tx2gene.tsv")

# Main output base (written under the D280 outputs directory)
OUT_BASE <- file.path(PROJ_OUT_DIR, "DESeq2", "DESeq2_pairs_hg38_tximport_shrinkLFC_gene_names")

# Outlier replicate to exclude for downstream analyses (but keep for "all samples" PCA)
EXCLUDE_SAMPLE_NAMES <- c("P1")

MIN_TOTAL_COUNT <- 10
ALPHA           <- 0.05
LFC_CUTOFF      <- 1

# Volcano axis consistency
VOLCANO_XLIM_MODE <- "max"   # "max" or "p99"
VOLCANO_XPAD      <- 0.25

dir.create(OUT_BASE, recursive = TRUE, showWarnings = FALSE)

## -----------------------------
## Helpers: ID handling + operators
## -----------------------------
strip_ens_version <- function(x) sub("\\..*$", "", as.character(x))
is_ensg_id <- function(x) grepl("^ENSG", strip_ens_version(x))
`%||%` <- function(a, b) if (!is.null(a) && nzchar(a)) a else b

## -----------------------------
## Locate quant.sf files
## -----------------------------
sample_dirs <- list.dirs(STAR_SALMON_DIR, recursive = FALSE, full.names = TRUE)

pick_quant <- function(d) {
  cands <- c(
    file.path(d, "quant.sf"),
    file.path(d, "salmon", "quant.sf"),
    file.path(d, "quant", "quant.sf"),
    file.path(d, "logs", "quant.sf")
  )
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

quant_files  <- vapply(sample_dirs, pick_quant, character(1))
sample_names <- basename(sample_dirs)

ok <- !is.na(quant_files) & nzchar(quant_files)
quant_files  <- quant_files[ok]
sample_names <- sample_names[ok]

infer_group <- function(s) {
  if (grepl("^A6_", s)) return("A6")
  if (grepl("^H",  s))  return("H")
  if (grepl("^P",  s))  return("P")
  NA_character_
}

groups <- vapply(sample_names, infer_group, character(1))
keep <- !is.na(groups)
quant_files  <- quant_files[keep]
sample_names <- sample_names[keep]
groups       <- groups[keep]

if (length(quant_files) == 0) stop("No quant.sf files found matching A6_/H*/P* under: ", STAR_SALMON_DIR)

files <- setNames(quant_files, sample_names)
meta  <- data.frame(
  sample = sample_names,
  group  = factor(groups, levels = c("A6","H","P")),
  row.names = sample_names,
  stringsAsFactors = FALSE
)

message("Found ", nrow(meta), " samples with quant.sf:")
print(table(meta$group))
write.csv(meta, file.path(OUT_BASE, "sample_metadata_detected.csv"), row.names = FALSE)

## -----------------------------
## Read tx2gene
## -----------------------------
if (!file.exists(TX2GENE_TSV)) stop("tx2gene.tsv not found: ", TX2GENE_TSV)
tx2gene <- read.delim(TX2GENE_TSV, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
if (ncol(tx2gene) < 2) stop("tx2gene has <2 columns, unexpected.")

tx2gene <- tx2gene[, 1:2]
colnames(tx2gene) <- c("TXNAME", "GENEID")
tx2gene$TXNAME <- strip_ens_version(tx2gene$TXNAME)
tx2gene$GENEID <- strip_ens_version(tx2gene$GENEID)
tx2gene <- tx2gene[!duplicated(tx2gene$TXNAME), , drop = FALSE]

## -----------------------------
## tximport (gene-level)
## -----------------------------
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE,
  ignoreAfterBar  = TRUE,
  countsFromAbundance = "no"
)

cat("Example txi gene IDs:\n")
print(head(rownames(txi$counts), 20))
cat("How many look like ENSG?: ", sum(grepl("^ENSG", rownames(txi$counts))), "\n")

## -----------------------------
## Gene labels + chromosome mapping
##   - Supports both SYMBOL and ENSG gene IDs
##   - If ENSG: maps to SYMBOL for label display, and maps ENSG -> CHR for chr calls
## -----------------------------
all_ids_raw <- rownames(txi$counts)
all_ids <- strip_ens_version(all_ids_raw)

if (any(is.na(all_ids) | !nzchar(all_ids))) stop("Some gene IDs are blank/NA after cleaning.")

USING_ENSG <- any(grepl("^ENSG", all_ids))
message("Detected gene ID type from tximport: ", if (USING_ENSG) "ENSEMBL (ENSG...)" else "SYMBOL")

sym_map_from_ens <- NULL
chr_map_from_ens <- NULL
chr_map_from_sym <- NULL

if (USING_ENSG) {
  sym_map_from_ens <- tryCatch({
    mapIds(
      org.Hs.eg.db,
      keys      = unique(all_ids),
      keytype   = "ENSEMBL",
      column    = "SYMBOL",
      multiVals = "first"
    )
  }, error = function(e) stop("Failed to map ENSEMBL -> SYMBOL. Error:\n", conditionMessage(e)))

  chr_map_from_ens <- tryCatch({
    mapIds(
      org.Hs.eg.db,
      keys      = unique(all_ids),
      keytype   = "ENSEMBL",
      column    = "CHR",
      multiVals = "first"
    )
  }, error = function(e) stop("Failed to map ENSEMBL -> CHR. Error:\n", conditionMessage(e)))
} else {
  chr_map_from_sym <- tryCatch({
    mapIds(
      org.Hs.eg.db,
      keys      = unique(all_ids),
      keytype   = "SYMBOL",
      column    = "CHR",
      multiVals = "first"
    )
  }, error = function(e) stop("Failed to map SYMBOL -> CHR. Error:\n", conditionMessage(e)))
}

make_labels <- function(gene_ids) {
  ids <- strip_ens_version(gene_ids)

  if (USING_ENSG) {
    sym <- unname(sym_map_from_ens[ids])
    sym <- as.character(sym)
    sym[is.na(sym) | !nzchar(sym)] <- ids[is.na(sym) | !nzchar(sym)]
    return(sym)
  } else {
    return(ids)
  }
}

chr_of_geneid <- function(gene_ids) {
  ids <- strip_ens_version(gene_ids)
  if (USING_ENSG) {
    chr <- unname(chr_map_from_ens[ids])
  } else {
    chr <- unname(chr_map_from_sym[ids])
  }
  chr <- as.character(chr)
  chr <- sub("^chr", "", chr, ignore.case = TRUE)
  chr
}

## -----------------------------
## PCA helper (bubbles + hulls)
## -----------------------------
make_pca_bubbles <- function(vsd_obj, intgroup_col = "group", title = NULL) {
  pca_df <- plotPCA(vsd_obj, intgroup = intgroup_col, returnData = TRUE)
  percentVar <- attr(pca_df, "percentVar")

  if (!("name" %in% colnames(pca_df))) pca_df$name <- rownames(pca_df)

  pca_df[[intgroup_col]] <- factor(pca_df[[intgroup_col]])
  groups2 <- levels(pca_df[[intgroup_col]])

  hull_list <- list()
  for (g in groups2) {
    df_g <- pca_df[pca_df[[intgroup_col]] == g, , drop = FALSE]
    hull_idx <- if (nrow(df_g) <= 2) seq_len(nrow(df_g)) else chull(df_g$PC1, df_g$PC2)
    hull_list[[g]] <- data.frame(PC1=df_g$PC1[hull_idx], PC2=df_g$PC2[hull_idx], group=g)
  }
  hulls <- do.call(rbind, hull_list)
  hulls$group <- factor(hulls$group, levels = levels(pca_df[[intgroup_col]]))

  ggplot() +
    geom_polygon(data=hulls, aes(x=PC1, y=PC2, group=group, fill=group), alpha=0.2, colour=NA) +
    geom_point(data=pca_df, aes(x=PC1, y=PC2, colour=.data[[intgroup_col]]), size=3) +
    geom_text_repel(
      data=pca_df,
      aes(x=PC1, y=PC2, label=name, colour=.data[[intgroup_col]]),
      size=2, box.padding=0.4, point.padding=0.6, max.overlaps=Inf,
      segment.alpha=0.6, show.legend=FALSE
    ) +
    xlab(paste0("PC1 (", round(percentVar[1] * 100, 1), "%)")) +
    ylab(paste0("PC2 (", round(percentVar[2] * 100, 1), "%)")) +
    ggtitle(title %||% "") +
    theme_bw() +
    theme(panel.grid=element_blank(), legend.title=element_blank(),
          plot.title=element_text(hjust=0.5))
}

subset_txi_to_samples <- function(txi_obj, sample_ids) {
  txi_sub <- txi_obj
  txi_sub$counts    <- txi_obj$counts[, sample_ids, drop = FALSE]
  txi_sub$abundance <- txi_obj$abundance[, sample_ids, drop = FALSE]
  txi_sub$length    <- txi_obj$length[, sample_ids, drop = FALSE]
  txi_sub
}

## -----------------------------
## PCA 1) All samples (includes P1)
## -----------------------------
dds_all <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ group)
dds_all <- dds_all[rowSums(counts(dds_all)) > MIN_TOTAL_COUNT, ]
dds_all <- estimateSizeFactors(dds_all, type = "poscounts")
vsd_all <- vst(dds_all, blind = FALSE)

p_all <- make_pca_bubbles(vsd_all, intgroup_col = "group", title = "All samples (A6 vs H vs P) [includes P1]")
ggsave(file.path(OUT_BASE, "PCA_all_samples_including_P1.png"), p_all, width = 6, height = 5, dpi = 600)
ggsave(file.path(OUT_BASE, "PCA_all_samples_including_P1.pdf"), p_all, width = 6, height = 5)

## -----------------------------
## Exclude P1 for downstream analyses (and make PCA 2)
## -----------------------------
present_excl <- intersect(EXCLUDE_SAMPLE_NAMES, meta$sample)
if (length(present_excl) == 0) {
  message("NOTE: None of EXCLUDE_SAMPLE_NAMES were found in detected samples. Downstream will include all samples.")
} else {
  message("Excluding samples for downstream analyses: ", paste(present_excl, collapse = ", "))
}
writeLines(present_excl, con = file.path(OUT_BASE, "excluded_samples_downstream.txt"))

meta_main <- meta[!(meta$sample %in% EXCLUDE_SAMPLE_NAMES), , drop = FALSE]
if (nrow(meta_main) < nrow(meta)) {
  txi_main <- subset_txi_to_samples(txi, rownames(meta_main))
} else {
  txi_main <- txi
}

dds_noP1 <- DESeqDataSetFromTximport(txi_main, colData = meta_main, design = ~ group)
dds_noP1 <- dds_noP1[rowSums(counts(dds_noP1)) > MIN_TOTAL_COUNT, ]
dds_noP1 <- estimateSizeFactors(dds_noP1, type = "poscounts")
vsd_noP1 <- vst(dds_noP1, blind = FALSE)

p_noP1 <- make_pca_bubbles(vsd_noP1, intgroup_col = "group", title = "All samples (A6 vs H vs P) [P1 excluded]")
ggsave(file.path(OUT_BASE, "PCA_all_samples_excluding_P1.png"), p_noP1, width = 6, height = 5, dpi = 600)
ggsave(file.path(OUT_BASE, "PCA_all_samples_excluding_P1.pdf"), p_noP1, width = 6, height = 5)

## -----------------------------
## Volcano helper (EnhancedVolcano + consistent x axis)
## -----------------------------
choose_tick_step <- function(xmax) {
  if (xmax <= 5) return(1)
  if (xmax <= 10) return(2)
  if (xmax <= 20) return(4)
  return(5)
}

make_volcano_plot <- function(res_df, title, xmax, breaks) {
  res_df$padj_plot <- res_df$padj
  res_df$padj_plot[is.na(res_df$padj_plot)] <- 1
  res_df$padj_plot <- pmax(res_df$padj_plot, .Machine$double.xmin)

  EnhancedVolcano(
    res_df,
    lab = res_df$label,
    x   = "log2FoldChange",
    y   = "padj_plot",
    pCutoff  = ALPHA,
    FCcutoff = LFC_CUTOFF,
    title = title,
    xlim = c(-xmax, xmax)
  ) + scale_x_continuous(breaks = breaks)
}

## -----------------------------
## DESeq2 per unique pair (DOWNSTREAM: P1 excluded)
## Also: additionally run a separate DESeq2 excluding chrX genes (noChrX)
## -----------------------------
pairs <- list(c("A6","H"), c("A6","P"), c("H","P"))

max_abs_lfc_seen_all  <- 0
max_abs_lfc_seen_noX  <- 0
contrast_dirs_all     <- character(0)
contrast_dirs_noX     <- character(0)
chrX_summary_rows     <- list()

run_pair_once <- function(g_den, g_num, txi_obj, meta_obj, mode = c("all", "noChrX")) {
  mode <- match.arg(mode)

  keep_samp <- meta_obj$group %in% c(g_den, g_num)
  meta_sub  <- droplevels(meta_obj[keep_samp, , drop = FALSE])

  if (nrow(meta_sub) < 2) stop("Not enough samples for contrast: ", g_num, " vs ", g_den)

  tab <- table(meta_sub$group)
  if (any(tab < 2)) {
    warning(
      "Contrast ", g_num, " vs ", g_den, " has <2 replicates in at least one group: ",
      paste(names(tab), tab, sep="=", collapse=", "),
      " (DESeq2 may be unstable / may error if no replication)."
    )
  }

  txi_sub <- subset_txi_to_samples(txi_obj, rownames(meta_sub))

  dds <- DESeqDataSetFromTximport(txi_sub, colData = meta_sub, design = ~ group)
  dds <- dds[rowSums(counts(dds)) > MIN_TOTAL_COUNT, ]

  if (mode == "noChrX") {
    chrs <- chr_of_geneid(rownames(dds))
    keep_rows <- is.na(chrs) | (chrs != "X")
    n_drop <- sum(!keep_rows, na.rm = TRUE)
    if (n_drop > 0) message("  [", g_num, " vs ", g_den, " | noChrX] Dropping ", n_drop, " chrX genes before DESeq2.")
    dds <- dds[keep_rows, ]
  }

  dds <- DESeq(dds, sfType = "poscounts")
  res <- results(dds, contrast = c("group", g_num, g_den), alpha = ALPHA)

  if (requireNamespace("apeglm", quietly = TRUE)) {
    rn <- resultsNames(dds)
    coef_forward <- paste0("group_", g_num, "_vs_", g_den)
    coef_reverse <- paste0("group_", g_den, "_vs_", g_num)

    if (coef_forward %in% rn) {
      res_shr <- lfcShrink(dds, coef = coef_forward, type = "apeglm", res = res)
    } else if (coef_reverse %in% rn) {
      res_tmp <- lfcShrink(dds, coef = coef_reverse, type = "apeglm")
      res_shr <- res
      res_shr$log2FoldChange <- -res_tmp$log2FoldChange
      if (!is.null(res_tmp$lfcSE)) res_shr$lfcSE <- res_tmp$lfcSE
      if (!is.null(res_tmp$stat))  res_shr$stat  <- -res_tmp$stat
    } else {
      res_shr <- res
      warning("Could not find coef for shrink; using unshrunken log2FoldChange for ", g_num, " vs ", g_den)
    }
  } else {
    res_shr <- res
    warning("apeglm not installed; using unshrunken log2FoldChange for ", g_num, " vs ", g_den)
  }

  res_df <- as.data.frame(res_shr)
  res_df$gene_id <- rownames(res_df)                # ENSG or SYMBOL, depending on tx2gene
  res_df$label   <- make_labels(res_df$gene_id)     # SYMBOL if possible, else fallback
  res_df$chr     <- chr_of_geneid(res_df$gene_id)   # CHR based on gene_id type
  res_df <- res_df[order(res_df$padj), ]

  list(dds = dds, res_df = res_df)
}

write_chrX_tables <- function(res_df, out_dir, dir_name) {
  res_x <- res_df[!is.na(res_df$chr) & res_df$chr == "X", , drop = FALSE]

  write.csv(res_x,
            file = file.path(out_dir, paste0(dir_name, "_DE_results_chrX.csv")),
            row.names = FALSE)

  sig_x <- res_x[!is.na(res_x$padj) &
                   (res_x$padj < ALPHA) &
                   (abs(res_x$log2FoldChange) >= LFC_CUTOFF), , drop = FALSE]

  write.csv(sig_x,
            file = file.path(out_dir, paste0(dir_name, "_DE_sig_padj_lfc_chrX.csv")),
            row.names = FALSE)

  up_x   <- sig_x[sig_x$log2FoldChange >=  LFC_CUTOFF, , drop = FALSE]
  down_x <- sig_x[sig_x$log2FoldChange <= -LFC_CUTOFF, , drop = FALSE]

  write.csv(up_x,
            file = file.path(out_dir, paste0(dir_name, "_DE_sig_UP_padj_lfc_chrX.csv")),
            row.names = FALSE)
  write.csv(down_x,
            file = file.path(out_dir, paste0(dir_name, "_DE_sig_DOWN_padj_lfc_chrX.csv")),
            row.names = FALSE)

  n_x_tested <- nrow(res_x)
  n_x_sig    <- nrow(sig_x)
  frac_x_sig <- if (n_x_tested > 0) n_x_sig / n_x_tested else NA_real_

  data.frame(
    comparison = dir_name,
    n_chrX_tested = n_x_tested,
    n_chrX_sig_padj_and_lfc = n_x_sig,
    frac_chrX_sig = frac_x_sig,
    stringsAsFactors = FALSE
  )
}

write_outputs_for_direction <- function(dir_name, dds, res_df, add_chrX = TRUE, add_pca = TRUE) {
  out_dir <- file.path(OUT_BASE, dir_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  norm_counts <- counts(dds, normalized = TRUE)
  write.csv(as.data.frame(norm_counts),
            file = file.path(out_dir, paste0(dir_name, "_normalized_counts.csv")))

  write.csv(res_df,
            file = file.path(out_dir, paste0(dir_name, "_DE_results.csv")),
            row.names = FALSE)

  sig_df <- res_df[!is.na(res_df$padj) &
                     (res_df$padj < ALPHA) &
                     (abs(res_df$log2FoldChange) >= LFC_CUTOFF), , drop = FALSE]
  write.csv(sig_df, file = file.path(out_dir, paste0(dir_name, "_DE_sig_padj_lfc.csv")), row.names = FALSE)

  up_df   <- sig_df[sig_df$log2FoldChange >=  LFC_CUTOFF, , drop = FALSE]
  down_df <- sig_df[sig_df$log2FoldChange <= -LFC_CUTOFF, , drop = FALSE]
  write.csv(up_df,   file = file.path(out_dir, paste0(dir_name, "_DE_sig_UP_padj_lfc.csv")), row.names = FALSE)
  write.csv(down_df, file = file.path(out_dir, paste0(dir_name, "_DE_sig_DOWN_padj_lfc.csv")), row.names = FALSE)

  if (add_pca) {
    vsd <- vst(dds, blind = TRUE)
    p_pca <- make_pca_bubbles(vsd, intgroup_col = "group", title = dir_name)
    ggsave(file.path(out_dir, paste0(dir_name, "_PCA.png")), p_pca, width = 6, height = 5, dpi = 600)
    ggsave(file.path(out_dir, paste0(dir_name, "_PCA.pdf")), p_pca, width = 6, height = 5)
  }

  if (add_chrX) {
    chrX_summary_rows[[dir_name]] <<- write_chrX_tables(res_df, out_dir, dir_name)
  }

  out_dir
}

## Run DE (P1 excluded) + write result tables first (ALL genes, then noChrX), tracking global x-ranges
for (pr in pairs) {
  g1 <- pr[1]
  g2 <- pr[2]

  message("\n=== Running unique pair (P1 excluded): ", g2, " vs ", g1, " (and reverse) ===")

  ## A) Standard: ALL genes
  obj <- run_pair_once(g_den = g1, g_num = g2, txi_obj = txi_main, meta_obj = meta_main, mode = "all")
  dds    <- obj$dds
  res_df <- obj$res_df

  lfc <- res_df$log2FoldChange
  lfc <- lfc[is.finite(lfc)]
  if (length(lfc) > 0) {
    max_here <- if (VOLCANO_XLIM_MODE == "p99") as.numeric(quantile(abs(lfc), 0.99, na.rm = TRUE)) else max(abs(lfc), na.rm = TRUE)
    max_abs_lfc_seen_all <- max(max_abs_lfc_seen_all, max_here)
  }

  dir_fwd <- paste0(g2, "_vs_", g1)
  out_fwd <- write_outputs_for_direction(dir_fwd, dds, res_df, add_chrX = TRUE, add_pca = TRUE)
  contrast_dirs_all <- c(contrast_dirs_all, out_fwd)

  res_rev <- res_df
  res_rev$log2FoldChange <- -res_rev$log2FoldChange
  if ("stat" %in% colnames(res_rev) && is.numeric(res_rev$stat)) res_rev$stat <- -res_rev$stat

  dir_rev <- paste0(g1, "_vs_", g2)
  out_rev <- write_outputs_for_direction(dir_rev, dds, res_rev, add_chrX = TRUE, add_pca = TRUE)
  contrast_dirs_all <- c(contrast_dirs_all, out_rev)

  ## B) Additional: DESeq2 excluding chrX genes (noChrX)
  message("  -> Also running DESeq2 excluding chrX genes (noChrX) for this pair.")
  obj_noX <- run_pair_once(g_den = g1, g_num = g2, txi_obj = txi_main, meta_obj = meta_main, mode = "noChrX")
  dds_noX    <- obj_noX$dds
  res_df_noX <- obj_noX$res_df

  lfc2 <- res_df_noX$log2FoldChange
  lfc2 <- lfc2[is.finite(lfc2)]
  if (length(lfc2) > 0) {
    max_here2 <- if (VOLCANO_XLIM_MODE == "p99") as.numeric(quantile(abs(lfc2), 0.99, na.rm = TRUE)) else max(abs(lfc2), na.rm = TRUE)
    max_abs_lfc_seen_noX <- max(max_abs_lfc_seen_noX, max_here2)
  }

  dir_fwd_noX <- paste0(g2, "_vs_", g1, "_noChrX")
  out_fwd_noX <- write_outputs_for_direction(dir_fwd_noX, dds_noX, res_df_noX, add_chrX = FALSE, add_pca = FALSE)
  contrast_dirs_noX <- c(contrast_dirs_noX, out_fwd_noX)

  res_rev_noX <- res_df_noX
  res_rev_noX$log2FoldChange <- -res_rev_noX$log2FoldChange
  if ("stat" %in% colnames(res_rev_noX) && is.numeric(res_rev_noX$stat)) res_rev_noX$stat <- -res_rev_noX$stat

  dir_rev_noX <- paste0(g1, "_vs_", g2, "_noChrX")
  out_rev_noX <- write_outputs_for_direction(dir_rev_noX, dds_noX, res_rev_noX, add_chrX = FALSE, add_pca = FALSE)
  contrast_dirs_noX <- c(contrast_dirs_noX, out_rev_noX)
}

## Decide consistent x-axes
xmax_all <- ceiling(max_abs_lfc_seen_all + VOLCANO_XPAD)
tick_all <- choose_tick_step(xmax_all)
xbreaks_all <- seq(-xmax_all, xmax_all, by = tick_all)
message("\n[ALL genes] Volcano x-axis: [-", xmax_all, ", ", xmax_all, "] with tick step ", tick_all)

xmax_noX <- ceiling(max_abs_lfc_seen_noX + VOLCANO_XPAD)
tick_noX <- choose_tick_step(xmax_noX)
xbreaks_noX <- seq(-xmax_noX, xmax_noX, by = tick_noX)
message("[noChrX]   Volcano x-axis: [-", xmax_noX, ", ", xmax_noX, "] with tick step ", tick_noX)

## Second pass: generate volcano plots
for (out_dir in contrast_dirs_all) {
  dir_name <- basename(out_dir)
  res_file <- file.path(out_dir, paste0(dir_name, "_DE_results.csv"))
  if (!file.exists(res_file)) next

  res_df <- read.csv(res_file, check.names = FALSE)
  if (!("gene_id" %in% colnames(res_df)) && ("gene" %in% colnames(res_df))) res_df$gene_id <- res_df$gene
  if (!("label" %in% colnames(res_df))) res_df$label <- make_labels(res_df$gene_id %||% res_df$gene_id)
  if (!("chr" %in% colnames(res_df)))   res_df$chr   <- chr_of_geneid(res_df$gene_id)

  p_all <- make_volcano_plot(res_df, title = dir_name, xmax = xmax_all, breaks = xbreaks_all)
  ggsave(file.path(out_dir, paste0(dir_name, "_volcano.png")), p_all, width = 12, height = 10, dpi = 300, limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(dir_name, "_volcano.pdf")), p_all, width = 12, height = 10, limitsize = FALSE)

  res_x <- res_df[!is.na(res_df$chr) & res_df$chr == "X", , drop = FALSE]
  if (nrow(res_x) > 0) {
    p_x <- make_volcano_plot(res_x, title = paste0(dir_name, " (chrX only)"), xmax = xmax_all, breaks = xbreaks_all)
    ggsave(file.path(out_dir, paste0(dir_name, "_volcano_chrX.png")), p_x, width = 12, height = 10, dpi = 300, limitsize = FALSE)
    ggsave(file.path(out_dir, paste0(dir_name, "_volcano_chrX.pdf")), p_x, width = 12, height = 10, limitsize = FALSE)
  } else {
    message("No chrX genes found/mapped for ", dir_name, " — skipping chrX volcano.")
  }
}

for (out_dir in contrast_dirs_noX) {
  dir_name <- basename(out_dir)
  res_file <- file.path(out_dir, paste0(dir_name, "_DE_results.csv"))
  if (!file.exists(res_file)) next

  res_df <- read.csv(res_file, check.names = FALSE)
  if (!("gene_id" %in% colnames(res_df)) && ("gene" %in% colnames(res_df))) res_df$gene_id <- res_df$gene
  if (!("label" %in% colnames(res_df))) res_df$label <- make_labels(res_df$gene_id %||% res_df$gene_id)
  if (!("chr" %in% colnames(res_df)))   res_df$chr   <- chr_of_geneid(res_df$gene_id)

  p_noX <- make_volcano_plot(res_df, title = paste0(dir_name, " (chrX excluded in DESeq2)"), xmax = xmax_noX, breaks = xbreaks_noX)
  ggsave(file.path(out_dir, paste0(dir_name, "_volcano.png")), p_noX, width = 12, height = 10, dpi = 300, limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(dir_name, "_volcano.pdf")), p_noX, width = 12, height = 10, limitsize = FALSE)
}

## -----------------------------
## Write chrX summary + optional Excel workbook (ALL-genes contrasts only)
## -----------------------------
chrX_summary_df <- do.call(rbind, chrX_summary_rows)
write.csv(chrX_summary_df,
          file = file.path(OUT_BASE, "chrX_significance_summary_all_contrasts.csv"),
          row.names = FALSE)

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "summary_chrX")
  openxlsx::writeData(wb, "summary_chrX", chrX_summary_df)

  for (out_dir in contrast_dirs_all) {
    dir_name <- basename(out_dir)
    f_x <- file.path(out_dir, paste0(dir_name, "_DE_results_chrX.csv"))
    if (!file.exists(f_x)) next
    df_x <- read.csv(f_x, check.names = FALSE)

    sheet <- substr(dir_name, 1, 31)
    if (sheet %in% names(wb)) sheet <- paste0(substr(sheet, 1, 28), "_x")
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, df_x)
  }

  openxlsx::saveWorkbook(wb,
                         file = file.path(OUT_BASE, "chrX_results_all_contrasts.xlsx"),
                         overwrite = TRUE)
} else {
  message("openxlsx not installed; skipping .xlsx workbook (CSVs were written).")
}

message("\nDone.")
message("Outputs in: ", OUT_BASE)
message("Downstream DE excludes samples: ", ifelse(length(present_excl) == 0, "(none found)", paste(present_excl, collapse = ", ")))
message("Extra DESeq2 runs written as *_noChrX directories (chrX excluded BEFORE DESeq2).")