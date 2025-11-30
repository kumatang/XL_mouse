
suppressPackageStartupMessages({
  library(tidyverse)       # dplyr/ggplot2/readr/tibble等
  library(rtracklayer)     # import(GTF)
  library(GenomicRanges)
  library(GenomicFeatures) # makeTxDbFromGFF, transcriptLengths
})

## ========= 参数 =========
gtf_file <- "decoy_rnaseq/Mus_musculus.GRCm39.112.chr.gtf.gz"
out_tsv  <- "gene_symbol_key.lengths.genomic_tx_cds.tsv"

## ========= 读取注释，构建 TxDb =========
message(">> Import GTF & build TxDb ...")
gtf  <- rtracklayer::import(gtf_file)
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")

## ========= 仅保留 protein_coding 基因的 genomic length（按 symbol 聚合）=========
message(">> Compute genomic length (protein_coding genes) ...")
genes_pc_gr <- gtf[gtf$type == "gene" & gtf$gene_biotype == "protein_coding"]

gene_len_by_symbol <- tibble::tibble(
  gene_symbol    = genes_pc_gr$gene_name,
  genomic_length = width(genes_pc_gr)
) |>
  dplyr::filter(!is.na(gene_symbol)) |>
  dplyr::group_by(gene_symbol) |>
  # 如果更偏好最大跨度，把 median(...) 改成 max(...)
  dplyr::summarise(
    genomic_length = median(genomic_length),
    n_gene_ids     = dplyr::n(),
    .groups = "drop"
  )

## ========= 构建 protein_coding transcript -> gene_symbol 映射 =========
message(">> Map protein_coding transcripts to gene_symbol ...")
has_tx_biotype <- "transcript_biotype" %in% colnames(mcols(gtf))
if (has_tx_biotype) {
  tx_rows <- gtf[
    gtf$type == "transcript" &
      !is.na(gtf$transcript_id) &
      gtf$gene_biotype == "protein_coding" &
      gtf$transcript_biotype == "protein_coding"
  ]
} else {
  # 某些 GTF 无 transcript_biotype 字段：退化为 gene_biotype 过滤
  tx_rows <- gtf[
    gtf$type == "transcript" &
      !is.na(gtf$transcript_id) &
      gtf$gene_biotype == "protein_coding"
  ]
}

tx2symbol_pc <- tibble::tibble(
  transcript_id = tx_rows$transcript_id,
  gene_symbol   = tx_rows$gene_name
) |>
  dplyr::filter(!is.na(transcript_id), !is.na(gene_symbol)) |>
  dplyr::distinct()

## ========= 每条转录本的长度（外显子并集）与 CDS 长度 =========
message(">> Compute per-transcript exon-union length and CDS length ...")
tl <- GenomicFeatures::transcriptLengths(
  txdb,
  with.cds_len  = TRUE,   # 提供 cds_len
  with.utr5_len = FALSE,
  with.utr3_len = FALSE
)

tl_tbl <- tibble::as_tibble(tl) |>
  dplyr::rename(transcript_id = tx_name) |>
  dplyr::select(transcript_id, gene_id, tx_len, cds_len)

## ========= 仅限 protein_coding 转录本，并合并 symbol =========
tl_pc <- tl_tbl |>
  dplyr::inner_join(tx2symbol_pc, by = "transcript_id")

## ========= 在 gene_symbol 层面聚合（Stoeger 口径：取中位数）=========
message(">> Aggregate to gene_symbol (median across transcripts) ...")
symbol_lengths <- tl_pc |>
  dplyr::group_by(gene_symbol) |>
  dplyr::summarise(
    tx_length_median  = median(tx_len,  na.rm = TRUE),
    cds_length_median = if (all(is.na(cds_len))) NA_real_ else median(cds_len, na.rm = TRUE),
    n_tx              = dplyr::n(),
    n_tx_with_cds     = sum(!is.na(cds_len)),
    .groups = "drop"
  )

## ========= 合并三类长度到一个表（以 gene_symbol 为主键）=========
message(">> Join genomic/transcript/CDS lengths ...")
gene_symbol_lengths <- gene_len_by_symbol |>
  dplyr::right_join(symbol_lengths, by = "gene_symbol") |>
  dplyr::relocate(gene_symbol)


## ========= 输出 =========
message(">> Preview:")
print(gene_symbol_lengths, n = 10)

message(">> Write TSV: ", out_tsv)
readr::write_tsv(gene_symbol_lengths, out_tsv)

message(">> Done.")



suppressPackageStartupMessages({
  library(tidyverse)  # dplyr/ggplot2/readr/stringr等
  library(readr)
  library(ggplot2)
  library(fgsea)
})

# 结果根目录（沿用你的差异分析输出结构）
if (!exists("out_root")) out_root <- "deseq2_by_subset"
stopifnot(dir.exists(out_root))

# 基因长度表（使用你前面生成的）
lengths_tsv <- "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/revised/XL/gene_symbol_key.lengths.genomic_tx_cds.tsv"
stopifnot(file.exists(lengths_tsv))
len_tbl <- readr::read_tsv(lengths_tsv, show_col_types = FALSE) %>%
  dplyr::select(gene_symbol, genomic_length) %>%
  dplyr::rename(Genomic = genomic_length)

# 仅分析 A2（Old_Heter_vs_WT）
a2_name <- "A2_Old_Heter_vs_WT.tsv"
subdirs <- list.dirs(out_root, full.names = TRUE, recursive = FALSE)
has_a2 <- file.exists(file.path(subdirs, a2_name))
subdirs <- subdirs[has_a2]
stopifnot(length(subdirs) > 0)

# 子集标签（Tissue_Gender）
lab_from_dir <- function(p) {
  nm <- basename(p)
  ps <- strsplit(nm, "_")[[1]]
  if (length(ps) >= 2) paste0(stringr::str_to_title(ps[1]), "_", stringr::str_to_title(ps[2])) else nm
}

# 阈值
thr_lfc  <- log2(1.5)
thr_padj <- 0.05


# 排序函数：优先用 DESeq2 的 stat，否则用 sign(log2FC)*-log10(padj)
make_rank <- function(df) {
  if ("stat" %in% names(df) && any(is.finite(df$stat))) {
    sc <- df %>%
      dplyr::transmute(gene, score = stat)
  } else {
    sc <- df %>%
      dplyr::mutate(
        score = sign(log2FoldChange) * (-log10(padj + 1e-300))
      ) %>%
      dplyr::transmute(gene, score)
  }
  sc %>%
    dplyr::filter(is.finite(score)) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(score))
}

# 输出目录
out_dir <- file.path(out_root, "length_correlation_genomic")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

wilcox_rows      <- list()  # Wilcoxon 结果（无BH）

for (sd in subdirs) {
  subset_label <- lab_from_dir(sd)
  message("Processing subset: ", subset_label)

  # 读 A2 结果（保留 stat 列以便排序）
  df <- read.delim(file.path(sd, a2_name), check.names = FALSE)
  if (!("gene" %in% colnames(df))) df <- tibble::rownames_to_column(df, var = "gene")
  df <- df %>%
    dplyr::select(
      gene,
      log2FoldChange,
      padj,
      dplyr::any_of("stat")   # 有 stat 就保留，没有就忽略
    )

  # 合并 Genomic 长度 & 计算 log10
  dat <- df %>%
    dplyr::left_join(len_tbl, by = c("gene" = "gene_symbol")) %>%
    dplyr::mutate(
      log10_Genomic = ifelse(is.finite(Genomic) & Genomic > 0, log10(Genomic), NA_real_),
      SigDir = dplyr::case_when(
        !is.na(padj) & padj <= thr_padj & !is.na(log2FoldChange) & log2FoldChange >  thr_lfc ~ "Up",
        !is.na(padj) & padj <= thr_padj & !is.na(log2FoldChange) & log2FoldChange < -thr_lfc ~ "Down",
        TRUE ~ "NS"
      ),
      SigDir = factor(SigDir, levels = c("Down","NS","Up"))
    )

  ## ---- Boxplot：显著 Up/Down 的 Genomic 长度分布 + Wilcoxon（两侧）----
  long_box <- dat %>%
    dplyr::filter(SigDir %in% c("Up","Down")) %>%
    dplyr::select(gene, SigDir, log10_Genomic) %>%
    dplyr::filter(is.finite(log10_Genomic))

  if (nrow(long_box) > 0 && length(unique(long_box$SigDir)) == 2) {
    # Wilcoxon（两侧）
    w_p <- tryCatch(
      wilcox.test(log10_Genomic ~ SigDir, data = long_box, exact = FALSE)$p.value,
      error = function(e) NA_real_
    )
    y_max <- max(long_box$log10_Genomic, na.rm = TRUE)
    y_min <- min(long_box$log10_Genomic, na.rm = TRUE)
    stats_box <- tibble::tibble(
      p = w_p,
      label = ifelse(is.finite(w_p),
                     paste0("Wilcoxon p=", formatC(w_p, format = "e", digits = 2)),
                     "Wilcoxon p=NA"),
      x = 1.5,
      y = y_max + (y_max - y_min) * 0.08,
      n_up = sum(long_box$SigDir == "Up"),
      n_down = sum(long_box$SigDir == "Down")
    )
    wilcox_rows[[length(wilcox_rows) + 1]] <- data.frame(
      Subset = subset_label,
      n_up = stats_box$n_up, n_down = stats_box$n_down, p = stats_box$p
    )

    p_box <- ggplot(long_box, aes(x = SigDir, y = log10_Genomic)) +
      geom_boxplot(fill = NA, color = "black", width = 0.55, outlier.shape = NA) +
      geom_jitter(aes(color = SigDir), width = 0.15, size = 0.2, alpha = 1) +
      scale_color_manual(values = c(Down = "#3B6DB3B2", Up = "#CC0C00B2"), name = "Direction") +
      labs(title = subset_label,
           x = NULL, y = "log10(genomic length, bp)") +
      theme_classic(base_size = 11) +
      theme(plot.title = element_text(face = "bold")) +
      geom_text(
        data = stats_box,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE, size = 3.2
      ) +
      scale_x_discrete(limits = c("Up","Down"))

    ggsave(filename = file.path(out_dir, paste0(subset_label, "_A2_boxplot_sig_log10Genomic.pdf")),
           plot = p_box, width = 2.8, height = 2.3, units = "in")
  } else {
    message("  [", subset_label, "] No both Up and Down significant genes -> skip boxplot & test.")
  }
}


