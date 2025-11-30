
suppressPackageStartupMessages({
  library(tidyverse)       # dplyr/ggplot2/readr/tibble等
  library(rtracklayer)     # import(GTF)
  library(GenomicRanges)
  library(GenomicFeatures) # makeTxDbFromGFF, transcriptLengths
})

## ========= 参数 =========
gtf_file <- "Mus_musculus.GRCm39.112.chr.gtf.gz"
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

