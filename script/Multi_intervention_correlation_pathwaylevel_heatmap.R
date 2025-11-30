## =====================================
## GSEA on combined gene sets (Hallmark + GO BP + Reactome)
## and correlation of NES between interventions
## Only pathways with padj < 0.1 in at least one signature are used.
## Heatmap: color = Spearman correlation, numbers = significance stars (from NES-correlation padj)
## =====================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(msigdbr)
  library(fgsea)
  library(pheatmap)
})

set.seed(123)

## --------- 0) 前提检查 ---------
if (!exists("de_all_filt") || !exists("sig_meta_filtered")) {
  stop("需要先按前面的脚本构建 de_all_filt 和 sig_meta_filtered 再运行本段代码。")
}

## 保证有 label
stopifnot(all(c("signature_id","label") %in% colnames(sig_meta_filtered)))

sig_ids   <- sig_meta_filtered$signature_id
sig_label <- sig_meta_filtered$label
names(sig_label) <- sig_ids   # signature_id -> label 的映射

## --------- 1) 构建 prerank 统计量：优先使用 DESeq2 Wald stat ---------

make_ranks <- function(df) {
  ## df 来自 de_all_filt 的某个 signature
  if (!("gene" %in% names(df))) stop("df 中缺少 gene 列")
  if (!("ranking_p" %in% names(df))) stop("df 中缺少 ranking_p 列；请确认用的是之前的 load_deseq_result。")
  
  if ("stat" %in% names(df) && any(is.finite(df$stat))) {
    r <- df$stat
  } else {
    r <- sign(df$log2FoldChange) * (-log10(df$ranking_p + 1e-300))
  }
  names(r) <- df$gene
  r <- r[is.finite(r)]
  r <- sort(r, decreasing = TRUE)
  r
}

de_split <- de_all_filt %>%
  filter(signature_id %in% sig_ids) %>%
  split(.$signature_id)

ranks_list <- lapply(de_split, make_ranks)

## --------- 2) 取 MSigDB 基因集：Hallmark + GO BP + Reactome ---------

## Hallmark (H)
msig_H <- msigdbr(
  species  = "Mus musculus",
  category = "H"
) %>%
  select(gs_name, gene_symbol)

## GO Biological Process (C5:BP)
msig_GOBP <- msigdbr(
  species     = "Mus musculus",
  category    = "C5",
  subcategory = "BP"
) %>%
  select(gs_name, gene_symbol)

## Reactome (C2:REACTOME)
msig_REACT <- msigdbr(
  species     = "Mus musculus",
  category    = "C2",
  subcategory = "CP:REACTOME"
) %>%
  select(gs_name, gene_symbol)

make_gs_list <- function(df) {
  split(df$gene_symbol, df$gs_name)
}

gs_H     <- make_gs_list(msig_H)
gs_GOBP  <- make_gs_list(msig_GOBP)
gs_REACT <- make_gs_list(msig_REACT)

## --------- 3) 对每个 signature 做 GSEA（三套基因集） ---------

run_fgsea_for_collection <- function(gs_list, ranks_list, collection_name,
                                     minSize = 15, maxSize = 500) {
  res_list <- list()
  for (sid in names(ranks_list)) {
    rnk <- ranks_list[[sid]]
    if (length(rnk) < minSize) next
    
    fg <- fgsea::fgsea(
      pathways = gs_list,
      stats    = rnk,
      minSize  = minSize,
      maxSize  = maxSize,
      eps      = 1e-50
    )
    if (!nrow(fg)) next
    
    fg2 <- fg %>%
      transmute(
        signature_id = sid,
        collection   = collection_name,
        pathway      = pathway,
        NES          = NES,
        pval         = pval,
        padj         = padj,
        size         = size
      )
    res_list[[sid]] <- fg2
  }
  if (length(res_list) == 0) {
    return(NULL)
  } else {
    return(bind_rows(res_list))
  }
}

gsea_H     <- run_fgsea_for_collection(gs_H,     ranks_list, "Hallmark")
gsea_GOBP  <- run_fgsea_for_collection(gs_GOBP,  ranks_list, "GOBP")
gsea_REACT <- run_fgsea_for_collection(gs_REACT, ranks_list, "Reactome")

gsea_all <- bind_rows(
  if (!is.null(gsea_H)) gsea_H else tibble(),
  if (!is.null(gsea_GOBP)) gsea_GOBP else tibble(),
  if (!is.null(gsea_REACT)) gsea_REACT else tibble()
)

if (!nrow(gsea_all)) stop("GSEA 没有得到任何结果，请检查 ranks 或基因集。")

## --------- 4) 只保留至少在一个 signature 中 padj < 0.1 的通路 ---------

gsea_sig <- gsea_all %>%
  group_by(collection, pathway) %>%
  filter(any(padj < 0.1, na.rm = TRUE)) %>%
  ungroup()

## 给每个 (collection, pathway) 一个唯一 ID
gsea_sig <- gsea_sig %>%
  mutate(
    path_id = paste(collection, pathway, sep = "::")
  )

## 加上 label（便于列名好看）
gsea_sig <- gsea_sig %>%
  left_join(
    sig_meta_filtered %>% select(signature_id, label),
    by = "signature_id"
  )

## --------- 5) 构建 NES 矩阵：行 = 通路，列 = signatures(label) ---------

nes_mat_long <- gsea_sig %>%
  select(path_id, label, NES) %>%
  distinct()

nes_mat <- nes_mat_long %>%
  tidyr::pivot_wider(
    id_cols     = path_id,
    names_from  = label,
    values_from = NES
  )

nes_mat <- nes_mat %>%
  column_to_rownames("path_id") %>%
  as.matrix()

## 至少有 2 个非 NA 的通路
nes_mat <- nes_mat[rowSums(!is.na(nes_mat)) >= 2, , drop = FALSE]

## 列顺序按 sig_meta_filtered$label 排序（方便和前面的 logFC 图对应）
common_labels <- intersect(sig_meta_filtered$label, colnames(nes_mat))
nes_mat <- nes_mat[, common_labels, drop = FALSE]

n_sig <- ncol(nes_mat)
if (n_sig < 2) stop("可用的 signature 数不足 2，无法计算相关性。")

sig_names <- colnames(nes_mat)

## --------- 6) 相关矩阵 & p 值（Spearman） ---------

cor_mat <- suppressWarnings(
  cor(nes_mat, method = "spearman", use = "pairwise.complete.obs")
)

raw_p <- matrix(NA_real_, n_sig, n_sig, dimnames = list(sig_names, sig_names))

for (i in seq_len(n_sig)) {
  for (j in i:n_sig) {
    if (i == j) {
      raw_p[i, j] <- 0
      next
    }
    x <- nes_mat[, i]
    y <- nes_mat[, j]
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 5) next
    ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
    raw_p[i, j] <- ct$p.value
    raw_p[j, i] <- ct$p.value
  }
}

## 对上三角做 BH 校正（得到 “exact adjusted p-values”）
upper_idx <- upper.tri(raw_p, diag = FALSE)
p_vec <- raw_p[upper_idx]

valid <- which(!is.na(p_vec))
padj_vec <- rep(NA_real_, length(p_vec))
padj_vec[valid] <- p.adjust(p_vec[valid], method = "BH")

padj_mat <- matrix(NA_real_, n_sig, n_sig, dimnames = list(sig_names, sig_names))
padj_mat[upper_idx] <- padj_vec
padj_mat[lower.tri(padj_mat)] <- t(padj_mat)[lower.tri(padj_mat)]
diag(padj_mat) <- 0  # 自己与自己设为 0

## --------- 7) 用 * 表示显著性 (来自 NES-correlation padj) ---------
##  padj <= 0.001: ***
##  padj <= 0.01 : **
##  padj <= 0.05 : *
##  其余为空

star_mat <- matrix("", n_sig, n_sig, dimnames = dimnames(padj_mat))

for (i in seq_len(n_sig)) {
  for (j in seq_len(n_sig)) {
    p <- padj_mat[i, j]
    if (is.na(p) || i == j) {
      star_mat[i, j] <- ""
    } else if (p <= 0.001) {
      star_mat[i, j] <- "***"
    } else if (p <= 0.01) {
      star_mat[i, j] <- "**"
    } else if (p <= 0.05) {
      star_mat[i, j] <- "*"
    } else {
      star_mat[i, j] <- ""
    }
  }
}

## --------- 8) 画 NES 相关性热图：颜色=相关系数，格子=星号 ---------

col_fun <- colorRampPalette(c("blue", "white", "red"))
range_val <- max(abs(cor_mat), na.rm = TRUE)
if (!is.finite(range_val) || range_val == 0) range_val <- 1
breaks <- seq(-range_val, range_val, length.out = 101)

out_dir_cor2 <- "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/decoy_rnaseq/all_quantification/gsea_all_collections_maxAge_with_17aE2"
dir.create(out_dir_cor2, showWarnings = FALSE, recursive = TRUE)

pdf_file <- file.path(out_dir_cor2, "NES_spearman_correlation_heatmap_Hallmark_GOBP_Reactome_with_stars.pdf")

pdf(pdf_file, width = 6, height = 6)
pheatmap(
  cor_mat,
  main             = "Spearman correlation of NES\n(Hallmark + GO BP + Reactome; pathways with padj<0.1 in ≥1 signature)",
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  treeheight_row   = 0,
  treeheight_col   = 0,
  color            = col_fun(100),
  breaks           = breaks,
  border_color     = NA,
  display_numbers  = star_mat,
  number_color     = "black",
  fontsize_number  = 6
)
dev.off()

## 保存矩阵（方便后面在 R 或 Python 里再玩）
nes_tsv   <- file.path(out_dir_cor2, "NES_matrix_Hallmark_GOBP_Reactome_used.tsv")
cor_tsv   <- file.path(out_dir_cor2, "NES_spearman_correlation_matrix_Hallmark_GOBP_Reactome.tsv")
padj_tsv  <- file.path(out_dir_cor2, "NES_spearman_correlation_padj_matrix_Hallmark_GOBP_Reactome.tsv")

write.table(nes_mat,  file = nes_tsv,  sep = "\t", quote = FALSE, row.names = TRUE)
write.table(cor_mat,  file = cor_tsv,  sep = "\t", quote = FALSE, row.names = TRUE)
write.table(padj_mat, file = padj_tsv, sep = "\t", quote = FALSE, row.names = TRUE)

message("GSEA NES 相关性热图(星号显著性)已输出: ", pdf_file)
message("NES 矩阵: ", nes_tsv)
message("相关性矩阵: ", cor_tsv)
message("相关性 padj 矩阵: ", padj_tsv)
