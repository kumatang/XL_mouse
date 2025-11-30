## =====================================
## 0) 加载 R 包
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
  library(ggplot2)
})

## （可选）强制把 dplyr::select 绑定到当前环境，避免被 AnnotationDbi::select 干扰
select <- dplyr::select

## =====================================
## 1) 定义需要做 GSEA 的 DESeq2 结果文件
##    只保留 XL 的 A2_Old_Heter_vs_WT：
##    - Liver:   male / female
##    - Muscle:  male / female
##    - Skin:    male / female
## =====================================

xl_root <- "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/revised/XL/deseq2_by_subset"

sig_defs_all <- tibble(
  dataset         = "XL",
  tissue          = c("Liver", "Liver", "Muscle", "Muscle", "Skin", "Skin"),
  sex             = c("male",  "female","male",   "female","male", "female"),
  age_label       = "Old",
  condition_short = "Heter"
) %>%
  dplyr::mutate(
    subdir = paste0(tolower(tissue), "_", sex),
    file   = file.path(xl_root, subdir, "A2_Old_Heter_vs_WT.tsv"),
    label  = paste(dataset, tissue, sex, age_label, condition_short, sep = "_")
  ) %>%
  dplyr::select(dataset, tissue, sex, age_label, condition_short, file, label)

## 检查文件是否存在
missing_files <- sig_defs_all$file[!file.exists(sig_defs_all$file)]
if (length(missing_files) > 0) {
  stop("以下 DESeq2 结果文件不存在，请检查路径或文件名:\n",
       paste(missing_files, collapse = "\n"))
}

## 输出 signature 定义，方便检查
out_dir <- file.path(
  "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/decoy_rnaseq/all_quantification",
  "GSEA_Hallmark_spearman_cor"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.table(
  sig_defs_all,
  file      = file.path(out_dir, "signature_defs_used.tsv"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

message("共发现 signature 数量: ", nrow(sig_defs_all))


## =====================================
## 2) 从 DESeq2 结果提取 prerank 统计量 (优先用 Wald stat)
## =====================================

load_deseq_for_gsea <- function(file) {
  df <- read.table(
    file,
    header      = TRUE,
    sep         = "\t",
    quote       = "",
    check.names = FALSE
  )
  ## 找 gene 列名
  gene_col <- if ("gene" %in% colnames(df)) {
    "gene"
  } else if ("gene_symbol" %in% colnames(df)) {
    "gene_symbol"
  } else {
    stop("找不到 gene/gene_symbol 列: ", file)
  }

  need_cols <- c(gene_col, "log2FoldChange")
  miss <- setdiff(need_cols, colnames(df))
  if (length(miss) > 0) {
    stop("文件缺少必要列 ", paste(miss, collapse = ","), " : ", file)
  }

  has_stat <- "stat" %in% colnames(df)

  df2 <- df %>%
    dplyr::rename(gene = !!gene_col) %>%
    dplyr::filter(!is.na(gene), !is.na(log2FoldChange))

  if (has_stat) {
    df2 <- df2 %>%
      dplyr::filter(!is.na(stat))
  } else {
    message("文件缺少 stat 列，使用 log2FoldChange 作为排序: ", file)
  }

  ## 按绝对值排序，保留每个基因一条记录
  if (has_stat) {
    df2 <- df2 %>%
      dplyr::arrange(dplyr::desc(abs(stat))) %>%
      dplyr::distinct(gene, .keep_all = TRUE)
    stats_vec <- df2$stat
  } else {
    df2 <- df2 %>%
      dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
      dplyr::distinct(gene, .keep_all = TRUE)
    stats_vec <- df2$log2FoldChange
  }

  names(stats_vec) <- df2$gene
  stats_vec <- stats_vec[!is.na(stats_vec)]
  stats_vec
}

## 加载所有 signature 的 stats 向量
stats_list <- list()
for (i in seq_len(nrow(sig_defs_all))) {
  lab  <- sig_defs_all$label[i]
  file <- sig_defs_all$file[i]
  message("加载 DESeq2 结果用于 GSEA: ", lab, " <- ", file)
  stats_vec <- load_deseq_for_gsea(file)
  stats_list[[lab]] <- stats_vec
}


## =====================================
## 3) 准备 MSigDB Hallmark 基因集（只用 Hallmark）
## =====================================

msig_h <- msigdbr(
  species    = "Mus musculus",
  collection = "H"
)

msig_to_list <- function(msig_df) {
  msig_df %>%
    split(.$gs_name) %>%
    lapply(function(df) unique(df$gene_symbol))
}

hallmark_pathways <- msig_to_list(msig_h)


## =====================================
## 4) 对每个 signature 做 Hallmark GSEA (fgseaMultilevel)
## =====================================

gsea_res_list <- list()

for (lab in names(stats_list)) {
  stats_vec <- stats_list[[lab]]
  message("GSEA (Hallmark): ", lab)

  gres <- tryCatch(
    {
      fgseaMultilevel(
        pathways = hallmark_pathways,
        stats    = stats_vec,
        minSize  = 15,
        maxSize  = 500
      )
    },
    error = function(e) {
      message("  fgseaMultilevel 出错: label = ", lab,
              ", error = ", e$message)
      NULL
    }
  )

  if (!is.null(gres) && nrow(gres) > 0) {
    gres2 <- tibble::as_tibble(gres) %>%
      dplyr::select(pathway, pval, padj, ES, NES) %>%
      dplyr::mutate(
        collection = "Hallmark",
        label      = lab
      )

    gsea_res_list[[length(gsea_res_list) + 1]] <- gres2
  }
}

gsea_all <- dplyr::bind_rows(gsea_res_list)

## 保存 Hallmark GSEA 结果
gsea_all_tsv <- file.path(out_dir, "GSEA_Hallmark_long.tsv")
write.table(
  gsea_all,
  file      = gsea_all_tsv,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
message("已写出 Hallmark GSEA 结果: ", gsea_all_tsv)



## =====================================
## 5) Hallmark bubble plot across signatures（cluster）
##    只画这 6 个 Old_Heter signatures
## =====================================

gsea_all_long <- readr::read_tsv(gsea_all_tsv, show_col_types = FALSE)


## 构造 signature_label_clean（这里基本等于 label）
gsea_all_long <- gsea_all_long %>%
  dplyr::mutate(
    signature_label_clean = label
  )

padj_thr <- 0.1

## ========= 用 NES 构建 Hallmark × Signature 矩阵，并做层次聚类 =========

hallmark_all <- gsea_all_long  ## 全部都是 Hallmark

nes_mat <- hallmark_all %>%
  dplyr::select(pathway, signature_label_clean, NES) %>%
  dplyr::distinct() %>%
  tidyr::pivot_wider(
    names_from  = signature_label_clean,
    values_from = NES
  )

path_vec <- nes_mat$pathway
nes_mat  <- as.data.frame(nes_mat)
rownames(nes_mat) <- path_vec
nes_mat$pathway <- NULL

## NA 用 0 填充，便于做距离
nes_mat[is.na(nes_mat)] <- 0

## 行聚类（Hallmark 通路）
if (nrow(nes_mat) >= 2) {
  row_ord    <- hclust(dist(nes_mat))$order
  row_levels <- rownames(nes_mat)[row_ord]
} else {
  row_levels <- rownames(nes_mat)
}

## 列聚类（6 个 signature）
if (ncol(nes_mat) >= 2) {
  col_ord    <- hclust(dist(t(nes_mat)))$order
  col_levels <- colnames(nes_mat)[col_ord]
} else {
  col_levels <- colnames(nes_mat)
}

## 通路名字清理（去掉 HALLMARK_，下划线变空格）
pathway_levels_clean <- row_levels %>%
  stringr::str_remove("^HALLMARK_") %>%
  stringr::str_replace_all("_", " ")

## ========= 构造显著条目的 long 表 =========

hallmark_long_sig <- hallmark_all %>%
  dplyr::filter(padj < padj_thr) %>%
  dplyr::mutate(
    pathway_clean = stringr::str_remove(pathway, "^HALLMARK_") %>%
      stringr::str_replace_all("_", " "),
    pathway_clean = factor(pathway_clean,
                           levels = rev(pathway_levels_clean)),
    signature_label_clean = factor(signature_label_clean,
                                   levels = col_levels)
  )

if (nrow(hallmark_long_sig) == 0) {
  warning("在 Hallmark 中没有 padj <", padj_thr, " 的显著通路，无法绘制 bubble plot。")
} else {

  pdf_file <- file.path(out_dir, "Hallmark_bubble_clustered_A2_Old_Heter_only_withPDF.pdf")

  pdf(pdf_file, width = 5.5, height = 9)
  print(
    ggplot(hallmark_long_sig,
           aes(x = signature_label_clean,
               y = pathway_clean,
               size = -log10(padj),
               color = NES)) +
      geom_point(alpha = 0.8) +
      scale_color_gradient2(
        low      = "blue",
        mid      = "white",
        high     = "red",
        midpoint = 0
      ) +
      scale_size_continuous(name = "-log10(adj. p)", range = c(1, 4.5)) +
      labs(
        x = "Interventions / signatures",
        y = "Hallmark gene sets",
        title = "Hallmark enrichment (padj < 0.1)\nXL A2 Old_Heter vs WT across tissues/sex"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
        axis.text.y  = element_text(color = "black"),
        axis.title.x = element_text(color = "black"),
        axis.title.y = element_text(color = "black"),
        plot.title   = element_text(hjust = 0.5)
      )
  )
  dev.off()

  message("已写出 Hallmark bubble plot (A2 Old_Heter only): ", pdf_file)
}
