## =====================================
## Fisher's exact test enrichment for all interventions
## 使用 DESeq2 结果 + Hallmark/Reactome/GO BP 基因集
## 分别对“每个 signature 的 top200 上调 / top200 下调基因”做 Fisher
## =====================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(msigdbr)
  library(ggplot2)
})

## ========== 0) 路径配置 ==========

## 你的 XL DESeq2 结果根目录
xl_root <- "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/revised/XL/deseq2_by_subset"

## nmrHas2（GSE234563）DESeq2 结果根目录
nmr_root <- "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/decoy_rnaseq/all_quantification/deseq2_GSE234563_by_subset"

## GSE131754 多干预 DESeq2 结果根目录
g131_root <- "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/decoy_rnaseq/all_quantification/deseq2_GSE131754_by_subset"

## 输出目录（Fisher 结果 + 图）
out_dir_fisher <- file.path(
  "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/decoy_rnaseq/all_quantification",
  "Fisher_all_collections"
)
dir.create(out_dir_fisher, showWarnings = FALSE, recursive = TRUE)

## ========== 1) 工具函数：读入 DESeq2 结果 ==========

load_deseq_result <- function(file, signature_id, dataset, extra_meta = list()) {
  if (!file.exists(file)) {
    stop("DESeq2 结果文件不存在: ", file)
  }
  df <- read.table(
    file,
    header      = TRUE,
    sep         = "\t",
    quote       = "",
    check.names = FALSE
  )
  need_cols <- c("gene", "log2FoldChange", "pvalue", "padj")
  miss <- setdiff(need_cols, colnames(df))
  if (length(miss) > 0) {
    stop("文件缺少必要列 [", paste(miss, collapse = ","), "]: ", file)
  }

  df <- df %>%
    mutate(
      ranking_p = case_when(
        !is.na(padj)   ~ padj,
        !is.na(pvalue) ~ pvalue,
        TRUE           ~ NA_real_
      )
    ) %>%
    filter(
      !is.na(gene),
      !is.na(log2FoldChange),
      !is.na(ranking_p)
    )

  df$signature_id <- signature_id
  df$dataset      <- dataset
  for (nm in names(extra_meta)) {
    df[[nm]] <- extra_meta[[nm]]
  }
  df
}

## 解析 GSE131754 文件名，例如:
##   Acarbose_12m_female.tsv
##   Acarbose_12m_female_Acarbose_vs_Control.tsv
parse_g131_meta_from_file <- function(file) {
  base <- tools::file_path_sans_ext(basename(file))
  parts <- strsplit(base, "_")[[1]]
  age_label       <- NA_character_
  condition_short <- NA_character_
  sex             <- NA_character_

  if (length(parts) >= 3) {
    condition_short <- parts[1]
    age_label       <- parts[2]   # "6m", "12m", "14m", "5m" ...
    sex             <- parts[3]   # "female" / "male"
  }
  list(
    age_label       = age_label,
    condition_short = condition_short,
    sex             = sex
  )
}

## ========== 2) 读取 XL + nmrHas2 的 DE 结果 ==========

## XL：Liver, male/female, Young/Old, Heter vs WT
xl_files <- tibble(
  dataset         = "XL",
  tissue          = "Liver",
  sex             = c("male", "male", "female", "female"),
  age_label       = c("Young", "Old", "Young", "Old"),
  condition_short = "Heter",
  file = c(
    file.path(xl_root, "liver_male",   "A1_Young_Heter_vs_WT.tsv"),
    file.path(xl_root, "liver_male",   "A2_Old_Heter_vs_WT.tsv"),
    file.path(xl_root, "liver_female", "A1_Young_Heter_vs_WT.tsv"),
    file.path(xl_root, "liver_female", "A2_Old_Heter_vs_WT.tsv")
  )
)

## nmrHas2：Liver, male/female, Young/Old, nmrHAS2 vs CreER
nmr_files <- tibble(
  dataset         = "nmrHas2",
  tissue          = "Liver",
  sex             = c("male", "male", "female", "female"),
  age_label       = c("Young", "Old", "Young", "Old"),
  condition_short = "nmrHAS2",
  file = c(
    file.path(nmr_root, "liver_male",   "A1_Young_nmrHAS2_vs_CreER.tsv"),
    file.path(nmr_root, "liver_male",   "A2_Old_nmrHAS2_vs_CreER.tsv"),
    file.path(nmr_root, "liver_female", "A1_Young_nmrHAS2_vs_CreER.tsv"),
    file.path(nmr_root, "liver_female", "A2_Old_nmrHAS2_vs_CreER.tsv")
  )
)

## 注意：signature_id 不包含 dataset
all_sig_defs <- bind_rows(xl_files, nmr_files) %>%
  mutate(
    signature_id = paste(tissue, sex, age_label, condition_short, sep = "_")
  )

de_list <- list()

for (i in seq_len(nrow(all_sig_defs))) {
  rowi <- all_sig_defs[i, ]
  message("读取: ", rowi$signature_id, " <- ", rowi$file)

  extra_meta <- list(
    tissue          = rowi$tissue,
    sex             = rowi$sex,
    age_label       = rowi$age_label,
    condition_short = rowi$condition_short
  )

  de_list[[rowi$signature_id]] <- load_deseq_result(
    file         = rowi$file,
    signature_id = rowi$signature_id,
    dataset      = rowi$dataset,
    extra_meta   = extra_meta
  )
}

## ========== 3) 读取 GSE131754 所有干预的 DE 结果 ==========

g131_files_all <- list.files(
  g131_root,
  pattern    = "\\.tsv$",
  recursive  = TRUE,
  full.names = TRUE
)

g131_files <- g131_files_all[
  basename(dirname(g131_files_all)) != basename(g131_root)
]

if (length(g131_files) == 0) {
  stop("没有在 ", g131_root, " 子目录中找到 .tsv 结果文件，请确认 GSE131754 的 DE 已经跑完。")
}

g131_defs <- tibble(
  dataset = "GSE131754",
  file    = g131_files
) %>%
  rowwise() %>%
  mutate(
    meta            = list(parse_g131_meta_from_file(file)),
    condition_short = meta$condition_short,
    age_label       = meta$age_label,
    sex             = meta$sex,
    tissue          = "Liver",
    ## signature_id 不包含 dataset
    signature_id    = paste(condition_short, age_label, sex, sep = "_")
  ) %>%
  ungroup() %>%
  select(-meta)

for (i in seq_len(nrow(g131_defs))) {
  rowi <- g131_defs[i, ]
  message("读取: ", rowi$signature_id, " <- ", rowi$file)

  extra_meta <- list(
    tissue          = rowi$tissue,
    sex             = rowi$sex,
    age_label       = rowi$age_label,
    condition_short = rowi$condition_short
  )

  sig_id <- rowi$signature_id
  if (!is.null(de_list[[sig_id]])) {
    sig_id <- paste0(rowi$signature_id, "_", i)
  }

  de_list[[sig_id]] <- load_deseq_result(
    file         = rowi$file,
    signature_id = sig_id,
    dataset      = rowi$dataset,
    extra_meta   = extra_meta
  )
}

## ========== 4) 合并所有 DE 结果并构建 meta 表 ==========

de_all <- bind_rows(de_list)

sig_meta <- de_all %>%
  distinct(
    signature_id, dataset, tissue, sex, age_label, condition_short
  ) %>%
  arrange(dataset, condition_short, age_label, sex)

## ========== 5) 按规则筛选签名 ==========
## - XL: 全部保留
## - nmrHas2: 只保留 Old
## - GSE131754: 每个干预只保留最大年龄

sig_meta_xl <- sig_meta %>%
  filter(dataset == "XL")

sig_meta_nmr <- sig_meta %>%
  filter(dataset == "nmrHas2", age_label == "Old")

sig_meta_g131 <- sig_meta %>%
  filter(dataset == "GSE131754") %>%
  mutate(
    age_numeric = as.numeric(stringr::str_extract(age_label, "[0-9]+"))
  ) %>%
  group_by(condition_short) %>%
  mutate(max_age = max(age_numeric, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(age_numeric), age_numeric == max_age) %>%
  select(-age_numeric, -max_age)

sig_meta_filtered <- bind_rows(sig_meta_xl, sig_meta_nmr, sig_meta_g131) %>%
  arrange(dataset, condition_short, age_label, sex) %>%
  mutate(
    signature_label = case_when(
      dataset == "XL" ~ paste("XL", sex, age_label, condition_short, sep = "_"),
      dataset == "nmrHas2" ~ paste("nmrHas2", sex, age_label, sep = "_"),
      dataset == "GSE131754" ~ paste(condition_short, age_label, sex, sep = "_"),
      TRUE ~ signature_id
    )
  )

sig_list_file <- file.path(out_dir_fisher, "signatures_used_for_Fisher.tsv")
write.table(
  sig_meta_filtered,
  file      = sig_list_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
message("本次 Fisher 分析使用的签名列表: ", sig_list_file)

de_all_filt <- de_all %>%
  filter(signature_id %in% sig_meta_filtered$signature_id)

## ========== 6) 准备 MSigDB 基因集（Hallmark / Reactome / GO BP） ==========

msig_h <- msigdbr(species = "Mus musculus", collection = "H")
msig_re <- msigdbr(
  species    = "Mus musculus",
  collection = "C2",
  subcategory = "CP:REACTOME"
)
msig_gobp <- msigdbr(
  species    = "Mus musculus",
  collection = "C5",
  subcategory = "BP"
)

msig_all <- bind_rows(
  msig_h %>% mutate(collection = "Hallmark"),
  msig_re %>% mutate(collection = "Reactome"),
  msig_gobp %>% mutate(collection = "GOBP")
)

## 每条通路对应的基因列表
pathway_list <- msig_all %>%
  group_by(collection, gs_name) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop")

## 通路使用的基因 universe：出现在 DE 中且在 MSigDB 里的基因
universe_genes <- intersect(
  unique(de_all_filt$gene),
  unique(msig_all$gene_symbol)
)

## ========== 7) 对每个 signature × 每条通路做 Fisher's exact test（Up / Down 分开）
##      使用“每个 signature 的 top200 上调 / top200 下调基因”
## ==========

padj_thr_sig <- 0.05
topN_per_direction <- 200

de_split <- split(de_all_filt, de_all_filt$signature_id)

fisher_res_list <- list()

for (i in seq_len(nrow(sig_meta_filtered))) {
  rowi <- sig_meta_filtered[i, ]
  sid  <- rowi$signature_id
  lab  <- rowi$signature_label

  df_sig <- de_split[[sid]]
  if (is.null(df_sig)) next

  ## 限制在 universe_genes，保证和背景一致
  df_sig <- df_sig %>%
    filter(gene %in% universe_genes) %>%
    filter(!is.na(log2FoldChange))

  ## ===== 选 top200 上调 / 下调基因 =====
  ## 如果你想“完全不看 padj，只按 log2FC 排序取前 200”，可以把下面 filter(!is.na(padj), padj < padj_thr_sig, ...) 改成只看 log2FoldChange 的条件。

  up_genes <- df_sig %>%
    filter(!is.na(padj),
           padj < padj_thr_sig,
           log2FoldChange > 0) %>%
    arrange(desc(log2FoldChange)) %>%
    pull(gene) %>%
    unique() %>%
    head(topN_per_direction)

  down_genes <- df_sig %>%
    filter(!is.na(padj),
           padj < padj_thr_sig,
           log2FoldChange < 0) %>%
    arrange(log2FoldChange) %>%   ## 越负越靠前
    pull(gene) %>%
    unique() %>%
    head(topN_per_direction)

  ## 两个方向分别做 Fisher
  for (direction in c("Up", "Down")) {
    sig_genes_dir <- if (direction == "Up") up_genes else down_genes

    if (length(sig_genes_dir) == 0) {
      next
    }

    is_sig_vec <- universe_genes %in% sig_genes_dir

    for (j in seq_len(nrow(pathway_list))) {
      rowp   <- pathway_list[j, ]
      coll   <- rowp$collection
      pname  <- rowp$gs_name
      genes_p <- rowp$genes[[1]]

      in_path <- universe_genes %in% genes_p

      a <- sum(is_sig_vec & in_path)
      b <- sum(is_sig_vec & !in_path)
      c <- sum(!is_sig_vec & in_path)
      d <- sum(!is_sig_vec & !in_path)

      if (a + b + c + d == 0 || a + c == 0 || b + d == 0) next

      mat <- matrix(c(a, c, b, d), nrow = 2)
      ft  <- fisher.test(mat)

      fisher_res_list[[length(fisher_res_list) + 1]] <- tibble(
        signature_id    = sid,
        signature_label = lab,
        dataset         = rowi$dataset,
        condition_short = rowi$condition_short,
        age_label       = rowi$age_label,
        sex             = rowi$sex,
        collection      = coll,
        pathway         = pname,
        direction       = direction,   ## Up / Down
        a = a, b = b, c = c, d = d,
        odds_ratio = unname(ft$estimate),
        pval       = ft$p.value
      )
    }
  }
}

fisher_res <- bind_rows(fisher_res_list)

## 组内做 BH 矫正：每个 signature × collection × direction 一次
fisher_res <- fisher_res %>%
  group_by(signature_label, collection, direction) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  mutate(
    log2OR = log2(odds_ratio)
  )

fisher_tsv <- file.path(out_dir_fisher, "Fisher_all_collections_all_signatures_UpDown_top200.tsv")
write.table(
  fisher_res,
  file      = fisher_tsv,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
message("已写出 Fisher 结果(Up/Down, top200): ", fisher_tsv)






## ========== 8) 工具函数：将 Up/Down 合并，并对处理 & 通路做 clustering ==========
## 规则：
##   - Fisher 结果仍然是 Up / Down 分开做；
##   - 对每个 signature × pathway，取 padj 最小的那条记录（Up/Down 二选一）；
##   - 如果来自 Down，则 signed_log2OR = -log2OR，
##     这样 Up 来源为正值，Down 来源为负值；
##   - 用 signed_log2OR 的矩阵对 pathway (行) 和 signature (列) 分别做层次聚类，
##     得到行列的排序，作为 bubble plot 的轴顺序。

prepare_combined_df <- function(collection_name,
                                fisher_res,
                                sig_meta_filtered,
                                padj_path_thr  = 0.1,
                                padj_point_thr = 0.05) {
  df <- fisher_res %>%
    dplyr::filter(collection == collection_name,
                  !is.na(padj))

  if (nrow(df) == 0) {
    message("【", collection_name, "】在 fisher_res 中没有结果。")
    return(NULL)
  }

  ## 对于每个 signature × pathway，选 padj 最小的那一条（Up / Down 二选一）
  combined <- df %>%
    dplyr::group_by(signature_id, signature_label, pathway) %>%
    dplyr::slice_min(padj, with_ties = FALSE) %>%
    dplyr::ungroup()

  ## 至少在一个 signature 中 padj < padj_path_thr 的通路才考虑
  sig_paths <- combined %>%
    dplyr::filter(padj < padj_path_thr) %>%
    dplyr::pull(pathway) %>%
    unique()

  plot_df <- combined %>%
    dplyr::filter(pathway %in% sig_paths,
                  padj <= padj_point_thr)

  if (nrow(plot_df) == 0) {
    message("【", collection_name, "】在合并后没有 padj <= ", padj_point_thr, " 的条目。")
    return(NULL)
  }

  ## 整理通路名字
  plot_df <- plot_df %>%
    dplyr::mutate(pathway_clean = pathway)

  if (collection_name == "Hallmark") {
    plot_df$pathway_clean <- gsub("^HALLMARK_", "", plot_df$pathway_clean)
  } else if (collection_name == "Reactome") {
    plot_df$pathway_clean <- gsub("^REACTOME_", "", plot_df$pathway_clean)
  } else if (collection_name == "GOBP") {
    plot_df$pathway_clean <- gsub("^GOBP_", "", plot_df$pathway_clean)
    plot_df$pathway_clean <- gsub("^GO_",    "", plot_df$pathway_clean)
  }
  plot_df$pathway_clean <- gsub("_", " ", plot_df$pathway_clean)

  ## 星号表示显著性
  plot_df <- plot_df %>%
    dplyr::mutate(
      signif_star    = dplyr::case_when(
        padj <= 0.001 ~ "***",
        padj <= 0.01  ~ "**",
        padj <= 0.05  ~ "*",
        TRUE          ~ ""
      ),
      neg_log10_padj = -log10(padj)
    )

  ## Up/Down -> 带符号的 log2(OR)
  plot_df <- plot_df %>%
    dplyr::mutate(
      direction     = factor(direction, levels = c("Up", "Down")),
      signed_log2OR = dplyr::case_when(
        direction == "Up"   ~ log2OR,
        direction == "Down" ~ -log2OR,
        TRUE                ~ NA_real_
      )
    )

  ## ========== clustering：构造 (pathway × signature) 矩阵 ==========
  ## 用 signed_log2OR 做聚类，NA 当成 0 处理

  mat_wide <- plot_df %>%
    dplyr::select(pathway_clean, signature_label, signed_log2OR) %>%
    tidyr::pivot_wider(
      names_from  = signature_label,
      values_from = signed_log2OR
    )

  rownames(mat_wide) <- mat_wide$pathway_clean
  mat_wide$pathway_clean <- NULL

  mat <- as.matrix(mat_wide)
  if (nrow(mat) > 0 && ncol(mat) > 0) {
    mat[is.na(mat)] <- 0
  }

  ## 行 clustering：pathway
  if (nrow(mat) > 1) {
    d_row <- stats::dist(mat)
    hc_row <- stats::hclust(d_row)
    pathway_order <- hc_row$labels[hc_row$order]
  } else {
    pathway_order <- rownames(mat)
  }

  ## 列 clustering：signature
  if (ncol(mat) > 1) {
    d_col <- stats::dist(t(mat))
    hc_col <- stats::hclust(d_col)
    signature_order <- hc_col$labels[hc_col$order]
  } else {
    signature_order <- colnames(mat)
  }

  ## 让 plot_df 中的 factor 顺序按 clustering 结果排
  plot_df$pathway_clean <- factor(
    plot_df$pathway_clean,
    levels = pathway_order
  )

  plot_df$signature_label_clean <- factor(
    plot_df$signature_label,
    levels = signature_order
  )

  plot_df
}

plot_combined_bubble <- function(plot_df,
                                 title_prefix,
                                 pdf_file,
                                 pdf_width  = 8,
                                 pdf_height = 10) {
  ## 控制点大小
  max_size_cap <- stats::quantile(plot_df$neg_log10_padj,
                                  0.95,
                                  na.rm = TRUE)
  plot_df <- plot_df %>%
    dplyr::mutate(
      neg_log10_padj_capped = pmin(neg_log10_padj, max_size_cap)
    )

  ## 控制颜色范围（按 signed_log2OR 的绝对值）
  finite_vals  <- plot_df$signed_log2OR[is.finite(plot_df$signed_log2OR)]
  range_log2OR <- max(abs(finite_vals), na.rm = TRUE)
  if (!is.finite(range_log2OR) || range_log2OR == 0) {
    range_log2OR <- 1
  }

  p <- ggplot(
    plot_df,
    aes(
      x = signature_label_clean,
      y = pathway_clean
    )
  ) +
    geom_point(
      aes(size = neg_log10_padj_capped, fill = signed_log2OR),
      shape  = 21,
      colour = "black",
      alpha  = 0.8
    ) +
    geom_text(
      aes(label = signif_star),
      size  = 2.5,
      vjust = 0.5,
      hjust = 0.5
    ) +
    scale_size_continuous(
      name  = expression(-log[10]("padj")),
      range = c(1, 4),
      guide = guide_legend(order = 2)
    ) +
    scale_fill_gradient2(
      name     = "signed log2(OR)\n(Up > 0, Down < 0)",
      low      = "blue",
      mid      = "white",
      high     = "red",
      midpoint = 0,
      limits   = c(-range_log2OR, range_log2OR),
      na.value = "white",
      guide    = guide_colourbar(order = 1)
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y      = element_text(colour = "black"),
      axis.title.x     = element_blank(),
      axis.title.y     = element_blank(),
      panel.grid       = element_blank(),
      legend.position  = "right",
      strip.background = element_rect(colour = "black", fill = "grey95")
    ) +
    ggtitle(paste0(
      title_prefix,
      "\nUp- and down-regulated genes combined (signed log2(OR), clustered)"
    ))

  pdf(pdf_file, width = pdf_width, height = pdf_height)
  print(p)
  dev.off()

  message("Bubble plot 已保存: ", pdf_file)
  invisible(p)
}





## ========== 9) Hallmark / Reactome / GO BP 可视化：合并 Up/Down + clustering ==========

## Hallmark
hallmark_plot_df <- prepare_combined_df(
  collection_name   = "Hallmark",
  fisher_res        = fisher_res,
  sig_meta_filtered = sig_meta_filtered,
  padj_path_thr     = 0.1,
  padj_point_thr    = 0.05
)

if (!is.null(hallmark_plot_df)) {
  pdf_file_hallmark <- file.path(
    out_dir_fisher,
    "Hallmark_Fisher_bubble_all_signatures_signed_top200_clustered.pdf"
  )
  plot_combined_bubble(
    plot_df      = hallmark_plot_df,
    title_prefix = "Fisher's exact test enrichment (Hallmark)",
    pdf_file     = pdf_file_hallmark,
    pdf_width    = 6.5,
    pdf_height   = 6
  )
}


