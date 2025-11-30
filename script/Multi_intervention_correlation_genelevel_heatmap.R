## =====================================
## logFC Spearman 相关矩阵（按“最大年龄”筛选）
## - XL: Liver Young/Old Heter vs WT（male/female，全保留）
## - nmrHas2: 只保留 Old（male/female）
## - GSE131754: 每个干预只保留最大年龄的结果（该年龄下的 male/female 全保留）
##   👉 包括 17-alpha-estradiol（目录/文件名形如 E17a_6m_female.tsv）
## =====================================
#region

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(pheatmap)
  library(tools)
})

## ========== 0) 配置区：路径 & 输出目录 ==========

## 1) 你的 XL 数据
xl_root <- "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/revised/XL/deseq2_by_subset"

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

## 2) nmrHas2（GSE234563）数据  ⭐ 修正根路径：deseq2_GSE234563_by_subset
nmr_root <- "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/decoy_rnaseq/all_quantification/deseq2_GSE234563_by_subset"

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

## 3) GSE131754 多干预 DE 结果根目录
g131_root <- "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/decoy_rnaseq/all_quantification/deseq2_GSE131754_by_subset"
if (!dir.exists(g131_root)) {
  stop("找不到 GSE131754 DE 结果目录: ", g131_root)
}

## 4) 输出目录（这里我保留原名，你也可以改成 *_with_17aE2 以示区分）
out_dir_cor <- "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/decoy_rnaseq/all_quantification/signature_correlation_maxAge_with_17aE2"
dir.create(out_dir_cor, showWarnings = FALSE, recursive = TRUE)

## ========== 1) 工具函数 ==========

## 读取一个 DESeq2 结果并标准化
load_deseq_result <- function(file, signature_id, dataset, extra_meta = list()) {
  if (!file.exists(file)) {
    stop("结果文件不存在: ", file)
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
    stop("文件缺少必要列 ", paste(miss, collapse = ","), " : ", file)
  }

  df <- df %>%
    mutate(
      ranking_p = dplyr::case_when(
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

## 解析 GSE131754 文件名中的 meta 信息
## 例如:
##   Acarbose_12m_female.tsv
##   Acarbose_12m_female_Acarbose_vs_Control.tsv
##   E17a_6m_male.tsv  （17-alpha-estradiol）
parse_g131_meta_from_file <- function(file) {
  base <- tools::file_path_sans_ext(basename(file))
  parts <- strsplit(base, "_")[[1]]
  age_label       <- NA_character_
  condition_short <- NA_character_
  sex             <- NA_character_

  if (length(parts) >= 3) {
    condition_short <- parts[1]   # Acarbose / CR / Rapamycin / Protandim / GHRKO / Snell / MR / E17a 等
    age_label       <- parts[2]   # "6m", "12m", "14m", "5m" ...
    sex             <- parts[3]   # "female" / "male"
  }
  list(
    age_label       = age_label,
    condition_short = condition_short,
    sex             = sex
  )
}

## ========== 2) 读取 XL & nmrHas2 的 DE 结果 ==========

all_sig_defs <- dplyr::bind_rows(xl_files, nmr_files) %>%
  mutate(
    signature_id = paste(dataset, tissue, sex, age_label, condition_short, sep = "_")
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

## ========== 3) 扫描并读取 GSE131754 的 DE 结果（包括 E17a） ==========

g131_files_all <- list.files(
  g131_root,
  pattern    = "\\.tsv$",
  recursive  = TRUE,
  full.names = TRUE
)

## 排除根目录下的汇总文件，只保留子目录里的单个对比结果
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
    condition_short = meta$condition_short,  # Acarbose / CR / Rapamycin / Protandim / MR / GHRKO / Snell / E17a ...
    age_label       = meta$age_label,        # 6m / 12m / 14m / 5m ...
    sex             = meta$sex,             # female / male
    tissue          = "Liver",
    signature_id    = paste(dataset, condition_short, age_label, sex, sep = "_")
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
    ## 极少数情况下 signature_id 重名，加一个序号避免覆盖
    sig_id <- paste0(rowi$signature_id, "_", i)
  }

  de_list[[sig_id]] <- load_deseq_result(
    file         = rowi$file,
    signature_id = sig_id,
    dataset      = rowi$dataset,
    extra_meta   = extra_meta
  )
}

## ========== 4) 合并所有 DE 结果，构建 meta 表 ==========

de_all <- dplyr::bind_rows(de_list)

sig_meta <- de_all %>%
  distinct(
    signature_id, dataset, tissue, sex, age_label, condition_short
  ) %>%
  arrange(dataset, condition_short, age_label, sex)

## ========== 5) 按规则筛选签名 ==========
##  - XL: 全部保留
##  - nmrHas2: 只保留 Old
##  - GSE131754: 每个 condition_short 只保留最大年龄（该年龄下的 male/female 全保留）

## XL
sig_meta_xl <- sig_meta %>%
  filter(dataset == "XL")

## nmrHas2: 只保留 Old
sig_meta_nmr <- sig_meta %>%
  filter(dataset == "nmrHas2", age_label == "Old")

## GSE131754: 每个干预只保留最大年龄
## E.g.
##  Acarbose: 6m & 12m -> 只保留 12m (male/female)
##  Rapamycin: 6m & 12m -> 只保留 12m
##  Protandim: 只有 6m -> 保留 6m
##  17-alpha-estradiol (E17a): 只有 6m -> 保留 6m
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
    label = case_when(
      dataset == "XL"        ~ paste("XL", sex, age_label, condition_short, sep = "_"),
      dataset == "nmrHas2"   ~ paste("nmrHas2", sex, age_label, sep = "_"),
      dataset == "GSE131754" ~ paste(condition_short, age_label, sex, sep = "_"),
      TRUE                   ~ signature_id
    )
  )

## 把这次实际使用的签名列表写出，便于你检查
sig_list_file <- file.path(out_dir_cor, "signatures_used_maxAge_with_17aE2.tsv")
write.table(
  sig_meta_filtered,
  file      = sig_list_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
message("本次相关性分析使用的签名列表: ", sig_list_file)

## 用筛选后的 signatures 过滤 de_all
de_all_filt <- de_all %>%
  filter(signature_id %in% sig_meta_filtered$signature_id)

message("de_all_filt 维度: genes x signatures = ", nrow(de_all_filt), " 行，",
        length(unique(de_all_filt$signature_id)), " 个 signature")

## ========== 6) 为每个 signature 准备 logFC & top400 基因集合 ==========

top_n <- 400L

de_split <- split(de_all_filt, de_all_filt$signature_id)

lfc_list <- lapply(de_split, function(df) {
  stats::setNames(df$log2FoldChange, df$gene)
})

top_genes_list <- lapply(de_split, function(df) {
  df2 <- df[order(df$ranking_p, na.last = NA), ]
  head(df2$gene, top_n)
})

## ========== 7) 计算 pairwise Spearman 相关矩阵 ==========

sig_ids <- sig_meta_filtered$signature_id
names(sig_ids) <- sig_meta_filtered$label

n_sig <- length(sig_ids)
cor_mat <- matrix(
  NA_real_, n_sig, n_sig,
  dimnames = list(names(sig_ids), names(sig_ids))
)
p_mat <- matrix(
  NA_real_, n_sig, n_sig,
  dimnames = list(names(sig_ids), names(sig_ids))
)

for (i in seq_len(n_sig)) {
  for (j in i:n_sig) {
    id_i <- sig_ids[i]
    id_j <- sig_ids[j]

    label_i <- names(sig_ids)[i]
    label_j <- names(sig_ids)[j]

    if (i == j) {
      cor_mat[label_i, label_j] <- 1
      p_mat[label_i, label_j]   <- 0
      next
    }

    gi <- top_genes_list[[id_i]]
    gj <- top_genes_list[[id_j]]

    union_genes <- union(gi, gj)
    if (length(union_genes) < 10) next

    v_i <- lfc_list[[id_i]][union_genes]
    v_j <- lfc_list[[id_j]][union_genes]

    keep <- !is.na(v_i) & !is.na(v_j)
    if (sum(keep) < 10) next

    ct <- suppressWarnings(cor.test(v_i[keep], v_j[keep], method = "spearman"))

    cor_mat[label_i, label_j] <- ct$estimate
    cor_mat[label_j, label_i] <- ct$estimate
    p_mat[label_i, label_j]   <- ct$p.value
    p_mat[label_j, label_i]   <- ct$p.value
  }
}

## ========== 8) 利用星号显著性 + 红白蓝色板画热图 ==========

## 8.1 构造星号矩阵（只显示显著性）
star_mat <- matrix(
  "",
  nrow = nrow(p_mat),
  ncol = ncol(p_mat),
  dimnames = dimnames(p_mat)
)

for (i in seq_len(nrow(p_mat))) {
  for (j in seq_len(ncol(p_mat))) {
    p <- p_mat[i, j]
    if (is.na(p) || i == j) {
      next
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

## 8.2 红–白–蓝配色：蓝(负) -> 白(0) -> 红(正)，白色严格对应 0
col_fun <- colorRampPalette(c("blue", "white", "red"))

range_val <- max(abs(cor_mat), na.rm = TRUE)
if (!is.finite(range_val) || range_val == 0) {
  range_val <- 1  # 兜底
}
breaks <- seq(-range_val, range_val, length.out = 101)  # 对称

## 屏幕预览（无 annotation）
pheatmap(
  cor_mat,
  main             = "Spearman correlation of log2FC signatures\n(top 400 genes per pair union, max-age per intervention, +17aE2)",
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  treeheight_row   = 0,
  treeheight_col   = 0,
  border_color     = NA,
  color            = col_fun(100),
  breaks           = breaks,
  display_numbers  = star_mat,
  number_color     = "black",
  fontsize_number  = 6,
  annotation_row   = NA,
  annotation_col   = NA
)

## 保存 PDF
pdf_file <- file.path(out_dir_cor, "logFC_signature_spearman_correlation_heatmap_maxAge_with_17aE2.pdf")
pdf(pdf_file, width = 5.5, height = 5.5)
pheatmap(
  cor_mat,
  main             = "Spearman correlation of log2FC signatures\n(top 400 genes per pair union, max-age per intervention, +17aE2)",
  cluster_rows     = TRUE,
  cluster_cols     = TRUE,
  treeheight_row   = 0,
  treeheight_col   = 0,
  border_color     = NA,
  color            = col_fun(100),
  breaks           = breaks,
  display_numbers  = star_mat,
  number_color     = "black",
  fontsize_number  = 6,
  annotation_row   = NA,
  annotation_col   = NA
)
dev.off()

## 矩阵保存
cor_tsv <- file.path(out_dir_cor, "logFC_signature_spearman_correlation_matrix_maxAge_with_17aE2.tsv")
p_tsv   <- file.path(out_dir_cor, "logFC_signature_spearman_pvalues_matrix_maxAge_with_17aE2.tsv")

write.table(
  cor_mat,
  file      = cor_tsv,
  sep       = "\t",
  quote     = FALSE,
  row.names = TRUE
)
write.table(
  p_mat,
  file      = p_tsv,
  sep       = "\t",
  quote     = FALSE,
  row.names = TRUE
)

message("相关性热图 PDF: ", pdf_file)
message("相关性矩阵 TSV: ", cor_tsv)
message("P 值矩阵 TSV: ", p_tsv)

