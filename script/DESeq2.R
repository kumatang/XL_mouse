suppressPackageStartupMessages({
  library(tximport);library(DESeq2);library(PCAtools);library(ggsci);library(limma);library(vegan);library(dplyr);library(tidyr);library(tibble);library(readr);library(ggplot2);library(variancePartition);library(edgeR);library(extrafont)
})

set.seed(1)



samplelist <- read.table("/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/XL_aging_metadata_final_sup.txt",header=TRUE, sep="\t", quote="", row.names=1)


samplelist <- subset(samplelist, rownames(samplelist) %in% selected_samples)
# samplelist <- subset(samplelist, Tissue == "Skin")

stopifnot("Batch" %in% colnames(samplelist))
samplelist$Batch     <- factor(samplelist$Batch)
samplelist$Condition <- factor(samplelist$Condition)
samplelist$Tissue    <- factor(samplelist$Tissue)
samplelist$Gender    <- factor(samplelist$Gender)
samplelist$Age       <- factor(samplelist$Age)
for (v in c("Date","Company")) if (v %in% colnames(samplelist)) samplelist[[v]] <- factor(samplelist[[v]])

# Salmon files
dir  <- "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/decoy_rnaseq/quantification/salmon_paired/"
files <- file.path(dir, samplelist$ID, "quant.sf"); names(files) <- rownames(samplelist)


# tx2gene
tx2gene_full <- read.table("/Users/biobear/Workspace/Coop/Wangxn/sirt2_sod1_mouse/tx2gene_GRCm39_full_modified_full.tsv",
                           header=FALSE, sep="\t", quote="")
colnames(tx2gene_full) <- c("tx_id","gene_id","gene_name","gene_type")
tx2gene <- tx2gene_full[, c("tx_id","gene_name")]
pc_gene <- unique(subset(tx2gene_full, gene_type == "protein_coding")$gene_name)

# 导入原始计数（不在此处全局过滤）
txi.salmon <- tximport(files, type="salmon", tx2gene=tx2gene)
exp <- txi.salmon$counts
exp <- exp[rownames(exp) != "G609G", , drop=FALSE]
exp_pc <- exp[rownames(exp) %in% pc_gene, , drop=FALSE]
exp_pc <- exp_pc[!rownames(exp_pc) %in% c("LIG4","XRCC4"), , drop=FALSE]

cts <- round(exp_pc)
coldata <- samplelist
stopifnot(all(rownames(coldata) == colnames(cts)))
coldata$LibSize    <- colSums(cts)
coldata$logLibSize <- log1p(coldata$LibSize)



#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

## ========== 用户指定路径（按你的要求） ==========
batch_tbl_file <- "/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/XL_aging_mouse/revised/XL/diagnostics_by_subset/batch_vars_to_consider.tsv"


suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(readr)
  library(dplyr)
  library(tidyr)
})

## ========== 前提检查 ==========
if (!exists("cts") || !exists("coldata")) {
  stop("缺少对象 cts / coldata。请先在本会话中运行你的前置脚本以构建它们。")
}
stopifnot(file.exists(batch_tbl_file))
batch_tbl <- readr::read_tsv(batch_tbl_file, show_col_types = FALSE)

out_root <- "deseq2_by_subset"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

## ========== 工具函数：满秩保障 & 过滤 ==========
# a 是否“决定” b：对 a 的每个水平，b 至多出现 1 个水平（嵌套/一对一）
det_by <- function(a, b) {
  a <- droplevels(factor(a)); b <- droplevels(factor(b))
  if (nlevels(a) < 2 || nlevels(b) < 2) return(FALSE)
  tb <- table(a, b)
  all(rowSums(tb > 0) <= 1)
}
is_confounded_either <- function(x, y) det_by(x, y) || det_by(y, x)

is_full_rank <- function(meta, rhs_terms) {
  X <- model.matrix(as.formula(paste("~", paste(rhs_terms, collapse = " + "))), data = meta)
  qr(X)$rank == ncol(X)
}

# 解析并去冗余批次因子，保证设计满秩；返回最终 terms 与 notes
resolve_batch_terms <- function(meta_sub, batch_terms,
                                keep_priority = c("Date","Company","Batch")) {
  notes <- c()
  terms <- batch_terms

  # 1) 去掉该子集内水平 < 2 的批次因子
  lv_ok <- sapply(terms, function(b) nlevels(droplevels(factor(meta_sub[[b]]))) >= 2)
  if (any(!lv_ok)) {
    notes <- c(notes, sprintf("drop{%s}: <2 levels", paste(terms[!lv_ok], collapse=",")))
    terms <- terms[lv_ok]
  }
  if (length(terms) == 0) return(list(terms=character(0), notes=notes))

  # 2) 与 Condition / Age / (Condition×Age) 完全共线则剔除
  CA <- interaction(meta_sub$Condition, meta_sub$Age, drop=TRUE)
  bad <- c()
  for (b in terms) {
    if (is_confounded_either(meta_sub[[b]], meta_sub$Condition)) bad <- c(bad, b)
    if (is_confounded_either(meta_sub[[b]], meta_sub$Age))       bad <- c(bad, b)
    if (is_confounded_either(meta_sub[[b]], CA))                 bad <- c(bad, b)
  }
  bad <- unique(bad)
  if (length(bad)) {
    notes <- c(notes, sprintf("drop{%s}: confounded with Condition/Age", paste(bad, collapse=",")))
    terms <- setdiff(terms, bad)
  }
  if (length(terms) == 0) return(list(terms=character(0), notes=notes))

  # 3) 批次因子之间的嵌套/决定：按优先级只保留一个
  if (length(terms) >= 2) {
    # 高优先级放前面
    ordered <- terms[order(match(terms, keep_priority, nomatch = length(keep_priority)+1))]
    keep <- c(); drop_nested <- c()
    for (bi in ordered) {
      if (bi %in% drop_nested) next
      # 丢掉被 bi 决定的其他批次
      for (bj in setdiff(terms, c(keep, drop_nested, bi))) {
        if (det_by(meta_sub[[bi]], meta_sub[[bj]])) drop_nested <- c(drop_nested, bj)
      }
      keep <- c(keep, bi)
    }
    if (length(drop_nested)) {
      notes <- c(notes, sprintf("drop{%s}: nested/redundant by higher-priority batch", paste(unique(drop_nested), collapse=",")))
      terms <- setdiff(terms, unique(drop_nested))
    }
  }

  # 4) 若仍不满秩，按优先级反向贪心删除直到满秩
  rhs0 <- c(terms, "Condition", "Age", "Condition:Age")
  if (!is_full_rank(meta_sub, rhs0)) {
    for (b in rev(keep_priority)) {
      if (!(b %in% terms)) next
      try_terms <- setdiff(terms, b)
      rhs_try <- c(try_terms, "Condition", "Age", "Condition:Age")
      if (is_full_rank(meta_sub, rhs_try)) {
        notes <- c(notes, sprintf("drop{%s}: ensure full rank", b))
        terms <- try_terms
      }
      if (is_full_rank(meta_sub, c(terms, "Condition","Age","Condition:Age"))) break
    }
  }
  list(terms=terms, notes=notes)
}

# 设计感知低表达过滤（广泛采用参数设置）
filter_by_expr_subset <- function(counts, meta) {
  y   <- edgeR::DGEList(counts = counts)
  grp <- interaction(meta$Condition, meta$Age, drop = TRUE)
  keep <- edgeR::filterByExpr(
    y, group = grp,
    min.count = 10,       # per-sample阈值（与CPM标尺/文库大小自适应）
    min.total.count = 15, # 全样本最小总计数
    large.n = 10,         # 组内样本数>large.n时允许用 min.prop
    min.prop = 0.7        # 大样本时至少在该比例样本中达到阈值
  )
  counts[keep, , drop = FALSE]
}

# 近似 CPM 阈值（报告用）：k ≈ 10 / 中位文库大小(百万)
calc_k_cpm <- function(counts_subset) {
  Lmed <- median(colSums(counts_subset))
  k    <- 10 / (Lmed / 1e6)
  c(Lmed = Lmed, k_cpm = k)
}

## ========== 主流程 ==========
summ_rows <- list()
batch_used_rows <- list()
filt_rows <- list()

tissues <- levels(coldata$Tissue); if (is.null(tissues)) tissues <- unique(coldata$Tissue)
genders <- levels(coldata$Gender); if (is.null(genders)) genders <- unique(coldata$Gender)

for (tissue in tissues) {
  for (gender in genders) {

    # 匹配 batch 表
    hit <- which(batch_tbl$Tissue == tissue & batch_tbl$Gender == gender)
    if (length(hit) != 1) {
      message(sprintf("跳过：无法唯一匹配 batch 记录 -> %s / %s", tissue, gender))
      next
    }
    used_batches_raw <- batch_tbl$batch_used_for_RBE[hit]
    used_batches <- if (is.na(used_batches_raw) || used_batches_raw=="" || used_batches_raw=="None") character(0) else strsplit(used_batches_raw, "\\+")[[1]] |> trimws()

    message(sprintf("\n===== DESeq2 | Tissue=%s | Gender=%s | batches={%s} =====",
                    tissue, gender, ifelse(length(used_batches)==0,"None", paste(used_batches, collapse=","))))

    # 子集
    idx <- (coldata$Tissue == tissue) & (coldata$Gender == gender)
    if (sum(idx) < 4) { message("  子集样本数过少，跳过"); next }

    meta_sub <- droplevels(coldata[idx, , drop = FALSE])
    cts_sub  <- cts[, rownames(meta_sub), drop = FALSE]

    # 仅保留 Age∈{Young,Old}, Condition∈{WT,Heter}
    keep_age  <- intersect(levels(droplevels(meta_sub$Age)),       c("Young","Old"))
    keep_cond <- intersect(levels(droplevels(meta_sub$Condition)), c("WT","Heter"))
    keep_mask <- meta_sub$Age %in% keep_age & meta_sub$Condition %in% keep_cond
    meta_sub  <- droplevels(meta_sub[keep_mask, , drop = FALSE])
    cts_sub   <- cts_sub[, rownames(meta_sub), drop = FALSE]

    if (nlevels(meta_sub$Age) < 2 || nlevels(meta_sub$Condition) < 2) {
      message("  Age 或 Condition 水平不足两类，跳过")
      next
    }

    # 固定参考水平，稳定系数命名
    meta_sub$Age       <- droplevels(factor(meta_sub$Age,       levels = c("Young","Old")))
    meta_sub$Condition <- droplevels(factor(meta_sub$Condition, levels = c("WT","Heter")))

    # 低表达过滤（广泛采用）+ 记录
    genes_before <- nrow(cts_sub)
    cts_sub_f <- filter_by_expr_subset(cts_sub, meta_sub)
    genes_after_fbe <- nrow(cts_sub_f)
    if (nrow(cts_sub_f) < 50) { message("  过滤后基因数 < 50，跳过"); next }

    kinfo <- calc_k_cpm(cts_sub)  # 报告用（以子集库中位数估近似阈值）

    # 解析批次因子，保证满秩
    res_bt <- resolve_batch_terms(meta_sub, used_batches,
                                  keep_priority = c("Date","Company","Batch"))
    batch_terms <- res_bt$terms
    if (length(res_bt$notes)) message("  batch resolution notes: ", paste(res_bt$notes, collapse = " | "))

    rhs <- c(batch_terms, "Condition", "Age", "Condition:Age")
    if (!is_full_rank(meta_sub, rhs)) stop("设计仍不满秩，请检查数据或规则。")

    design_fml <- as.formula(paste("~", paste(rhs, collapse = " + ")))
    message("  design(final): ", deparse(design_fml))

    # 拟合
    dds <- DESeqDataSetFromMatrix(countData = cts_sub_f,
                                  colData   = meta_sub,
                                  design    = design_fml)
    dds <- DESeq(dds)

    rn <- resultsNames(dds)
    message("  coefficients: ", paste(rn, collapse = ", "))

    coef_cond <- "Condition_Heter_vs_WT"
    coef_age  <- "Age_Old_vs_Young"
    coef_int  <- "ConditionHeter.AgeOld"
    need <- c(coef_cond, coef_age, coef_int)
    if (!all(need %in% rn)) { message("  系数命名与预期不一致，跳过"); next }

    # 四个比较（未经收缩；alpha=0.05 仅用于独立过滤优化与显著阈值提示，不会删行）
    res_A1 <- results(dds, name = coef_cond, alpha = 0.05)                               # Young 内 Heter vs WT
    res_A2 <- results(dds, contrast = list(c(coef_cond, coef_int)), alpha = 0.05)        # Old   内 Heter vs WT
    res_B1 <- results(dds, name = coef_age,  alpha = 0.05)                               # WT    内 Old vs Young
    res_B2 <- results(dds, contrast = list(c(coef_age,  coef_int)), alpha = 0.05)        # Heter 内 Old vs Young


    # 写结果
    subdir <- file.path(out_root, paste0(tolower(tissue), "_", tolower(gender)))
    dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
    write.table(as.data.frame(res_A1) |> tibble::rownames_to_column("gene"),
                file = file.path(subdir, "A1_Young_Heter_vs_WT.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(as.data.frame(res_A2) |> tibble::rownames_to_column("gene"),
                file = file.path(subdir, "A2_Old_Heter_vs_WT.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(as.data.frame(res_B1) |> tibble::rownames_to_column("gene"),
                file = file.path(subdir, "B1_WT_Old_vs_Young.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(as.data.frame(res_B2) |> tibble::rownames_to_column("gene"),
                file = file.path(subdir, "B2_Heter_Old_vs_Young.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    # 汇总统计
    count_tested <- function(res) sum(!is.na(res$padj))
    count_sig    <- function(res) sum(!is.na(res$padj) & res$padj < 0.05)

    summ_rows[[paste(tissue, gender, sep="_")]] <- tibble(
      Tissue  = tissue, Gender = gender,
      genes_before_filter      = genes_before,
      genes_after_filterByExpr = genes_after_fbe,
      A1_tested = count_tested(res_A1), A1_sig = count_sig(res_A1),
      A2_tested = count_tested(res_A2), A2_sig = count_sig(res_A2),
      B1_tested = count_tested(res_B1), B1_sig = count_sig(res_B1),
      B2_tested = count_tested(res_B2), B2_sig = count_sig(res_B2)
    )

    batch_used_rows[[paste(tissue, gender, sep="_")]] <- tibble(
      Tissue = tissue, Gender = gender,
      batches_used_in_design = ifelse(length(batch_terms)==0, "None", paste(batch_terms, collapse="+")),
      notes = ifelse(length(res_bt$notes)==0, "", paste(res_bt$notes, collapse="; "))
    )

    filt_rows[[paste(tissue, gender, sep="_")]] <- tibble(
      Tissue = tissue, Gender = gender,
      Lmed_reads = round(kinfo["Lmed"]),            # 文库中位数（reads）
      approx_k_cpm = signif(kinfo["k_cpm"], 3)      # 近似 CPM 阈值（报告用）
    )

    message(sprintf("  DONE | A1_sig=%d  A2_sig=%d  B1_sig=%d  B2_sig=%d (padj<0.05)",
                    count_sig(res_A1), count_sig(res_A2), count_sig(res_B1), count_sig(res_B2)))
  }
}

## ========== 写总表 ==========
if (length(summ_rows)) {
  summary_all <- dplyr::bind_rows(summ_rows)
  write.table(summary_all,
              file = file.path(out_root, "DE_summary_counts.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  print(summary_all)
}
if (length(batch_used_rows)) {
  batch_used_all <- dplyr::bind_rows(batch_used_rows)
  write.table(batch_used_all,
              file = file.path(out_root, "DE_batch_factors_used.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  print(batch_used_all)
}
if (length(filt_rows)) {
  filt_all <- dplyr::bind_rows(filt_rows)
  write.table(filt_all,
              file = file.path(out_root, "DE_filtering_qc.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  print(filt_all)
}

message("\nDESeq2 差异分析完成；结果见 deseq2_by_subset/*/")


# A1（Young 内 Heter vs WT）= name="Condition_Heter_vs_WT"
# A2（Old  内 Heter vs WT）  = contrast = list(c("Condition_Heter_vs_WT","ConditionHeter.AgeOld"))
# B1（WT   内 Old vs Young）= name="Age_Old_vs_Young"
# B2（Heter 内 Old vs Young）= contrast = list(c("Age_Old_vs_Young","ConditionHeter.AgeOld"))


## ====== 增量：统计每组上/下调（|log2FC|>1 & padj<=0.05） ======
thr_lfc  <- log2(1.5)
thr_padj <- 0.05

subdirs <- list.dirs(out_root, full.names = TRUE, recursive = FALSE)
stats_rows <- list()

for (sd in subdirs) {
  nm <- basename(sd)                # e.g., liver_female
  parts <- strsplit(nm, "_")[[1]]
  tissue <- parts[1]; gender <- parts[2]

  fls <- list.files(sd, pattern="\\.tsv$", full.names = TRUE)
  for (f in fls) {
    base <- tools::file_path_sans_ext(basename(f))  # e.g., A1_Young_Heter_vs_WT
    df <- tryCatch(read.delim(f, check.names = FALSE), error = function(e) NULL)
    if (is.null(df)) next
    if (!all(c("log2FoldChange","padj") %in% names(df))) next

    n_up   <- sum(!is.na(df$padj) & df$padj <= thr_padj & !is.na(df$log2FoldChange) & df$log2FoldChange >  thr_lfc)
    n_down <- sum(!is.na(df$padj) & df$padj <= thr_padj & !is.na(df$log2FoldChange) & df$log2FoldChange < -thr_lfc)

    # 统一对比标签
    if (grepl("^A1", base)) contrast <- "A1_Young_Heter_vs_WT"
    else if (grepl("^A2", base)) contrast <- "A2_Old_Heter_vs_WT"
    else if (grepl("^B1", base)) contrast <- "B1_WT_Old_vs_Young"
    else if (grepl("^B2", base)) contrast <- "B2_Heter_Old_vs_Young"
    else contrast <- base

    stats_rows[[length(stats_rows) + 1]] <- data.frame(
      Tissue   = tissue,
      Gender   = gender,
      Contrast = contrast,
      Up       = n_up,
      Down     = n_down,
      stringsAsFactors = FALSE
    )
  }
}

if (length(stats_rows)) {
  updown_counts <- dplyr::bind_rows(stats_rows) %>%
    dplyr::arrange(Tissue, Gender, Contrast)

  out_long <- file.path(out_root, "DE_updown_counts_LFC1_padj0.05.tsv")
  write.table(updown_counts, out_long, sep = "\t", quote = FALSE, row.names = FALSE)
  print(updown_counts)

  # 宽表（每个对比拆成 Up/Down 两列）
  updown_wide <- updown_counts %>%
    tidyr::pivot_longer(cols = c("Up","Down"), names_to = "Direction", values_to = "N") %>%
    tidyr::unite("Contrast_Dir", Contrast, Direction, sep = ":") %>%
    tidyr::pivot_wider(names_from = Contrast_Dir, values_from = N, values_fill = 0) %>%
    dplyr::arrange(Tissue, Gender)

  out_wide <- file.path(out_root, "DE_updown_counts_LFC1_padj0.05_wide.tsv")
  write.table(updown_wide, out_wide, sep = "\t", quote = FALSE, row.names = FALSE)
  print(updown_wide)
} else {
  message("未找到任何对比结果文件；请确认 deseq2_by_subset/*/*.tsv 是否已生成。")
}













