

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
})

## ===== 你已读入并拆分好的数据帧 =====
## XL, WT 已在环境中；若没有，请取消注释下两行并改路径
XL <- read.table("/path/to/C03_Old_XL_vs_Old_WT_sig_q0.05.tsv", header=TRUE, sep="\t", check.names=FALSE)

XL_gain <- subset(XL, Fold > 0)
XL_loss <- subset(XL, Fold < 0)


## 背景 BED（mm39；无表头更稳，取前三列）
bg_df <- read.table("/Users/biobear/Library/CloudStorage/OneDrive-Personal/Transfer/cuttag_peaks/DiffBind_background_windows_mm39.bed",
                    header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(bg_df)[1:6] <- c("seq", "start0", "end", "name", "score", "strand")  # 可能只有前三列有用

## ===== 小工具：把 data.frame → GRanges（无 chr → UCSC +chr） =====
to_ucsc_chr <- function(x) {
  x <- as.character(x)
  x <- ifelse(grepl("^chr", x, ignore.case = FALSE), x, paste0("chr", x))
  ## 常见别名修正：chrMT → chrM
  x[x == "chrMT"] <- "chrM"
  x
}
df_to_gr_ucsc <- function(df, seq_col="seqnames", start_col="start", end_col="end") {
  stopifnot(all(c(seq_col, start_col, end_col) %in% colnames(df)))
  chr <- to_ucsc_chr(df[[seq_col]])
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = as.integer(df[[start_col]]),
                                                 end   = as.integer(df[[end_col]])))
  ## 只留标准染色体（chr1-19, chrX, chrY, chrM）
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr
}
bed3_to_gr_ucsc <- function(bed_df) {
  ## BED: 0-based, half-open; start = V2+1, end = V3
  chr <- to_ucsc_chr(bed_df[[1]])
  st  <- as.integer(bed_df[[2]]) + 1L
  en  <- as.integer(bed_df[[3]])
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = st, end = en))
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr
}

## ===== 1) 四个集合转 GRanges（UCSC 风格）=====
gr_XL_gain <- df_to_gr_ucsc(XL_gain)
gr_XL_loss <- df_to_gr_ucsc(XL_loss)

## ===== 2) 背景 BED → GRanges（UCSC 风格；仅前三列）=====
gr_bg <- bed3_to_gr_ucsc(bg_df)



res <- great(gr_XL_gain, "GO:BP", "mm39",background = gr_bg)
tb = getEnrichmentTable(res)
head(tb)
png("XL_gain_GO_TSS.png",width = 10, height = 4,units = "in",res = 1200)
plotRegionGeneAssociations(res)
dev.off()

## 1) 先选出显著条目（以 GREAT 的 binomial FDR 为主）
sig_tb <- subset(tb, fold_enrichment >= 2 & p_adjust_hyper <= 0.05 & p_adjust <= 0.05 & observed_gene_hits >= 3)

## 2) 计算 GO:BP 语义相似矩阵（用 "Rel" 或 "Wang" 都行）
sim <- rrvgo::calculateSimMatrix(
  sig_tb$id, orgdb = "org.Mm.eg.db",
  ont = "BP", method = "Rel"   # 或 method="Wang"
)

## 3) 用 -log10(FDR) 作为权重，聚合相似 term（threshold 越低越严格）
scores <- setNames(-log10(sig_tb$p_adjust + .Machine$double.xmin), sig_tb$id)
red <- rrvgo::reduceSimMatrix(
  simMatrix = sim,
  scores    = scores,
  threshold = 0.7,        # 0.7~0.9 常用；越低保留越少
  orgdb     = "org.Mm.eg.db"
)

## 4) 得到代表性条目表；保留你关心的列并导出
keep_ids <- red$parent      # 代表 term 的 GO ID
tb_rep <- sig_tb %>%
  filter(id %in% keep_ids) %>%
  arrange(p_adjust)

tb_use <- tb_rep %>%
  mutate(
    mlog10      = -log10(p_value  + .Machine$double.xmin),  # -log10(p)
    mlog10_adj  = -log10(p_adjust + .Machine$double.xmin),  # -log10(FDR)
    desc_clean  = str_replace_all(description, "_", " "),
    desc_title  = str_to_sentence(desc_clean),              # 首字母大写（句首）
    desc_wrap   = str_wrap(desc_title, width = 40)          # 自动换行
  ) %>%
  arrange(p_adjust, desc(fold_enrichment))

topN   <- 15
tb_top <- tb_use %>% slice_head(n = topN)

## 画 GO dotplot
pdf("XL_gain_GO_GREAT_dotplot_top15_rrvgo.pdf", width = 5, height = 5)
ggplot(
  tb_top,
  aes(
    x = fold_enrichment,
    y = reorder(desc_wrap, fold_enrichment)
  )
) +
  geom_point(
    aes(color = mlog10_adj),  # 颜色表示 -log10(FDR)
    size  = 1.8,
    alpha = 0.9
  ) +
  scale_color_gradient(
    name = "-log10(FDR)",
    low  = "steelblue3",
    high = "firebrick3"
  ) +
  labs(
    x     = "Fold enrichment (binomial)",
    y     = NULL,
    title = "Top enriched GO:BP terms (GREAT + rrvgo)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(color = "black", size = 9),
    axis.text.x        = element_text(color = "black"),
    axis.title.x       = element_text(color = "black"),
    axis.title.y       = element_text(color = "black"),
    plot.title         = element_text(color = "black", face = "bold"),
    panel.border       = element_rect(color = "black", linewidth = 1),
    axis.ticks         = element_line(color = "black", linewidth = 1)
  )
dev.off()





res <- great(gr_XL_loss, "GO:BP", "mm39",background = gr_bg)
tb = getEnrichmentTable(res)
head(tb)
png("XL_loss_GO_TSS.png",width = 10, height = 4,units = "in",res = 1200)
plotRegionGeneAssociations(res)
dev.off()

## 1) 先选出显著条目（以 GREAT 的 binomial FDR 为主）
sig_tb <- subset(tb, fold_enrichment >= 2 & p_adjust_hyper <= 0.05 & p_adjust <= 0.05 & observed_gene_hits >= 3)


sim <- rrvgo::calculateSimMatrix(
  sig_tb$id, orgdb = "org.Mm.eg.db",
  ont = "BP", method = "Rel"   
)

## 3) 用 -log10(FDR) 作为权重，聚合相似 term（threshold 越低越严格）
scores <- setNames(-log10(sig_tb$p_adjust + .Machine$double.xmin), sig_tb$id)
red <- rrvgo::reduceSimMatrix(
  simMatrix = sim,
  scores    = scores,
  threshold = 0.7,        # 0.7~0.9 常用；越低保留越少
  orgdb     = "org.Mm.eg.db"
)

## 4) 得到代表性条目表；保留你关心的列并导出
keep_ids <- red$parent      # 代表 term 的 GO ID
tb_rep <- sig_tb %>%
  filter(id %in% keep_ids) %>%
  arrange(p_adjust)


tb_use <- tb_rep %>%
  mutate(
    mlog10      = -log10(p_value  + .Machine$double.xmin),  # -log10(p)
    mlog10_adj  = -log10(p_adjust + .Machine$double.xmin),  # -log10(FDR)
    desc_clean  = str_replace_all(description, "_", " "),
    desc_title  = str_to_sentence(desc_clean),              # 首字母大写（句首）
    desc_wrap   = str_wrap(desc_title, width = 40)          # 自动换行
  ) %>%
  arrange(p_adjust, desc(fold_enrichment))

topN   <- 15
tb_top <- tb_use %>% slice_head(n = topN)

## 画 GO dotplot
pdf("XL_loss_GO_GREAT_dotplot_top15_rrvgo.pdf", width = 5, height = 5)
ggplot(
  tb_top,
  aes(
    x = fold_enrichment,
    y = reorder(desc_wrap, mlog10_adj)
  )
) +
  geom_point(
    aes(color = mlog10_adj),  # 颜色表示 -log10(FDR)
    size  = 1.8,
    alpha = 0.9
  ) +
  scale_color_gradient(
    name = "-log10(FDR)",
    low  = "steelblue3",
    high = "firebrick3"
  ) +
  labs(
    x     = "Fold enrichment (binomial)",
    y     = NULL,
    title = "Top enriched GO:BP terms (GREAT + rrvgo)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(color = "black", size = 9),
    axis.text.x        = element_text(color = "black"),
    axis.title.x       = element_text(color = "black"),
    axis.title.y       = element_text(color = "black"),
    plot.title         = element_text(color = "black", face = "bold"),
    panel.border       = element_rect(color = "black", linewidth = 1),
    axis.ticks         = element_line(color = "black", linewidth = 1)
  )
dev.off()
