
## 读入与预处理（去掉重复读入）
df <- read.table("/Users/biobear/Workspace/Coop/Zhangwn/mouse_model_rmaseq/rebuttal_letter/rrbs/rrbs_liver.txt",
                 header=TRUE, sep="\t", quote="", check.names=FALSE, stringsAsFactors=FALSE)
df$`Chronological age` <- as.numeric(df$`Chronological age`)
df$Age       <- factor(df$Age, levels=c("Young","Old"))
df$Tissue    <- ifelse(grepl("Liver", df$Sample), "Liver", "Muscle")
df$Condition <- factor(df$Condition, levels=c("WT","XL"))

## 仅保留三种主分析时钟（删除 Petkovich）
clocks <- list(
  Meer      = c(y="DNAm age (Meer multi-tissue)",        w="OverlapProp Meer multi-tissue"),
  Stubbs    = c(y="DNAm age (Stubbs multi-tissue)",      w="OverlapProp Stubbs multi-tissue"),
  Thompson  = c(y="DNAm age (Thompson multi-tissue EN)", w="OverlapProp Thompson multi-tissue EN")
)

## 1) 组织内回归得到 EAA（加权）
for(nm in names(clocks)){
  y <- clocks[[nm]]["y"]; w <- clocks[[nm]]["w"]; eaa_col <- paste0("EAA_", nm)
  df[[eaa_col]] <- NA_real_
  for(tt in c("Liver","Muscle")){
    idx <- df$Tissue==tt & is.finite(df[[y]]) & is.finite(df$`Chronological age`) & is.finite(df[[w]])
    if(sum(idx) >= 3){
      fit <- lm(as.formula(paste0("`", y, "` ~ `Chronological age`")), data=df[idx,], weights=df[[w]][idx])
      df[[eaa_col]][idx] <- resid(fit)
    }
  }
}

## Old + Liver：双侧 Welch t + 加权 Stouffer（双侧）
old_liver <- subset(df, Age=="Old" & Tissue=="Liver")

pvals_two <- c(); diffs <- c(); ws <- c()
for(nm in names(clocks)){
  eaa <- old_liver[[paste0("EAA_", nm)]]
  ok  <- is.finite(eaa) & !is.na(old_liver$Condition)
  if(sum(ok) >= 3){
    tt <- t.test(eaa[ok] ~ old_liver$Condition[ok], var.equal=TRUE)  # 双侧
    p  <- tt$p.value
    dWT <- mean(eaa[ok & old_liver$Condition=="WT"])
    dXL <- mean(eaa[ok & old_liver$Condition=="XL"])
    del <- dWT - dXL
    cat(sprintf("Old-Liver %s: WT-XL = %.2f, two-sided p = %.4f\n", nm, del, p))
    pvals_two <- c(pvals_two, p)
    diffs     <- c(diffs, del)
    ws        <- c(ws, mean(old_liver[[clocks[[nm]]["w"]]][ok], na.rm=TRUE))
  }
}

# 带符号 z 合并（双侧）
z_signed <- sign(diffs) * qnorm(1 - pvals_two/2)
Z_w <- sum(ws * z_signed) / sqrt(sum(ws^2))
p_stouffer_two <- 2 * (1 - pnorm(abs(Z_w)))
cat(sprintf("Old-Liver (3 clocks) combined (weighted Stouffer): two-sided p = %.4g\n", p_stouffer_two))

## 3) 斜率交互：DNAm age ~ ChronAge * Condition（双侧检验）
suppressPackageStartupMessages({library(sandwich); library(lmtest)})
get_stat_col <- function(ct) intersect(c("t value","z value"), colnames(ct))[1]
for(nm in names(clocks)){
  y <- clocks[[nm]]["y"]; sub <- subset(df, Tissue=="Liver" & is.finite(`Chronological age`) & is.finite(get(y)))
  fit <- lm(as.formula(paste0("`", y, "` ~ `Chronological age` * Condition")), data=sub)
  vc  <- vcovHC(fit, type="HC3")
  ct  <- coeftest(fit, vcov.=vc)
  term <- if("ConditionXL:`Chronological age`"%in%rownames(ct)) "ConditionXL:`Chronological age`" else "`Chronological age`:ConditionXL"
  stat_col <- get_stat_col(ct); stat <- ct[term, stat_col]
  ## 双侧 p（正态近似；如想用常规 t，自行改成 pt）
  p_two <- 2*pnorm(-abs(stat))
  cat(sprintf("Liver slope interaction (%s): beta=%.3f, two-sided p=%.4g\n", nm, coef(fit)[term], p_two))
}

## 4) 可视化
suppressPackageStartupMessages({library(ggplot2); library(tidyr)})

# 4a) Old-Liver 的 EAA（仅三钟）
png("rrbs_liver_eaa.png", width=12, height=6, units="in", res=1200)
dl <- subset(df, Tissue=="Liver" & Age=="Old")
keep_cols <- paste0("EAA_", names(clocks))
dl_long <- pivot_longer(dl, cols=all_of(keep_cols), names_to="Clock", values_to="EAA")
dl_long$Clock <- factor(gsub("^EAA_","",dl_long$Clock), levels=names(clocks))
ggplot(dl_long, aes(x=Condition, y=EAA, fill=Condition)) +
  geom_violin(width=0.9, alpha=0.15, color=NA, trim=FALSE) +
  geom_boxplot(width=0.4, outlier.shape=NA, alpha=0.6) +
  geom_jitter(width=0.08, size=2, alpha=0.8) +
  stat_summary(fun=mean, geom="point", shape=23, size=3, fill="white") +
  facet_wrap(~Clock, scales="free_y") +
  labs(title="Old-Liver: EAA by condition (Meer / Stubbs / Thompson)", y="EAA (residuals)", x=NULL) +
  theme_classic() + theme(legend.position="none")
dev.off()


