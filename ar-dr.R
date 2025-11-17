rm(list = ls())
library(dplyr)
library(data.table)
library(vroom)
library(gwasvcf)
library(gwasglue)
library(ieugwasr)
library(plinkbinr)
library(TwoSampleMR)
library(MRPRESSO)
library(LDlinkR)
library(MREILLS)


#DR
AD_exposure_finn = fread("GCST007320_GRCh37.tsv")
head(AD_exposure_finn)
colnames(AD_exposure_finn)
AD_exposure_finn <- subset(AD_exposure_finn,p_value< 5e-08)
dim(AD_exposure_finn)#[1] 826  13
AD_exposure_finn<-data.frame(AD_exposure_finn)
AD_exposure_finn = format_data(dat = AD_exposure_finn,
                               type = "exposure",
                               snp_col = "variant_id",
                               beta_col = "beta",
                               pval_col = "p_value",
                               se_col = "standard_error",
                               eaf_col = "effect_allele_frequency",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele")

write.csv(AD_exposure_finn,file = "Results/AD_exposure_finn_ALL.csv")
AD_exposure_dat <- AD_exposure_finn
#修改列名
#去除连锁不平衡
AD_exposure_dat_clumped = ieugwasr::ld_clump(dplyr::tibble(rsid = AD_exposure_dat$SNP,
                                                           pval = AD_exposure_dat$pval.exposure),
                                             clump_kb = 10000,
                                             clump_r2 = 0.001,
                                             clump_p = 1,
                                             bfile = "1kg.v3/EUR",
                                             plink_bin = plinkbinr::get_plink_exe(),
                                             pop = "EUR")
AD_exposure_dat_clumped$id.exposure = "AD"
AD_exposure_dat_clumped = AD_exposure_dat[AD_exposure_dat$SNP %in% AD_exposure_dat_clumped$rsid,]
#保留F值>10

# Calculate R2 and F stat for exposure data
# Liberal hypertension F stat
AD_exposure_dat_clumped$samplesize.exposure = 455258
AD_exposure_dat_clumped$r2 <- (2 * (AD_exposure_dat_clumped$beta.exposure^2) * AD_exposure_dat_clumped$eaf.exposure * (1 - AD_exposure_dat_clumped$eaf.exposure)) /
  (2 * (AD_exposure_dat_clumped$beta.exposure^2) * AD_exposure_dat_clumped$eaf.exposure * (1 - AD_exposure_dat_clumped$eaf.exposure) +
     2 * AD_exposure_dat_clumped$samplesize.exposure * AD_exposure_dat_clumped$eaf.exposure * (1 - AD_exposure_dat_clumped$eaf.exposure) * AD_exposure_dat_clumped$se.exposure^2)

AD_exposure_dat_clumped$F <- AD_exposure_dat_clumped$r2 * (AD_exposure_dat_clumped$samplesize.exposure - 2) / (1 - AD_exposure_dat_clumped$r2)

AD_exposure_dat_f <-AD_exposure_dat_clumped[AD_exposure_dat_clumped$F>10,]
dim(AD_exposure_dat_clumped)#24 14
write.csv(AD_exposure_dat_clumped,file = "Results/AD_exposure_dat_clumped_ALL.csv")
save(AD_exposure_dat_clumped,file = "Results/Rdata/AD.Rdata")


#结局----
outcome_finn = fread("GCST90479891.tsv.gz")
colnames(outcome_finn)
outcome_finn_dat = subset(outcome_finn, outcome_finn$rsid %in% AD_exposure_dat_clumped$SNP)
dim(outcome_finn_dat)# 24 13
outcome_finn_dat<-data.frame(outcome_finn_dat)
head(outcome_finn_dat)
outcome_finn_dat$standard_error = abs(log(outcome_finn_dat$odds_ratio) / qnorm(outcome_finn_dat$p_value / 2))
outcome_finn_dat$Z =  log(outcome_finn_dat$odds_ratio) / outcome_finn_dat$standard_error
outcome_finn_dat$beta <- log(outcome_finn_dat$odds_ratio)

outcome_finn_dat = format_data(dat = outcome_finn_dat,
                               type = "outcome",
                               snp_col = "rsid",
                               beta_col = "beta",
                               se_col = "standard_error",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               eaf_col = "effect_allele_frequency",
                               pval_col = "p_value",
                               samplesize_col = "n" )
outcome_finn_dat = subset(outcome_finn_dat, pval.outcome>5e-08)
dim(outcome_finn_dat)#[1] 24 12

#修改暴露和结局的名称
AD_exposure_dat_clumped$id.exposure = "AD"
outcome_finn_dat$id.outcome = "DR"

mr_data_finn = harmonise_data(exposure_dat = AD_exposure_dat_clumped,
                              outcome_dat = outcome_finn_dat,
                              action= 2)   #去除回文序列，2去除，1不去除
table(mr_data_finn$mr_keep)#1 23

mr_data_finn<-mr_data_finn[mr_data_finn$mr_keep==TRUE,]

write.csv(mr_data_finn, file = "Results/mr_data_ALL.csv")
View(mr_data_finn)
#mr_data_finn <- data.table::fread("data/AD/mr_DR_data_PGC.csv",data.table = F)

#MR----
#MR分析
res = mr(mr_data_finn)
res

#异质性检验
het = mr_heterogeneity(mr_data_finn)
het
#异质性可视化
sin = mr_singlesnp(mr_data_finn)
# pdf(file = "data/AD/yizhixing.pdf", width=10, heigh=8)
mr_funnel_plot(singlesnp_results = sin)
# dev.off()

#MR-PRESSO异常值检测
presso = mr_presso(BetaOutcome ="beta.outcome",
                   BetaExposure = "beta.exposure",
                   SdOutcome ="se.outcome",
                   SdExposure = "se.exposure",
                   OUTLIERtest = TRUE,
                   DISTORTIONtest = TRUE,
                   data = mr_data_finn,
                   NbDistribution = 3000,
                   SignifThreshold = 0.05)
presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices` #如果显示NULL，则表示不存在异常值

het_cleDR = mr_heterogeneity(mr_data_cleDR) #异质性检验
het_cleDR
write.csv(het_cleDR,file = "Results/het_cleDR_ALL.csv")
#异质性可视化
sin_cleDR = mr_singlesnp(mr_data_cleDR)
sin_cleDR
pdf(file = "Results/yizhixing.pdf", width=10, heigh=8)
mr_funnel_plot(singlesnp_results = sin_cleDR)
dev.off()


#多效性检验
pleio = mr_pleiotropy_test(mr_data_cleDR)
pleio
write.csv(pleio,file = "Results/pleio_ALL.csv")
#重新进行mr分析
res_cleDR = mr(mr_data_cleDR)
OR_cleDR = generate_odds_ratios(res_cleDR)

OR_cleDR$pval.fdr <- p.adjust(OR_cleDR$pval, method = "fdr")
OR_cleDR$pval.bonferroni <- p.adjust(OR_cleDR$pval, method = "bonferroni")
OR_cleDR
write.csv(OR_cleDR, "Results/MRresult_AD_DR.csv")

#MR分析结果可视化
pdf(file = "Results/plot/MR_ALL.pdf", width=10, heigh=8)
mr_scatter_plot(res_cleDR, mr_data_cleDR)
dev.off()

#MR分析结果森林图
col = c("id.exposure", "id.outcome", "nsnp", "method", "or", "pval", "or_lci95", "or_uci95")
OR = OR_cleDR
OR = OR[,col]
OR$`OR (95% CI)` = sprintf("%.2f (%.2f to %.2f)", OR$or, OR$or_lci95, OR$or_uci95)
OR$pval = signif(OR$pval, 3)
#write.csv(OR, "OR.csv", row.names = F)
#OR = read.csv("OR.csv", header = T, check.names = F)
#OR$nsnp = ifelse(is.na(OR$nsnp), "", OR$nsnp)
OR$` ` = paste(rep(" ", 40), collapse = " ")
#绘制
library(forestploter)
library(grid)
pdf("Results/plot/MRresult_eQTLGen_ALL.pdf", width = 11, height = 7, onefile = FALSE)
tm = forest_theme(base_size = 10,   # 文本的大小
                  ci_pch = 16,      # 可信区间点的形状
                  ci_col = "black", # CI的颜色
                  ci_fill = "red",  # CI中se点的颜色填充
                  ci_alpha = 0.8,   # CI透明度
                  ci_lty = 1,       # CI的线型
                  ci_lwd = 1.5,     # CI的线宽
                  ci_Theight = 0.2, # CI的高度，默认是NULL
                  #参考线默认的参数，中间的竖的虚线
                  refline_lwd = 1,  # 中间的竖的虚线
                  refline_lty = "dashed",
                  refline_col = "grey20")
p = forest(OR[,c(1:4,10,6,9)],
           est = OR$or,         #效应值
           lower = OR$or_lci95, #置信区间下限
           upper = OR$or_uci95, #置信区间上限
           ci_column = 5,       #在哪一列画森林图，选空的那一列
           ref_line = 1,        #参考线位置
           xlim = c(-1, 3),      #设置轴范围
           ticks_at = c(0, 1, 2, 3),#设置刻度
           theme = tm)
edit_plot(p,
          which = "background",
          row = c(3),
          gp = gpar(fill = "lightsteelblue1"))
dev.off()

#SNP森林图
res_single = mr_singlesnp(mr_data_cleDR)
#绘制SNP森林图
pdf(file = "Results/plot/森林图_ALL.pdf", width = 9, heigh = 11)
mr_forest_plot(res_single)
dev.off()

#leave-one-out DRalysis(留一法分析)
pdf(file = "Results/plot/leave-one-out_ALL.pdf", width = 9, heigh = 11)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(mr_data_cleDR))
dev.off()
