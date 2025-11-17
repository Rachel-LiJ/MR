rm(list = ls())
setwd("/Users/lijingli/data/MR")
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

source("function.R")

#DR
DR_exposure_finn = fread("~/data/MR/GWAS/DR/GCST90479891.tsv.gz")
DR_exposure_finn$standard_error = abs(log(DR_exposure_finn$odds_ratio) / qnorm(DR_exposure_finn$p_value / 2))
DR_exposure_finn$Z =  log(DR_exposure_finn$odds_ratio) / DR_exposure_finn$standard_error
DR_exposure_finn$beta <- log(DR_exposure_finn$odds_ratio)


head(DR_exposure_finn)
colnames(DR_exposure_finn)
DR_exposure_finn <- subset(DR_exposure_finn,p_value< 5e-08)
dim(DR_exposure_finn)#[1] 826  13
DR_exposure_finn<-data.frame(DR_exposure_finn)
DR_exposure_finn = format_data(dat = DR_exposure_finn,
                               type = "exposure",
                               snp_col = "rsid",
                               beta_col = "beta",
                               se_col = "standard_error",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele",
                               eaf_col = "effect_allele_frequency",
                               pval_col = "p_value",
                               samplesize_col = "n" )

write.csv(DR_exposure_finn,file = "Results/DR_exposure_finn_ALL.csv")
DR_exposure_dat <- DR_exposure_finn
#修改列名
#去除连锁不平衡
DR_exposure_dat_clumped = ieugwasr::ld_clump(dplyr::tibble(rsid = DR_exposure_dat$SNP,
                                                           pval = DR_exposure_dat$pval.exposure),
                                             clump_kb = 10000,
                                             clump_r2 = 0.001,
                                             clump_p = 1,
                                             bfile = "1kg.v3/EUR",
                                             plink_bin = plinkbinr::get_plink_exe(),
                                             pop = "EUR")
DR_exposure_dat_clumped$id.exposure = "DR"
DR_exposure_dat_clumped = DR_exposure_dat[DR_exposure_dat$SNP %in% DR_exposure_dat_clumped$rsid,]
#保留F值>10

# Calculate R2 and F stat for exposure data
# Liberal hypertension F stat
#DR_exposure_dat_clumped$N<-7049
DR_exposure_dat_clumped$r2 <- (2 * (DR_exposure_dat_clumped$beta.exposure^2) * DR_exposure_dat_clumped$eaf.exposure * (1 - DR_exposure_dat_clumped$eaf.exposure)) /
  (2 * (DR_exposure_dat_clumped$beta.exposure^2) * DR_exposure_dat_clumped$eaf.exposure * (1 - DR_exposure_dat_clumped$eaf.exposure) +
     2 * DR_exposure_dat_clumped$samplesize.exposure * DR_exposure_dat_clumped$eaf.exposure * (1 - DR_exposure_dat_clumped$eaf.exposure) * DR_exposure_dat_clumped$se.exposure^2)

DR_exposure_dat_clumped$F <- DR_exposure_dat_clumped$r2 * (DR_exposure_dat_clumped$samplesize.exposure - 2) / (1 - DR_exposure_dat_clumped$r2)

DR_exposure_dat_f <-DR_exposure_dat_clumped[DR_exposure_dat_clumped$F>10,]
dim(DR_exposure_dat_clumped)#24 14
write.csv(DR_exposure_dat_clumped,file = "Results/DR_exposure_dat_clumped_ALL.csv")
save(DR_exposure_dat_clumped,file = "Results/Rdata/DR_ALL.Rdata")

data<-LDtrait(
  DR_exposure_dat_f$SNP,#	
  #between 1 - 50 variants, using an rsID or chromosome coordinate (e.g. "chr7:24966446"). All input variants must match a bi-allelic variant.
  pop = "EUR",#人群，可以使用list_pop()查看支持人群缩写
  r2d = "r2",#	
  #use "r2" to filter desired output from a threshold based on estimated LD R2 (R squared) or "d" for LD D' (D-prime), default = "r2".
  r2d_threshold = 0.1,#筛选阈值，小的将被筛除
  win_size = 5e+05,
  token = "76d9a4a9589a",#输入邮箱获得token
  file = FALSE,
  genome_build = "grch38",#在三个选项中任选其一…'grch37'用于基因组构建grch37 (hg19)， 'grch38'用于grch38 (hg38)， 'grch38_high_coverage'用于grch38 High Coverage (hg38) 1000基因组计划数据集。默认为GRCh37 (hg19)。
  api_root = "https://ldlink.nih.gov/LDlinkRest"
)

trait <- c(
  # Alzheimer’s disease & dementia
  "Alzheimer", "Dementia", "Cognitive", "Memory", "Neurodegenerative", 
  
  # Diabetic retinopathy & eye diseases
  "Diabetic retinopathy", "Macular", "Retina", "Glaucoma", "Ophthalm", "Optic",
  
  # Lipid metabolism
  "LDL", "HDL", "Triglyceride", "Cholesterol", "Lipid", "Apolipoprotein", 
  
  # Inflammation & immune response
  "Inflammation", "C-reactive protein", "Cytokine", "IL-6", "IL6", "TNF", "Immune",
  
  # Cardiometabolic risk
  "Type 2 diabetes", "Diabetes", "HbA1c", "Insulin", "Obesity", "BMI", "Blood pressure",
  
  # Stroke & vascular
  "Stroke", "Myocardial infarction", "Coronary artery", "Atherosclerosis",
  
  # Other possible confounders
  "Smoking", "Alcohol", "Education", "Depression", "Sleep", "Physical activity"
)

tcL<-function(snps,trait){
  as=LDtrait(
    snps=snps,
    pop = "ALL",
    r2d = "r2",
    r2d_threshold = 0.1,
    win_size = 5e+05,
    token = "76d9a4a9589a",
    file = FALSE,
    genome_build = "grch38",
    api_root = "https://ldlink.nih.gov/LDlinkRest")
  #去除含有trait特征的SNP
  # 直接匹配多个关键词（忽略大小写）
  
  pattern <- paste(trait, collapse = "|")
  match_index <- str_detect(as[,"GWAS_Trait"], regex(pattern, ignore_case = TRUE))
  
  if (any(match_index)) {
    matched_snps <- unique(as[match_index, "Query"])
    filtered_snps <- snps[!(snps %in% matched_snps)]
  } else {
    filtered_snps <- snps
  }
  
  return(filtered_snps)
}

snp<-tcL(DR_exposure_dat_clumped$SNP,trait)

#结局----
outcome_finn = fread("~/data/MR/GWAS/AD/GCST007320_GRCh37.tsv")
colnames(outcome_finn)
outcome_finn_dat = subset(outcome_finn, outcome_finn$variant_id %in% DR_exposure_dat_clumped$SNP)
dim(outcome_finn_dat)# 24 13
outcome_finn_dat<-data.frame(outcome_finn_dat)
outcome_finn_dat = format_data(dat = outcome_finn_dat,
                               type = "outcome",
                               snp_col = "variant_id",
                               beta_col = "beta",
                               pval_col = "p_value",
                               se_col = "standard_error",
                               eaf_col = "effect_allele_frequency",
                               effect_allele_col = "effect_allele",
                               other_allele_col = "other_allele")

outcome_finn_dat = subset(outcome_finn_dat, pval.outcome>5e-08)
dim(outcome_finn_dat)#[1] 24 12

#修改暴露和结局的名称
DR_exposure_dat_clumped$id.exposure = "DR"
outcome_finn_dat$id.outcome = "AD"

mr_data_finn = harmonise_data(exposure_dat = DR_exposure_dat_clumped,
                              outcome_dat = outcome_finn_dat,
                              action= 2)   #去除回文序列，2去除，1不去除
table(mr_data_finn$mr_keep)#1 23

mr_data_finn<-mr_data_finn[mr_data_finn$mr_keep==TRUE,]

write.csv(mr_data_finn, file = "Results/mr_data_ALL.csv")
View(mr_data_finn)
mr_data_finn <- data.table::fread("data/AD/mr_DR_data_PGC.csv",data.table = F)

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
mr_data_cleDR= mr_data_finn
mr_data_cleDR = mr_data_finn[-c(44),]
mr_data_cleDR=mr_data_cleDR[!mr_data_cleDR$SNP %in% c("rs7903146"),]
#mr_data_cleDR = mr_data_finn
#去除离群值
#rs = sin$SNP[sin$b < -3|sin$b > 3] #找到离群的SNP
#mr_data_cleDR = mr_data_finn[!mr_data_finn$SNP %in% rs,] #从原始数据框中移除这些SNP

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
write.csv(OR_cleDR, "Results/MRresult_DR_AD_ALL.csv")

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
