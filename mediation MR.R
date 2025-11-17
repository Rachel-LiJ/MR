rm(list = ls())
library(dplyr)
library(data.table)
library(vroom)
library(VariantAnnotation)
library(gwasvcf)
library(gwasglue)
library(ieugwasr)
library(plinkbinr)
library(TwoSampleMR)
library(MRPRESSO)
library(vcfR)

setwd("/Users/lijingli/data/MR")


########################################
#暴露—中介两样本MR
########################################
exposure_dat_clumped<-read.csv("DR_AD_Results/DR_exposure_dat_clumped.csv")
exposure_dat_clumped$id.exposure="DR"
outcome_finn = fread("~/data/MR/GWAS/AD/GCST007320_GRCh37.tsv")

#ieu-b-110 LDL 440546	
#ieu-b-109 HDL 403943 ebi-a-GCST90025956 	400754
#ieu-b-111 triglycerides 441016
#####ieu-b-118 HOMA-IR 37037
#ebi-a-GCST000568 Fasting blood glucose 46186
#/ieu-a-1 Adiponectin 39883
#met-d-SFA Saturated fatty acids 114999
#met-d-PUFA Polyunsaturated fatty acids
#met-d-Omega_6 Omega-6 fatty acids
#met-d-Omega_3 Omega-3 fatty acids
#met-d-MUFA Monounsaturated fatty acids
#met-a-466 Glutamate

#ukb-b-19953 Body mass index (BMI)
#ebi-a-GCST001212 Proinsulin levels
#ukb-b-15445 Treatment/medication code: insulin product

#####prot-a-2087 Oxysterols receptor LXR-beta
#ieu-b-108 apolipoprotein B
#####ieu-a-1012 Plasma cortisol

#prot-c-4337_49_2 CRP
#prot-c-4673_13_2 IL-6
#prot-c-2597_8_3 VEGF
#ebi-a-GCST90014006 HbA1c

#prot-a-131 Apolipoprotein E (isoform E3)

data <- readVcf("mediation/prot-a-132.vcf.gz")

Total_cholesterol_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")

LDL_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
HDL_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
triglycerides_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
glucose_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
HOMA_IR_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
Adiponectin_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
SFA_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
PUFA_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
Omega_6_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
Omega_3_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
MUFA_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
BMI_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
Proinsulin_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
apolipoproteinB_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
LXR_beta_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
CRP_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
IL_6_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
VEGF_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
cortisol_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
Glutamate_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
HbA1c_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")
AOPE_data = gwasvcf_to_TwoSampleMR(vcf = data, type="outcome")

colnames(LDL_data)

#修改暴露和中介的名称
Total_cholesterol_data$id.outcome = "Total_cholesterol"
LDL_data$id.outcome = "LDL" 
HDL_data$id.outcome = "HDL" 
triglycerides_data$id.outcome = "Triglycerides"
glucose_data$id.outcome = "Fasting_blood_glucose"
HOMA_IR_data$id.outcome = "HOMA_IR"
Adiponectin_data$id.outcome = "Adiponectin"
SFA_data$id.outcome = "SFA"
PUFA_data$id.outcome = "PUFA"
Omega_6_data$id.outcome = "Omega_6"
Omega_3_data$id.outcome = "Omega_3"
MUFA_data$id.outcome = "MUFA"
BMI_data$id.outcome = "BMI"
Proinsulin_data$id.outcome = "Proinsulin"
apolipoproteinB_data$id.outcome = "ApolipoproteinB"
LXR_beta_data$id.outcome = "LXR_beta"
CRP_data$id.outcome = "CRP"
IL_6_data$id.outcome = "IL_6"
VEGF_data$id.outcome = "VEGF"
cortisol_data$id.outcome = "Cortisol"
Glutamate_data$id.outcome = "Glutamate"
HbA1c_data$id.outcome = "HbA1c"
AOPE_data$id.outcome = "APOE"

dat_XM <- harmonise_data(exposure_dat_clumped, AOPE_data)

dat_XM<-dat_XM[dat_XM$mr_keep==TRUE,]

res_mi_XM <- generate_odds_ratios(mr(dat_XM,method_list="mr_ivw"))
res_mi_XM


write.csv(res_mi_XM, "mediation/result/MRresult_XM_APOE.csv")
#res_mi_XM<-read.csv("mediation/result/MRresult_XM_Total_cholesterol.csv")

#####################################
#中介—结局两样本MR
#####################################
mediation_dat_exposure<-HbA1c_data
mediation_dat_exposure<-format_data(dat = as.data.frame(mediation_dat_exposure),                            
                                    type = "exposure",
                                    snp_col = "SNP",
                                    beta_col = "beta.outcome",
                                    se_col = "se.outcome",
                                    effect_allele_col = "effect_allele.outcome",
                                    other_allele_col = "other_allele.outcome",
                                    pval_col = "pval.outcome",
                                    eaf_col = "eaf.outcome",
                                    samplesize_col = "samplesize.outcome",
                                    id_col = "id.outcome"
)

mediation_dat_exposure<- subset(mediation_dat_exposure, pval.exposure<=5e-08) #筛选P<5e-08

mediation_dat_exposure_clumped = ieugwasr::ld_clump(dplyr::tibble(rsid = mediation_dat_exposure$SNP,
                                                          pval = mediation_dat_exposure$pval.exposure),
                                            clump_kb = 10000,
                                            clump_r2 = 0.001,
                                            clump_p = 1,
                                            bfile = "1kg.v3/EUR",
                                            plink_bin = plinkbinr::get_plink_exe(),
                                            pop = "EUR")
mediation_dat_exposure_clumped = mediation_dat_exposure[mediation_dat_exposure$SNP %in% mediation_dat_exposure_clumped$rsid,]
dim(mediation_dat_exposure_clumped)#24 14

outcome_finn_dat = subset(outcome_finn, outcome_finn$variant_id %in% mediation_dat_exposure_clumped$SNP)
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
                               other_allele_col = "other_allele",
                               id_col = "AD")

#outcome_finn_dat = subset(outcome_finn_dat, pval.outcome>5e-08)
dim(outcome_finn_dat)#[1] 24 12
outcome_finn_dat$id.outcome = "AD" #修改中介和结局的名称

dat_MY <- harmonise_data(mediation_dat_exposure_clumped, outcome_finn_dat)
dat_MY<-dat_MY[dat_MY$mr_keep==TRUE,]

res_mi_MY <- generate_odds_ratios(mr(dat_MY,method_list="mr_ivw"))
res_mi_MY

#write.csv(dat_MY, "mediation/result/MRdata_MY_Total_cholesterol.csv")
write.csv(res_mi_MY, "mediation/result/MRresult_MY_HbA1c.csv")

##########################################
#计算中介效应
##########################################
res_XY<-read.csv("~/data/MR/Results/MRresult_DR_AD_E.csv")
all_beta = res_XY[3,8]  #提取暴露与结局的总效应（均为IVW方法结果的beta值）
all_beta
## [1] -0.3785459
XM_beta = res_mi_XM$b  #提取暴露与中介的效应XM_beta
XM_beta
## [1] -0.2112891
MY_beta = res_mi_MY$b  #提取中介与结局的效应MY_beta
MY_beta
## [1] 0.8643837
aa = XM_beta * MY_beta  #计算中介效应aa
aa
## [1] -0.1826348
bb = aa/all_beta  #计算中介占比bb
bb
## [1] 0.4824641
cc = all_beta - aa  #计算直接效应cc
cc
## [1] -0.1959111

##########################################
#系数乘积检验以及中介效应的置信区间
##########################################
XM_se = res_mi_XM$se  #提取暴露与中介的se值（均为IVW方法结果的se值）
XM_se
## [1] 0.03785756
MY_se = res_mi_MY$se   #提取中介与结局的se值MY_se
MY_se
## [1] 0.05008443
S <- sqrt( (MY_beta^2) * (XM_se^2) + (XM_beta^2) * (MY_se^2) )
S
## [1] 0.04402492
Z = aa/S  #计算统计量Z
Z
## [1] -4.148442
Pval = 2*pnorm(q = abs(Z), lower.tail=FALSE)  #计算中介效应的P值
Pval

## [1] 3.347454e-05
#中介效应的置信区间
a1 = aa + S         ## [1] -0.1386099
a1
a2 = aa - S          ## [1] -0.2266597
a2
b1 = a1/all_beta  ## [1] 0.366164
b1
b2 = a2/all_beta  ## [1] 0.5987642
b2

mediation_result<-data.frame(
  beta_XM=XM_beta,
  beta_MY=MY_beta,
  beta_mediation=aa,
  beta_percent=bb,
  se_XM=XM_se,
  se_MY=MY_se,
  S_standard=S,
  Z=Z,
  pval=Pval,
  beta_a1=a1,
  beta_a2=a2,
  beta_b1=b1,
  beta_b2=b2
)

write.csv(mediation_result, "mediation/result/mediation_result_HbA1c.csv")

rm(data)

