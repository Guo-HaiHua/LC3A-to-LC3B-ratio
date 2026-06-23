###
path_data <- "G:\\Project\\ADNI_LC3AB\\Data\\"
path_adni <- "G:\\Project\\Cohort\\ADNI\\Original\\"
path_model <- "G:\\Project\\Model_R\\HHG\\"
path_output <- "G:\\Project\\ADNI_LC3AB\\Results\\"
###
library(data.table)
library(dplyr)
library(tidyr)
library(writexl)
###
#
LC3 <- read.csv(paste0(path_data,"LC3AB_p62.csv"))
#
LC3 <- subset(LC3, !is.na(LC3$LC3AB))
LC3$LC3AB <- scale(LC3$LC3AB)
#
PGRN_STREM2 <- read.csv(paste0(path_adni,"PGRN_STREM2_used.csv"))
PHC <- read.csv(paste0(path_adni,"ADSP_PHC_COGN_10_05_22_29Mar2025.csv"))
PET <- read.csv(paste0(path_adni,"PET.csv"))
demo <- read.csv(paste0(path_adni,"demo_fo_blnoAD.csv"))
PET_TAU_processed_tobl <- read.csv(paste0(path_adni,"PET_TAU_processed_tobl.csv"))
plasma_nfl <- read.csv(paste0(path_adni,"ADNI_BLENNOWPLASMANFLLONG_10_03_18_15Mar2026.csv"))
###
LC3_bl_demo <- merge(LC3, demo, by = c("RID"),all.x=T)
names(LC3_bl_demo)[14] <- "VISCODE2"
merge_df0 <- merge(LC3_bl_demo, PHC, by = c("RID","VISCODE2"),all.x=T)
merge_df1 <- merge(merge_df0, PET, by = c("RID","VISCODE2"),all.x=T)
merge_df2 <- merge(merge_df1, PGRN_STREM2, by = c("RID","VISCODE2"),all.x=T)
merge_df3 <- merge(merge_df2, PET_TAU_processed_tobl, by = c("RID","VISCODE2"),all.x=T)
merge_df <- merge(merge_df3, plasma_nfl, by = c("RID","VISCODE2"),all.x=T)
###
source(paste0(path_model,"ghhlmm.R"))
#Cog
lmm_Nestedlist <- NULL
for (i in names(merge_df)[c(10)]){
  for (j in names(merge_df)[c(36,38,40,49,142:145)]){
    lmm_cov <- names(merge_df)[c(7,18,20,21,26)]
    lmm_Nestedlist[[i]][[j]] <- ghhlmm(merge_df[i],merge_df[j],merge_df["Years_bl"],merge_df["RID"],merge_df)
  }
}
lmm_Flatlist <- NULL
for (i in names(merge_df)[c(10)]){
  lmm_Flatlist[[i]] <- Reduce(rbind.data.frame,lmm_Nestedlist[[i]],accumulate = F)
}
lmm_df_results <- Reduce(rbind.data.frame,lmm_Flatlist,accumulate = F)
#
lmm_df_results$P_fdr <- ave(lmm_df_results$P, lmm_df_results$Factor, 
                            FUN = function(x) p.adjust(x, method = "fdr"))
 ## output
# write_xlsx(lmm_df_results, path = paste0(path_output, "lmm_LC3AB_Cog_results_fo_M1_SOMA.xlsx"))
write_xlsx(lmm_df_results, path = paste0(path_output, "lmm_LC3AB_Cog_results_fo_M1_SOMA_covp62.xlsx"))

# d <- merge_df
# d <- d[!is.na(d$LC3B) & !is.na(d$CDRSB) & !is.na(d$Years_bl),]
# d <- d %>%
#   group_by(RID) %>%
#   filter(n() > 1)
# 
# fit <- lmer(scale(CDRSB) ~ LC3B*Years_bl + AGE + PTGENDER_C + PTEDUCAT + APOE4_C+(1+Years_bl | RID), data=d)
# summary(fit)


27,29,31:35,148,150,162,169

#ADbio
lmm_Nestedlist <- NULL
for (i in names(merge_df)[c(10)]){
  for (j in names(merge_df)[c(27,29,31:33,176,180,195)]){
    lmm_cov <- names(merge_df)[c(7,18,20,21,26)]
    lmm_Nestedlist[[i]][[j]] <- ghhlmm(merge_df[i],merge_df[j],merge_df["Years_bl"],merge_df["RID"],merge_df)
  }
}
lmm_Flatlist <- NULL
for (i in names(merge_df)[c(10)]){
  lmm_Flatlist[[i]] <- Reduce(rbind.data.frame,lmm_Nestedlist[[i]],accumulate = F)
}
lmm_df_results1 <- Reduce(rbind.data.frame,lmm_Flatlist,accumulate = F)
#


#PET_TAU
lmm_Nestedlist <- NULL
for (i in names(merge_df)[c(10)]){
  for (j in names(merge_df)[c(187)]){
    lmm_cov <- names(merge_df)[c(7,18,20,21,26,189)]
    lmm_Nestedlist[[i]][[j]] <- ghhlmm(merge_df[i],merge_df[j],merge_df["Years_bl"],merge_df["RID"],merge_df)
  }
}
lmm_Flatlist <- NULL
for (i in names(merge_df)[c(10)]){
  lmm_Flatlist[[i]] <- Reduce(rbind.data.frame,lmm_Nestedlist[[i]],accumulate = F)
}
lmm_df_results2 <- Reduce(rbind.data.frame,lmm_Flatlist,accumulate = F)
#
lmm_df_results <- rbind(lmm_df_results1,lmm_df_results2)
lmm_df_results$P_fdr <- ave(lmm_df_results$P, lmm_df_results$Factor, 
                            FUN = function(x) p.adjust(x, method = "fdr"))
## output
# write_xlsx(lmm_df_results, path = paste0(path_output, "lmm_LC3AB_ADbio_results_fo_M1_SOMA.xlsx"))
write_xlsx(lmm_df_results, path = paste0(path_output, "lmm_LC3AB_ADbio_results_fo_M1_SOMA_covp62.xlsx"))



#Bra
lmm_Nestedlist <- NULL
for (i in names(merge_df)[c(10)]){
  for (j in names(merge_df)[c(67:72)]){
    lmm_cov <- names(merge_df)[c(7,18,20,21,26,73)]
    lmm_Nestedlist[[i]][[j]] <- ghhlmm(merge_df[i],merge_df[j],merge_df["Years_bl"],merge_df["RID"],merge_df)
  }
}
lmm_Flatlist <- NULL
for (i in names(merge_df)[c(10)]){
  lmm_Flatlist[[i]] <- Reduce(rbind.data.frame,lmm_Nestedlist[[i]],accumulate = F)
}
lmm_df_results <- Reduce(rbind.data.frame,lmm_Flatlist,accumulate = F)
#
lmm_df_results$P_fdr <- ave(lmm_df_results$P, lmm_df_results$Factor, 
                            FUN = function(x) p.adjust(x, method = "fdr"))
## output
# write_xlsx(lmm_df_results, path = paste0(path_output, "lmm_LC3AB_Bra_results_fo_M1_SOMA.xlsx"))
write_xlsx(lmm_df_results, path = paste0(path_output, "lmm_LC3AB_Bra_results_fo_M1_SOMA_covp62.xlsx"))

# #PGRN
# lmm_Nestedlist <- NULL
# for (i in names(merge_df)[c(10)]){
#   for (j in names(merge_df)[c(176,180)]){
#     lmm_cov <- names(merge_df)[c(18,20,21,26)]
#     lmm_Nestedlist[[i]][[j]] <- ghhlmm(merge_df[i],merge_df[j],merge_df["Years_bl"],merge_df["RID"],merge_df)
#   }
# }
# lmm_Flatlist <- NULL
# for (i in names(merge_df)[c(10)]){
#   lmm_Flatlist[[i]] <- Reduce(rbind.data.frame,lmm_Nestedlist[[i]],accumulate = F)
# }
# lmm_df_results <- Reduce(rbind.data.frame,lmm_Flatlist,accumulate = F)
# #
# lmm_df_results$P_fdr <- ave(lmm_df_results$P, lmm_df_results$Factor, 
#                             FUN = function(x) p.adjust(x, method = "fdr"))
# ## output
# write_xlsx(lmm_df_results, path = paste0(path_output, "lmm_LC3AB_PGRN_results_fo_SOMA.xlsx"))

# 
# merge_df <- merge_df[!is.na(merge_df$LC3B) & !is.na(merge_df$MSD_PGRNCORRECTED) & !is.na(merge_df$Years_bl),]
# # Only duplicate ID rows are retained
# merge_df <- merge_df %>% group_by(across(all_of("RID"))) %>%  filter(n() > 1)
# a <- lmer(scale(MSD_PGRNCORRECTED) ~ LC3B*Years_bl + AGE + PTGENDER_C + PTEDUCAT + APOE4_C + (1+Years_bl | RID), data = merge_df)
# summary(a)
