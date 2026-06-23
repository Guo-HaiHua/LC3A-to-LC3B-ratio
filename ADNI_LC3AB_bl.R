###
path_data <- "G:\\Project\\ADNI_LC3AB\\Data\\"
path_adni <- "G:\\Project\\Cohort\\ADNI\\Original\\"
path_model <- "G:\\Project\\Model_R\\HHG\\"
path_output <- "G:\\Project\\ADNI_LC3AB\\Results\\"
###
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(writexl)
###
#
LC3 <- read.csv(paste0(path_data,"LC3AB_p62.csv"))
LC3[, 7:9] <- scale(LC3[, 7:9])
a <- lm(LC3A ~ LC3B, data = LC3)
summary(a)
#
LC3 <- subset(LC3, !is.na(LC3$LC3AB))
###
PGRN_STREM2 <- read.csv(paste0(path_adni,"PGRN_STREM2_used.csv"))
PHC <- read.csv(paste0(path_adni,"ADSP_PHC_COGN_10_05_22_29Mar2025.csv"))
PET <- read.csv(paste0(path_adni,"PET.csv"))
demo_bl <- read.csv(paste0(path_adni,"demo_bl_blnoAD.csv"))
PET_TAU_processed_tobl <- read.csv(paste0(path_adni,"PET_TAU_processed_tobl.csv"))
plasma_nfl_bl <- subset(read.csv(paste0(path_adni,"ADNI_BLENNOWPLASMANFL_10_03_18_19Apr2026.csv")),VISCODE2 == "bl")
###
LC3_bl_demo <- merge(LC3, demo_bl, by = c("RID","VISCODE2"))
merge_df0 <- merge(LC3_bl_demo, PHC, by = c("RID","VISCODE2"),all.x=T)
merge_df1 <- merge(merge_df0, PET, by = c("RID","VISCODE2"),all.x=T)
merge_df2 <- merge(merge_df1, PGRN_STREM2, by = c("RID","VISCODE2"),all.x=T)
merge_df3 <- merge(merge_df2, PET_TAU_processed_tobl, by = c("RID","VISCODE2"),all.x=T)
merge_df <- merge(merge_df3, plasma_nfl_bl, by = c("RID","VISCODE2"),all.x=T)


###
source(paste0(path_model,"ghhlm.R"))
#Cog
lm_Nestedlist <- NULL
for (i in names(merge_df)[c(7:9)]){
  for (j in names(merge_df)[c(35,37,39,48,141:144)]){
    lm_cov <- names(merge_df)[c(6,17,19,20,25)]
    lm_Nestedlist[[i]][[j]] <- ghhlm(merge_df[i],merge_df[j],merge_df)
  }
}
lm_Flatlist <- NULL
for (i in names(merge_df)[c(7:9)]){
  lm_Flatlist[[i]] <- Reduce(rbind.data.frame,lm_Nestedlist[[i]],accumulate = F)
}
lm_df_results <- Reduce(rbind.data.frame,lm_Flatlist,accumulate = F)
##
lm_df_results$P_fdr <- ave(lm_df_results$P, lm_df_results$Factor,
                            FUN = function(x) p.adjust(x, method = "fdr"))
##
# write_xlsx(lm_df_results, path = paste0(path_output, "lm_LC3AB_Cog_results_bl_M1_SOMA.xlsx"))
write_xlsx(lm_df_results, path = paste0(path_output, "lm_LC3AB_Cog_results_bl_M1_SOMA_covp62.xlsx"))

26,28,30:32,148,150,162
#ADbio
lm_Nestedlist <- NULL
for (i in names(merge_df)[c(7:9)]){
  for (j in names(merge_df)[c(26,28,30:32,175,179,193)]){
    lm_cov <- names(merge_df)[c(6,17,19,20,25)]
    lm_Nestedlist[[i]][[j]] <- ghhlm(merge_df[i],merge_df[j],merge_df)
  }
}
lm_Flatlist <- NULL
for (i in names(merge_df)[c(7:9)]){
  lm_Flatlist[[i]] <- Reduce(rbind.data.frame,lm_Nestedlist[[i]],accumulate = F)
}
lm_df_results <- Reduce(rbind.data.frame,lm_Flatlist,accumulate = F)
##



#PET_TAU
lm_Nestedlist <- NULL
for (i in names(merge_df)[c(9)]){
  for (j in names(merge_df)[c(186)]){
    lm_cov <- names(merge_df)[c(6,17,19,20,25,188)]
    lm_Nestedlist[[i]][[j]] <- ghhlm(merge_df[i],merge_df[j],merge_df)
  }
}
lm_Flatlist <- NULL
for (i in names(merge_df)[c(9)]){
  lm_Flatlist[[i]] <- Reduce(rbind.data.frame,lm_Nestedlist[[i]],accumulate = F)
}
lm_df_results2 <- Reduce(rbind.data.frame,lm_Flatlist,accumulate = F)
##
lm_df_results <- rbind(lm_df_results1,lm_df_results2)
##
lm_df_results$P_fdr <- ave(lm_df_results$P, lm_df_results$Factor,
                           FUN = function(x) p.adjust(x, method = "fdr"))
##
# write_xlsx(lm_df_results, path = paste0(path_output, "lm_LC3AB_ADbio_results_bl_M1_SOMA.xlsx"))
write_xlsx(lm_df_results, path = paste0(path_output, "lm_LC3AB_ADbio_results_bl_M1_SOMA_covp62.xlsx"))


#Bra
lm_Nestedlist <- NULL
for (i in names(merge_df)[c(9)]){
  for (j in names(merge_df)[c(66:71)]){
    lm_cov <- names(merge_df)[c(6,17,19,20,25,72)]
    lm_Nestedlist[[i]][[j]] <- ghhlm(merge_df[i],merge_df[j],merge_df)
  }
}
lm_Flatlist <- NULL
for (i in names(merge_df)[c(9)]){
  lm_Flatlist[[i]] <- Reduce(rbind.data.frame,lm_Nestedlist[[i]],accumulate = F)
}
lm_df_results <- Reduce(rbind.data.frame,lm_Flatlist,accumulate = F)
#
##
lm_df_results$P_fdr <- ave(lm_df_results$P, lm_df_results$Factor,
                           FUN = function(x) p.adjust(x, method = "fdr"))
##
# write_xlsx(lm_df_results, path = paste0(path_output, "lm_LC3AB_Bra_results_bl_M1_SOMA.xlsx"))
write_xlsx(lm_df_results, path = paste0(path_output, "lm_LC3AB_Bra_results_bl_M1_SOMA_covp62.xlsx"))

# #PGRN
# lm_Nestedlist <- NULL
# for (i in names(merge_df)[c(9)]){
#   for (j in names(merge_df)[c(175,179)]){
#     lm_cov <- names(merge_df)[c(17,19,20,25)]
#     lm_Nestedlist[[i]][[j]] <- ghhlm(merge_df[i],merge_df[j],merge_df)
#   }
# }
# lm_Flatlist <- NULL
# for (i in names(merge_df)[c(9)]){
#   lm_Flatlist[[i]] <- Reduce(rbind.data.frame,lm_Nestedlist[[i]],accumulate = F)
# }
# lm_df_results <- Reduce(rbind.data.frame,lm_Flatlist,accumulate = F)
# #
# ##
# lm_df_results$P_fdr <- ave(lm_df_results$P, lm_df_results$Factor,
#                            FUN = function(x) p.adjust(x, method = "fdr"))
# # write.csv(lm_df_results,file = paste0(path_output,"LC3_PGRN_results_bl_M1.csv"),row.names = F)
# write_xlsx(lm_df_results, path = paste0(path_output, "lm_LC3AB_PGRN_results_bl_M1_SOMA.xlsx"))


# ### source ghhlm_inter function----
# source(paste0(path_model,"ghhlm_inter.R"))
# ### ghhlm_inter
# lm_Nestedlist_inter <- NULL
# for (i in names(merge_df)[c(9)]){
#   for (j in names(merge_df)[c(26,28,30:32,148,150,162)]){
#     for (k in names(merge_df)[c(35,37,39,48,141:144)]){
#       lm_intercov <- names(merge_df)[c(17,19,20,25)]
#       lm_Nestedlist_inter[[i]][[j]][[k]] <- ghhlm_inter(merge_df[i],merge_df[j],merge_df[k],merge_df)
#     }
#   }
# }
# lm_Flatlist_inter <- NULL
# for (i in names(merge_df)[c(9)]){
#   for (j in names(merge_df)[c(26,28,30:32,148,150,162)]){
#     lm_Flatlist_inter[[i]][[j]] <- Reduce(rbind.data.frame,lm_Nestedlist_inter[[i]][[j]],accumulate = F)
#     }
#   lm_Flatlist_inter[[i]] <- Reduce(rbind.data.frame,lm_Flatlist_inter[[i]],accumulate = F)
# }
# lm_df_inter_results <- Reduce(rbind.data.frame,lm_Flatlist_inter,accumulate = F)
# lm_df_inter_results$P_fdr <- ave(lm_df_inter_results$PI,
#                            FUN = function(x) p.adjust(x, method = "fdr"))
# ## output
# # write.csv(lm_df_inter_results,file = paste0(path_output,"LC3_bio_Cog_inter_results_bl_M1.csv"),row.names = F)
# write_xlsx(lm_df_inter_results, path = paste0(path_output, "LC3_bio_Cog_inter_results_bl_M1.xlsx"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #################################
# path_adni <- "G:\\Project\\Cohort\\ADNI\\Original\\"
# lifestyle <- read.csv(paste0(path_adni,"lifestyle.csv"))
# med <- read.csv(paste0(path_adni,"med_history_Hamilton.csv"))
# lifestyle_med <- merge(lifestyle, med, by = c("RID","VISCODE2"), all.x=T)
# colSums(is.na(lifestyle_med[, c(5:6,9:10)]))
# # 
# # colSums(is.na(med[, 5:6]))
# write.csv(lifestyle_med,file = paste0(path_adni,"lifestyle_med.csv"),row.names = F)
# 
# 
# ###
# demo <- read.csv(paste0(path_adni,"ADNIMERGE_26Jun2024_bl.csv"))
# demo$Stage[demo$ABETA > 976.6 & demo$PTAU < 21.8 & demo$TAU < 245] <- "Stage 0" #A-T-N-
# demo$Stage[demo$ABETA > 976.6 & (demo$PTAU >= 21.8 | demo$TAU >= 245)] <- "SNAP" #A-T+(N)-, A-T-(N)+, or A-T+(N)+
# demo$Stage[demo$ABETA <= 976.6 & demo$PTAU < 21.8 & demo$TAU < 245] <- "Stage 1" #A+T-N-
# demo$Stage[demo$ABETA <= 976.6 & (demo$PTAU >= 21.8 | demo$TAU >= 245)] <- "Stage 2" #A+T+N- #A+T+N+ #A+T-N+
# write.csv(demo,file = paste0(path_adni,"demo_bl.csv"),row.names = F)
# demo <- demo %>% filter(DX_bl_C != "AD", DX_bl_C != "")
# write.csv(demo,file = paste0(path_adni,"demo_bl_noAD.csv"),row.names = F)
# #
# #A+  ABETA <= 976.6
# #T+ PTAU >= 21.8
# #N+ TAU >= 245
# 
# 
# demo <- read.csv(paste0(path_adni,"ADNIMERGE_26Jun2024.csv"))
# demo$Stage[demo$ABETA > 976.6 & demo$PTAU < 21.8 & demo$TAU < 245] <- "Stage 0" #A-T-N-
# demo$Stage[demo$ABETA > 976.6 & (demo$PTAU >= 21.8 | demo$TAU >= 245)] <- "SNAP" #A-T+(N)-, A-T-(N)+, or A-T+(N)+
# demo$Stage[demo$ABETA <= 976.6 & demo$PTAU < 21.8 & demo$TAU < 245] <- "Stage 1" #A+T-N-
# demo$Stage[demo$ABETA <= 976.6 & (demo$PTAU >= 21.8 | demo$TAU >= 245)] <- "Stage 2" #A+T+N- #A+T+N+ #A+T-N+
# write.csv(demo,file = paste0(path_adni,"demo_fo.csv"),row.names = F)
# demo <- demo %>% filter(DX_bl_C != "AD", DX_bl_C != "")
# write.csv(demo,file = paste0(path_adni,"demo_fo_noAD.csv"),row.names = F)
# 
# 
# ###
# PET_FDG <- read.csv(paste0(path_adni,"PET_FDG.csv"))
# MetaROI <- subset(PET_FDG,ROINAME=="MetaROI")
# Top50PonsVermis <- subset(PET_FDG,ROINAME=="Top50PonsVermis")
# PET_FDG_2 <- merge(MetaROI,Top50PonsVermis,by=c("RID","VISCODE2","EXAMDATE"))
# # write.csv(PET_FDG_2,file = paste0(path_adni,"PET_FDG_2.csv"),row.names = F)
# 
# PET_ABETA <- read.csv(paste0(path_adni,"PET_ABETA.csv"))
# 
# PET_TAU <- read.csv(paste0(path_adni,"PET_TAU.csv"))
# 
# FDG_ABETA <- merge(PET_FDG_2, PET_ABETA, by=c("RID","VISCODE2"), all=T)
# PET <- merge(FDG_ABETA, PET_TAU, by=c("RID","VISCODE2"), all=T)
# write.csv(PET,file = paste0(path_adni,"PET.csv"),row.names = F)
# 
# PET <- read.csv(paste0(path_adni,"PET.csv"))
# PET$ABETA_Stage[PET$TRACER.ABETA == "FBB" & PET$SUMMARY_SUVR <= 1.08] <- "A-"
# PET$ABETA_Stage[PET$TRACER.ABETA == "FBB" & PET$SUMMARY_SUVR > 1.08] <- "A+"
# PET$ABETA_Stage[PET$TRACER.ABETA == "FBP" & PET$SUMMARY_SUVR <= 1.32] <- "A-"
# PET$ABETA_Stage[PET$TRACER.ABETA == "FBP" & PET$SUMMARY_SUVR > 1.32] <- "A+"
# write.csv(PET,file = paste0(path_adni,"PET_Stage.csv"),row.names = F)
# 
# 
