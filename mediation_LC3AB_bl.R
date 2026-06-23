###
path_data <- "G:\\Project\\ADNI_LC3AB\\Data\\"
path_adni <- "G:\\Project\\Cohort\\ADNI\\Original\\"
path_model <- "G:\\Project\\Model_R\\HHG\\"
path_output <- "G:\\Project\\ADNI_LC3AB\\Results\\"
###
library(lavaan)
library(dplyr)
library(data.table)
library(writexl)
library(progress)
###
###
#
LC3 <- read.csv(paste0(path_data,"LC3AB_p62.csv"))
LC3 <- subset(LC3, !is.na(LC3$LC3AB))
#
PGRN_STREM2 <- read.csv(paste0(path_adni,"PGRN_STREM2_used.csv"))
PHC <- read.csv(paste0(path_adni,"ADSP_PHC_COGN_10_05_22_29Mar2025.csv"))
PET <- read.csv(paste0(path_adni,"PET.csv"))
demo_bl <- read.csv(paste0(path_adni,"demo_bl_blnoAD.csv"))
###
LC3_bl_demo <- merge(LC3, demo_bl, by = c("RID","VISCODE2"),all.x=T)
merge_df <- merge(LC3_bl_demo, PHC, by = c("RID","VISCODE2"),all.x=T)
# merge_df1 <- merge(merge_df0, PET, by = c("RID","VISCODE2"),all.x=T)
# merge_df <- merge(merge_df1, PGRN_STREM2, by = c("RID","VISCODE2"),all.x=T)
###
merge_df[, c(66:71)] <- scale(merge_df[, c(66:71)])
###
source(paste0(path_model, "ghhmediate.R"))
### ghhmed
pb <- progress_bar$new(
  format = "  Running [:bar] :percent  | elapsed: :elapsed  | eta: :eta",
  total = 1 * 6* 8,
  clear = FALSE, width = 60, complete = "=", incomplete = "-", current = ">")
###
rs_list <- NULL
for (i in names(merge_df)[c(9)]) {
  for (j in names(merge_df)[c(66:71)]) {
    for (k in names(merge_df)[c(35,37,39,48,141:144)]) {
      #定义协变量
      mediate_covariates <- names(merge_df)[c(17,19,20,25,72)]
      # 定义模型公式
      mod_m <- as.formula(paste0(j, "~", i, "+", paste(mediate_covariates, collapse = "+")))
      mod_y <- as.formula(paste0(k, "~", i, "+", j, "+", paste(mediate_covariates, collapse = "+")))
      # 构建因变量的模型
      rs_list[[i]][[j]][[k]] <- ghhmediate(merge_df[i], merge_df[j], merge_df[k], merge_df, mod_m, mod_y)
      pb$tick()
    }
  }
}
rs <- NULL
for (i in names(merge_df)[c(9)]) {
  for (j in names(merge_df)[c(66:71)]) {
    rs[[i]][[j]] <- Reduce(rbind.data.frame, rs_list[[i]][[j]], accumulate = FALSE)
  }
  rs[[i]] <- Reduce(rbind.data.frame, rs[[i]], accumulate = FALSE)
}
med_rs <- Reduce(rbind.data.frame, rs, accumulate = FALSE)
###
write_xlsx(med_rs, path = paste0(path_output, "med_LC3AB_Bra_Cog_results_bl.xlsx"))
