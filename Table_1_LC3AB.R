###
path_data <- "G:\\Project\\ADNI_LC3AB\\Data\\"
path_adni <- "G:\\Project\\Cohort\\ADNI\\Original\\"
path_model <- "G:\\Project\\Model_R\\HHG\\"
path_output <- "G:\\Project\\ADNI_LC3AB\\Results\\"
#
library(data.table)
library(dplyr)
library(ztable)
library(moonBook)
library(magrittr)
library(officer)
library(dplyr)
library(flextable)
library(data.table)
###
LC3 <- read.csv(paste0(path_data,"LC3AB_p62.csv"))
LC3 <- subset(LC3, !is.na(LC3$LC3AB))
LC3$LC3AB <- scale(LC3$LC3AB)
###
LC3AB_median <- median(LC3$LC3AB, na.rm = TRUE)
LC3$LC3AB_group <- ifelse(LC3$LC3AB >= LC3AB_median, 1, 0)
table(LC3$LC3AB_group)
###
PGRN_STREM2 <- read.csv(paste0(path_adni,"PGRN_STREM2_used.csv"))
PHC <- read.csv(paste0(path_adni,"ADSP_PHC_COGN_10_05_22_29Mar2025.csv"))
PET <- read.csv(paste0(path_adni,"PET_Stage.csv"))
demo_bl <- read.csv(paste0(path_adni,"demo_bl_blnoAD.csv"))
lifestyle_med <- fread(paste0(path_adni,"lifestyle_med.csv"))
PET_TAU_processed_tobl <- read.csv(paste0(path_adni,"PET_TAU_processed_tobl.csv"))
plasma_nfl_bl <- subset(read.csv(paste0(path_adni,"ADNI_BLENNOWPLASMANFLLONG_10_03_18_19Apr2026.csv")),VISCODE2 == "bl")
###
LC3_bl_demo <- merge(LC3, demo_bl, by = c("RID","VISCODE2"))
merge_df0 <- merge(LC3_bl_demo, PHC, by = c("RID","VISCODE2"),all.x=T)
merge_df1 <- merge(merge_df0, PET, by = c("RID","VISCODE2"),all.x=T)
merge_df2 <- merge(merge_df1, PGRN_STREM2, by = c("RID","VISCODE2"),all.x=T)
merge_df3 <- merge(merge_df2, lifestyle_med, by = c("RID","VISCODE2"),all.x=T)
merge_df4 <- merge(merge_df3, PET_TAU_processed_tobl, by = c("RID","VISCODE2"),all.x=T)
merge_df5 <- merge(merge_df4, plasma_nfl_bl, by = c("RID","VISCODE2"),all.x=T)
###
###
merge_df4 <- subset(merge_df3,
                    !is.na(merge_df3$AGE) & 
                      !is.na(merge_df3$PTGENDER_C) & 
                      !is.na(merge_df3$PTEDUCAT) & 
                      !is.na(merge_df3$APOE4_C) &
                      !is.na(merge_df3$ABETA) &
                      !is.na(merge_df3$LC3B) &
                      !is.na(merge_df3$LC3A) &
                      !is.na(merge_df3$PTAU) &
                      !is.na(merge_df3$TAU) &
                      !is.na(merge_df3$Years_bl)
)
###
Characteristics <- merge_df5[,c("LC3AB_group","AGE","PTGENDER_C","PTEDUCAT","APOE4_C",
                               "ABETA","PTAU","TAU","FDG","AV45","META_TEMPORAL_SUVR_tobl","MSD_PGRNCORRECTED","MSD_STREM2CORRECTED","PLASMA_NFL",
                               "ADAS13","CDRSB","MMSE","MOCA","PHC_EXF","PHC_LAN","PHC_MEM","PHC_VSP",
                               "Entorhinal","Fusiform","Hippocampus","MidTemp","Ventricles","WholeBrain")]
sapply(Characteristics[, c(1:28)], function(x) sum(is.na(x) | x == ""))
###
Characteristics <- Characteristics %>% mutate_at(vars(2,4,6:28), as.numeric)
Characteristics <- Characteristics %>% mutate_at(vars(1,3,5), as.character)
### Group compute
my_table <- mytable(LC3AB_group~AGE+PTGENDER_C+PTEDUCAT+APOE4_C+
                        ABETA+PTAU+TAU+FDG+AV45+META_TEMPORAL_SUVR_tobl+MSD_PGRNCORRECTED+MSD_STREM2CORRECTED+PLASMA_NFL+
                      ADAS13+CDRSB+MMSE+MOCA+PHC_EXF+PHC_LAN+PHC_MEM+PHC_VSP+
                      Entorhinal+Fusiform+Hippocampus+MidTemp+Ventricles+WholeBrain,data=Characteristics,digits=2)
my_table_df <- mytable2df(my_table)
### overall
Characteristics$Overall[!is.na(Characteristics$LC3AB_group)] <- 1
merge1 <- Characteristics
Characteristics_all <- merge1[,c("Overall","AGE","PTGENDER_C","PTEDUCAT","APOE4_C",
                                            "ABETA","PTAU","TAU","FDG","AV45","META_TEMPORAL_SUVR_tobl","MSD_PGRNCORRECTED","MSD_STREM2CORRECTED","PLASMA_NFL",
                                            "ADAS13","CDRSB","MMSE","MOCA","PHC_EXF","PHC_LAN","PHC_MEM","PHC_VSP",
                                            "Entorhinal","Fusiform","Hippocampus","MidTemp","Ventricles","WholeBrain")]
###
Characteristics_all <- Characteristics_all %>%
  mutate_at(vars(2,4,6:28), as.numeric)
###
Characteristics_all <- Characteristics_all %>%
  mutate_at(vars(1,3,5), as.character)
###
my_table_all <- mytable(Overall~AGE+PTGENDER_C+PTEDUCAT+APOE4_C+
                          ABETA+PTAU+TAU+FDG+AV45+META_TEMPORAL_SUVR_tobl+MSD_PGRNCORRECTED+MSD_STREM2CORRECTED+PLASMA_NFL+
                          ADAS13+CDRSB+MMSE+MOCA+PHC_EXF+PHC_LAN+PHC_MEM+PHC_VSP+
                          Entorhinal+Fusiform+Hippocampus+MidTemp+Ventricles+WholeBrain,data=Characteristics,digits=2)
my_table_df_all <- mytable2df(my_table_all)
###
my_table_df$Overall <- my_table_df_all[c(2)]
my_table_df <- my_table_df[c(1,5,2:4)]
colnames(my_table_df[,2]) <- "Overall"
### output
my_ztable <- ztable_sub(my_table_df,include.rownames = FALSE)
options(ztable.type="html")
my_ztable_ft<-ztable2flextable(my_ztable)
read_docx() %>%
  body_add_flextable(my_ztable_ft) %>%
  print(target = paste0(path_output,"Table_11.docx"))




