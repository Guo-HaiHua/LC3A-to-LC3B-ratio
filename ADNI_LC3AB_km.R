###
path_data <- "G:\\Project\\ADNI_LC3AB\\Data\\"
path_adni <- "G:\\Project\\Cohort\\ADNI\\Original\\"
path_model <- "G:\\Project\\Paper_new\\CodeBank\\Model\\"
path_output <- "G:\\Project\\ADNI_LC3AB\\Results\\"
###
library(data.table)
library(dplyr)
library(tidyr)
library(survminer)
library(survival)
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
lifestyle_med <- read.csv(paste0(path_adni,"lifestyle_med.csv"))
###
LC3_bl_demo <- merge(LC3, demo, by = c("RID"),all.x=T)
names(LC3_bl_demo)[14] <- "VISCODE2"
merge_df0 <- merge(LC3_bl_demo, PHC, by = c("RID","VISCODE2"),all.x=T)
merge_df1 <- merge(merge_df0, PET, by = c("RID","VISCODE2"),all.x=T)
merge_df2 <- merge(merge_df1, PGRN_STREM2, by = c("RID","VISCODE2"),all.x=T)
merge_df3 <- merge(merge_df2, lifestyle_med, by = c("RID","VISCODE2"),all.x=T)
###
merge_df4 <- subset(merge_df3,!is.na(merge_df3$AGE) & !is.na(merge_df3$PTGENDER_C) & !is.na(merge_df3$PTEDUCAT) & !is.na(merge_df3$APOE4_C) & 
                      merge_df3$DX != "" & !is.na(merge_df3$Years_bl))
###
case <- merge_df4 %>%
  filter(DX == "Dementia") %>%
  group_by(RID) %>%
  slice(which.min(Years_bl)) %>%
  ungroup()

#
control <- merge_df4 %>%
  # 去掉曾经出现过 Dementia 的 RID
  group_by(RID) %>%
  filter(!any(DX == "Dementia")) %>%
  # 取每个 RID 最大 Years_bl
  slice(which.max(Years_bl)) %>%
  ungroup()

merge_df <-as.data.frame(rbind(control, case))
# merge_df$DX[merge_df$DX %in% c("CN", "MCI")] <- "noDementia"
###
LC3AB_median <- median(merge_df$LC3AB, na.rm = TRUE)
merge_df$LC3AB_group <- ifelse(merge_df$LC3AB >= LC3AB_median, 1, 0)
table(merge_df$LC3AB_group)
###
p1 <- suppressWarnings(
  ggsurvplot(
    survfit(Surv(Years_bl, DX_C) ~ LC3AB_group, data = merge_df),
    palette = c("#005b9c", "#d62728"),
    risk.table = TRUE,
    conf.int = T,
    linewidth = 1.3,
    pval = TRUE,
    pval.method = TRUE, 
    risk.table.fontsize = 5,
    xlab = "Time (years, model 1)",
    ylab = "Dementia-free (%)",
    legend.labs = c("Low", "High"),
    legend.title = "LC3A/LC3B ratio",
    ggtheme = theme_bw() + 
      theme(plot.title =element_text(size = 16),
            axis.text.x=element_text(size=16,face = "bold"),
            axis.text.y=element_text(size=18,face = "bold"),
            axis.title.x =element_text(size=16),
            axis.title.y=element_text(size=16),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.border = element_blank(),#去掉图片边框
            axis.line = element_line(linewidth = 0.7)) # 去掉次网格线
  )
)
p1

###
merge_dff <- subset(merge_df, Years_bl >= 2)
p2 <- suppressWarnings(
  ggsurvplot(
    survfit(Surv(Years_bl, DX_C) ~ LC3AB_group, data = merge_dff),
    palette = c("#005b9c", "#d62728"),
    risk.table = TRUE,
    conf.int = T,
    linewidth = 1.3,
    pval = TRUE,
    pval.method = TRUE, 
    risk.table.fontsize = 5,
    xlab = "Time (years, model 1)",
    ylab = "Dementia-free (%)",
    legend.labs = c("Low", "High"),
    legend.title = "LC3A/LC3B ratio",
    ggtheme = theme_bw() + 
      theme(plot.title =element_text(size = 16),
            axis.text.x=element_text(size=16,face = "bold"),
            axis.text.y=element_text(size=18,face = "bold"),
            axis.title.x =element_text(size=16),
            axis.title.y=element_text(size=16),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.border = element_blank(),#去掉图片边框
            axis.line = element_line(linewidth = 0.7)) # 去掉次网格线
  )
)
p2
###
###
merge_df9 <- subset(merge_df, !is.na(SQSTM1_S2))
p3 <- suppressWarnings(
  ggsurvplot(
    survfit(Surv(Years_bl, DX_C) ~ LC3AB_group, data = merge_df9),
    palette = c("#005b9c", "#d62728"),
    risk.table = TRUE,
    conf.int = T,
    linewidth = 1.3,
    pval = TRUE,
    pval.method = TRUE, 
    risk.table.fontsize = 5,
    xlab = "Time (years, model 2)",
    ylab = "Dementia-free (%)",
    legend.labs = c("Low", "High"),
    legend.title = "LC3A/LC3B ratio",
    ggtheme = theme_bw() + 
      theme(plot.title =element_text(size = 16),
            axis.text.x=element_text(size=16,face = "bold"),
            axis.text.y=element_text(size=18,face = "bold"),
            axis.title.x =element_text(size=16),
            axis.title.y=element_text(size=16),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.border = element_blank(),#去掉图片边框
            axis.line = element_line(linewidth = 0.7)) # 去掉次网格线
  )
)
p3

###
merge_dff9 <- subset(merge_df9, Years_bl >= 2)
p4 <- suppressWarnings(
  ggsurvplot(
    survfit(Surv(Years_bl, DX_C) ~ LC3AB_group, data = merge_dff9),
    palette = c("#005b9c", "#d62728"),
    risk.table = TRUE,
    conf.int = T,
    linewidth = 1.3,
    pval = TRUE,
    pval.method = TRUE, 
    risk.table.fontsize = 5,
    xlab = "Time (years, model 2)",
    ylab = "Dementia-free (%)",
    legend.labs = c("Low", "High"),
    legend.title = "LC3A/LC3B ratio",
    ggtheme = theme_bw() + 
      theme(plot.title =element_text(size = 16),
            axis.text.x=element_text(size=16,face = "bold"),
            axis.text.y=element_text(size=18,face = "bold"),
            axis.title.x =element_text(size=16),
            axis.title.y=element_text(size=16),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 16),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            # panel.border = element_blank(),#去掉图片边框
            axis.line = element_line(linewidth = 0.7)) # 去掉次网格线
  )
)
p4
###
p1_all <- ggarrange(p1$plot, p1$table, ncol = 1, heights = c(0.8, 0.2))
p2_all <- ggarrange(p2$plot, p2$table, ncol = 1, heights = c(0.8, 0.2))
p3_all <- ggarrange(p3$plot, p3$table, ncol = 1, heights = c(0.8, 0.2))
p4_all <- ggarrange(p4$plot, p4$table, ncol = 1, heights = c(0.8, 0.2))
###
# p1$plot <- p1$plot + labs(color = "LC3A")
# p2$plot <- p2$plot + labs(color = "LC3B")
# p3$plot <- p3$plot + labs(color = "LC3A")
# p4$table <- p4$table + labs(y = "LC3B")
# p4


# library(gridExtra)
# all <- grid.arrange(p1, p2$plot, p3$plot, p8$plot, p4$plot, p5$plot, p6$plot, p7$plot,ncol = 4,nrow=2)
# 
# 
# ggsave(all,file = paste0(path_out,"KM_new_RT.pdf"),
#         width = 18, height = 10, dpi = 800)


# 提取ggplot对象
plot_list <- list(p1$plot, p2$plot, p3$plot, p4$plot)
table_list <- list(p1$table, p2$table, p3$table, p4$table)


# 组合生存曲线和风险表为单个列表
library(cowplot)
combined_plots <- list()
for(i in 1:4) {
  combined_plots[[i]] <- plot_grid(plot_list[[i]], table_list[[i]], ncol = 1, rel_heights = c(2.5, 1))
}

# 组合图表为2行4列
final_plot <- plot_grid(plotlist = combined_plots, ncol = 2, nrow = 2,
                        labels = LETTERS[1:4],
                        label_size = 20)
final_plot
library(ggpubr)
ggsave(final_plot,file = paste0(path_output,"km_LC3AB1.pdf"),
        width = 10, height = 12, dpi = 300)



