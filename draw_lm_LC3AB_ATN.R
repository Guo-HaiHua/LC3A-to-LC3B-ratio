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

a <- lm(LC3A ~ LC3B, data = LC3)
summary(a)
#
LC3 <- subset(LC3, !is.na(LC3$LC3AB))
LC3$LC3AB <- scale(LC3$LC3AB)
###
PGRN_STREM2 <- read.csv(paste0(path_adni,"PGRN_STREM2_used.csv"))
PHC <- read.csv(paste0(path_adni,"ADSP_PHC_COGN_10_05_22_29Mar2025.csv"))
PET <- read.csv(paste0(path_adni,"PET.csv"))
demo_bl <- read.csv(paste0(path_adni,"demo_bl_blnoAD.csv"))
PET_TAU_processed_tobl <- read.csv(paste0(path_adni,"PET_TAU_processed_tobl.csv"))
plasma_nfl_bl <- subset(read.csv(paste0(path_adni,"ADNI_BLENNOWPLASMANFLLONG_10_03_18_15Mar2026.csv")),VISCODE2 == "bl")
###
LC3_bl_demo <- merge(LC3, demo_bl, by = c("RID","VISCODE2"))
merge_df0 <- merge(LC3_bl_demo, PHC, by = c("RID","VISCODE2"),all.x=T)
merge_df1 <- merge(merge_df0, PET, by = c("RID","VISCODE2"),all.x=T)
merge_df2 <- merge(merge_df1, PGRN_STREM2, by = c("RID","VISCODE2"),all.x=T)
merge_df3 <- merge(merge_df2, PET_TAU_processed_tobl, by = c("RID","VISCODE2"),all.x=T)
merge_df <- merge(merge_df3, plasma_nfl_bl, by = c("RID","VISCODE2"),all.x=T)

##############################################################################################
###
merge_df$ATN_Stage[(merge_df$TRACER.ABETA=="FBB" & merge_df$SUMMARY_SUVR <1.08 | merge_df$TRACER.ABETA=="FBP" & merge_df$SUMMARY_SUVR <1.32) & merge_df$PTAU < 21.8 & merge_df$TAU < 245] <- "Stage 0" #A-T-N-
merge_df$ATN_Stage[(merge_df$TRACER.ABETA=="FBB" & merge_df$SUMMARY_SUVR <1.08 | merge_df$TRACER.ABETA=="FBP" & merge_df$SUMMARY_SUVR <1.32) & (merge_df$PTAU >= 21.8 | merge_df$TAU >= 245)] <- "SNAP" #A-T+(N)-, A-T-(N)+, or A-T+(N)+
merge_df$ATN_Stage[(merge_df$TRACER.ABETA=="FBB" & merge_df$SUMMARY_SUVR >=1.08 | merge_df$TRACER.ABETA=="FBP" & merge_df$SUMMARY_SUVR >=1.32) & merge_df$PTAU < 21.8 & merge_df$TAU < 245] <- "Stage 1" #A+T-N-
merge_df$ATN_Stage[(merge_df$TRACER.ABETA=="FBB" & merge_df$SUMMARY_SUVR >=1.08 | merge_df$TRACER.ABETA=="FBP" & merge_df$SUMMARY_SUVR >=1.32) & (merge_df$PTAU >= 21.8 | merge_df$TAU >= 245)] <- "Stage 2" #A+T+N- #A+T+N+ #A+T-N+
###
merge_df$CSNAP_S0[merge_df$Stage=="SNAP"] <- 0
merge_df$CSNAP_S0[merge_df$Stage=="Stage 0"] <- 1
merge_df$CSNAP_S1[merge_df$Stage=="SNAP"] <- 0
merge_df$CSNAP_S1[merge_df$Stage=="Stage 1"] <- 1
merge_df$CSNAP_S2[merge_df$Stage=="SNAP"] <- 0
merge_df$CSNAP_S2[merge_df$Stage=="Stage 2"] <- 1
merge_df$CS0_S1[merge_df$Stage=="Stage 0"] <- 0
merge_df$CS0_S1[merge_df$Stage=="Stage 1"] <- 1
merge_df$CS0_S2[merge_df$Stage=="Stage 0"] <- 0
merge_df$CS0_S2[merge_df$Stage=="Stage 2"] <- 1
merge_df$CS1_S2[merge_df$Stage=="Stage 1"] <- 0
merge_df$CS1_S2[merge_df$Stage=="Stage 2"] <- 1
###
merge_df$SNAP_S0[merge_df$ATN_Stage=="SNAP"] <- 0
merge_df$SNAP_S0[merge_df$ATN_Stage=="Stage 0"] <- 1
merge_df$SNAP_S1[merge_df$ATN_Stage=="SNAP"] <- 0
merge_df$SNAP_S1[merge_df$ATN_Stage=="Stage 1"] <- 1
merge_df$SNAP_S2[merge_df$ATN_Stage=="SNAP"] <- 0
merge_df$SNAP_S2[merge_df$ATN_Stage=="Stage 2"] <- 1
merge_df$S0_S1[merge_df$ATN_Stage=="Stage 0"] <- 0
merge_df$S0_S1[merge_df$ATN_Stage=="Stage 1"] <- 1
merge_df$S0_S2[merge_df$ATN_Stage=="Stage 0"] <- 0
merge_df$S0_S2[merge_df$ATN_Stage=="Stage 2"] <- 1
merge_df$S1_S2[merge_df$ATN_Stage=="Stage 1"] <- 0
merge_df$S1_S2[merge_df$ATN_Stage=="Stage 2"] <- 1
###
merge_df$A_Stage[merge_df$Stage=="Stage 0" | merge_df$Stage=="SNAP"] <- "A- group"
merge_df$A_Stage[merge_df$Stage=="Stage 1" | merge_df$Stage=="Stage 2"] <- "A+ group"
merge_df$T_Stage[merge_df$Stage=="Stage 0" | merge_df$Stage=="Stage 1"] <- "T- group"
merge_df$T_Stage[merge_df$Stage=="SNAP" | merge_df$Stage=="Stage 2"] <- "T+ group"
###
merge_df9 <- subset(merge_df, !is.na(SQSTM1_S2))
table(merge_df9$A_Stage)
# ---- 准备包 ----
library(ggplot2)
library(dplyr)
library(gghalves)
library(ggpubr)
# ###
o1 <- ggplot(merge_df) +
  geom_boxplot(aes(x = A_Stage, y = LC3AB,color = A_Stage), width = 0.3, alpha = 1, outlier.shape = NA) +
  geom_half_violin(aes(x = A_Stage, y = LC3AB, fill = A_Stage),
                   position = position_nudge(x=0.2), #箱式图和小提琴图间隔
                   width = 0.6, side = "r", alpha = 0.8, trim = F, color = NA) +
  geom_jitter(aes(x = A_Stage, y = LC3AB, color = A_Stage),
              width = 0.1, alpha = 0.8, size = 2,
              # shape = 15 #实菱形
              shape = 16 #实正方形
              # shape = 17 #实三角形
              # shape = 18 #实菱形
  ) +
  scale_fill_manual(values = c("A- group" = "#354e55",  "A+ group" = "#d95725")) +
  scale_color_manual(values = c("A- group" = "#354e55",  "A+ group" = "#d95725")) +
  # coord_cartesian(xlim = c(-5, 6)) + #设定X轴范围
  # coord_cartesian(ylim = c(-10, -2)) + #设定Y轴范围
  geom_signif(aes(x = A_Stage, y = LC3AB),
              comparisons=list(c("A- group","A+ group")),
              annotations =c("***"),#手动添加显著性*
              vjust=0,# 星号与比较线上下距离
              textsize=6, # 星号大小
              color="black",
              size = 0.3, #比较线粗细
              step_increase = 0.1, # 两条比较线上下距离
              tip_length = 0 # 比较横线的两端短线长度为0
              # tip_length = rep(0.03, 5) # 4条横线的两端短线长度为0.03
  )+
  geom_segment(data = data.frame(xstart = c(0.6,1.6,  0.6,1.6,  1.4,2.4), # X轴坐标分组化
                                 xend   = c(1.4,2.4,  0.6,1.6,  1.4,2.4),
                                 ystart = c(-5.6,-5.6,  -5.6,-5.6,  -5.6,-5.6),
                                 yend   = c(-5.6,-5.6,  -5.8,-5.8,  -5.8,-5.8)),
               aes(x=xstart, xend=xend, y=ystart, yend=yend), size=0.6) +
  # coord_cartesian(ylim = c(-3, 3.5)) + # 限制显示范围，但不裁剪图中数据
  scale_x_discrete(
    limits = c("A- group", "A+ group"),  # 原始分类顺序
    labels = c("A- group \nn = 169", "A+ group \nn = 297")) + # 新显示名称
  labs(
    # title = paste0("  β = 0.073 *  n = 690"), 
    x = "CSF Aβ42 (model 1)", y = "LC3A/LC3B ratio", fill = "A_Stage",color = "A_Stage") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),   # 去掉网格
    panel.border = element_blank(), # 去掉框线
    # axis.line.x = element_line(color = "black", linewidth = 0.6), # 隐藏x 轴
    axis.ticks.x = element_blank(),  # 隐藏 X 轴所有刻度线
    axis.line.y = element_line(color = "black", linewidth = 0.6), # 只画 y 轴
    # axis.line = element_line(color = "black", linewidth = 0.6), # 加粗坐标轴
    axis.ticks = element_line(color = "black", linewidth = 0.6), # 坐标轴刻度线
    axis.title.y = element_text(margin = margin(r = 0)), #Y轴标签和轴距离
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(margin = margin(t = -10)) ,  # X轴分组标签向上移动
    strip.text = element_text(size = 16)
  )
o1

###
o2 <- ggplot(merge_df9) +
  geom_boxplot(aes(x = A_Stage, y = LC3AB,color = A_Stage), width = 0.3, alpha = 1, outlier.shape = NA) +
  geom_half_violin(aes(x = A_Stage, y = LC3AB, fill = A_Stage),
                   position = position_nudge(x=0.2), #箱式图和小提琴图间隔
                   width = 0.6, side = "r", alpha = 0.8, trim = F, color = NA) +
  geom_jitter(aes(x = A_Stage, y = LC3AB, color = A_Stage),
              width = 0.1, alpha = 0.8, size = 2,
              # shape = 15 #实菱形
              shape = 17 #实正方形
              # shape = 17 #实三角形
              # shape = 18 #实菱形
  ) +
  scale_fill_manual(values = c("A- group" = "#354e55",  "A+ group" = "#d95725")) +
  scale_color_manual(values = c("A- group" = "#354e55",  "A+ group" = "#d95725")) +
  # coord_cartesian(xlim = c(-5, 6)) + #设定X轴范围
  # coord_cartesian(ylim = c(-10, -2)) + #设定Y轴范围
  geom_signif(aes(x = A_Stage, y = LC3AB),
              comparisons=list(c("A- group","A+ group")),
              annotations =c("***"),#手动添加显著性*
              vjust=0,# 星号与比较线上下距离
              textsize=6, # 星号大小
              color="black",
              size = 0.3, #比较线粗细
              step_increase = 0.1, # 两条比较线上下距离
              tip_length = 0 # 比较横线的两端短线长度为0
              # tip_length = rep(0.03, 5) # 4条横线的两端短线长度为0.03
  )+
  geom_segment(data = data.frame(xstart = c(0.6,1.6,  0.6,1.6,  1.4,2.4), # X轴坐标分组化
                                 xend   = c(1.4,2.4,  0.6,1.6,  1.4,2.4),
                                 ystart = c(-5.6,-5.6,  -5.6,-5.6,  -5.6,-5.6),
                                 yend   = c(-5.6,-5.6,  -5.8,-5.8,  -5.8,-5.8)),
               aes(x=xstart, xend=xend, y=ystart, yend=yend), size=0.6) +
  # coord_cartesian(ylim = c(-3, 3.5)) + # 限制显示范围，但不裁剪图中数据
  scale_x_discrete(
    limits = c("A- group", "A+ group"),  # 原始分类顺序
    labels = c("A- group \nn = 166", "A+ group \nn = 287")) + # 新显示名称
  labs(
    # title = paste0("  β = 0.073 *  n = 690"), 
    x = "CSF Aβ42 (model 2)", y = "LC3A/LC3B ratio", fill = "A_Stage",color = "A_Stage") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),   # 去掉网格
    panel.border = element_blank(), # 去掉框线
    # axis.line.x = element_line(color = "black", linewidth = 0.6), # 隐藏x 轴
    axis.ticks.x = element_blank(),  # 隐藏 X 轴所有刻度线
    axis.line.y = element_line(color = "black", linewidth = 0.6), # 只画 y 轴
    # axis.line = element_line(color = "black", linewidth = 0.6), # 加粗坐标轴
    axis.ticks = element_line(color = "black", linewidth = 0.6), # 坐标轴刻度线
    axis.title.y = element_text(margin = margin(r = 0)), #Y轴标签和轴距离
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(margin = margin(t = -10)) ,  # X轴分组标签向上移动
    strip.text = element_text(size = 16)
  )
o2


###
o3 <- ggplot(merge_df) +
  geom_boxplot(aes(x = T_Stage, y = LC3AB,color = T_Stage), width = 0.3, alpha = 1, outlier.shape = NA) +
  geom_half_violin(aes(x = T_Stage, y = LC3AB, fill = T_Stage),
                   position = position_nudge(x=0.2), #箱式图和小提琴图间隔
                   width = 0.6, side = "r", alpha = 0.8, trim = F, color = NA) +
  geom_jitter(aes(x = T_Stage, y = LC3AB, color = T_Stage),
              width = 0.1, alpha = 0.8, size = 2,
              # shape = 15 #实菱形
              shape = 16 #实正方形
              # shape = 17 #实三角形
              # shape = 18 #实菱形
  ) +
  scale_fill_manual(values = c("T- group" = "#225081",  "T+ group" = "#6e1a00")) +
  scale_color_manual(values = c("T- group" = "#225081",  "T+ group" = "#6e1a00")) +
  # coord_cartesian(xlim = c(-5, 6)) + #设定X轴范围
  # coord_cartesian(ylim = c(-10, -2)) + #设定Y轴范围
  geom_signif(aes(x = T_Stage, y = LC3AB),
              comparisons=list(c("T- group","T+ group")),
              annotations =c("n.s."),#手动添加显著性*
              vjust=0,# 星号与比较线上下距离
              textsize=6, # 星号大小
              color="black",
              size = 0.3, #比较线粗细
              step_increase = 0.1, # 两条比较线上下距离
              tip_length = 0 # 比较横线的两端短线长度为0
              # tip_length = rep(0.03, 5) # 4条横线的两端短线长度为0.03
  )+
  geom_segment(data = data.frame(xstart = c(0.6,1.6,  0.6,1.6,  1.4,2.4), # X轴坐标分组化
                                 xend   = c(1.4,2.4,  0.6,1.6,  1.4,2.4),
                                 ystart = c(-5.6,-5.6,  -5.6,-5.6,  -5.6,-5.6),
                                 yend   = c(-5.6,-5.6,  -5.8,-5.8,  -5.8,-5.8)),
               aes(x=xstart, xend=xend, y=ystart, yend=yend), size=0.6) +
  # coord_cartesian(ylim = c(-3, 3.5)) + # 限制显示范围，但不裁剪图中数据
  scale_x_discrete(
    limits = c("T- group", "T+ group"),  # 原始分类顺序
    labels = c("T- group \nn = 169", "T+ group \nn = 297")) + # 新显示名称
  labs(
    # title = paste0("  β = 0.073 *  n = 690"), 
    x = "CSF t-tau, p-tau (model 1)", y = "LC3A/LC3B ratio", fill = "T_Stage",color = "T_Stage") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),   # 去掉网格
    panel.border = element_blank(), # 去掉框线
    # axis.line.x = element_line(color = "black", linewidth = 0.6), # 隐藏x 轴
    axis.ticks.x = element_blank(),  # 隐藏 X 轴所有刻度线
    axis.line.y = element_line(color = "black", linewidth = 0.6), # 只画 y 轴
    # axis.line = element_line(color = "black", linewidth = 0.6), # 加粗坐标轴
    axis.ticks = element_line(color = "black", linewidth = 0.6), # 坐标轴刻度线
    axis.title.y = element_text(margin = margin(r = 0)), #Y轴标签和轴距离
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(margin = margin(t = -10)) ,  # X轴分组标签向上移动
    strip.text = element_text(size = 16)
  )
o3

###
o4 <- ggplot(merge_df9) +
  geom_boxplot(aes(x = T_Stage, y = LC3AB,color = T_Stage), width = 0.3, alpha = 1, outlier.shape = NA) +
  geom_half_violin(aes(x = T_Stage, y = LC3AB, fill = T_Stage),
                   position = position_nudge(x=0.2), #箱式图和小提琴图间隔
                   width = 0.6, side = "r", alpha = 0.8, trim = F, color = NA) +
  geom_jitter(aes(x = T_Stage, y = LC3AB, color = T_Stage),
              width = 0.1, alpha = 0.8, size = 2,
              # shape = 15 #实菱形
              shape = 17 #实正方形
              # shape = 17 #实三角形
              # shape = 18 #实菱形
  ) +
  scale_fill_manual(values = c("T- group" = "#225081",  "T+ group" = "#6e1a00")) +
  scale_color_manual(values = c("T- group" = "#225081",  "T+ group" = "#6e1a00")) +
  # coord_cartesian(xlim = c(-5, 6)) + #设定X轴范围
  # coord_cartesian(ylim = c(-10, -2)) + #设定Y轴范围
  geom_signif(aes(x = T_Stage, y = LC3AB),
              comparisons=list(c("T- group","T+ group")),
              annotations =c("n.s."),#手动添加显著性*
              vjust=0,# 星号与比较线上下距离
              textsize=6, # 星号大小
              color="black",
              size = 0.3, #比较线粗细
              step_increase = 0.1, # 两条比较线上下距离
              tip_length = 0 # 比较横线的两端短线长度为0
              # tip_length = rep(0.03, 5) # 4条横线的两端短线长度为0.03
  )+
  geom_segment(data = data.frame(xstart = c(0.6,1.6,  0.6,1.6,  1.4,2.4), # X轴坐标分组化
                                 xend   = c(1.4,2.4,  0.6,1.6,  1.4,2.4),
                                 ystart = c(-5.6,-5.6,  -5.6,-5.6,  -5.6,-5.6),
                                 yend   = c(-5.6,-5.6,  -5.8,-5.8,  -5.8,-5.8)),
               aes(x=xstart, xend=xend, y=ystart, yend=yend), size=0.6) +
  # coord_cartesian(ylim = c(-3, 3.5)) + # 限制显示范围，但不裁剪图中数据
  scale_x_discrete(
    limits = c("T- group", "T+ group"),  # 原始分类顺序
    labels = c("T- group \nn = 166", "T+ group \nn = 287")) + # 新显示名称
  labs(
    # title = paste0("  β = 0.073 *  n = 690"), 
    x = "CSF t-tau, p-tau (model 2)", y = "LC3A/LC3B ratio", fill = "T_Stage",color = "T_Stage") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),   # 去掉网格
    panel.border = element_blank(), # 去掉框线
    # axis.line.x = element_line(color = "black", linewidth = 0.6), # 隐藏x 轴
    axis.ticks.x = element_blank(),  # 隐藏 X 轴所有刻度线
    axis.line.y = element_line(color = "black", linewidth = 0.6), # 只画 y 轴
    # axis.line = element_line(color = "black", linewidth = 0.6), # 加粗坐标轴
    axis.ticks = element_line(color = "black", linewidth = 0.6), # 坐标轴刻度线
    axis.title.y = element_text(margin = margin(r = 0)), #Y轴标签和轴距离
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(margin = margin(t = -10)) ,  # X轴分组标签向上移动
    strip.text = element_text(size = 16)
  )
o4

#### 作图
a1 <- ggplot(merge_df) +
  geom_boxplot(aes(x = ATN_Stage, y = LC3AB,color = ATN_Stage), width = 0.3, alpha = 1, outlier.shape = NA) +
  geom_half_violin(aes(x = ATN_Stage, y = LC3AB, fill = ATN_Stage),
                   position = position_nudge(x=0.2), #箱式图和小提琴图间隔
                   width = 0.6, side = "r", alpha = 0.8, trim = F, color = NA) +
  geom_jitter(aes(x = ATN_Stage, y = LC3AB, color = ATN_Stage),
              width = 0.1, alpha = 0.8, size = 2,
              # shape = 15 #实菱形
              shape = 16 #实正方形
              # shape = 17 #实三角形
              # shape = 18 #实菱形
  ) +
  scale_fill_manual(values = c("SNAP" = "#354e55", "Stage 0" = "#225081", "Stage 1" = "#d95725", "Stage 2" = "#6e1a00")) +
  scale_color_manual(values = c("SNAP" = "#354e55", "Stage 0" = "#225081", "Stage 1" = "#d95725", "Stage 2" = "#6e1a00")) +
  # coord_cartesian(xlim = c(-5, 6)) + #设定X轴范围
  # coord_cartesian(ylim = c(-10, -2)) + #设定Y轴范围
  geom_signif(aes(x = ATN_Stage, y = LC3AB),
              comparisons=list(c("Stage 0","Stage 1"), c("Stage 0","Stage 2"), c("Stage 1","Stage 2"),
                               c("SNAP","Stage 0"), c("SNAP","Stage 1"), c("SNAP","Stage 2")),
              annotations =c("n.s.","n.s.","n.s.","n.s.","n.s.","*"),#手动添加显著性*
              vjust=0,# 星号与比较线上下距离
              textsize=6, # 星号大小
              color="black",
              size = 0.3, #比较线粗细
              step_increase = 0.1, # 两条比较线上下距离
              tip_length = 0 # 比较横线的两端短线长度为0
              # tip_length = rep(0.03, 5) # 4条横线的两端短线长度为0.03
  )+
  geom_segment(data = data.frame(xstart = c(0.6,1.6,2.6,3.6,  0.6,1.6,2.6,3.6,  1.4,2.4,3.4,4.4), # X轴坐标分组化
                                 xend   = c(1.4,2.4,3.4,4.4,  0.6,1.6,2.6,3.6,  1.4,2.4,3.4,4.4),
                                 ystart = c(-5.6,-5.6,-5.6,-5.6,  -5.6,-5.6,-5.6,-5.6,  -5.6,-5.6,-5.6,-5.6),
                                 yend   = c(-5.6,-5.6,-5.6,-5.6,  -5.8,-5.8,-5.8,-5.8,  -5.8,-5.8,-5.8,-5.8)),
               aes(x=xstart, xend=xend, y=ystart, yend=yend), size=0.6) +
  # coord_cartesian(ylim = c(-3, 3.5)) + # 限制显示范围，但不裁剪图中数据
  scale_x_discrete(
    limits = c("SNAP", "Stage 0", "Stage 1", "Stage 2"),  # 原始分类顺序
    labels = c("SNAP \nn = 101", "Stage 0 \nn = 163", "Stage 1 \nn = 25", "Stage 2 \nn = 101")) + # 新显示名称
  labs(
    # title = paste0("  β = 0.073 *  n = 690"), 
    x = "Amyloid PET, CSF t-tau, p-tau (model 1)", y = "LC3A/LC3B ratio", fill = "ATN_Stage",color = "ATN_Stage") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),   # 去掉网格
    panel.border = element_blank(), # 去掉框线
    # axis.line.x = element_line(color = "black", linewidth = 0.6), # 隐藏x 轴
    axis.ticks.x = element_blank(),  # 隐藏 X 轴所有刻度线
    axis.line.y = element_line(color = "black", linewidth = 0.6), # 只画 y 轴
    # axis.line = element_line(color = "black", linewidth = 0.6), # 加粗坐标轴
    axis.ticks = element_line(color = "black", linewidth = 0.6), # 坐标轴刻度线
    axis.title.y = element_text(margin = margin(r = 0)), #Y轴标签和轴距离
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(margin = margin(t = -10)) ,  # X轴分组标签向上移动
    strip.text = element_text(size = 16)
  )
a1 

###
a2 <- ggplot(merge_df) +
  geom_boxplot(aes(x = ATN_Stage, y = LC3AB,color = ATN_Stage), width = 0.3, alpha = 1, outlier.shape = NA) +
  geom_half_violin(aes(x = ATN_Stage, y = LC3AB, fill = ATN_Stage),
                   position = position_nudge(x=0.2), #箱式图和小提琴图间隔
                   width = 0.6, side = "r", alpha = 0.8, trim = FALSE, color = NA) +
  geom_jitter(aes(x = ATN_Stage, y = LC3AB, color = ATN_Stage),
              width = 0.1, alpha = 0.8, size = 2,
              # shape = 15 #实菱形
              # shape = 16 #实正方形
              shape = 17 #实三角形
              # shape = 18 #实菱形
  ) +
  scale_fill_manual(values = c("SNAP" = "#354e55", "Stage 0" = "#225081", "Stage 1" = "#d95725", "Stage 2" = "#6e1a00")) +
  scale_color_manual(values = c("SNAP" = "#354e55", "Stage 0" = "#225081", "Stage 1" = "#d95725", "Stage 2" = "#6e1a00")) +
  # coord_cartesian(xlim = c(-5, 6)) + #设定X轴范围
  # coord_cartesian(ylim = c(-10, -2)) + #设定Y轴范围
  geom_signif(aes(x = ATN_Stage, y = LC3AB),
              comparisons=list(c("Stage 0","Stage 1"), c("Stage 0","Stage 2"), c("Stage 1","Stage 2"),
                               c("SNAP","Stage 0"), c("SNAP","Stage 1"), c("SNAP","Stage 2")),
              annotations =c("n.s.","n.s.","n.s.","n.s.","n.s.","*"),#手动添加显著性*
              vjust=0,# 星号与比较线上下距离
              textsize=6, # 星号大小
              color="black",
              size = 0.3, #比较线粗细
              step_increase = 0.1, # 两条比较线上下距离
              tip_length = 0 # 比较横线的两端短线长度为0
              # tip_length = rep(0.03, 5) # 4条横线的两端短线长度为0.03
  )+
  geom_segment(data = data.frame(xstart = c(0.6,1.6,2.6,3.6,  0.6,1.6,2.6,3.6,  1.4,2.4,3.4,4.4), # X轴坐标分组化
                                 xend   = c(1.4,2.4,3.4,4.4,  0.6,1.6,2.6,3.6,  1.4,2.4,3.4,4.4),
                                 ystart = c(-5.6,-5.6,-5.6,-5.6,  -5.6,-5.6,-5.6,-5.6,  -5.6,-5.6,-5.6,-5.6),
                                 yend   = c(-5.6,-5.6,-5.6,-5.6,  -5.8,-5.8,-5.8,-5.8,  -5.8,-5.8,-5.8,-5.8)),
               aes(x=xstart, xend=xend, y=ystart, yend=yend), size=0.6) +
  # coord_cartesian(ylim = c(-3, 3.5)) + # 限制显示范围，但不裁剪图中数据
  scale_x_discrete(
    limits = c("SNAP", "Stage 0", "Stage 1", "Stage 2"),  # 原始分类顺序
    labels = c("SNAP \nn = 101", "Stage 0 \nn = 163", "Stage 1 \nn = 25", "Stage 2 \nn = 101")) + # 新显示名称
    labs(
    # title = paste0("  β = 0.073 *  n = 690"), 
       x = "Amyloid PET, CSF t-tau, p-tau (model 2)", y = "LC3A/LC3B ratio", fill = "ATN_Stage",color = "ATN_Stage") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),   # 去掉网格
    panel.border = element_blank(), # 去掉框线
    # axis.line.x = element_line(color = "black", linewidth = 0.6), # 隐藏x 轴
    axis.ticks.x = element_blank(),  # 隐藏 X 轴所有刻度线
    axis.line.y = element_line(color = "black", linewidth = 0.6), # 只画 y 轴
    # axis.line = element_line(color = "black", linewidth = 0.6), # 加粗坐标轴
    axis.ticks = element_line(color = "black", linewidth = 0.6), # 坐标轴刻度线
    axis.title.y = element_text(margin = margin(r = 0)), #Y轴标签和轴距离
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(margin = margin(t = -10)) ,  # X轴分组标签向上移动
    strip.text = element_text(size = 12)
  )
a2 

###
a3 <- ggplot(merge_df) +
  geom_boxplot(aes(x = Stage, y = LC3AB,color = Stage), width = 0.3, alpha = 1, outlier.shape = NA) +
  geom_half_violin(aes(x = Stage, y = LC3AB, fill = Stage),
                   position = position_nudge(x=0.2), #箱式图和小提琴图间隔
                   width = 0.6, side = "r", alpha = 0.8, trim = FALSE, color = NA) +
  geom_jitter(aes(x = Stage, y = LC3AB, color = Stage),
              width = 0.1, alpha = 0.8, size = 2,
              # shape = 15 #实菱形
              shape = 16 #实正方形
              # shape = 17 #实三角形
              # shape = 18 #实菱形
  ) +
  scale_fill_manual(values = c("SNAP" = "#354e55", "Stage 0" = "#225081", "Stage 1" = "#d95725", "Stage 2" = "#6e1a00")) +
  scale_color_manual(values = c("SNAP" = "#354e55", "Stage 0" = "#225081", "Stage 1" = "#d95725", "Stage 2" = "#6e1a00")) +
  # coord_cartesian(xlim = c(-5, 6)) + #设定X轴范围
  # coord_cartesian(ylim = c(-10, -2)) + #设定Y轴范围
  geom_signif(aes(x = ATN_Stage, y = LC3AB),
              comparisons=list(c("Stage 0","Stage 1"), c("Stage 0","Stage 2"), c("Stage 1","Stage 2"),
                               c("SNAP","Stage 0"), c("SNAP","Stage 1"), c("SNAP","Stage 2")),
              annotations =c("**","***","n.s.","n.s.","*","***"),#手动添加显著性*
              vjust=0,# 星号与比较线上下距离
              textsize=6, # 星号大小
              color="black",
              size = 0.3, #比较线粗细
              step_increase = 0.1, # 两条比较线上下距离
              tip_length = 0 # 比较横线的两端短线长度为0
              # tip_length = rep(0.03, 5) # 4条横线的两端短线长度为0.03
  )+
  geom_segment(data = data.frame(xstart = c(0.6,1.6,2.6,3.6,  0.6,1.6,2.6,3.6,  1.4,2.4,3.4,4.4), # X轴坐标分组化
                                 xend   = c(1.4,2.4,3.4,4.4,  0.6,1.6,2.6,3.6,  1.4,2.4,3.4,4.4),
                                 ystart = c(-5.6,-5.6,-5.6,-5.6,  -5.6,-5.6,-5.6,-5.6,  -5.6,-5.6,-5.6,-5.6),
                                 yend   = c(-5.6,-5.6,-5.6,-5.6,  -5.8,-5.8,-5.8,-5.8,  -5.8,-5.8,-5.8,-5.8)),
               aes(x=xstart, xend=xend, y=ystart, yend=yend), size=0.6) +
  # coord_cartesian(ylim = c(-3, 3.5)) + # 限制显示范围，但不裁剪图中数据
  scale_x_discrete(
    limits = c("SNAP", "Stage 0", "Stage 1", "Stage 2"),  # 原始分类顺序
    labels = c("SNAP \nn = 49", "Stage 0 \nn = 120", "Stage 1 \nn = 97", "Stage 2 \nn = 200")) + # 新显示名称
  labs(
    # title = paste0("  β = 0.073 *  n = 690"), 
    x = "CSF Aβ42, t-tau, p-tau (model 1)", y = "LC3A/LC3B ratio", fill = "Stage",color = "Stage") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),   # 去掉网格
    panel.border = element_blank(), # 去掉框线
    # axis.line.x = element_line(color = "black", linewidth = 0.6), # 隐藏x 轴
    axis.ticks.x = element_blank(),  # 隐藏 X 轴所有刻度线
    axis.line.y = element_line(color = "black", linewidth = 0.6), # 只画 y 轴
    # axis.line = element_line(color = "black", linewidth = 0.6), # 加粗坐标轴
    axis.ticks = element_line(color = "black", linewidth = 0.6), # 坐标轴刻度线
    axis.title.y = element_text(margin = margin(r = 0)), #Y轴标签和轴距离
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(margin = margin(t = -10)) ,  # X轴分组标签向上移动
    strip.text = element_text( size = 12)
  )
a3 

###
a4 <- ggplot(merge_df) +
  geom_boxplot(aes(x = Stage, y = LC3AB,color = Stage), width = 0.3, alpha = 1, outlier.shape = NA) +
  geom_half_violin(aes(x = Stage, y = LC3AB, fill = Stage),
                   position = position_nudge(x=0.2), #箱式图和小提琴图间隔
                   width = 0.6, side = "r", alpha = 0.8, trim = FALSE, color = NA) +
  geom_jitter(aes(x = Stage, y = LC3AB, color = Stage),
              width = 0.1, alpha = 0.8, size = 2,
              # shape = 15 #实菱形
              # shape = 16 #实正方形
              shape = 17 #实三角形
              # shape = 18 #实菱形
  ) +
  scale_fill_manual(values = c("SNAP" = "#354e55", "Stage 0" = "#225081", "Stage 1" = "#d95725", "Stage 2" = "#6e1a00")) +
  scale_color_manual(values = c("SNAP" = "#354e55", "Stage 0" = "#225081", "Stage 1" = "#d95725", "Stage 2" = "#6e1a00")) +
  # coord_cartesian(xlim = c(-5, 6)) + #设定X轴范围
  # coord_cartesian(ylim = c(-10, -2)) + #设定Y轴范围
  geom_signif(aes(x = ATN_Stage, y = LC3AB),
              comparisons=list(c("Stage 0","Stage 1"), c("Stage 0","Stage 2"), c("Stage 1","Stage 2"),
                               c("SNAP","Stage 0"), c("SNAP","Stage 1"), c("SNAP","Stage 2")),
              annotations =c("**","***","n.s.","n.s.","*","**"),#手动添加显著性*
              vjust=0,# 星号与比较线上下距离
              textsize=6, # 星号大小
              color="black",
              size = 0.3, #比较线粗细
              step_increase = 0.1, # 两条比较线上下距离
              tip_length = 0 # 比较横线的两端短线长度为0
              # tip_length = rep(0.03, 5) # 4条横线的两端短线长度为0.03
  )+
  geom_segment(data = data.frame(xstart = c(0.6,1.6,2.6,3.6,  0.6,1.6,2.6,3.6,  1.4,2.4,3.4,4.4), # X轴坐标分组化
                                 xend   = c(1.4,2.4,3.4,4.4,  0.6,1.6,2.6,3.6,  1.4,2.4,3.4,4.4),
                                 ystart = c(-5.6,-5.6,-5.6,-5.6,  -5.6,-5.6,-5.6,-5.6,  -5.6,-5.6,-5.6,-5.6),
                                 yend   = c(-5.6,-5.6,-5.6,-5.6,  -5.8,-5.8,-5.8,-5.8,  -5.8,-5.8,-5.8,-5.8)),
               aes(x=xstart, xend=xend, y=ystart, yend=yend), size=0.6) +
  # coord_cartesian(ylim = c(-3, 3.5)) + # 限制显示范围，但不裁剪图中数据
  scale_x_discrete(
    limits = c("SNAP", "Stage 0", "Stage 1", "Stage 2"),  # 原始分类顺序
    labels = c("SNAP \nn = 49", "Stage 0 \nn = 120", "Stage 1 \nn = 97", "Stage 2 \nn = 200")) + # 新显示名称
  labs(
    # title = paste0("  β = 0.073 *  n = 690"), 
    x = "CSF Aβ42, t-tau, p-tau (model 2)", y = "LC3A/LC3B ratio", fill = "Stage",color = "Stage") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),   # 去掉网格
    panel.border = element_blank(), # 去掉框线
    # axis.line.x = element_line(color = "black", linewidth = 0.6), # 隐藏x 轴
    axis.ticks.x = element_blank(),  # 隐藏 X 轴所有刻度线
    axis.line.y = element_line(color = "black", linewidth = 0.6), # 只画 y 轴
    # axis.line = element_line(color = "black", linewidth = 0.6), # 加粗坐标轴
    axis.ticks = element_line(color = "black", linewidth = 0.6), # 坐标轴刻度线
    axis.title.y = element_text(margin = margin(r = 0)), #Y轴标签和轴距离
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(margin = margin(t = -10)) ,  # X轴分组标签向上移动
    strip.text = element_text(size = 12)
  )
a4
# ###
# 第二部分：p1_all–p4_all
bottom_part <- ggarrange(
  o1, o3, a3, a1, o2, o4, a4, a2,
  ncol = 4, nrow = 2,
  labels = LETTERS[1:8],
  # heights = c(5, 1.1),  
  widths = c(1,1,1.5,1.5),
  font.label = list(size = 16)
)
bottom_part

ggsave(bottom_part, file = paste0(path_output,"lm_LC3_ATN.pdf"),
       height = 10, width = 16, dpi = 300)


# ###
# ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1,
#           labels = c("M","N","O","P"),
#           font.label = list(size = 16),
#           widths = c(0.5, 0.5, 0.5),   # 控制每个图片间宽度
#           heights = c(0.5, 0.5, 0.5),    # 控制每个图片间高度
#           align = "hv" # 控制图之间坐标轴是否严格对齐
#           # common.legend = TRUE, #控制是否合并图例，减少重复
#           # legend = "right" #合并图例位置
# )
# ggsave(file = paste0(path_output,"lm_LC3_ATN_41.pdf"),height = 4.5, width = 20, dpi = 300)

