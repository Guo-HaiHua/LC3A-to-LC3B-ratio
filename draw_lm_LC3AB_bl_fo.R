###
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(grid)
library(ggpubr)
library(ggstance)
### 设置路径
path_data <- "G:\\Project\\ADNI_LC3AB\\results\\"
### 定义颜色模板
# mycolor <- c("#e56254","#cc5292","#ec7e04","#42a1d3","#9a9e10","#9f032d","#024f7d","#016f56","#5a3e83")
mycolor <- c(
  "#E53935",  # 更鲜红
  "#D81B60",  # 更鲜玫红
  "#FB8C00",  # 更亮橙
  "#1E88E5",  # 更亮蓝
  "#C0CA33",  # 更亮黄绿
  "#B71C1C",  # 深红（增强对比）
  "#01579B",  # 深蓝（增强对比）
  "#00796B",  # 深青绿
  "#4527A0"   # 更深紫
)
####################################
# 读取数据
data1 <- read_xlsx(paste0(path_data,"lm_LC3AB_Cog_results_bl_M1_SOMA.xlsx")) %>% mutate(Model = "Model 1")
data2 <- read_xlsx(paste0(path_data,"lm_LC3AB_Cog_results_bl_M1_SOMA_covp62.xlsx")) %>% mutate(Model = "Model 2")
data <- rbind(data1, data2)
data <- rbind(data1, data2)
data[data == "CDRSB"] <- "CDRSB"
data[data == "ADAS13"] <- "ADAS13"
data[data == "MMSE"] <- "MMSE"
data[data == "MOCA"] <- "MoCA"
data[data == "PHC_MEM"] <- "PHC_MEM"
data[data == "PHC_EXF"] <- "PHC_EXF"
data[data == "PHC_LAN"] <- "PHC_LAN"
data[data == "PHC_VSP"] <- "PHC_VSP"
# 分割 CI 并转为数值
data <- data %>%
  mutate(CI = gsub("[()]", "", `95% CI`)) %>%
  separate(CI, into = c("LCI", "UCI"), sep = ",") %>%
  mutate(
    LCI = as.numeric(LCI),
    UCI = as.numeric(UCI),
    `ES value` = as.numeric(`ES value`)
  ) %>%
  group_by(Model) %>%
  arrange(desc(Outcome), .by_group = TRUE) %>%  # 倒序
  mutate(
    Outcome = factor(Outcome, levels = unique(Outcome)),
    mycolor = rep(mycolor, length.out = n())
  ) %>%
  ungroup()

### 绘图：Model dodge + mycolor
a1 <- ggplot(data, aes(x = `ES value`, y = Outcome)) +
  geom_pointrangeh(
    aes(
      xmin = LCI,
      xmax = UCI,
      color = mycolor,
      shape = Model,
      group = Model
    ),
    size = 1.3,
    fatten = 6,   # 控制点大小
    position = position_dodge(width = 0.8)
  )+
  scale_shape_manual(values = c("Model 1" = 16, "Model 2" = 17)) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, color = mycolor, group = Model), 
                  width = 0, 
                 position = position_dodge(width = 0.8)) +
  geom_text(
    aes(
      x = UCI,
      label = paste0(
        round(`ES value`, 3), 
        " (", round(LCI, 3), ", ", round(UCI, 3), ") ",
        case_when(
          as.numeric(P_fdr) < 0.001 ~ "***",
          as.numeric(P_fdr) < 0.01  ~ "**",
          as.numeric(P_fdr) < 0.05  ~ "*",
          TRUE ~ ""
        )
      ),
      group = Model
    ),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 6,
    color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.3, color = "#8e8c8c") +
  scale_x_continuous(limits = c(-0.25, 0.7)) +
  scale_color_identity() +  # 使用数据框中的颜色
  theme_minimal() +
  labs(x = "β value (95% CI)", y = "") +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7, color = "black"),
    axis.ticks.length = unit(0.3, "cm")
  )
a1

##########################################
# 读取数据
data1 <- read_xlsx(paste0(path_data,"lmm_LC3AB_Cog_results_fo_M1_SOMA.xlsx")) %>% mutate(Model = "Model 1")
data2 <- read_xlsx(paste0(path_data,"lmm_LC3AB_Cog_results_fo_M1_SOMA_covp62.xlsx")) %>% mutate(Model = "Model 2")
data <- rbind(data1, data2)
data[data == "CDRSB"] <- "CDRSB change"
data[data == "ADAS13"] <- "ADAS13 change"
data[data == "MMSE"] <- "MMSE change"
data[data == "MOCA"] <- "MoCA change"
data[data == "PHC_MEM"] <- "PHC_MEM change"
data[data == "PHC_EXF"] <- "PHC_EXF change"
data[data == "PHC_LAN"] <- "PHC_LAN change"
data[data == "PHC_VSP"] <- "PHC_VSP change"
# 分割 CI 并转为数值
data <- data %>%
  mutate(CI = gsub("[()]", "", `95% CI_xt`)) %>%
  separate(CI, into = c("LCI", "UCI"), sep = ",") %>%
  mutate(
    LCI = as.numeric(LCI),
    UCI = as.numeric(UCI),
    `ES_xt value` = as.numeric(`ES_xt value`)
  ) %>%
  group_by(Model) %>%
  arrange(desc(Outcome), .by_group = TRUE) %>%  # 倒序
  mutate(
    Outcome = factor(Outcome, levels = unique(Outcome)),
    mycolor = rep(mycolor, length.out = n())
  ) %>%
  ungroup()
### 绘图：Model dodge + mycolor
a11 <- ggplot(data, aes(x = `ES_xt value`, y = Outcome)) +
  geom_pointrangeh(
    aes(
      xmin = LCI,
      xmax = UCI,
      color = mycolor,
      shape = Model,
      group = Model
    ),
    size = 1.3,
    fatten = 6,   # 控制点大小
    position = position_dodge(width = 0.8)
  )+
  scale_shape_manual(values = c("Model 1" = 16, "Model 2" = 17)) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, color = mycolor, group = Model), 
                 width = 0, 
                 position = position_dodge(width = 0.8)) +
  geom_text(
    aes(
      x = UCI,
      label = paste0(
        round(`ES_xt value`, 3), 
        " (", round(LCI, 3), ", ", round(UCI, 3), ") ",
        case_when(
          as.numeric(P_fdr) < 0.001 ~ "***",
          as.numeric(P_fdr) < 0.01  ~ "**",
          as.numeric(P_fdr) < 0.05  ~ "*",
          TRUE ~ ""
        )
      ),
      group = Model
    ),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 6,
    color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.3, color = "#8e8c8c") +
  scale_x_continuous(limits = c(-0.1, 0.30)) +
  scale_color_identity() +  # 使用数据框中的颜色
  theme_minimal() +
  labs(x = "β value (95% CI)", y = "") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7, color = "black"),
    axis.ticks.length = unit(0.3, "cm")
  )
a11
###
part1 <- ggarrange(
  a1, a11,
  ncol = 2, nrow = 1,
  labels = LETTERS[1:2],
  font.label = list(size = 20),
  align = "h",          # 水平对齐坐标轴
  widths = c(1,1.2)     # 强制每个 panel 宽度相等
)
part1

# 保存PDF
ggsave(
  filename = paste0(path_data,"LC3_Cog.pdf"),
  plot = part1,
  height = 12, width = 20, dpi = 300
)
##########################################
# 第二份数据：生物指标
data1 <- read_xlsx(paste0(path_data,"lm_LC3AB_ADbio_results_bl_M1_SOMA.xlsx")) %>% mutate(Model = "Model 1")
data2 <- read_xlsx(paste0(path_data,"lm_LC3AB_ADbio_results_bl_M1_SOMA_covp62.xlsx")) %>% mutate(Model = "Model 2")

data <- rbind(data1, data2)
data[data == "MSD_PGRNCORRECTED"] <- "CSF PGRN"
data[data == "MSD_STREM2CORRECTED"] <- "CSF sTREM2"
data[data == "FDG"] <- "FDG PET"
data[data == "AV45"] <- "Amyloid PET"
data[data == "ABETA"] <- "CSF Aβ42"
data[data == "TAU"] <- "CSF t-tau"
data[data == "PTAU"] <- "CSF p-tau"
data[data == "PLASMA_NFL"] <- "Plasma NFL"
data[data == "META_TEMPORAL_SUVR_tobl"] <- "Tau PET"

# 分割 CI 并转为数值
data <- data %>%
  mutate(CI = gsub("[()]", "", `95% CI`)) %>%
  separate(CI, into = c("LCI", "UCI"), sep = ",") %>%
  mutate(
    LCI = as.numeric(LCI),
    UCI = as.numeric(UCI),
    `ES value` = as.numeric(`ES value`)
  ) %>%
  group_by(Model) %>%
  arrange(desc(Outcome), .by_group = TRUE) %>%  # 倒序
  mutate(
    Outcome = factor(Outcome, levels = unique(Outcome)),
    mycolor = rep(mycolor, length.out = n())
  ) %>%
  ungroup()

### 绘图：Model dodge + mycolor
a2 <- ggplot(data, aes(x = `ES value`, y = Outcome)) +
  geom_pointrangeh(
    aes(
      xmin = LCI,
      xmax = UCI,
      color = mycolor,
      shape = Model,
      group = Model
    ),
    size = 1.3,
    fatten = 6,   # 控制点大小
    position = position_dodge(width = 0.8)
  )+
  scale_shape_manual(values = c("Model 1" = 16, "Model 2" = 17)) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, color = mycolor, group = Model), 
                 width = 0, 
                 position = position_dodge(width = 0.8)) +
  geom_text(
    aes(
      x = UCI,
      label = paste0(
        round(`ES value`, 3), 
        " (", round(LCI, 3), ", ", round(UCI, 3), ") ",
        case_when(
          as.numeric(P_fdr) < 0.001 ~ "***",
          as.numeric(P_fdr) < 0.01  ~ "**",
          as.numeric(P_fdr) < 0.05  ~ "*",
          TRUE ~ ""
        )
      ),
      group = Model
    ),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 6,
    color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.3, color = "#8e8c8c") +
  scale_x_continuous(limits = c(-0.32, 0.80)) +
  scale_color_identity() +  # 使用数据框中的颜色
  theme_minimal() +
  labs(x = "β value (95% CI)", y = "") +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7, color = "black"),
    axis.ticks.length = unit(0.3, "cm")
  )
a2

##########################################
# 第二份数据：生物指标
data1 <- read_xlsx(paste0(path_data,"lmm_LC3AB_ADbio_results_fo_M1_SOMA.xlsx")) %>% mutate(Model = "Model 1")
data2 <- read_xlsx(paste0(path_data,"lmm_LC3AB_ADbio_results_fo_M1_SOMA_covp62.xlsx")) %>% mutate(Model = "Model 2")

data <- rbind(data1, data2)
data[data == "MSD_PGRNCORRECTED"] <- "CSF PGRN change"
data[data == "MSD_STREM2CORRECTED"] <- "CSF sTREM2 change"
data[data == "FDG"] <- "FDG PET change"
data[data == "AV45"] <- "Amyloid PET change"
data[data == "ABETA"] <- "CSF Aβ42 change"
data[data == "TAU"] <- "CSF t-tau change"
data[data == "PTAU"] <- "CSF p-tau change"
data[data == "PLASMA_NFL"] <- "Plasma NFL change"
data[data == "META_TEMPORAL_SUVR_tobl"] <- "Tau PET change"

# 分割 CI 并转为数值
data <- data %>%
  mutate(CI = gsub("[()]", "", `95% CI_xt`)) %>%
  separate(CI, into = c("LCI", "UCI"), sep = ",") %>%
  mutate(
    LCI = as.numeric(LCI),
    UCI = as.numeric(UCI),
    `ES_xt value` = as.numeric(`ES_xt value`)
  ) %>%
  group_by(Model) %>%
  arrange(desc(Outcome), .by_group = TRUE) %>%  # 倒序
  mutate(
    Outcome = factor(Outcome, levels = unique(Outcome)),
    mycolor = rep(mycolor, length.out = n())
  ) %>%
  ungroup()

### 绘图：Model dodge + mycolor
a22 <- ggplot(data, aes(x = `ES_xt value`, y = Outcome)) +
  geom_pointrangeh(
    aes(
      xmin = LCI,
      xmax = UCI,
      color = mycolor,
      shape = Model,
      group = Model
    ),
    size = 1.3,
    fatten = 6,   # 控制点大小
    position = position_dodge(width = 0.8)
  )+
  scale_shape_manual(values = c("Model 1" = 16, "Model 2" = 17)) +
  geom_errorbarh(aes(xmin = LCI, xmax = UCI, color = mycolor, group = Model), 
                 width = 0, 
                 position = position_dodge(width = 0.8)) +
  geom_text(
    aes(
      x = UCI ,
      label = paste0(
        round(`ES_xt value`, 3), 
        " (", round(LCI, 3), ", ", round(UCI, 3), ") ",
        case_when(
          as.numeric(P_fdr) < 0.001 ~ "***",
          as.numeric(P_fdr) < 0.01  ~ "**",
          as.numeric(P_fdr) < 0.05  ~ "*",
          TRUE ~ ""
        )
      ),
      group = Model
    ),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 6,
    color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.3, color = "#8e8c8c") +
  scale_x_continuous(limits = c(-0.05, 0.15)) +
  scale_color_identity() +  # 使用数据框中的颜色
  theme_minimal() +
  labs(x = "β value (95% CI)", y = "") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7, color = "black"),
    axis.ticks.length = unit(0.3, "cm")
  )
a22

###
part2 <- ggarrange(
  a2, a22,
  ncol = 2, nrow = 1,
  labels = LETTERS[1:2],
  font.label = list(size = 20),
  align = "h",          # 水平对齐坐标轴
  widths = c(1,1.2)     # 强制每个 panel 宽度相等
)
part2

# 保存PDF
ggsave(
  filename = paste0(path_data,"LC3_ADbio.pdf"),
  plot = part2,
  height = 12, width = 20, dpi = 300
)
##########################################
# 第三份数据：生物指标
data1 <- read_xlsx(paste0(path_data,"lm_LC3AB_Bra_results_bl_M1_SOMA.xlsx")) %>% mutate(Model = "Model 1")
data2 <- read_xlsx(paste0(path_data,"lm_LC3AB_Bra_results_bl_M1_SOMA_covp62.xlsx")) %>% mutate(Model = "Model 2")

data <- rbind(data1, data2)
data[data == "Entorhinal"] <- "Entorhinal"
data[data == "Fusiform"] <- "Fusiform"
data[data == "Hippocampus"] <- "Hippocampus"
data[data == "MidTemp"] <- "MidTemp"
data[data == "Ventricles"] <- "Ventricles"
data[data == "WholeBrain"] <- "WholeBrain"


# 分割 CI 并转为数值
data <- data %>%
  mutate(CI = gsub("[()]", "", `95% CI`)) %>%
  separate(CI, into = c("LCI", "UCI"), sep = ",") %>%
  mutate(
    LCI = as.numeric(LCI),
    UCI = as.numeric(UCI),
    `ES value` = as.numeric(`ES value`)
  ) %>%
  group_by(Model) %>%
  arrange(desc(Outcome), .by_group = TRUE) %>%  # 倒序
  mutate(
    Outcome = factor(Outcome, levels = unique(Outcome)),
    mycolor = rep(mycolor, length.out = n())
  ) %>%
  ungroup()

### 绘图：Model dodge + mycolor
a3 <- ggplot(data, aes(x = `ES value`, y = Outcome)) +
  geom_pointrangeh(
    aes(
      xmin = LCI,
      xmax = UCI,
      color = mycolor,
      shape = Model,
      group = Model
    ),
    size = 1.3,
    fatten = 6,   # 控制点大小
    position = position_dodge(width = 0.8)
  )+
  scale_shape_manual(values = c("Model 1" = 16, "Model 2" = 17)) +
  geom_text(
    aes(
      x = UCI + 0.001,
      label = paste0(
        round(`ES value`, 3), 
        " (", round(LCI, 3), ", ", round(UCI, 3), ") ",
        case_when(
          as.numeric(P_fdr) < 0.001 ~ "***",
          as.numeric(P_fdr) < 0.01  ~ "**",
          as.numeric(P_fdr) < 0.05  ~ "*",
          TRUE ~ ""
        )
      ),
      group = Model
    ),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 6,
    color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.3, color = "#8e8c8c") +
  scale_x_continuous(limits = c(-0.292, 0.40)) +
  scale_color_identity() +  # 使用数据框中的颜色
  theme_minimal() +
  labs(x = "β value (95% CI)", y = "") +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7, color = "black"),
    axis.ticks.length = unit(0.3, "cm")
  )
a3

##########################################
# 第二份数据：生物指标
data1 <- read_xlsx(paste0(path_data,"lmm_LC3AB_Bra_results_fo_M1_SOMA.xlsx")) %>% mutate(Model = "Model 1")
data2 <- read_xlsx(paste0(path_data,"lmm_LC3AB_Bra_results_fo_M1_SOMA_covp62.xlsx")) %>% mutate(Model = "Model 2")

data <- rbind(data1, data2)
data[data == "Entorhinal"] <- "Entorhinal change"
data[data == "Fusiform"] <- "Fusiform change"
data[data == "Hippocampus"] <- "Hippocampus change"
data[data == "MidTemp"] <- "MidTemp change"
data[data == "Ventricles"] <- "Ventricles change"
data[data == "WholeBrain"] <- "WholeBrain change"

# 分割 CI 并转为数值
data <- data %>%
  mutate(CI = gsub("[()]", "", `95% CI_xt`)) %>%
  separate(CI, into = c("LCI", "UCI"), sep = ",") %>%
  mutate(
    LCI = as.numeric(LCI),
    UCI = as.numeric(UCI),
    `ES_xt value` = as.numeric(`ES_xt value`)
  ) %>%
  group_by(Model) %>%
  arrange(desc(Outcome), .by_group = TRUE) %>%  # 倒序
  mutate(
    Outcome = factor(Outcome, levels = unique(Outcome)),
    mycolor = rep(mycolor, length.out = n())
  ) %>%
  ungroup()

### 绘图：Model dodge + mycolor
a33 <- ggplot(data, aes(x = `ES_xt value`, y = Outcome)) +
  geom_pointrangeh(
    aes(
      xmin = LCI,
      xmax = UCI,
      color = mycolor,
      shape = Model,
      group = Model
    ),
    size = 1.3,
    fatten = 6,   # 控制点大小
    position = position_dodge(width = 0.8)
  )+
  scale_shape_manual(values = c("Model 1" = 16, "Model 2" = 17)) +
  geom_text(
    aes(
      x = UCI,
      label = paste0(
        round(`ES_xt value`, 3), 
        " (", round(LCI, 3), ", ", round(UCI, 3), ") ",
        case_when(
          as.numeric(P_fdr) < 0.001 ~ "***",
          as.numeric(P_fdr) < 0.01  ~ "**",
          as.numeric(P_fdr) < 0.05  ~ "*",
          TRUE ~ ""
        )
      ),
      group = Model
    ),
    position = position_dodge(width = 0.8),
    hjust = -0.2,
    size = 6,
    color = "black"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.3, color = "#8e8c8c") +
  scale_x_continuous(limits = c(-0.035, 0.075)) +
  scale_color_identity() +  # 使用数据框中的颜色
  theme_minimal() +
  labs(x = "β value (95% CI)", y = "") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18),
    axis.text.x = element_text(size = 18, color = "black"),
    axis.text.y = element_text(size = 18, color = "black", angle = 45, hjust = 1),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7, color = "black"),
    axis.ticks.length = unit(0.3, "cm")
  )
a33

###
part3 <- ggarrange(
  a3, a33,
  ncol = 2, nrow = 1,
  labels = LETTERS[1:3],
  font.label = list(size = 20),
  align = "h",          # 水平对齐坐标轴
  widths = c(1,1)     # 强制每个 panel 宽度相等
)
part3
# 保存PDF
ggsave(
  filename = paste0(path_data,"LC3_Bra.pdf"),
  plot = part3,
  height = 8, width = 20, dpi = 300
)

