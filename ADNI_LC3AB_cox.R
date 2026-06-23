###
path_data <- "G:\\Project\\ADNI_LC3AB\\Data\\"
path_adni <- "G:\\Project\\Cohort\\ADNI\\Original\\"
path_model <- "G:\\Project\\Model_R\\HHG\\"
path_output <- "G:\\Project\\ADNI_LC3AB\\Results\\"
###
library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(pROC)
library(readxl)
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
#
merge_df5 <-rbind(control, case)







### extra
hist(merge_df5$Years_bl)
table(merge_df5$DX_C)
quantile(merge_df5$Years_bl, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
max(merge_df5$Years_bl)
merge_df5 %>%
  group_by(DX_bl_C) %>%
  summarise(
    median_Years = median(Years_bl, na.rm = TRUE),
    Q1 = quantile(Years_bl, 0.25, na.rm = TRUE),
    Q3 = quantile(Years_bl, 0.75, na.rm = TRUE),
    IQR = IQR(Years_bl, na.rm = TRUE)
  )

t.test(merge_df5$Years_bl, merge_df5$DX_bl_C)
###
# merge_df5 <- subset(merge_df5, DX !="MCI")
# merge_df5$DX[merge_df5$DX %in% c("CN", "MCI")] <- "noDementia"
merge_df <- subset(merge_df5, Years_bl>=0)
# merge_df <- subset(merge_df5, Years_bl>=2)







###
source(paste0(path_model,"ghhcox.R"))
#
cox_Nestedlist <- NULL
for (i in names(merge_df)[c(10)]){
  for (j in names(merge_df)[c(75)]){
    cox_cov <- names(merge_df)[c(18,20,21,26)]
    cox_Nestedlist[[i]][[j]] <- ghhcox(merge_df[i],merge_df[j],merge_df["Years_bl"],merge_df["RID"],merge_df)
  }
}
cox_Flatlist <- NULL
for (i in names(merge_df)[10]){
  cox_Flatlist[[i]] <- Reduce(rbind.data.frame,cox_Nestedlist[[i]],accumulate = F)}
cox_df_results <- Reduce(rbind.data.frame,cox_Flatlist,accumulate = F)
#
write_xlsx(cox_df_results, path = paste0(path_output, "cox_LC3AB_results_noexclude.xlsx"))
# write_xlsx(cox_df_results, path = paste0(path_output, "cox_LC3AB_results_2y.xlsx"))
write_xlsx(cox_df_results, path = paste0(path_output, "cox_LC3AB_results_noexclude_covp62.xlsx"))
write_xlsx(cox_df_results, path = paste0(path_output, "cox_LC3AB_results_2y_covp62.xlsx"))


###------------
n <- nrow(merge_df)
predicted_risks <- numeric(n)
true_labels <- merge_df$DX_C
for (i in 1:n) {
  train_data <- merge_df[-i, ]
  test_data <- merge_df[i, , drop = FALSE]
  if (sum(train_data$DX_C) < 1) {
    predicted_risks[i] <- NA
    next
  }
  # Cox模型
  model <- coxph(Surv(Years_bl, DX_C) ~ LC3A + LC3B + AV45 + TAU,
                 data = train_data)
  
  # 使用线性预测值 (lp)，而不是risk
  predicted_risks[i] <- predict(model, newdata = test_data, type = "lp")
}
# 计算ROC
valid_idx <- !is.na(predicted_risks)
roc_result <- roc(true_labels[valid_idx], predicted_risks[valid_idx], quiet = TRUE)
auc_val <- auc(roc_result)
# 95% CI via bootstrap
ci_val <- ci.auc(roc_result, method = "bootstrap", boot.n = 2000)
# 输出结果
print(paste0("AUC: ", round(auc_val, 3)))
print(ci_val)
plot(roc_result, main = "ROC Curve for LC3B in Predicting AD Risk", col = "#1f77b4")













##########################
merge_df4$P_Stage[merge_df4$PTAU <= 21.8] <- 0
merge_df4$P_Stage[merge_df4$PTAU > 21.8] <- 1

merge_df4$T_Stage[merge_df4$TAU <= 245] <- 0
merge_df4$T_Stage[merge_df4$TAU > 245] <- 1
###
case <- merge_df4 %>%
  filter(P_Stage == "1") %>%
  group_by(RID) %>%
  slice(which.min(Years_bl)) %>%
  ungroup()

#
control <- merge_df4 %>%
  # 去掉曾经出现过 Dementia 的 RID
  group_by(RID) %>%
  filter(!any(P_Stage == "1")) %>%
  # 取每个 RID 最大 Years_bl
  slice(which.max(Years_bl)) %>%
  ungroup()

merge_df <-rbind(control, case)
###
source(paste0(path_model,"ghhcox.R"))
#
cox_Nestedlist <- NULL
for (i in names(merge_df)[c(9:10)]){
  for (j in names(merge_df)[c(185)]){
    cox_cov <- names(merge_df)[c(18,20,21,26)]
    cox_Nestedlist[[i]][[j]] <- ghhcox(merge_df[i],merge_df[j],merge_df["Years_bl"],merge_df["RID"],merge_df)
  }
}
cox_Flatlist <- NULL
for (i in names(merge_df)[9:10]){
  cox_Flatlist[[i]] <- Reduce(rbind.data.frame,cox_Nestedlist[[i]],accumulate = F)}
cox_df_results <- Reduce(rbind.data.frame,cox_Flatlist,accumulate = F)
#
# write.csv(cox_df_results,file = paste0(path_output,"LC3_COX_results_M1.csv"),row.names = F)
write_xlsx(cox_df_results, path = paste0(path_output, "LC3_PTAU_COX_results_M1.xlsx"))


# ########## ROC
# # 设置随机种子以便重复实验
# set.seed(123)
# # 随机划分训练集和测试集（80% 训练集，20% 测试集）
# n <- nrow(merge_df)
# train_idx <- sample(1:n, size = floor(0.7 * n))  # 随机选择 80% 的索引作为训练集
# train_data <- merge_df[train_idx, ]
# test_data <- merge_df[-train_idx, ]
# # 构建Cox回归模型
# cox_model <- coxph(Surv(Years_bl, P_Stage) ~ LC3A, data = train_data)
# summary(cox_model)
# # 预测测试集的风险值
# predicted_risks <- predict(cox_model, newdata = test_data, type = "risk")
# # 检查是否有足够的数据计算AUC
# if (length(unique(test_data$P_Stage)) < 2) {
#   stop("Not enough variation in DX_C to compute AUC (need both 0 and 1)")
# }
# # 计算整体AUC
# roc_result <- roc(test_data$P_Stage, predicted_risks, quiet = TRUE)
# auc_value <- auc(roc_result)
# # 输出结果
# cat("AUC from 80-20 split:", auc_value, "\n")
# # 可视化ROC曲线（可选）
# plot(roc_result, main = "ROC Curve for LC3B in Predicting AD Risk", col = "#1f77b4")


# 初始化LOOCV存储预测结果
n <- nrow(merge_df)
predicted_risks <- numeric(n)  # 存储预测风险值
true_labels <- merge_df$P_Stage      # 存储真实标签
# LOOCV循环
for (i in 1:n) {
  # 划分训练集和测试集
  train_data <- merge_df[-i, ]
  test_data <- merge_df[i, ]
  # 检查训练集中是否有足够的事件（DX_C=1）
  if (sum(train_data$P_Stage) < 1) {
    warning(paste("No events (DX_C=1) in training set for iteration", i, "- skipping"))
    predicted_risks[i] <- NA
    next
  }
  # 构建Cox回归模型
  cox_model <- coxph(Surv(Years_bl, P_Stage) ~ LC3B, data = train_data)
  summary(cox_model)
  # 预测测试集的风险值
  predicted_risks[i] <- predict(cox_model, newdata = test_data, type = "risk")
}
# 移除预测失败的行（如果有NA）
valid_idx <- !is.na(predicted_risks)
predicted_risks <- predicted_risks[valid_idx]
true_labels <- true_labels[valid_idx]
# 检查是否有足够的数据计算AUC
if (length(unique(true_labels)) < 2) {
  stop("Not enough variation in DX_C to compute AUC (need both 0 and 1)")
}
# 计算整体AUC
roc_result <- roc(true_labels, predicted_risks, quiet = TRUE)
auc_value <- auc(roc_result)
# 输出结果
cat("AUC from LOOCV:", auc_value, "\n")
# 可视化ROC曲线（可选）
plot(roc_result, main = "ROC Curve for LC3B in Predicting AD Risk", col = "#1f77b4")









###
case <- merge_df4 %>%
  filter(T_Stage == "1") %>%
  group_by(RID) %>%
  slice(which.min(Years_bl)) %>%
  ungroup()

#
control <- merge_df4 %>%
  # 去掉曾经出现过 Dementia 的 RID
  group_by(RID) %>%
  filter(!any(T_Stage == "1")) %>%
  # 取每个 RID 最大 Years_bl
  slice(which.max(Years_bl)) %>%
  ungroup()

merge_df <-rbind(control, case)
###
source(paste0(path_model,"ghhcox.R"))
#
cox_Nestedlist <- NULL
for (i in names(merge_df)[c(9:10)]){
  for (j in names(merge_df)[c(186)]){
    cox_cov <- names(merge_df)[c(18,20,21,26)]
    cox_Nestedlist[[i]][[j]] <- ghhcox(merge_df[i],merge_df[j],merge_df["Years_bl"],merge_df["RID"],merge_df)
  }
}
cox_Flatlist <- NULL
for (i in names(merge_df)[9:10]){
  cox_Flatlist[[i]] <- Reduce(rbind.data.frame,cox_Nestedlist[[i]],accumulate = F)}
cox_df_results <- Reduce(rbind.data.frame,cox_Flatlist,accumulate = F)
#
# write.csv(cox_df_results,file = paste0(path_output,"LC3_COX_results_M1.csv"),row.names = F)
write_xlsx(cox_df_results, path = paste0(path_output, "LC3_TAU_COX_results_M1.xlsx"))


########## ROC
# # 设置随机种子以便重复实验
# set.seed(123)
# # 随机划分训练集和测试集（80% 训练集，20% 测试集）
# n <- nrow(merge_df)
# train_idx <- sample(1:n, size = floor(0.8 * n))  # 随机选择 80% 的索引作为训练集
# train_data <- merge_df[train_idx, ]
# test_data <- merge_df[-train_idx, ]
# # 构建Cox回归模型
# cox_model <- coxph(Surv(Years_bl, T_Stage) ~ LC3A, data = train_data)
# summary(cox_model)
# # 预测测试集的风险值
# predicted_risks <- predict(cox_model, newdata = test_data, type = "risk")
# # 检查是否有足够的数据计算AUC
# if (length(unique(test_data$T_Stage)) < 2) {
#   stop("Not enough variation in DX_C to compute AUC (need both 0 and 1)")
# }
# # 计算整体AUC
# roc_result <- roc(test_data$T_Stage, predicted_risks, quiet = TRUE)
# auc_value <- auc(roc_result)
# # 输出结果
# cat("AUC from 80-20 split:", auc_value, "\n")
# # 可视化ROC曲线（可选）
# plot(roc_result, main = "ROC Curve for LC3B in Predicting AD Risk", col = "#1f77b4")




# 初始化LOOCV存储预测结果
n <- nrow(merge_df)
predicted_risks <- numeric(n)  # 存储预测风险值
true_labels <- merge_df$T_Stage      # 存储真实标签
# LOOCV循环
for (i in 1:n) {
  # 划分训练集和测试集
  train_data <- merge_df[-i, ]
  test_data <- merge_df[i, ]
  # 检查训练集中是否有足够的事件（DX_C=1）
  if (sum(train_data$T_Stage) < 1) {
    warning(paste("No events (DX_C=1) in training set for iteration", i, "- skipping"))
    predicted_risks[i] <- NA
    next
  }
  # 构建Cox回归模型
  cox_model <- coxph(Surv(Years_bl, T_Stage) ~ LC3B, data = train_data)
  summary(cox_model)
  # 预测测试集的风险值
  predicted_risks[i] <- predict(cox_model, newdata = test_data, type = "risk")
}
# 移除预测失败的行（如果有NA）
valid_idx <- !is.na(predicted_risks)
predicted_risks <- predicted_risks[valid_idx]
true_labels <- true_labels[valid_idx]
# 检查是否有足够的数据计算AUC
if (length(unique(true_labels)) < 2) {
  stop("Not enough variation in DX_C to compute AUC (need both 0 and 1)")
}
# 计算整体AUC
roc_result <- roc(true_labels, predicted_risks, quiet = TRUE)
auc_value <- auc(roc_result)
# 输出结果
cat("AUC from LOOCV:", auc_value, "\n")
# 可视化ROC曲线（可选）
plot(roc_result, main = "ROC Curve for LC3B in Predicting AD Risk", col = "#1f77b4")



























########################################draw
roc1 <- roc(merge_df$P_Stage, merge_df$LC3A) 
roc2 <- roc(merge_df$P_Stage, merge_df$LC3B)
roc1;roc2
# 颜色设置
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "black")
# 按特征名的首字母的大小写排序
features_sorted <- names(auc_values)[order(names(auc_values),decreasing = T)]  # 按首字母的大小写排序

# 绘制 ROC 曲线
ggplot(merge_df, aes(x = 1 - specificity, y = sensitivity, color = feature)) +
  geom_step(size = 1.2) +  # 使用阶梯样式绘制
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # 参考线（y=x）
  scale_color_manual(values = colors) +  # 设置颜色
  labs(title = "ROC Curves for Each Method", x = "1 - Specificity", y = "Sensitivity", color = "Feature") +
  ggprism::theme_prism(border = T)+
  # 添加 AUC 值，左对齐，按首字母的大小写排序
  annotate("text", x = 0.55, y = seq(0.05, 0.25, length.out = length(features_sorted)), 
           label = paste0(features_sorted, " AUC: ", round(auc_values[features_sorted], 2)), 
           size = 5, hjust = 0,color = rev(colors))+  # 左对齐
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = 1)





# 加载必要的R包
library(pROC)
library(ggplot2)
library(ggprism)  # 如果未安装ggprism，可先运行 install.packages("ggprism")

# 假设数据框merge_df已经存在，且包含P_Stage, LC3A, LC3B列
# 计算ROC曲线
roc1 <- roc(merge_df$P_Stage, merge_df$LC3A, quiet = TRUE)
roc2 <- roc(merge_df$P_Stage, merge_df$LC3B, quiet = TRUE)

# 提取ROC曲线数据
roc_data1 <- data.frame(
  specificity = roc1$specificities,
  sensitivity = roc1$sensitivities,
  feature = "LC3A"
)

roc_data2 <- data.frame(
  specificity = roc2$specificities,
  sensitivity = roc2$sensitivities,
  feature = "LC3B"
)

# 合并ROC数据
roc_data <- rbind(roc_data1, roc_data2)

# 提取AUC值
auc_values <- c(LC3A = auc(roc1), LC3B = auc(roc2))

# 按特征名的首字母大小写排序
features_sorted <- names(auc_values)[order(names(auc_values), decreasing = TRUE)]

# 颜色设置（为两个特征分配颜色）
colors <- c("#E41A1C", "#377EB8")  # 为LC3A和LC3B分配两种颜色

# 绘制ROC曲线
ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, color = feature)) +
  geom_step(size = 1.2) +  # 使用阶梯样式绘制ROC曲线
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # 参考线（y=x）
  scale_color_manual(values = colors) +  # 设置自定义颜色
  # labs(
  #   title = "各特征的ROC曲线",
  #   x = "1 - 特异度",
  #   y = "灵敏度",
  #   color = "特征"
  # ) +
  ggprism::theme_prism(border = TRUE) +  # 使用ggprism主题
  annotate(
    "text",
    x = 0.55,
    y = seq(0.05, 0.15, length.out = length(features_sorted)),
    label = paste0(features_sorted, " AUC: ", round(auc_values[features_sorted], 2)),
    size = 5,
    hjust = 0,  # 左对齐
    color = colors  # 使用对应颜色
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +  # 标题居中
  coord_fixed(ratio = 1)  # 固定坐标轴比例

coords(roc1, "best", ret = c("threshold", "sensitivity", "specificity"))




