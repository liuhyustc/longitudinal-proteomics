# ===================================================================
# Section 0: 加载必要的包
# ===================================================================
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")

library(limma)
library(dplyr)
library(ggplot2)
library(patchwork)

# ===================================================================
# Section 1: 全局模拟参数
# ===================================================================
N_simulations <- 100  # 模拟次数
set.seed(3476)

# 样本量场景
#sample_size_scenarios <- 2:4
sample_size_scenarios <- 2:10

# 固定参数
N_total <- 3000
N_de <- 30
N_null <- N_total - N_de
time_points <- 1:10
n_time_points <- length(time_points)
prior_d0 <- 100
prior_s0_2 <- 0.09
fdr_threshold <- 0.05

performance_summary_list <- list()

cat(sprintf("开始多样本量模拟分析 (Simulations per scenario: %d)...\n", N_simulations))

# ===================================================================
# Section 2: 主循环
# ===================================================================

for (n_per_group in sample_size_scenarios) {
  
  n_case <- n_per_group
  n_control <- n_per_group
  n_subjects <- n_case + n_control
  n_samples <- n_subjects * n_time_points
  
  cat(sprintf("\n>>> 正在运行场景: 每组样本量 = %d (总 Subject = %d) <<<\n", n_per_group, n_subjects))
  
  # --- 构建矩阵 ---
  groups <- c(rep("Case", n_case), rep("Control", n_control))
  subject_ids <- paste0("S", 1:n_subjects)
  
  col_info <- expand.grid(Time = time_points, SubjectID = subject_ids)
  col_info <- col_info[order(match(col_info$SubjectID, subject_ids), col_info$Time), ]
  col_info$Group <- rep(groups, each = n_time_points)
  col_info$Time <- as.factor(col_info$Time)
  col_info$Group <- as.factor(col_info$Group)
  
  design_factor <- factor(paste(col_info$Group, col_info$Time, sep="."))
  design_matrix <- model.matrix(~0 + design_factor)
  colnames(design_matrix) <- levels(design_factor)
  
  contrast_string <- paste(paste("Case.", time_points, sep=""), "-", paste("Control.", time_points, sep=""))
  contrast_matrix <- makeContrasts(contrasts = contrast_string, levels = design_matrix)
  
  current_sim_results <- list()
  
  # --- 内部模拟 ---
  for (sim in 1:N_simulations) {
    if (sim %% 10 == 0) cat(sprintf("   Running sim %d/%d...\n", sim, N_simulations))
    
    # 1. 生成真实状态
    true_status <- factor(c(rep("DE", N_de), rep("Null", N_null)), levels = c("Null", "DE"))
    peak_time <- 5.5
    signal_shape <- dnorm(time_points, mean = peak_time, sd = 2.5)
    signal_shape <- signal_shape / max(signal_shape)
    peak_log2fc <- runif(N_de, min = 0.8, max = 1.2)
    peak_log2fc[1:(N_de/2)] <- -peak_log2fc[1:(N_de/2)]
    true_log2fc_matrix <- matrix(0, nrow = N_total, ncol = n_time_points)
    true_log2fc_matrix[1:N_de, ] <- peak_log2fc %o% signal_shape
    
    # 2. 向量化生成数据
    base_expression <- rnorm(N_total, mean = 8, sd = 2)
    subject_eff_matrix <- matrix(rnorm(N_total * n_subjects, 0, 0.2), nrow = N_total)
    subject_component <- subject_eff_matrix[, rep(1:n_subjects, each = n_time_points)]
    
    group_component <- matrix(0, nrow = N_total, ncol = n_samples)
    case_indices <- which(col_info$Group == "Case")
    time_indices <- as.numeric(col_info$Time[case_indices])
    group_component[, case_indices] <- true_log2fc_matrix[, time_indices]
    
    sigma2_g <- (prior_d0 * prior_s0_2) / rchisq(N_total, df = prior_d0)
    sd_g <- sqrt(sigma2_g)
    error_component <- matrix(rnorm(N_total * n_samples), nrow = N_total) * sd_g
    
    LogAbundance <- base_expression + subject_component + group_component + error_component
    
    # 3. Limma 分析
    corfit <- duplicateCorrelation(LogAbundance, design = design_matrix, block = col_info$SubjectID)
    fit <- lmFit(LogAbundance, design_matrix, block = col_info$SubjectID, correlation = corfit$consensus)
    fit_cont <- contrasts.fit(fit, contrast_matrix)
    eb_fit <- eBayes(fit_cont, robust = TRUE, trend = TRUE)
    
    is_de_ground_truth <- (true_status == "DE")
    
    # --- Method 1: Moderated F-test ---
    # 关键修正：指定 coef=1:10 确保是对所有时间点的联合检验
    # 这是一个 P 值 per Gene，然后做 BH 校正
    f_test_res <- topTable(eb_fit, coef = 1:n_time_points, number = Inf, sort.by = "none")
    pred_f <- f_test_res$adj.P.Val < fdr_threshold
    
    # --- Method 2: Moderated t-test (Min-P -> BH) ---
    min_mod_p <- apply(eb_fit$p.value, 1, min)
    pred_mod_t <- p.adjust(min_mod_p, method = "BH") < fdr_threshold
    
    # --- Method 3: Ordinary t-test (Min-P -> BH) ---
    ord_t_stat <- fit_cont$coefficients / (fit_cont$stdev.unscaled * fit$sigma)
    ord_t_pvals <- 2 * pt(-abs(ord_t_stat), df = fit$df.residual)
    min_ord_p <- apply(ord_t_pvals, 1, min)
    pred_ord_t <- p.adjust(min_ord_p, method = "BH") < fdr_threshold
    
    # 4. 统计
    calc_stats <- function(pred, method_name) {
      TP <- sum(is_de_ground_truth & pred)
      FP <- sum(!is_de_ground_truth & pred)
      data.frame(Method = method_name, TP = TP, FP = FP)
    }
    
    current_sim_results[[sim]] <- bind_rows(
      calc_stats(pred_f, "Moderated F-test"),
      calc_stats(pred_mod_t, "Moderated t-test"),
      calc_stats(pred_ord_t, "Ordinary t-test")
    )
  }
  
  # 汇总
  sim_df <- bind_rows(current_sim_results)
  summary_df <- sim_df %>%
    group_by(Method) %>%
    summarise(
      Mean_TP = mean(TP),
      Mean_FP = mean(FP),
      .groups = 'drop'
    ) %>%
    mutate(
      FDR = Mean_FP / (Mean_TP + Mean_FP),
      Power = Mean_TP / N_de,
      Sample_Size = n_per_group
    ) %>%
    mutate(FDR = ifelse(is.nan(FDR), 0, FDR))
  
  performance_summary_list[[as.character(n_per_group)]] <- summary_df
  
}

final_results <- bind_rows(performance_summary_list)

cat("\n--- 模拟全部完成，准备绘图 ---\n")
print(final_results, n =27)


# ===================================================================
# Section 3: 绘图 (定制版：保留网格线调整，原始标注风格)
# ===================================================================

# 1. 定义颜色和形状
method_colors <- c("Ordinary t-test" = "#d95f02", 
                   "Moderated t-test" = "#1b9e77", 
                   "Moderated F-test" = "#7570b3")

# 形状: 16=圆, 17=三角, 15=方块
method_shapes <- c("Ordinary t-test" = 16, 
                   "Moderated t-test" = 17, 
                   "Moderated F-test" = 15)
save(final_results,file = "simulation/final_results.RData")
load("~/Desktop/horse/horse/simulation/final_results.RData")
# 确保因子顺序
final_results$Method <- factor(final_results$Method, 
                               levels = c("Moderated F-test", "Moderated t-test", "Ordinary t-test"))

# 2. 定义主题 (保留网格线调整，去除坐标轴增强，去除加粗)
clean_theme <- theme_minimal() +
  theme(
    # 网格线优化：移除次网格，调淡主网格 (保留您的要求)
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.5),
    
    # 图例放底部，无标题
    legend.position = "bottom",
    legend.title = element_blank()
  )

# --- 图 1: FDR ---
p_fdr <- ggplot(final_results, aes(x = Sample_Size, y = FDR, color = Method, group = Method, shape = Method)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 3) +
  geom_hline(yintercept = fdr_threshold, linetype = "dashed", color = "red") +
  
  # 恢复您最原始的标注风格：红色文字，简单位置
  annotate("text", x = max(sample_size_scenarios), y = fdr_threshold, 
           label = "Target FDR = 0.05", 
           vjust = -0.5, hjust = 1, color = "red", size=3.5) +
  
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  
  scale_x_continuous(breaks = sample_size_scenarios) +
  scale_y_continuous(limits = c(0, max(0.15, max(final_results$FDR) * 1.1))) +
  
  labs(title = NULL, x = "Sample Size per Group", y = "Empirical FDR") +
  clean_theme + 
  theme(legend.position = "none") # 在组合图中隐藏图例

# --- 图 2: Power ---
p_power <- ggplot(final_results, aes(x = Sample_Size, y = Power, color = Method, group = Method, shape = Method)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 3) +
  
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  
  scale_x_continuous(breaks = sample_size_scenarios) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.05)) +
  
  labs(title = NULL, x = "Sample Size per Group", y = "Power (TPR)") +
  clean_theme

# --- 组合图形 ---
combined_plot <- p_fdr + p_power + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",
        legend.title = element_blank())

# --- 打印与保存 ---
print(combined_plot)



## 根据 N_total 和 prior_d0 自动生成文件名
file_name <- sprintf(
  "simulation/results/longitudinal_performance_plot_final_N%d_d%d.png",
  N_total, prior_d0
)

ggsave(
  file_name,
  plot   = combined_plot,
  width  = 10,
  height = 5,
  dpi    = 300,
  bg     = "white"
)

cat(paste("图表已保存为:", file_name, "\n"))

