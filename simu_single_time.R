# 用于生成模拟结果
# ===================================================================
# Step 0: 加载必要的包
# ===================================================================
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork") # 用于组合图形

library(limma)
library(dplyr)
library(ggplot2)
library(patchwork)

# ===================================================================
# Step 1: 核心模拟函数 (此部分无需修改)
# ===================================================================
run_simulation <- function(params) {
  sigma2_g <- (params$d0 * params$s0_2) / rchisq(params$G, df = params$d0)
  beta_g <- rep(0, params$G)
  de_indices <- sample(1:params$G, params$n_DE)
  beta_g[de_indices] <- rnorm(params$n_DE, mean = 0, sd = sqrt(params$v0 * sigma2_g[de_indices]))
  
  Y <- matrix(0, nrow = params$G, ncol = params$n_samples)
  for (g in 1:params$G) {
    true_means <- params$design %*% c(8, beta_g[g])
    Y[g, ] <- rnorm(params$n_samples, mean = true_means, sd = sqrt(sigma2_g[g]))
  }
  
  fit <- lmFit(Y, params$design)
  fit_contrast <- contrasts.fit(fit, contrasts = params$contrast_matrix)
  fit_mod <- eBayes(fit_contrast)
  
  unmoderated_t <- fit_contrast$coefficients[, 1] / (fit_contrast$stdev.unscaled[, 1] * fit_contrast$sigma)
  df_residual <- fit_contrast$df.residual[1]
  p_val_ord_t <- 2 * pt(abs(unmoderated_t), df = df_residual, lower.tail = FALSE)
  
  is_de <- rep("Non-DE", params$G)
  is_de[de_indices] <- "DE"
  is_de <- factor(is_de, levels = c("Non-DE", "DE"))
  
  results <- data.frame(
    status = is_de,
    p_ord_t = p_val_ord_t,
    p_mod_t = fit_mod$p.value[, 1]
  )
  return(results)
}

# ===================================================================
# Step 2: 定义模拟参数
# ===================================================================
# --- 循环参数 ---
sample_size_scenarios <- 2:10 # 每组的样本量范围
d0_scenarios <- c(1, 4, 1000) # d0场景

# --- 固定参数 ---
n_reps <- 100               # 重复模拟次数 (可设为100以获得更平滑曲线)
fdr_level <- 0.05          # FDR阈值
set.seed(2025)             # 保证结果可复现

base_params <- list(
  G = 10000,               # 总基因数
  n_DE = 100,              # 真实差异基因数
  s0_2 = 1,                # 方差的先验尺度 (设为1，挑战性更高)
  v0 = 10                  # 效应大小的方差
)

# ===================================================================
# Step 3: 运行完整的嵌套模拟循环
# ===================================================================
# 创建一个空的data.frame来存储所有结果
performance_summary <- data.frame()

cat("--- 开始运行模拟... ---\n")

# 外层循环：遍历样本量
for (n_samp in sample_size_scenarios) {
  
  # 内层循环：遍历d0场景
  for (d0_val in d0_scenarios) {
    
    cat(paste0("正在运行: d0 = ", d0_val, ", 样本量 = ", n_samp, " vs ", n_samp, "\n"))
    
    # 为当前样本量设置参数
    current_params <- c(base_params, list(
      d0 = d0_val,
      n_samples = n_samp * 2,
      group = factor(c(rep("Control", n_samp), rep("Case", n_samp))),
      design = model.matrix(~factor(c(rep("Control", n_samp), rep("Case", n_samp)))),
      contrast_matrix = matrix(c(0, 1), ncol = 1)
    ))
    
    # 用于存储当前场景下n_reps次重复的结果
    ord_t_conf_list <- list()
    mod_t_conf_list <- list()
    
    # 重复n_reps次
    for (i in 1:n_reps) {
      sim_data <- run_simulation(current_params)
      
      adj_p_ord_t <- p.adjust(sim_data$p_ord_t, method = "BH")
      adj_p_mod_t <- p.adjust(sim_data$p_mod_t, method = "BH")
      
      is_de <- sim_data$status == "DE"
      
      pred_ord_t <- adj_p_ord_t < fdr_level
      ord_t_conf_list[[i]] <- data.frame(
        TP = sum(is_de & pred_ord_t), FP = sum(!is_de & pred_ord_t)
      )
      
      pred_mod_t <- adj_p_mod_t < fdr_level
      mod_t_conf_list[[i]] <- data.frame(
        TP = sum(is_de & pred_mod_t), FP = sum(!is_de & pred_mod_t)
      )
    }
    
    # 计算平均TP和FP
    avg_ord_t_counts <- bind_rows(ord_t_conf_list) %>% summarise(across(everything(), mean))
    avg_mod_t_counts <- bind_rows(mod_t_conf_list) %>% summarise(across(everything(), mean))
    
    # 计算FDR和TPR (Power)
    calc_metrics <- function(avg_counts) {
      FDR <- avg_counts$FP / (avg_counts$FP + avg_counts$TP)
      FDR[is.nan(FDR)] <- 0
      TPR <- avg_counts$TP / base_params$n_DE
      return(list(FDR = FDR, TPR = TPR))
    }
    
    metrics_ord <- calc_metrics(avg_ord_t_counts)
    metrics_mod <- calc_metrics(avg_mod_t_counts)
    
    # 将结果添加到汇总表中
    performance_summary <- bind_rows(performance_summary,
                                     data.frame(d0 = d0_val, n_samples_per_group = n_samp, Method = "Ordinary t-test", FDR = metrics_ord$FDR, TPR = metrics_ord$TPR),
                                     data.frame(d0 = d0_val, n_samples_per_group = n_samp, Method = "Moderated t-test", FDR = metrics_mod$FDR, TPR = metrics_mod$TPR)
    )
  }
}
cat("--- 模拟完成！正在生成图形... ---\n")


# # ===================================================================
# # Step 4: 为每个d0场景生成并排图并保存
# # ===================================================================
# 
# for (d0_val in d0_scenarios) {
#   
#   # 筛选出当前d0场景的数据
#   plot_data <- performance_summary %>% 
#     filter(d0 == d0_val) %>%
#     mutate(Method = factor(Method, levels = c("Moderated t-test", "Ordinary t-test"))) # 调整顺序
#   
#   # 定义颜色
#   method_colors <- c("Ordinary t-test" = "#d95f02", "Moderated t-test" = "#1b9e77")
#   
#   # 图1: FDR vs. 样本量
#   p_fdr <- ggplot(plot_data, aes(x = n_samples_per_group, y = FDR, color = Method, group = Method)) +
#     geom_line(linewidth = 1.2) +
#     geom_point(size = 3) +
#     geom_hline(yintercept = fdr_level, linetype = "dashed", color = "red") +
#     annotate("text", x = max(sample_size_scenarios), y = fdr_level, label = "FDR = 0.05", vjust = -0.5, hjust = 1, color = "red") +
#     scale_color_manual(values = method_colors) +
#     scale_y_continuous(limits = c(0, max(0.1, max(plot_data$FDR, na.rm = TRUE) * 1.1))) +
#     labs(
#       title = "False Discovery Rate (FDR)",
#       x = "Sample Size per Group",
#       y = "Observed FDR"
#     ) +
#     theme_minimal() +
#     theme(legend.position = "none")
#   
#   # 图2: Power (TPR) vs. 样本量
#   p_tpr <- ggplot(plot_data, aes(x = n_samples_per_group, y = TPR, color = Method, group = Method)) +
#     geom_line(linewidth = 1.2) +
#     geom_point(size = 3) +
#     scale_color_manual(values = method_colors) +
#     scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
#     labs(
#       title = "Power (True Positive Rate)",
#       x = "Sample Size per Group",
#       y = "Power / TPR"
#     ) +
#     theme_minimal()
#   
#   # 使用patchwork组合两个图
#   combined_plot <- p_fdr + p_tpr + 
#     plot_layout(guides = "collect") + # 收集图例并放在一起
#     plot_annotation(
#       # title = paste0("Performance Comparison for d0 = ", d0_val),
#       # subtitle = paste0("d0=", d0_val, " represents ", 
#       #                   ifelse(d0_val < 10, "highly variable", "stable"), 
#       #                   " gene variances across samples."),
#       theme = theme(plot.title = element_text(size = 16, face = "bold"))
#     )
#   
#   # *** 新增部分：保存图形 ***
#   # 1. 定义动态文件名
#   file_name <- paste0("simulation/performance_plot_d0_", d0_val, ".png")
#   
#   # 2. 使用ggsave保存
#   ggsave(
#     file_name,
#     plot = combined_plot,
#     width = 10,          # 宽度 (英寸)
#     height = 5,          # 高度 (英寸)
#     dpi = 300,           # 分辨率
#     bg = "white"         # 设置背景为白色，避免透明背景
#   )
#   
#   # 打印组合图到R的绘图窗口
#   print(combined_plot)
#   
#   # 在控制台输出提示信息
#   cat(paste("图已保存为:", file_name, "\n"))
# }
# 



# combine plot 
# ===================================================================
# Step 4: 生成带形状区分的网格组合图
# ===================================================================
cat("--- 模拟完成！正在生成组合图形... ---\n")

performance_summary <- performance_summary %>%
  mutate(
    d0_label = factor(paste0("d0 = ", d0), levels = paste0("d0 = ", d0_scenarios)),
    Method = factor(Method, levels = c("Moderated t-test", "Ordinary t-test"))
  )
#save(performance_summary,file = "simulation/performance_summary.RData")

# --- 1. 定义颜色和形状 (与Longitudinal一致) ---
method_colors <- c("Ordinary t-test" = "#d95f02", 
                   "Moderated t-test" = "#1b9e77")

# 形状: 16=实心圆, 17=实心三角
method_shapes <- c("Ordinary t-test" = 16, 
                   "Moderated t-test" = 17)

# --- 2. 创建绘图函数 ---
create_plot <- function(metric, d0_val, show_y_axis = TRUE, show_x_axis = TRUE) {
  
  plot_data <- performance_summary %>% filter(d0 == d0_val)
  
  # 在 aes 中添加 shape = Method
  p <- ggplot(plot_data, aes(x = n_samples_per_group, color = Method, group = Method, shape = Method)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.5) + # 稍微调大点的大小以看清形状
    
    # 应用颜色和形状
    scale_color_manual(values = method_colors) +
    scale_shape_manual(values = method_shapes) +
    
    labs(x = if(show_x_axis) "Sample Size per Group" else NULL) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none", 
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      # 保持之前的美化调整 (去除次网格，调淡主网格)
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.5)
    )
  
  if (metric == "FDR") {
    p <- p + 
      aes(y = FDR) +
      geom_hline(yintercept = fdr_level, linetype = "dashed", color = "red") +
      scale_y_continuous(limits = c(0, max(0.12, max(plot_data$FDR, na.rm = TRUE) * 1.1))) +
      labs(
        title = bquote(d[0] == .(d0_val)), 
        y = if(show_y_axis) "Empirical FDR" else NULL
      )
  } else { 
    p <- p + 
      aes(y = TPR) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(
        title = NULL, 
        y = if(show_y_axis) "Power (TPR)" else NULL
      )
  }
  
  if (!show_y_axis) {
    p <- p + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  
  if (!show_x_axis) {
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  return(p)
}

# --- 3. 生成子图 ---
p_fdr1 <- create_plot("FDR", 1,    TRUE,  FALSE)
p_fdr2 <- create_plot("FDR", 4,    FALSE, FALSE)
p_fdr3 <- create_plot("FDR", 1000, FALSE, FALSE)

p_tpr1 <- create_plot("TPR", 1,    TRUE,  TRUE)
p_tpr2 <- create_plot("TPR", 4,    FALSE, TRUE)
p_tpr3 <- create_plot("TPR", 1000, FALSE, TRUE)

# --- 4. 组合图形 ---
final_plot <- (p_fdr1 + p_fdr2 + p_fdr3) / (p_tpr1 + p_tpr2 + p_tpr3) +
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom", 
    legend.title = element_blank(), 
    legend.text = element_text(size = 12)
  )

# --- 5. 保存 ---
print(final_plot)

ggsave(
  "simulation/results/performance_grid_plot.png",
  plot = final_plot,
  width = 14,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("组合图已保存为: performance_grid_plot.png\n")

