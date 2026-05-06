# ===================================================================
# Step 0: 加载必要的包
# ===================================================================
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("e1071", quietly = TRUE)) install.packages("e1071") 
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")  
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")

library(limma)
library(e1071)
library(pROC)
library(dplyr)
library(ggplot2)
library(tidyr)

# ===================================================================
# Step 1: 数据生成函数 (N=3 vs N=3, 更具挑战性的参数)
# ===================================================================
generate_data_small_n <- function(G=3000, time_points=1:10, 
                                  n_DE=30, 
                                  d0=100,        # <--- 修改1: 降低d0 (增加方差异质性，模拟真实的高噪声环境)
                                  s0_2=0.09) {
  
  n_case <- 3
  n_control <- 3
  n_subjects <- n_case + n_control
  n_time_points <- length(time_points)
  n_samples <- n_subjects * n_time_points
  
  # 构建元数据
  groups <- c(rep("Case", n_case), rep("Control", n_control))
  subject_ids <- paste0("S", 1:n_subjects)
  
  col_info <- expand.grid(Time = time_points, SubjectID = subject_ids)
  col_info <- col_info[order(match(col_info$SubjectID, subject_ids), col_info$Time), ]
  col_info$Group <- rep(groups, each = n_time_points)
  
  # --- 信号生成 ---
  peak_time <- 5.5
  signal_shape <- dnorm(time_points, mean = peak_time, sd = 2.5)
  signal_shape <- signal_shape / max(signal_shape)
  
  true_log2fc_matrix <- matrix(0, nrow = G, ncol = n_time_points)
  
  # <--- 修改2: 稍微降低 LogFC (0.5-1.0)，增加检测难度
  peak_log2fc <- runif(n_DE, min = 0.5, max = 1.0) 
  peak_log2fc[1:(n_DE/2)] <- -peak_log2fc[1:(n_DE/2)]
  true_log2fc_matrix[1:n_DE, ] <- peak_log2fc %o% signal_shape
  
  # 1. 基线
  base_expression <- rnorm(G, mean = 8, sd = 2)
  
  # 2. 个体效应
  subject_eff_matrix <- matrix(rnorm(G * n_subjects, 0, 0.2), nrow = G)
  subject_component <- subject_eff_matrix[, rep(1:n_subjects, each = n_time_points)]
  
  # 3. 组别效应
  group_component <- matrix(0, nrow = G, ncol = n_samples)
  case_indices <- which(col_info$Group == "Case")
  group_component[, case_indices] <- true_log2fc_matrix[, col_info$Time[case_indices]]
  
  # 4. 误差 (Inverse Chi-square)
  # 当 d0 很小时，这里会产生一些极小和极大的方差
  sigma2_g <- (d0 * s0_2) / rchisq(G, df = d0)
  sd_g <- sqrt(sigma2_g)
  error_component <- matrix(rnorm(G * n_samples), nrow = G) * sd_g
  
  LogAbundance <- base_expression + subject_component + group_component + error_component
  
  list(data = LogAbundance, col_info = col_info)
}

# ===================================================================
# Step 2: 模拟参数设置
# ===================================================================
set.seed(3476)
N_simulations <- 100       
# <--- 修改3: 扩大 Top K 范围，直到 50
top_k_seq <- 1:10

cv_results_list <- list()

cat(sprintf("开始 9-Fold CV 预测模拟 (N=6, Reps=%d, TopK max=50)...\n", N_simulations))

# ===================================================================
# Step 3: 主模拟循环 (保持原有逻辑)
# ===================================================================
for (sim in 1:N_simulations) {
  if (sim %% 5 == 0) cat(sprintf("  Simulation %d / %d\n", sim, N_simulations))
  
  dataset <- generate_data_small_n()
  Y_full <- dataset$data
  info_full <- dataset$col_info
  
  # Case: 1,2,3; Control: 4,5,6
  case_combos <- combn(1:3, 2, simplify = FALSE) 
  ctrl_combos <- combn(4:6, 2, simplify = FALSE) 
  
  folds <- list()
  fold_idx <- 1
  for (cc in case_combos) {
    for (ct in ctrl_combos) {
      train_subjects_idx <- c(cc, ct)
      test_subjects_idx <- setdiff(1:6, train_subjects_idx)
      
      get_col_indices <- function(subj_indices) {
        which(info_full$SubjectID %in% paste0("S", subj_indices))
      }
      
      folds[[fold_idx]] <- list(
        train_cols = get_col_indices(train_subjects_idx),
        test_cols = get_col_indices(test_subjects_idx)
      )
      fold_idx <- fold_idx + 1
    }
  }
  
  for (f_i in 1:length(folds)) {
    fold <- folds[[f_i]]
    Y_train <- Y_full[, fold$train_cols]
    info_train <- info_full[fold$train_cols, ]
    Y_test <- Y_full[, fold$test_cols]
    info_test <- info_full[fold$test_cols, ]
    
    # --- 特征选择 ---
    design_factor <- factor(paste(info_train$Group, info_train$Time, sep="."))
    design_mat <- model.matrix(~0 + design_factor)
    colnames(design_mat) <- levels(design_factor)
    
    corfit <- tryCatch({
      duplicateCorrelation(Y_train, design = design_mat, block = info_train$SubjectID)
    }, error = function(e) list(consensus = 0.5)) 
    
    fit <- lmFit(Y_train, design_mat, block = info_train$SubjectID, correlation = corfit$consensus)
    
    contrast_str <- paste(paste("Case.", 1:10, sep=""), "-", paste("Control.", 1:10, sep=""))
    cont_mat <- makeContrasts(contrasts = contrast_str, levels = design_mat)
    fit_cont <- contrasts.fit(fit, cont_mat)
    eb_fit <- eBayes(fit_cont, robust = TRUE, trend = TRUE)
    
    # Method 1: Mod F-test
    f_res <- topTable(eb_fit, coef = 1:10, number = Inf, sort.by = "none")
    rank_f <- as.numeric(rownames(f_res)) 
    
    # Method 2: Mod t-test (Min-P)
    min_p_mod <- apply(eb_fit$p.value, 1, min)
    rank_mod_t <- order(min_p_mod)
    
    # Method 3: Ord t-test (Min-P)
    ord_t_stat <- fit_cont$coefficients / (fit_cont$stdev.unscaled * fit$sigma)
    ord_t_pvals <- 2 * pt(-abs(ord_t_stat), df = fit$df.residual)
    min_p_ord <- apply(ord_t_pvals, 1, min)
    rank_ord_t <- order(min_p_ord)
    
    # --- 预测 ---
    Y_train_t <- t(Y_train)
    Y_test_t <- t(Y_test)
    train_labels <- factor(info_train$Group, levels = c("Control", "Case"))
    test_labels <- factor(info_test$Group, levels = c("Control", "Case"))
    
    evaluate_method <- function(ranked_indices, method_name) {
      aucs <- numeric(length(top_k_seq))
      for (k_idx in seq_along(top_k_seq)) {
        k <- top_k_seq[k_idx]
        top_k_idx <- ranked_indices[1:k]
        
        train_df <- as.data.frame(Y_train_t[, top_k_idx, drop=FALSE])
        colnames(train_df) <- paste0("Gene", 1:k)
        train_df$Label <- train_labels
        
        test_df <- as.data.frame(Y_test_t[, top_k_idx, drop=FALSE])
        colnames(test_df) <- paste0("Gene", 1:k)
        
        # 捕获 Naive Bayes 可能的错误
        tryCatch({
          model <- naiveBayes(Label ~ ., data = train_df)
          preds <- predict(model, test_df, type = "raw")[, "Case"]
          
          if (length(unique(test_labels)) < 2) {
            aucs[k_idx] <- NA 
          } else {
            roc_obj <- roc(test_labels, preds, levels = c("Control", "Case"), direction = "<", quiet = TRUE)
            aucs[k_idx] <- as.numeric(roc_obj$auc)
          }
        }, error = function(e) { aucs[k_idx] <- NA })
      }
      return(data.frame(Sim = sim, Fold = f_i, Method = method_name, k = top_k_seq, AUC = aucs))
    }
    
    res_f <- evaluate_method(rank_f, "Moderated F-test")
    res_mod <- evaluate_method(rank_mod_t, "Moderated t-test")
    res_ord <- evaluate_method(rank_ord_t, "Ordinary t-test")
    
    cv_results_list[[length(cv_results_list) + 1]] <- bind_rows(res_f, res_mod, res_ord)
  }
}

# ===================================================================
# Step 4: 汇总与绘图
# ===================================================================
all_res <- bind_rows(cv_results_list)

final_summary <- all_res %>%
  group_by(Method, k) %>%
  summarise(
    Mean_AUC = mean(AUC, na.rm = TRUE),
    SE_AUC = sd(AUC, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

cat("模拟完成，正在绘图...\n")

method_colors <- c("Ordinary t-test" = "#d95f02", 
                   "Moderated t-test" = "#1b9e77", 
                   "Moderated F-test" = "#7570b3")
method_shapes <- c("Ordinary t-test" = 16, 
                   "Moderated t-test" = 17, 
                   "Moderated F-test" = 15)

final_summary$Method <- factor(final_summary$Method, 
                               levels = c("Moderated F-test", "Moderated t-test", "Ordinary t-test"))

p_auc <- ggplot(final_summary, aes(x = k, y = Mean_AUC, color = Method, group = Method, shape = Method)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.0) +
  geom_errorbar(aes(ymin = Mean_AUC - SE_AUC, ymax = Mean_AUC + SE_AUC), width = 0.5, alpha = 0.5) +
  
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  
  scale_x_continuous(breaks = seq(0, 10, by = 1)) + # X轴刻度调整
  scale_y_continuous(limits = c(0.4, 1.0), n.breaks = 6) +
  
  labs(
    title = NULL,
    x = "Number of Top Proteins Included in Model",
    y = "AUC"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray92", linewidth = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  )

print(p_auc)
## 把当前使用的 G 和 d0 读出来（从函数默认值）
G_val  <- formals(generate_data_small_n)$G
d0_val <- formals(generate_data_small_n)$d0

## 生成带 G 和 d0 的文件名，例如：auc_cv_G2000_d100.png
out_png <- sprintf("simulation/results/auc_cv_G%d_d%d.png", G_val, d0_val)

#print(p_auc)
ggsave(out_png, plot = p_auc, width = 8, height = 6, bg = "white")
cat("图表已保存为:", out_png, "\n")