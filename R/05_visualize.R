# 可视化：暴露-结局分布、趋势与森林图

source("R/00_packages.R")

analysis_path <- "data/analysis.rds"
if (!file.exists(analysis_path)) {
  stop("Missing data/analysis.rds. Run R/03_egfr.R first.")
}

analysis <- readRDS(analysis_path)

out_dir <- "report/_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!"pb_q5" %in% names(analysis)) {
  analysis <- analysis %>%
    mutate(
      pb_q5 = ntile(log_pb, 5),
      cd_q5 = ntile(log_cd, 5)
    )
}

# 密度热力图 + 线性趋势
plot_pb_hex <- analysis %>%
  ggplot(aes(log_pb, egfr)) +
  stat_bin2d(bins = 35) +
  scale_fill_gradient(low = "#edf2f7", high = "#2b6cb0") +
  geom_smooth(method = "lm", se = FALSE, color = "#2b6cb0", linewidth = 1) +
  labs(x = "log(Lead, ug/dL)", y = "eGFR", title = "Lead and eGFR (density)")

ggsave(filename = file.path(out_dir, "hex_log_pb_egfr.png"), plot = plot_pb_hex,
       width = 6, height = 4, dpi = 300)

plot_cd_hex <- analysis %>%
  ggplot(aes(log_cd, egfr)) +
  stat_bin2d(bins = 35) +
  scale_fill_gradient(low = "#fff7ed", high = "#c05621") +
  geom_smooth(method = "lm", se = FALSE, color = "#c05621", linewidth = 1) +
  labs(x = "log(Cadmium, ug/L)", y = "eGFR", title = "Cadmium and eGFR (density)")

ggsave(filename = file.path(out_dir, "hex_log_cd_egfr.png"), plot = plot_cd_hex,
       width = 6, height = 4, dpi = 300)

# 暴露分布（log 变换）
plot_log_pb <- analysis %>%
  ggplot(aes(log_pb)) +
  geom_histogram(bins = 30, fill = "#2b6cb0", color = "white") +
  labs(x = "log(Lead, ug/dL)", y = "Count", title = "Distribution of log blood lead")

ggsave(filename = file.path(out_dir, "dist_log_pb.png"), plot = plot_log_pb,
       width = 6, height = 4, dpi = 300)

plot_log_cd <- analysis %>%
  ggplot(aes(log_cd)) +
  geom_histogram(bins = 30, fill = "#c05621", color = "white") +
  labs(x = "log(Cadmium, ug/L)", y = "Count", title = "Distribution of log blood cadmium")

ggsave(filename = file.path(out_dir, "dist_log_cd.png"), plot = plot_log_cd,
       width = 6, height = 4, dpi = 300)

# 散点 + 线性拟合（直观看趋势）
plot_pb_scatter <- analysis %>%
  ggplot(aes(log_pb, egfr)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "#2b6cb0", linewidth = 1) +
  labs(x = "log(Lead, ug/dL)", y = "eGFR", title = "Lead and eGFR (linear fit)")

ggsave(filename = file.path(out_dir, "scatter_log_pb_egfr.png"), plot = plot_pb_scatter,
       width = 6, height = 4, dpi = 300)

plot_cd_scatter <- analysis %>%
  ggplot(aes(log_cd, egfr)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "#c05621", linewidth = 1) +
  labs(x = "log(Cadmium, ug/L)", y = "eGFR", title = "Cadmium and eGFR (linear fit)")

ggsave(filename = file.path(out_dir, "scatter_log_cd_egfr.png"), plot = plot_cd_scatter,
       width = 6, height = 4, dpi = 300)

# 镉的样条曲线（非线性形状）
# 取类别型变量的众数
mode_value <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}

# 变量是否可用于建模（避免单一水平）
is_usable_var <- function(data, var) {
  if (!var %in% names(data)) return(FALSE)
  x <- data[[var]]
  if (is.factor(x)) return(nlevels(droplevels(x)) >= 2)
  any(!is.na(x))
}

# 样条模型协变量（会自动剔除不可用项）
cd_covars <- c("log_pb", "RIDAGEYR", "sex", "race", "education", "INDFMPIR", "bmi", "smoking")
cd_covars <- cd_covars[vapply(cd_covars, is_usable_var, logical(1), data = analysis)]
rhs_terms <- c(paste0("splines::ns(log_cd, df = 3)"), cd_covars)
formula_cd_spline <- as.formula(paste("egfr ~", paste(rhs_terms, collapse = " + ")))

# 固定其他协变量，画镉的非线性曲线
fit_cd_spline <- lm(formula_cd_spline, data = analysis)
grid_cd <- seq(min(analysis$log_cd, na.rm = TRUE), max(analysis$log_cd, na.rm = TRUE), length.out = 100)

newdata_cd <- data.frame(log_cd = grid_cd)
for (v in cd_covars) {
  x <- analysis[[v]]
  if (is.numeric(x)) {
    newdata_cd[[v]] <- median(x, na.rm = TRUE)
  } else {
    mv <- mode_value(x)
    newdata_cd[[v]] <- factor(mv, levels = levels(x))
  }
}

newdata_cd$spline_pred <- predict(fit_cd_spline, newdata = newdata_cd)

# 背景用密度，曲线用样条预测
plot_cd_spline <- ggplot(analysis, aes(log_cd, egfr)) +
  stat_bin2d(bins = 35) +
  scale_fill_gradient(low = "#fff7ed", high = "#c05621") +
  geom_line(data = newdata_cd, aes(y = spline_pred), color = "#2f855a", linewidth = 1) +
  labs(x = "log(Cadmium, ug/L)", y = "eGFR",
       title = "Cadmium and eGFR (spline curve)") +
  theme(legend.position = "right")

ggsave(filename = file.path(out_dir, "cadmium_spline_gam.png"), plot = plot_cd_spline,
       width = 6, height = 4, dpi = 300)

# 分位数箱线图（稳健比较分布差异）
plot_pb_box <- analysis %>%
  mutate(pb_q5 = factor(pb_q5, levels = 1:5, labels = paste0("Q", 1:5))) %>%
  ggplot(aes(pb_q5, egfr)) +
  geom_boxplot(outlier.alpha = 0.3, fill = "#bee3f8") +
  labs(x = "Lead quintile (Q1-Q5)", y = "eGFR", title = "eGFR by lead quintiles")

ggsave(filename = file.path(out_dir, "box_egfr_by_lead_q5.png"), plot = plot_pb_box,
       width = 6, height = 4, dpi = 300)

plot_cd_box <- analysis %>%
  mutate(cd_q5 = factor(cd_q5, levels = 1:5, labels = paste0("Q", 1:5))) %>%
  ggplot(aes(cd_q5, egfr)) +
  geom_boxplot(outlier.alpha = 0.3, fill = "#fed7aa") +
  labs(x = "Cadmium quintile (Q1-Q5)", y = "eGFR", title = "eGFR by cadmium quintiles")

ggsave(filename = file.path(out_dir, "box_egfr_by_cadmium_q5.png"), plot = plot_cd_box,
       width = 6, height = 4, dpi = 300)

# 四分位趋势（均值）
plot_pb <- analysis %>%
  mutate(pb_q = ntile(LBXBPB, 4)) %>%
  group_by(pb_q) %>%
  summarize(mean_egfr = mean(egfr, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(pb_q, mean_egfr)) +
  geom_line() +
  geom_point() +
  labs(x = "Blood lead quartile (Q1-Q4)", y = "Mean eGFR")

ggsave(filename = file.path(out_dir, "lead_quartile_egfr.png"), plot = plot_pb,
       width = 6, height = 4, dpi = 300)

# 四分位趋势（均值）
plot_cd <- analysis %>%
  mutate(cd_q = ntile(LBXBCD, 4)) %>%
  group_by(cd_q) %>%
  summarize(mean_egfr = mean(egfr, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(cd_q, mean_egfr)) +
  geom_line() +
  geom_point() +
  labs(x = "Blood cadmium quartile (Q1-Q4)", y = "Mean eGFR")

ggsave(filename = file.path(out_dir, "cadmium_quartile_egfr.png"), plot = plot_cd,
       width = 6, height = 4, dpi = 300)

# 五分位趋势（未加权）
plot_pb_q5 <- analysis %>%
  group_by(pb_q5) %>%
  summarize(mean_egfr = mean(egfr, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(pb_q5, mean_egfr)) +
  geom_line() + geom_point() +
  labs(x = "Lead quintile (Q1-Q5)", y = "Mean eGFR",
       title = "Mean eGFR by lead quintiles")

ggsave(filename = file.path(out_dir, "trend_egfr_lead_q5.png"), plot = plot_pb_q5,
       width = 6, height = 4, dpi = 300)

plot_cd_q5 <- analysis %>%
  group_by(cd_q5) %>%
  summarize(mean_egfr = mean(egfr, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(cd_q5, mean_egfr)) +
  geom_line() + geom_point() +
  labs(x = "Cadmium quintile (Q1-Q5)", y = "Mean eGFR",
       title = "Mean eGFR by cadmium quintiles")

ggsave(filename = file.path(out_dir, "trend_egfr_cadmium_q5.png"), plot = plot_cd_q5,
       width = 6, height = 4, dpi = 300)

# 五分位趋势（加权）
design <- survey::svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  nest = TRUE,
  data = analysis
)

pb_q5_svy <- survey::svyby(~egfr, ~factor(pb_q5), design, survey::svymean, na.rm = TRUE)
pb_q5_svy <- pb_q5_svy %>% rename(pb_q5 = `factor(pb_q5)`)
plot_pb_q5_svy <- ggplot(pb_q5_svy, aes(x = pb_q5, y = egfr)) +
  geom_line(group = 1) + geom_point() +
  labs(x = "Lead quintile (Q1-Q5)", y = "Weighted mean eGFR",
       title = "Weighted mean eGFR by lead quintiles")

ggsave(filename = file.path(out_dir, "trend_egfr_lead_q5_weighted.png"), plot = plot_pb_q5_svy,
       width = 6, height = 4, dpi = 300)

cd_q5_svy <- survey::svyby(~egfr, ~factor(cd_q5), design, survey::svymean, na.rm = TRUE)
cd_q5_svy <- cd_q5_svy %>% rename(cd_q5 = `factor(cd_q5)`)
plot_cd_q5_svy <- ggplot(cd_q5_svy, aes(x = cd_q5, y = egfr)) +
  geom_line(group = 1) + geom_point() +
  labs(x = "Cadmium quintile (Q1-Q5)", y = "Weighted mean eGFR",
       title = "Weighted mean eGFR by cadmium quintiles")

ggsave(filename = file.path(out_dir, "trend_egfr_cadmium_q5_weighted.png"), plot = plot_cd_q5_svy,
       width = 6, height = 4, dpi = 300)

# CKD 分位数患病率
ckd_by_pb <- analysis %>%
  mutate(pb_q = ntile(LBXBPB, 4)) %>%
  group_by(pb_q) %>%
  summarize(ckd_rate = mean(ckd, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(pb_q, ckd_rate)) +
  geom_col(fill = "#2b6cb0", width = 0.6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Blood lead quartile (Q1-Q4)", y = "CKD prevalence",
       title = "CKD prevalence by lead quartiles")

ggsave(filename = file.path(out_dir, "ckd_by_lead_quartile.png"), plot = ckd_by_pb,
       width = 6, height = 4, dpi = 300)

ckd_by_cd <- analysis %>%
  mutate(cd_q = ntile(LBXBCD, 4)) %>%
  group_by(cd_q) %>%
  summarize(ckd_rate = mean(ckd, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(cd_q, ckd_rate)) +
  geom_col(fill = "#c05621", width = 0.6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Blood cadmium quartile (Q1-Q4)", y = "CKD prevalence",
       title = "CKD prevalence by cadmium quartiles")

ggsave(filename = file.path(out_dir, "ckd_by_cadmium_quartile.png"), plot = ckd_by_cd,
       width = 6, height = 4, dpi = 300)

# CKD OR 森林图（未加权）
models_tidy <- readRDS("data/models_tidy.rds")

or_unweighted <- models_tidy$fit_glm_full %>%
  filter(term %in% c("log_pb", "log_cd")) %>%
  mutate(term = recode(term, log_pb = "Lead (log)", log_cd = "Cadmium (log)"))

plot_or_unweighted <- ggplot(or_unweighted, aes(x = estimate, y = term)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_log10() +
  labs(x = "Odds Ratio (log scale)", y = "", title = "CKD ORs (unweighted, full model)")

ggsave(filename = file.path(out_dir, "ckd_or_forest_unweighted.png"), plot = plot_or_unweighted,
       width = 6, height = 3, dpi = 300)

# CKD OR 森林图（加权）
or_weighted <- models_tidy$fit_svy_ckd_full %>%
  filter(term %in% c("log_pb", "log_cd")) %>%
  mutate(term = recode(term, log_pb = "Lead (log)", log_cd = "Cadmium (log)"))

plot_or_weighted <- ggplot(or_weighted, aes(x = estimate, y = term)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_log10() +
  labs(x = "Odds Ratio (log scale)", y = "", title = "CKD ORs (survey-weighted, full model)")

ggsave(filename = file.path(out_dir, "ckd_or_forest_weighted.png"), plot = plot_or_weighted,
       width = 6, height = 3, dpi = 300)

invisible(TRUE)
