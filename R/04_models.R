# 建模：eGFR 线性回归 + CKD Logistic + survey 加权

source("R/00_packages.R")

analysis_path <- "data/analysis.rds"
if (!file.exists(analysis_path)) {
  stop("Missing data/analysis.rds. Run R/03_egfr.R first.")
}

analysis <- readRDS(analysis_path)

# 判断变量是否可用于回归（因子>=2水平、数值有变异）
is_usable_covariate <- function(data, var) {
  if (!var %in% names(data)) return(FALSE)
  x <- data[[var]]
  if (is.factor(x)) {
    return(nlevels(droplevels(x)) >= 2)
  }
  uniq <- unique(x[!is.na(x)])
  length(uniq) >= 2
}

# 自动生成公式，剔除不合格变量
make_formula <- function(outcome, base_vars, data) {
  vars_ok <- base_vars[vapply(base_vars, is_usable_covariate, logical(1), data = data)]
  if (length(vars_ok) == 0) {
    stop("No usable covariates available for model: ", outcome)
  }
  as.formula(paste(outcome, "~", paste(vars_ok, collapse = " + ")))
}

# 生成样条模型公式（用于非线性检验）
make_spline_formula <- function(outcome, spline_var, base_vars, data, df = 3) {
  vars_ok <- base_vars[vapply(base_vars, is_usable_covariate, logical(1), data = data)]
  vars_ok <- unique(vars_ok)
  vars_ok <- setdiff(vars_ok, spline_var)
  rhs <- c(vars_ok, paste0("splines::ns(", spline_var, ", df = ", df, ")"))
  as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
}
# 分层回归（每个层单独拟合）
fit_by_strata <- function(data, strata_var, outcome, vars, family = NULL) {
  strata_vals <- na.omit(unique(data[[strata_var]]))
  res <- list()
  for (val in strata_vals) {
    d_sub <- data[data[[strata_var]] == val, , drop = FALSE]
    formula_sub <- make_formula(outcome, vars, d_sub)
    fit <- if (is.null(family)) {
      lm(formula_sub, data = d_sub)
    } else {
      glm(formula_sub, data = d_sub, family = family)
    }
    res[[as.character(val)]] <- fit
  }
  res
}

# 分层时去掉分层变量本身
fit_by_strata_drop <- function(data, strata_var, outcome, vars, family = NULL) {
  vars_use <- setdiff(vars, strata_var)
  fit_by_strata(data, strata_var, outcome, vars_use, family = family)
}

# 分层协变量：最小 -> 人口学 -> 生活方式
base_vars <- c("log_pb", "log_cd", "RIDAGEYR", "sex")
demo_vars <- c("race", "education", "INDFMPIR")
life_vars <- c("bmi", "smoking")

formula_min <- make_formula("egfr", base_vars, analysis)
formula_demo <- make_formula("egfr", c(base_vars, demo_vars), analysis)
formula_full <- make_formula("egfr", c(base_vars, demo_vars, life_vars), analysis)

formula_ckd_min <- make_formula("ckd", base_vars, analysis)
formula_ckd_demo <- make_formula("ckd", c(base_vars, demo_vars), analysis)
formula_ckd_full <- make_formula("ckd", c(base_vars, demo_vars, life_vars), analysis)

# 非加权线性模型（eGFR）
fit_lm_min <- lm(formula_min, data = analysis)
fit_lm_demo <- lm(formula_demo, data = analysis)
fit_lm_full <- lm(formula_full, data = analysis)

# 非加权 Logistic（CKD）
fit_glm_min <- glm(formula_ckd_min, data = analysis, family = binomial())
fit_glm_demo <- glm(formula_ckd_demo, data = analysis, family = binomial())
fit_glm_full <- glm(formula_ckd_full, data = analysis, family = binomial())

# winsorize 敏感性分析（截尾后的暴露）
base_vars_w <- c("log_pb_w", "log_cd_w", "RIDAGEYR", "sex")
formula_w_min <- make_formula("egfr", base_vars_w, analysis)
formula_w_demo <- make_formula("egfr", c(base_vars_w, demo_vars), analysis)
formula_w_full <- make_formula("egfr", c(base_vars_w, demo_vars, life_vars), analysis)

formula_ckd_w_min <- make_formula("ckd", base_vars_w, analysis)
formula_ckd_w_demo <- make_formula("ckd", c(base_vars_w, demo_vars), analysis)
formula_ckd_w_full <- make_formula("ckd", c(base_vars_w, demo_vars, life_vars), analysis)

fit_lm_w_min <- lm(formula_w_min, data = analysis)
fit_lm_w_demo <- lm(formula_w_demo, data = analysis)
fit_lm_w_full <- lm(formula_w_full, data = analysis)

fit_glm_w_min <- glm(formula_ckd_w_min, data = analysis, family = binomial())
fit_glm_w_demo <- glm(formula_ckd_w_demo, data = analysis, family = binomial())
fit_glm_w_full <- glm(formula_ckd_w_full, data = analysis, family = binomial())

# 镉非线性检验（样条）
vars_spline <- c("log_pb", "log_cd", "RIDAGEYR", "sex", "race", "education", "INDFMPIR", "bmi", "smoking")
formula_cd_spline <- make_spline_formula("egfr", "log_cd", vars_spline, analysis)
formula_cd_spline_ckd <- make_spline_formula("ckd", "log_cd", vars_spline, analysis)

fit_lm_cd_linear <- fit_lm_full
fit_lm_cd_spline <- lm(formula_cd_spline, data = analysis)

fit_glm_cd_linear <- fit_glm_full
fit_glm_cd_spline <- glm(formula_cd_spline_ckd, data = analysis, family = binomial())

# Survey 加权模型（NHANES 复杂抽样）
options(survey.lonely.psu = "adjust")
# 选择可用的 survey 设计（自动降阶避免样本为空或自由度为负）
pick_svy_design <- function(data, var_sets, label, weight_var = "WTMEC2YR") {
  for (vars in var_sets) {
    svy_vars <- unique(c(vars, "SDMVPSU", "SDMVSTRA", weight_var))
    d <- data %>%
      select(any_of(svy_vars)) %>%
      tidyr::drop_na() %>%
      mutate(across(where(is.factor), droplevels))

    if (nrow(d) == 0) {
      next
    }

    dsgn <- survey::svydesign(
      ids = ~SDMVPSU,
      strata = ~SDMVSTRA,
      weights = as.formula(paste0("~", weight_var)),
      nest = TRUE,
      data = d
    )

    deg <- survey::degf(dsgn)
    if (is.finite(deg) && deg > 0) {
      if (!identical(vars, var_sets[[1]])) {
        message("Survey fallback for ", label, ": using reduced covariates.")
      }
      return(list(vars = vars, design = dsgn))
    }
    message("Survey fallback for ", label, ": degf ", deg, " -> trying reduced covariates.")
  }
  stop("No valid survey design available after dropping missing values: ", label)
}

weight_var <- if ("WTMEC4YR" %in% names(analysis)) "WTMEC4YR" else "WTMEC2YR"

svy_egfr_min <- pick_svy_design(analysis, list(c("egfr", base_vars)), "egfr_min", weight_var = weight_var)
svy_egfr_demo <- pick_svy_design(analysis, list(c("egfr", base_vars, demo_vars), c("egfr", base_vars)), "egfr_demo", weight_var = weight_var)
svy_egfr_full <- pick_svy_design(analysis, list(c("egfr", base_vars, demo_vars, life_vars),
                                               c("egfr", base_vars, demo_vars),
                                               c("egfr", base_vars)), "egfr_full", weight_var = weight_var)

svy_ckd_min <- pick_svy_design(analysis, list(c("ckd", base_vars)), "ckd_min", weight_var = weight_var)
svy_ckd_demo <- pick_svy_design(analysis, list(c("ckd", base_vars, demo_vars), c("ckd", base_vars)), "ckd_demo", weight_var = weight_var)
svy_ckd_full <- pick_svy_design(analysis, list(c("ckd", base_vars, demo_vars, life_vars),
                                              c("ckd", base_vars, demo_vars),
                                              c("ckd", base_vars)), "ckd_full", weight_var = weight_var)

formula_svy_min <- make_formula("egfr", setdiff(svy_egfr_min$vars, "egfr"), analysis)
formula_svy_demo <- make_formula("egfr", setdiff(svy_egfr_demo$vars, "egfr"), analysis)
formula_svy_full <- make_formula("egfr", setdiff(svy_egfr_full$vars, "egfr"), analysis)

formula_svy_ckd_min <- make_formula("ckd", setdiff(svy_ckd_min$vars, "ckd"), analysis)
formula_svy_ckd_demo <- make_formula("ckd", setdiff(svy_ckd_demo$vars, "ckd"), analysis)
formula_svy_ckd_full <- make_formula("ckd", setdiff(svy_ckd_full$vars, "ckd"), analysis)

# 加权 eGFR 模型
fit_svy_min <- survey::svyglm(formula_svy_min, design = svy_egfr_min$design)
fit_svy_demo <- survey::svyglm(formula_svy_demo, design = svy_egfr_demo$design)
fit_svy_full <- survey::svyglm(formula_svy_full, design = svy_egfr_full$design)

# 加权 CKD 模型
fit_svy_ckd_min <- survey::svyglm(formula_svy_ckd_min, design = svy_ckd_min$design, family = quasibinomial())
fit_svy_ckd_demo <- survey::svyglm(formula_svy_ckd_demo, design = svy_ckd_demo$design, family = quasibinomial())
fit_svy_ckd_full <- survey::svyglm(formula_svy_ckd_full, design = svy_ckd_full$design, family = quasibinomial())

svy_egfr_w_min <- pick_svy_design(analysis, list(c("egfr", base_vars_w)), "egfr_w_min", weight_var = weight_var)
svy_egfr_w_demo <- pick_svy_design(analysis, list(c("egfr", base_vars_w, demo_vars), c("egfr", base_vars_w)), "egfr_w_demo", weight_var = weight_var)
svy_egfr_w_full <- pick_svy_design(analysis, list(c("egfr", base_vars_w, demo_vars, life_vars),
                                                 c("egfr", base_vars_w, demo_vars),
                                                 c("egfr", base_vars_w)), "egfr_w_full", weight_var = weight_var)

svy_ckd_w_min <- pick_svy_design(analysis, list(c("ckd", base_vars_w)), "ckd_w_min", weight_var = weight_var)
svy_ckd_w_demo <- pick_svy_design(analysis, list(c("ckd", base_vars_w, demo_vars), c("ckd", base_vars_w)), "ckd_w_demo", weight_var = weight_var)
svy_ckd_w_full <- pick_svy_design(analysis, list(c("ckd", base_vars_w, demo_vars, life_vars),
                                                c("ckd", base_vars_w, demo_vars),
                                                c("ckd", base_vars_w)), "ckd_w_full", weight_var = weight_var)

formula_svy_w_min <- make_formula("egfr", setdiff(svy_egfr_w_min$vars, "egfr"), analysis)
formula_svy_w_demo <- make_formula("egfr", setdiff(svy_egfr_w_demo$vars, "egfr"), analysis)
formula_svy_w_full <- make_formula("egfr", setdiff(svy_egfr_w_full$vars, "egfr"), analysis)

formula_svy_ckd_w_min <- make_formula("ckd", setdiff(svy_ckd_w_min$vars, "ckd"), analysis)
formula_svy_ckd_w_demo <- make_formula("ckd", setdiff(svy_ckd_w_demo$vars, "ckd"), analysis)
formula_svy_ckd_w_full <- make_formula("ckd", setdiff(svy_ckd_w_full$vars, "ckd"), analysis)

fit_svy_w_min <- survey::svyglm(formula_svy_w_min, design = svy_egfr_w_min$design)
fit_svy_w_demo <- survey::svyglm(formula_svy_w_demo, design = svy_egfr_w_demo$design)
fit_svy_w_full <- survey::svyglm(formula_svy_w_full, design = svy_egfr_w_full$design)

fit_svy_ckd_w_min <- survey::svyglm(formula_svy_ckd_w_min, design = svy_ckd_w_min$design, family = quasibinomial())
fit_svy_ckd_w_demo <- survey::svyglm(formula_svy_ckd_w_demo, design = svy_ckd_w_demo$design, family = quasibinomial())
fit_svy_ckd_w_full <- survey::svyglm(formula_svy_ckd_w_full, design = svy_ckd_w_full$design, family = quasibinomial())

fit_svy_cd_linear <- fit_svy_full
svy_cd_spline <- pick_svy_design(analysis, list(c("egfr", vars_spline),
                                               c("egfr", base_vars, demo_vars),
                                               c("egfr", base_vars)), "egfr_cd_spline", weight_var = weight_var)
svy_ckd_cd_spline <- pick_svy_design(analysis, list(c("ckd", vars_spline),
                                                   c("ckd", base_vars, demo_vars),
                                                   c("ckd", base_vars)), "ckd_cd_spline", weight_var = weight_var)

formula_svy_cd_spline <- make_spline_formula("egfr", "log_cd", setdiff(svy_cd_spline$vars, "egfr"), analysis)
formula_svy_ckd_cd_spline <- make_spline_formula("ckd", "log_cd", setdiff(svy_ckd_cd_spline$vars, "ckd"), analysis)
fit_svy_cd_spline <- survey::svyglm(formula_svy_cd_spline, design = svy_cd_spline$design)
fit_svy_ckd_cd_linear <- fit_svy_ckd_full
fit_svy_ckd_cd_spline <- survey::svyglm(formula_svy_ckd_cd_spline, design = svy_ckd_cd_spline$design, family = quasibinomial())

# 分层模型（探索性）
base_vars_sex <- setdiff(c("log_pb", "log_cd", "RIDAGEYR", "sex"), "sex")
full_vars_sex <- c(base_vars_sex, demo_vars, life_vars)
strata_sex_lm <- fit_by_strata(analysis, "sex", "egfr", full_vars_sex)
strata_sex_glm <- fit_by_strata(analysis, "sex", "ckd", full_vars_sex, family = binomial())

base_vars_age <- c("log_pb", "log_cd", "RIDAGEYR", "sex")
full_vars_age <- c(base_vars_age, demo_vars, life_vars)
strata_age_lm <- fit_by_strata_drop(analysis, "age_group", "egfr", full_vars_age)
strata_age_glm <- fit_by_strata_drop(analysis, "age_group", "ckd", full_vars_age, family = binomial())

full_vars_smoke <- c(base_vars, demo_vars, life_vars)
strata_smoke_lm <- fit_by_strata_drop(analysis, "smoking", "egfr", full_vars_smoke)
strata_smoke_glm <- fit_by_strata_drop(analysis, "smoking", "ckd", full_vars_smoke, family = binomial())

models <- list(
  fit_lm_min = fit_lm_min,
  fit_lm_demo = fit_lm_demo,
  fit_lm_full = fit_lm_full,
  fit_glm_min = fit_glm_min,
  fit_glm_demo = fit_glm_demo,
  fit_glm_full = fit_glm_full,
  fit_lm_w_min = fit_lm_w_min,
  fit_lm_w_demo = fit_lm_w_demo,
  fit_lm_w_full = fit_lm_w_full,
  fit_glm_w_min = fit_glm_w_min,
  fit_glm_w_demo = fit_glm_w_demo,
  fit_glm_w_full = fit_glm_w_full,
  fit_svy_min = fit_svy_min,
  fit_svy_demo = fit_svy_demo,
  fit_svy_full = fit_svy_full,
  fit_svy_ckd_min = fit_svy_ckd_min,
  fit_svy_ckd_demo = fit_svy_ckd_demo,
  fit_svy_ckd_full = fit_svy_ckd_full,
  fit_svy_w_min = fit_svy_w_min,
  fit_svy_w_demo = fit_svy_w_demo,
  fit_svy_w_full = fit_svy_w_full,
  fit_svy_ckd_w_min = fit_svy_ckd_w_min,
  fit_svy_ckd_w_demo = fit_svy_ckd_w_demo,
  fit_svy_ckd_w_full = fit_svy_ckd_w_full,
  fit_lm_cd_spline = fit_lm_cd_spline,
  fit_glm_cd_spline = fit_glm_cd_spline,
  fit_svy_cd_spline = fit_svy_cd_spline,
  fit_svy_ckd_cd_spline = fit_svy_ckd_cd_spline,
  strata_sex_lm = strata_sex_lm,
  strata_sex_glm = strata_sex_glm,
  strata_age_lm = strata_age_lm,
  strata_age_glm = strata_age_glm,
  strata_smoke_lm = strata_smoke_lm,
  strata_smoke_glm = strata_smoke_glm
)

saveRDS(models, file = "data/models.rds")
message("Saved data/models.rds")

# Tidy 输出（便于汇总与作图）
model_tidy <- list(
  fit_lm_min = broom::tidy(fit_lm_min),
  fit_lm_demo = broom::tidy(fit_lm_demo),
  fit_lm_full = broom::tidy(fit_lm_full),
  fit_glm_min = broom::tidy(fit_glm_min, conf.int = TRUE, exponentiate = TRUE),
  fit_glm_demo = broom::tidy(fit_glm_demo, conf.int = TRUE, exponentiate = TRUE),
  fit_glm_full = broom::tidy(fit_glm_full, conf.int = TRUE, exponentiate = TRUE),
  fit_lm_w_min = broom::tidy(fit_lm_w_min),
  fit_lm_w_demo = broom::tidy(fit_lm_w_demo),
  fit_lm_w_full = broom::tidy(fit_lm_w_full),
  fit_glm_w_min = broom::tidy(fit_glm_w_min, conf.int = TRUE, exponentiate = TRUE),
  fit_glm_w_demo = broom::tidy(fit_glm_w_demo, conf.int = TRUE, exponentiate = TRUE),
  fit_glm_w_full = broom::tidy(fit_glm_w_full, conf.int = TRUE, exponentiate = TRUE),
  fit_svy_min = broom::tidy(fit_svy_min),
  fit_svy_demo = broom::tidy(fit_svy_demo),
  fit_svy_full = broom::tidy(fit_svy_full),
  fit_svy_ckd_min = broom::tidy(fit_svy_ckd_min, conf.int = TRUE, exponentiate = TRUE),
  fit_svy_ckd_demo = broom::tidy(fit_svy_ckd_demo, conf.int = TRUE, exponentiate = TRUE),
  fit_svy_ckd_full = broom::tidy(fit_svy_ckd_full, conf.int = TRUE, exponentiate = TRUE),
  fit_svy_w_min = broom::tidy(fit_svy_w_min),
  fit_svy_w_demo = broom::tidy(fit_svy_w_demo),
  fit_svy_w_full = broom::tidy(fit_svy_w_full),
  fit_svy_ckd_w_min = broom::tidy(fit_svy_ckd_w_min, conf.int = TRUE, exponentiate = TRUE),
  fit_svy_ckd_w_demo = broom::tidy(fit_svy_ckd_w_demo, conf.int = TRUE, exponentiate = TRUE),
  fit_svy_ckd_w_full = broom::tidy(fit_svy_ckd_w_full, conf.int = TRUE, exponentiate = TRUE),
  fit_lm_cd_spline = broom::tidy(fit_lm_cd_spline),
  fit_glm_cd_spline = broom::tidy(fit_glm_cd_spline, conf.int = TRUE, exponentiate = TRUE),
  fit_svy_cd_spline = broom::tidy(fit_svy_cd_spline),
  fit_svy_ckd_cd_spline = broom::tidy(fit_svy_ckd_cd_spline, conf.int = TRUE, exponentiate = TRUE)
)

saveRDS(model_tidy, file = "data/models_tidy.rds")
message("Saved data/models_tidy.rds")

# 样条检验结果保存
spline_tests <- list(
  lm_cd = anova(fit_lm_cd_linear, fit_lm_cd_spline),
  glm_cd = anova(fit_glm_cd_linear, fit_glm_cd_spline, test = "Chisq"),
  svy_cd = survey::regTermTest(fit_svy_cd_spline, ~splines::ns(log_cd, df = 3)),
  svy_ckd_cd = survey::regTermTest(fit_svy_ckd_cd_spline, ~splines::ns(log_cd, df = 3))
)

saveRDS(spline_tests, file = "data/spline_tests.rds")
message("Saved data/spline_tests.rds")

# 分层模型结果保存
strata_tidy <- list(
  sex_lm = lapply(strata_sex_lm, broom::tidy),
  sex_glm = lapply(strata_sex_glm, broom::tidy, conf.int = TRUE, exponentiate = TRUE),
  age_lm = lapply(strata_age_lm, broom::tidy),
  age_glm = lapply(strata_age_glm, broom::tidy, conf.int = TRUE, exponentiate = TRUE),
  smoke_lm = lapply(strata_smoke_lm, broom::tidy),
  smoke_glm = lapply(strata_smoke_glm, broom::tidy, conf.int = TRUE, exponentiate = TRUE)
)

saveRDS(strata_tidy, file = "data/strata_models_tidy.rds")
message("Saved data/strata_models_tidy.rds")

invisible(TRUE)
