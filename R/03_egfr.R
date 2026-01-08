# eGFR 计算与 CKD 构造（CKD-EPI 2021）

source("R/00_packages.R")

# CKD-EPI 2021（肌酐）公式
calc_egfr_2021 <- function(scr_mg_dl, age, sex) {
  k <- ifelse(sex == "Female", 0.7, 0.9)
  a <- ifelse(sex == "Female", -0.241, -0.302)
  scr_k <- scr_mg_dl / k
  egfr <- 142 * (pmin(scr_k, 1) ^ a) * (pmax(scr_k, 1) ^ -1.200) * (0.9938 ^ age)
  egfr <- ifelse(sex == "Female", egfr * 1.012, egfr)
  egfr
}

base_path <- "data/clean_base.rds"
if (!file.exists(base_path)) {
  stop("Missing data/clean_base.rds. Run R/02_clean_construct.R first.")
}

# 读取清洗后的数据
base <- readRDS(base_path)

# 计算 eGFR 与 CKD 指标
analysis <- base %>%
  mutate(
    egfr = calc_egfr_2021(LBXSCR, RIDAGEYR, sex),
    ckd = if_else(!is.na(egfr) & egfr < 60, 1L, 0L)
  )

saveRDS(analysis, file = "data/analysis.rds")
message("Saved data/analysis.rds")

invisible(TRUE)
