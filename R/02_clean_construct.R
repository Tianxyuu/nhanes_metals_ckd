# 数据清洗与变量构造（暴露、协变量、分位数）

source("R/00_packages.R")

raw_path <- "data/merged_raw.rds"
if (!file.exists(raw_path)) {
  stop("Missing data/merged_raw.rds. Run R/01_download_merge.R first.")
}

raw <- readRDS(raw_path)

# 性别编码兼容：既支持数值(1/2)也支持字符(Male/Female)
sex_from_riagendr <- function(x) {
  if (is.numeric(x)) {
    return(factor(x, levels = c(1, 2), labels = c("Male", "Female")))
  }
  x_chr <- as.character(x)
  x_chr[x_chr %in% c("1", "2")] <- ifelse(x_chr[x_chr %in% c("1", "2")] == "1", "Male", "Female")
  factor(x_chr, levels = c("Male", "Female"))
}

# 吸烟状态：Never/Former/Current
smoking_status <- function(smq020, smq040) {
  # SMQ020: smoked 100 cigarettes (1=Yes, 2=No)
  # SMQ040: do you now smoke (1=Every day, 2=Some days, 3=Not at all)
  status <- rep(NA_character_, length(smq020))
  status[smq020 == 2] <- "Never"
  status[smq020 == 1 & smq040 %in% c(1, 2)] <- "Current"
  status[smq020 == 1 & smq040 == 3] <- "Former"
  factor(status, levels = c("Never", "Former", "Current"))
}

# 极端值截尾（1%/99%）
winsorize <- function(x, probs = c(0.01, 0.99)) {
  qs <- quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  x <- pmin(pmax(x, qs[1]), qs[2])
  x
}

# 去除关键变量缺失与未成年人
clean <- raw %>%
  filter(!is.na(LBXBPB), !is.na(LBXBCD), !is.na(LBXSCR), RIDAGEYR >= 20) %>%
  mutate(
    sex = sex_from_riagendr(RIAGENDR),
    race = factor(RIDRETH3),
    education = factor(DMDEDUC2),
    log_pb = log(LBXBPB),
    log_cd = log(LBXBCD),
    bmi = BMXBMI,
    smoking = smoking_status(SMQ020, SMQ040),
    log_pb_w = winsorize(log_pb),
    log_cd_w = winsorize(log_cd),
    pb_q5 = ntile(log_pb, 5),
    cd_q5 = ntile(log_cd, 5),
    # 年龄分层用于亚组分析
    age_group = cut(RIDAGEYR, breaks = c(20, 40, 60, Inf),
                    labels = c("20-39", "40-59", "60+"), right = FALSE)
  )

saveRDS(clean, file = "data/clean_base.rds")
message("Saved data/clean_base.rds")

invisible(TRUE)
