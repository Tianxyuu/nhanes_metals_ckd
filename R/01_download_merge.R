# 下载并合并 NHANES 2015-2016 + 2017-2018 表格

source("R/00_packages.R")

# 缓存目录（可减少重复下载）
options(nhanesA.cache = "data")

get_creatinine_table <- function(suffix) {
  # Try standard biochemistry file first
  biochem <- nhanesA::nhanes(paste0("BIOPRO_", suffix))
  if ("LBXSCR" %in% names(biochem)) {
    return(biochem)
  }

  # Fallback: search for LBXSCR in the specified cycle
  yr <- if (suffix == "I") "2015-2016" else "2017-2018"
  srch <- nhanesA::nhanesSearch("LBXSCR", yr)
  if (nrow(srch) == 0) {
    srch <- nhanesA::nhanesSearch("creatinine", yr)
  }

  if (nrow(srch) == 0) {
  stop("Could not find a table with serum creatinine (LBXSCR) for ", yr, ".")
  }

  for (file in unique(srch$Data_File)) {
    tbl <- nhanesA::nhanes(file)
    if ("LBXSCR" %in% names(tbl)) {
      return(tbl)
    }
  }

  stop("Searched candidate tables but none contained LBXSCR for ", yr, ".")
}

message("Downloading tables...")

load_cycle <- function(suffix, cycle_label) {
  # DEMO: 人口学；PBCD: 血铅/血镉；BMX: 体测(BMI)；SMQ: 吸烟
  demo <- nhanesA::nhanes(paste0("DEMO_", suffix))
  pbcd <- nhanesA::nhanes(paste0("PBCD_", suffix))
  bmx <- nhanesA::nhanes(paste0("BMX_", suffix))
  smq <- nhanesA::nhanes(paste0("SMQ_", suffix))
  biochem <- get_creatinine_table(suffix)

  dat <- demo %>%
    select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH3, DMDEDUC2, INDFMPIR,
           SDMVPSU, SDMVSTRA, WTMEC2YR) %>%
    left_join(pbcd %>% select(SEQN, LBXBPB, LBXBCD), by = "SEQN") %>%
    left_join(biochem %>% select(SEQN, LBXSCR), by = "SEQN") %>%
    left_join(bmx %>% select(SEQN, BMXBMI), by = "SEQN") %>%
    left_join(smq %>% select(SEQN, SMQ020, SMQ040), by = "SEQN") %>%
    mutate(cycle = cycle_label)

  dat
}

message("Merging tables...")

dat_i <- load_cycle("I", "2015-2016")
dat_j <- load_cycle("J", "2017-2018")
dat <- bind_rows(dat_i, dat_j)

# 4年合并权重（2个周期）
dat <- dat %>% mutate(WTMEC4YR = WTMEC2YR / 2)

saveRDS(dat, file = "data/merged_raw.rds")
message("Saved data/merged_raw.rds")

invisible(TRUE)
