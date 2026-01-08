# 包管理：自动设置 CRAN 镜像与 conda 库路径

# Ensure packages install into the active conda environment
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (nzchar(conda_prefix)) {
  conda_lib <- file.path(conda_prefix, "lib", "R", "library")
  Sys.setenv(R_LIBS_USER = conda_lib)
  .libPaths(unique(c(conda_lib, .libPaths())))
}

# 非交互模式下设置默认 CRAN 镜像
if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# 核心依赖
required_pkgs <- c("nhanesA", "tidyverse", "survey", "broom", "here", "scales", "mgcv")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs)
}

lapply(required_pkgs, library, character.only = TRUE)

invisible(TRUE)
