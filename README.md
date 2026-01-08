# nhanes_metals_ckd
My R language assignment

## 运行说明

项目链接：GitHub：

1) **安装/加载包**  
   `R/00_packages.R` 会自动设置 CRAN 镜像，并安装缺失包。

2) **下载与合并数据**  
   `R/01_download_merge.R` 拉取 DEMO/PBCD/BMX/SMQ 及生化肌酐表，按 `SEQN` 合并。

3) **清洗与变量构造**  
   `R/02_clean_construct.R`  
   - 生成 `sex`、`race`、`education`、`bmi`、`smoking`  
   - 计算 `log_pb`、`log_cd`  
   - 生成分位数（Q5）与 winsorize 版本

4) **eGFR 与 CKD**  
   `R/03_egfr.R` 使用 CKD-EPI 2021 公式计算 eGFR，并生成 CKD（<60）。

5) **模型**  
   `R/04_models.R`  
   - 线性/Logistic（非加权）  
   - survey 加权（NHANES 复杂抽样）  
   - winsorize 敏感性分析  
   - 分层模型（性别/年龄/吸烟）  
   - 镉的非线性样条检验

6) **可视化**  
   `R/05_visualize.R`  
   - 密度热力图、散点趋势、分位数趋势  
   - CKD OR 森林图  
   - 镉 spline 曲线
