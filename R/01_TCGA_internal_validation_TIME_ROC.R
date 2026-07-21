# ==============================================================================
# 01_TCGA_internal_validation.R
# ==============================================================================
# Title:
#   Development and internal validation of prognostic Cox models in TCGA-OV
#
# Models:
#   M1 — Clinical model:  age + stage
#   M2 — Molecular model: ZNF429
#   M3 — Combined model:  age + stage + ZNF429
#
# Endpoint:
#   Overall survival (OS), with time expressed in years.
#
# Main analyses:
#   - Cox proportional hazards models
#   - AIC comparison
#   - Likelihood ratio test
#   - Bootstrap-corrected C-index
#   - Bootstrap-corrected calibration slope
#   - Time-dependent AUC at 2, 3 and 4 years
#   - Time-dependent Brier score at 2, 3 and 4 years
#   - Calibration curves
#   - Kaplan–Meier curves based on the combined-model risk score
#   - Hazard-ratio table for the combined model
#
# Input file:
#   dados_ZNF429_TCGA_INCA.RData
#
# Required object in the RData file:
#   tcga_ov_final
#
# Required columns:
#   bcr_patient_barcode
#   Stage
#   OS
#   OS.time
#   age_at_initial_pathologic_diagnosis
#   ZNF429_survminer
#
# Reproducibility:
#   Bootstrap analyses use seed 123 and 1,000 resamples.
#
# Important:
#   This script preserves the original model definitions, filters, time points,
#   seeds and risk-group definition.
# ==============================================================================


# ==============================================================================
# 0. LOAD DATA AND PACKAGES
# ==============================================================================

load("data/dados_ZNF429_TCGA_INCA.RData")

library(survival)
library(rms)
library(dplyr)
library(timeROC)
library(ggplot2)
library(patchwork)
library(survminer)
library(riskRegression)


# ==============================================================================
# 1. PREPARE TCGA DATA
# ==============================================================================

cox_data_znf <- tcga_ov_final[, c(
  "bcr_patient_barcode",
  "Stage",
  "OS",
  "OS.time",
  "age_at_initial_pathologic_diagnosis",
  "ZNF429_survminer"
)]

cox_data_znf <- cox_data_znf %>%
  mutate(
    OS = as.numeric(OS),
    OS.time = as.numeric(OS.time),

    age_at_initial_pathologic_diagnosis =
      as.numeric(age_at_initial_pathologic_diagnosis),

    Stage = as.factor(Stage),

    ZNF429 = as.factor(ZNF429_survminer)
  ) %>%
  filter(
    !is.na(age_at_initial_pathologic_diagnosis),
    !is.na(ZNF429)
  )


# ==============================================================================
# 2. RMS SETTINGS
# ==============================================================================

# datadist() stores variable distributions used internally by rms functions.
dd <- datadist(cox_data_znf)
options(datadist = "dd")


# ==============================================================================
# 3. FIT RMS COX MODELS
# ==============================================================================

# x = TRUE and y = TRUE retain the design matrix and response.
# surv = TRUE allows survival-probability calculations.

model_clinical <- cph(
  Surv(OS.time, OS) ~
    age_at_initial_pathologic_diagnosis +
    Stage,
  data = cox_data_znf,
  x = TRUE,
  y = TRUE,
  surv = TRUE
)

model_molecular <- cph(
  Surv(OS.time, OS) ~
    ZNF429,
  data = cox_data_znf,
  x = TRUE,
  y = TRUE,
  surv = TRUE
)

model_combined <- cph(
  Surv(OS.time, OS) ~
    age_at_initial_pathologic_diagnosis +
    Stage +
    ZNF429,
  data = cox_data_znf,
  x = TRUE,
  y = TRUE,
  surv = TRUE
)

print(model_clinical)
print(model_molecular)
print(model_combined)


# ==============================================================================
# 4. FIT SURVIVAL::COXPH MODELS
# ==============================================================================

# These parallel coxph objects are used for AIC, LRT, the forest table and
# riskRegression::Score(). Setting x = TRUE stores the design matrix and does
# not alter coefficients or fitted values.

model_clinical_coxph <- coxph(
  Surv(OS.time, OS) ~
    age_at_initial_pathologic_diagnosis +
    Stage,
  data = cox_data_znf,
  x = TRUE
)

model_molecular_coxph <- coxph(
  Surv(OS.time, OS) ~
    ZNF429,
  data = cox_data_znf,
  x = TRUE
)

model_combined_coxph <- coxph(
  Surv(OS.time, OS) ~
    age_at_initial_pathologic_diagnosis +
    Stage +
    ZNF429,
  data = cox_data_znf,
  x = TRUE
)

print(summary(model_clinical_coxph))
print(summary(model_molecular_coxph))
print(summary(model_combined_coxph))


# ==============================================================================
# 5. MODEL COMPARISON: AIC
# ==============================================================================

# Lower AIC indicates a better balance between model fit and complexity.

AIC_results <- AIC(
  model_clinical_coxph,
  model_molecular_coxph,
  model_combined_coxph
)

print(AIC_results)


# ==============================================================================
# 6. LIKELIHOOD RATIO TEST
# ==============================================================================

# The clinical model is nested within the combined model.
# This test evaluates whether adding ZNF429 improves model fit.

lrt_results <- anova(
  model_clinical_coxph,
  model_combined_coxph,
  test = "LRT"
)

print(lrt_results)


# ==============================================================================
# 7. INTERNAL VALIDATION BY BOOTSTRAP
# ==============================================================================

set.seed(123)

val_clinical <- validate(
  model_clinical,
  method = "boot",
  B = 1000,
  dxy = TRUE
)

set.seed(123)

val_molecular <- validate(
  model_molecular,
  method = "boot",
  B = 1000,
  dxy = TRUE
)

set.seed(123)

val_combined <- validate(
  model_combined,
  method = "boot",
  B = 1000,
  dxy = TRUE
)

print(val_clinical)
print(val_molecular)
print(val_combined)


# ==============================================================================
# 8. BOOTSTRAP-CORRECTED C-INDEX
# ==============================================================================

# rms reports Somers' Dxy.
# The corresponding C-index is calculated as: C = (Dxy + 1) / 2.

cindex_table <- data.frame(
  Model = c(
    "Clinical",
    "Molecular",
    "Combined"
  ),

  Dxy_corrected = c(
    val_clinical["Dxy", "index.corrected"],
    val_molecular["Dxy", "index.corrected"],
    val_combined["Dxy", "index.corrected"]
  )
)

cindex_table$C_index_corrected <- (
  cindex_table$Dxy_corrected + 1
) / 2

cindex_table <- cindex_table %>%
  mutate(
    Dxy_corrected = round(Dxy_corrected, 3),
    C_index_corrected = round(C_index_corrected, 3)
  )

print(cindex_table)


# ==============================================================================
# 9. BOOTSTRAP-CORRECTED CALIBRATION SLOPE
# ==============================================================================

# Interpretation:
#   slope = 1: ideal global calibration
#   slope < 1: predictions may be too extreme
#   slope > 1: predictions may be insufficiently extreme

calibration_slope <- data.frame(
  Model = c(
    "Clinical",
    "ZNF429",
    "Combined"
  ),

  Slope_corrected = c(
    val_clinical["Slope", "index.corrected"],
    val_molecular["Slope", "index.corrected"],
    val_combined["Slope", "index.corrected"]
  )
)

calibration_slope$Slope_corrected <- round(
  calibration_slope$Slope_corrected,
  3
)

print(calibration_slope)


# ==============================================================================
# 10. LINEAR PREDICTORS / RISK SCORES
# ==============================================================================

# No newdata argument is used because predictions correspond to the same
# observations used to fit each rms::cph model.

cox_data_znf$risk_clinical <- predict(
  model_clinical,
  type = "lp"
)

cox_data_znf$risk_molecular <- predict(
  model_molecular,
  type = "lp"
)

cox_data_znf$risk_combined <- predict(
  model_combined,
  type = "lp"
)


# ==============================================================================
# 11. TIME-DEPENDENT ROC AND AUC
# ==============================================================================

roc_clinical <- timeROC(
  T = cox_data_znf$OS.time,
  delta = cox_data_znf$OS,
  marker = cox_data_znf$risk_clinical,
  cause = 1,
  times = c(2, 3, 4),
  iid = TRUE
)

roc_molecular <- timeROC(
  T = cox_data_znf$OS.time,
  delta = cox_data_znf$OS,
  marker = cox_data_znf$risk_molecular,
  cause = 1,
  times = c(2, 3, 4),
  iid = TRUE
)

roc_combined <- timeROC(
  T = cox_data_znf$OS.time,
  delta = cox_data_znf$OS,
  marker = cox_data_znf$risk_combined,
  cause = 1,
  times = c(2, 3, 4),
  iid = TRUE
)

auc_table <- data.frame(
  Time = c(
    "2 years",
    "3 years",
    "4 years"
  ),

  Clinical = roc_clinical$AUC,
  ZNF429 = roc_molecular$AUC,
  Combined = roc_combined$AUC
) %>%
  mutate(
    Clinical = round(Clinical, 3),
    ZNF429 = round(ZNF429, 3),
    Combined = round(Combined, 3)
  )

print(auc_table)


# ==============================================================================
# 12. COMPLETE-CASE DATASET USED ONLY BY SCORE()
# ==============================================================================

# riskRegression::Score() does not accept missing values in the response.
# A separate complete-case object is therefore created only for Brier-score and
# calibration calculations. The original analysis object is not overwritten.

cox_data_znf_brier <- cox_data_znf %>%
  filter(
    complete.cases(
      OS,
      OS.time,
      age_at_initial_pathologic_diagnosis,
      Stage,
      ZNF429
    )
  )

cat(
  "\nNumber of observations in cox_data_znf:",
  nrow(cox_data_znf),
  "\n"
)

cat(
  "Number of observations used for Brier score:",
  nrow(cox_data_znf_brier),
  "\n"
)


# ==============================================================================
# 13. TIME-DEPENDENT BRIER SCORE
# ==============================================================================

# Lower Brier scores indicate lower overall prediction error.
# IPCW is used to account for censoring.

score_calibration <- riskRegression::Score(
  object = list(
    Clinical = model_clinical_coxph,
    ZNF429 = model_molecular_coxph,
    Combined = model_combined_coxph
  ),

  formula = Surv(OS.time, OS) ~ 1,

  data = cox_data_znf_brier,

  metrics = "Brier",

  plots = "calibration",

  times = c(
    2,
    3,
    4
  ),

  cens.method = "ipcw",

  cens.model = "km",

  null.model = TRUE,

  conf.int = TRUE
)

brier_table_complete <- score_calibration$Brier$score %>%
  as.data.frame() %>%
  filter(
    model %in% c(
      "Clinical",
      "ZNF429",
      "Combined"
    )
  ) %>%
  arrange(
    model,
    times
  )

print(brier_table_complete)

brier_table <- brier_table_complete %>%
  select(
    model,
    times,
    Brier
  ) %>%
  mutate(
    Brier = round(Brier, 3)
  )

print(brier_table)


# ==============================================================================
# 14. NAMED BRIER-SCORE VECTORS
# ==============================================================================

# These named vectors simplify metric extraction for the ROC annotations and
# are also used by 03_TCGA_INCA_final_panel.R.

brier_clinical <- setNames(
  brier_table$Brier[
    brier_table$model == "Clinical"
  ],
  brier_table$times[
    brier_table$model == "Clinical"
  ]
)

brier_znf429 <- setNames(
  brier_table$Brier[
    brier_table$model == "ZNF429"
  ],
  brier_table$times[
    brier_table$model == "ZNF429"
  ]
)

brier_combined <- setNames(
  brier_table$Brier[
    brier_table$model == "Combined"
  ],
  brier_table$times[
    brier_table$model == "Combined"
  ]
)

print(brier_clinical)
print(brier_znf429)
print(brier_combined)


# ==============================================================================
# 15. TIME-DEPENDENT ROC FIGURES
# ==============================================================================

extract_roc_df <- function(
    roc_obj,
    timepoint,
    model_name
) {

  idx <- which(
    roc_obj$times == timepoint
  )

  data.frame(
    FPR = roc_obj$FP[, idx],
    TPR = roc_obj$TP[, idx],
    Model = model_name
  )
}

make_roc_plot <- function(
    timepoint,
    auc_index,
    title_text
) {

  roc_df <- bind_rows(
    extract_roc_df(
      roc_clinical,
      timepoint,
      "Clinical"
    ),

    extract_roc_df(
      roc_molecular,
      timepoint,
      "ZNF429"
    ),

    extract_roc_df(
      roc_combined,
      timepoint,
      "Combined"
    )
  )

  metrics_label <- paste0(
    "Clinical: AUC = ",
    round(roc_clinical$AUC[auc_index], 3),
    " | Brier = ",
    round(brier_clinical[as.character(timepoint)], 3),

    "\nZNF429: AUC = ",
    round(roc_molecular$AUC[auc_index], 3),
    " | Brier = ",
    round(brier_znf429[as.character(timepoint)], 3),

    "\nCombined: AUC = ",
    round(roc_combined$AUC[auc_index], 3),
    " | Brier = ",
    round(brier_combined[as.character(timepoint)], 3)
  )

  ggplot(
    roc_df,
    aes(
      x = FPR,
      y = TPR,
      color = Model
    )
  ) +
    geom_line(
      linewidth = 1.3
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray60"
    ) +
    annotate(
      "text",
      x = 0.32,
      y = 0.22,
      hjust = 0,
      size = 3.5,
      label = metrics_label
    ) +
    labs(
      title = title_text,
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = NULL
    ) +
    coord_equal() +
    theme_classic(
      base_size = 14
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "bottom",
      legend.text = element_text(
        size = 11
      )
    )
}

p2 <- make_roc_plot(
  timepoint = 2,
  auc_index = 1,
  title_text = "2-Year ROC Curve"
)

p3 <- make_roc_plot(
  timepoint = 3,
  auc_index = 2,
  title_text = "3-Year ROC Curve"
)

p4 <- make_roc_plot(
  timepoint = 4,
  auc_index = 3,
  title_text = "4-Year ROC Curve"
)

p_roc_all <-
  p2 |
  p3 |
  p4

print(p_roc_all)



# ==============================================================================
# 17. KAPLAN–MEIER CURVE USING COMBINED-MODEL RISK
# ==============================================================================

# The original analysis uses the cohort median linear predictor to define
# high- and low-risk groups.

cox_data_znf$risk_group_combined <- ifelse(
  cox_data_znf$risk_combined >= median(
    cox_data_znf$risk_combined,
    na.rm = TRUE
  ),
  "High risk",
  "Low risk"
)

cox_data_znf$risk_group_combined <- factor(
  cox_data_znf$risk_group_combined,
  levels = c(
    "Low risk",
    "High risk"
  )
)

km_fit <- survfit(
  Surv(OS.time, OS) ~
    risk_group_combined,
  data = cox_data_znf
)

p_km <- ggsurvplot(
  km_fit,
  data = cox_data_znf,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c(
    "#3B82F6",
    "#EF4444"
  ),
  legend.title = "",
  legend.labs = c(
    "Low risk",
    "High risk"
  ),
  xlab = "Time (years)",
  ylab = "Overall Survival Probability",
  ggtheme = theme_classic(
    base_size = 14
  ),
  risk.table.height = 0.25,
  risk.table.y.text = FALSE
)

print(p_km)

km_logrank <- survdiff(
  Surv(OS.time, OS) ~
    risk_group_combined,
  data = cox_data_znf
)

print(km_logrank)


# ==============================================================================
# 18. COMBINED-MODEL HAZARD-RATIO TABLE
# ==============================================================================

cox_summary <- summary(
  model_combined_coxph
)

forest_table <- data.frame(
  Variable = rownames(
    cox_summary$coefficients
  ),

  HR = cox_summary$conf.int[
    ,
    "exp(coef)"
  ],

  Lower95 = cox_summary$conf.int[
    ,
    "lower .95"
  ],

  Upper95 = cox_summary$conf.int[
    ,
    "upper .95"
  ],

  Pvalue = cox_summary$coefficients[
    ,
    "Pr(>|z|)"
  ]
)

forest_table <- forest_table %>%
  mutate(
    HR = round(HR, 3),
    Lower95 = round(Lower95, 3),
    Upper95 = round(Upper95, 3),
    Pvalue = signif(Pvalue, 3)
  )

print(forest_table)


# ==============================================================================
# 19. FINAL CONSOLE SUMMARY
# ==============================================================================

cat(
  "\n============================================================\n",
  "AIC RESULTS\n",
  "============================================================\n",
  sep = ""
)
print(AIC_results)

cat(
  "\n============================================================\n",
  "LIKELIHOOD RATIO TEST\n",
  "============================================================\n",
  sep = ""
)
print(lrt_results)

cat(
  "\n============================================================\n",
  "BOOTSTRAP-CORRECTED C-INDEX\n",
  "============================================================\n",
  sep = ""
)
print(cindex_table)

cat(
  "\n============================================================\n",
  "BOOTSTRAP-CORRECTED CALIBRATION SLOPE\n",
  "============================================================\n",
  sep = ""
)
print(calibration_slope)

cat(
  "\n============================================================\n",
  "TIME-DEPENDENT AUC\n",
  "============================================================\n",
  sep = ""
)
print(auc_table)

cat(
  "\n============================================================\n",
  "TIME-DEPENDENT BRIER SCORE\n",
  "============================================================\n",
  sep = ""
)
print(brier_table)

cat(
  "\n============================================================\n",
  "COMBINED MODEL HAZARD-RATIO TABLE\n",
  "============================================================\n",
  sep = ""
)
print(forest_table)


# ==============================================================================
# 20. OPTIONAL EXPORTS
# ==============================================================================

# Uncomment the commands below to export results.
#
# write.csv(
#   cindex_table,
#   "results/TCGA_bootstrap_corrected_cindex.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   calibration_slope,
#   "results/TCGA_bootstrap_corrected_calibration_slope.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   auc_table,
#   "results/TCGA_time_dependent_AUC.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   brier_table,
#   "results/TCGA_time_dependent_Brier.csv",
#   row.names = FALSE
# )
#
# write.csv(
#   forest_table,
#   "results/TCGA_combined_model_HR.csv",
#   row.names = FALSE
# )
#
# ggsave(
#   filename = "figures/TCGA_time_dependent_ROC.png",
#   plot = p_roc_all,
#   width = 15,
#   height = 5,
#   dpi = 300
# )
