# ==============================================================================
# 02_INCA_external_validation.R
# ==============================================================================
# Title:
#   External validation of the TCGA-derived prognostic model in the INCA cohort
#
# Development model:
#   age + stage + ZNF429
#
# External validation outcomes:
#   - TCGA-derived linear predictor applied to INCA
#   - Time-dependent AUC at 2, 3 and 4 years
#   - Kaplan–Meier curves based on the INCA median predicted risk
#
# Input file:
#   dados_ZNF429_TCGA_INCA.RData
#
# Required objects in the RData file:
#   tcga_ov_final
#   znf_inca
#
# Required TCGA columns:
#   bcr_patient_barcode
#   Stage
#   OS
#   OS.time
#   age_at_initial_pathologic_diagnosis
#   ZNF429_survminer
#
# Required INCA columns:
#   Sample
#   ZNF429
#   OS_time5
#   OS_status5
#   Stage
#   Idade_do_diagnostico
#
# Important:
#   The TCGA model is fitted without using INCA outcomes.
#   INCA observations are used only for external performance evaluation.
# ==============================================================================


# ==============================================================================
# 0. LOAD DATA AND PACKAGES
# ==============================================================================

load("dados_ZNF429_TCGA_INCA.RData")

library(survival)
library(dplyr)
library(timeROC)
library(ggplot2)
library(patchwork)
library(survminer)


# ==============================================================================
# 1. RECREATE THE TCGA ANALYSIS DATA
# ==============================================================================

# This section allows the external-validation script to be run independently
# from 01_TCGA_internal_validation.R.

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
    # These filters remain commented to preserve the original analysis:
    # !is.na(OS),
    # !is.na(OS.time),
    # !is.na(Stage),

    !is.na(age_at_initial_pathologic_diagnosis),
    !is.na(ZNF429)
  )


# ==============================================================================
# 2. HARMONIZE TCGA FACTOR LEVELS
# ==============================================================================

# The validation data must use the same factor levels as the development data.

cox_data_znf$Stage <- factor(
  as.character(cox_data_znf$Stage),
  levels = c(
    "1",
    "2",
    "3",
    "4"
  )
)

cox_data_znf$ZNF429 <- factor(
  as.character(cox_data_znf$ZNF429),
  levels = levels(cox_data_znf$ZNF429)
)


# ==============================================================================
# 3. FIT THE TCGA COMBINED MODEL
# ==============================================================================

# This is the model transported to the INCA cohort.
# x = TRUE stores the model matrix but does not change the coefficients.

model_combined_tcga <- coxph(
  Surv(OS.time, OS) ~
    age_at_initial_pathologic_diagnosis +
    Stage +
    ZNF429,
  data = cox_data_znf,
  x = TRUE
)

print(summary(model_combined_tcga))


# ==============================================================================
# 4. PREPARE THE INCA EXTERNAL-VALIDATION DATA
# ==============================================================================

cox_data_znf_inca <- znf_inca[, c(
  "Sample",
  "ZNF429",
  "OS_time5",
  "OS_status5",
  "Stage",
  "Idade_do_diagnostico"
)]

cox_data_znf_inca <- cox_data_znf_inca %>%
  rename(
    OS.time = OS_time5,
    OS = OS_status5,

    age_at_initial_pathologic_diagnosis =
      Idade_do_diagnostico
  ) %>%
  mutate(
    OS = as.numeric(OS),

    # INCA survival time is converted from months to years.
    OS.time = as.numeric(OS.time) / 12,

    age_at_initial_pathologic_diagnosis =
      as.numeric(age_at_initial_pathologic_diagnosis),

    # Keep the original stage variable for quality-control tables.
    Stage_original = Stage,

    # Harmonize INCA stage categories with the TCGA model coding.
    Stage = case_when(
      Stage == "I/II" ~ "2",
      Stage == "III"  ~ "3",
      Stage == "IV"   ~ "4",
      TRUE ~ NA_character_
    ),

    Stage = factor(
      Stage,
      levels = c(
        "1",
        "2",
        "3",
        "4"
      )
    ),

    # Apply the same ZNF429 factor levels used in TCGA.
    ZNF429 = factor(
      as.character(ZNF429),
      levels = levels(cox_data_znf$ZNF429)
    )
  ) %>%
  filter(
    # These filters remain commented to preserve the original analysis:
    # !is.na(OS),
    # !is.na(OS.time),
    # !is.na(Stage),

    !is.na(age_at_initial_pathologic_diagnosis),
    !is.na(ZNF429)
  )

# Quality-control summaries.
print(
  table(
    cox_data_znf_inca$Stage_original,
    cox_data_znf_inca$Stage
  )
)

print(
  table(
    cox_data_znf_inca$ZNF429
  )
)


# ==============================================================================
# 5. APPLY THE TCGA MODEL TO INCA
# ==============================================================================

# type = "lp" returns the linear predictor from the fixed TCGA coefficients.

cox_data_znf_inca$risk_score <- predict(
  model_combined_tcga,
  newdata = cox_data_znf_inca,
  type = "lp"
)

print(
  summary(
    cox_data_znf_inca$risk_score
  )
)


# ==============================================================================
# 6. TIME-DEPENDENT ROC AND AUC IN INCA
# ==============================================================================

roc_inca <- timeROC(
  T = cox_data_znf_inca$OS.time,
  delta = cox_data_znf_inca$OS,
  marker = cox_data_znf_inca$risk_score,
  cause = 1,
  times = c(
    2,
    3,
    4
  ),
  iid = TRUE
)

auc_inca <- data.frame(
  Time = c(
    "2 years",
    "3 years",
    "4 years"
  ),
  AUC = roc_inca$AUC
)

print(auc_inca)


# ==============================================================================
# 7. INCA ROC FIGURES
# ==============================================================================

extract_roc <- function(
    roc_obj,
    timepoint
) {

  idx <- which(
    roc_obj$times == timepoint
  )

  data.frame(
    FPR = roc_obj$FP[, idx],
    TPR = roc_obj$TP[, idx]
  )
}

make_roc_plot <- function(
    timepoint,
    auc_index,
    title_text
) {

  roc_df <- extract_roc(
    roc_inca,
    timepoint
  )

  ggplot(
    roc_df,
    aes(
      x = FPR,
      y = TPR
    )
  ) +
    geom_line(
      linewidth = 1.4,
      color = "#DC2626"
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray60"
    ) +
    annotate(
      "text",
      x = 0.62,
      y = 0.15,
      size = 5,
      fontface = "bold",
      label = paste0(
        "AUC = ",
        round(
          roc_inca$AUC[auc_index],
          3
        )
      )
    ) +
    labs(
      title = title_text,
      x = "False positive rate",
      y = "True positive rate"
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
      axis.title = element_text(
        face = "bold"
      ),
      axis.text = element_text(
        color = "black"
      )
    )
}

p2_inca <- make_roc_plot(
  timepoint = 2,
  auc_index = 1,
  title_text = "INCA validation: 2-Year ROC"
)

p3_inca <- make_roc_plot(
  timepoint = 3,
  auc_index = 2,
  title_text = "INCA validation: 3-Year ROC"
)

p4_inca <- make_roc_plot(
  timepoint = 4,
  auc_index = 3,
  title_text = "INCA validation: 4-Year ROC"
)

inca_roc_panel <-
  p2_inca |
  p3_inca |
  p4_inca

print(inca_roc_panel)


# ==============================================================================
# 8. KAPLAN–MEIER CURVES IN THE INCA COHORT
# ==============================================================================

# The original analysis defines external high- and low-risk groups using the
# median INCA risk score predicted by the fixed TCGA model.

cox_data_znf_inca$risk_group <- ifelse(
  cox_data_znf_inca$risk_score >= median(
    cox_data_znf_inca$risk_score,
    na.rm = TRUE
  ),
  "High risk",
  "Low risk"
)

cox_data_znf_inca$risk_group <- factor(
  cox_data_znf_inca$risk_group,
  levels = c(
    "Low risk",
    "High risk"
  )
)

km_fit_inca <- survfit(
  Surv(OS.time, OS) ~
    risk_group,
  data = cox_data_znf_inca
)

p_km_inca <- ggsurvplot(
  km_fit_inca,
  data = cox_data_znf_inca,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c(
    "#2563EB",
    "#DC2626"
  ),
  legend.title = "",
  legend.labs = c(
    "Low risk",
    "High risk"
  ),
  xlab = "Time (years)",
  ylab = "Overall survival probability",
  title = "INCA validation using TCGA model",
  ggtheme = theme_classic(
    base_size = 14
  ),
  risk.table.height = 0.25,
  risk.table.y.text = FALSE
)

print(p_km_inca)


# ==============================================================================
# 9. OPTIONAL EXPORTS
# ==============================================================================

# Uncomment the commands below to export results.
#
# write.csv(
#   auc_inca,
#   "results/INCA_time_dependent_AUC.csv",
#   row.names = FALSE
# )
#
# ggsave(
#   filename = "figures/INCA_time_dependent_ROC.png",
#   plot = inca_roc_panel,
#   width = 15,
#   height = 5,
#   dpi = 300
# )
#
# ggsave(
#   filename = "figures/INCA_Kaplan_Meier.png",
#   plot = p_km_inca$plot,
#   width = 7,
#   height = 6,
#   dpi = 300
# )
