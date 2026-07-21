# ==============================================================================
# ZNF429 MULTIVARIATE COX ANALYSIS
# TCGA-OV, INCA ADJUVANT AND INCA NEOADJUVANT
# ==============================================================================
#
# Models:
#   TCGA-OV:
#     OS ~ ZNF429 group + stage group (I/II vs III/IV) + age
#
#   INCA adjuvant:
#     OS ~ ZNF429 group + stage group (I/II vs III/IV) + age
#
#   INCA neoadjuvant:
#     OS ~ ZNF429 group + stage (III vs IV) + age
#
# Input:
#   dados_ZNF429_cox.RData
#
#
# Important:
#   - TCGA OS.time is assumed to already be expressed in years.
#   - INCA OS_time5 is converted from months to years.
#   - "low" is the reference category for ZNF429.
#   - I/II is the reference category for grouped stage.
#   - Stage III is the reference category in the neoadjuvant cohort.
# ==============================================================================


# ==============================================================================
# 0. LOAD DATA AND PACKAGES
# ==============================================================================

load("data/dados_ZNF429_cox.RData")

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(patchwork)


# ==============================================================================
# 1. TCGA-OV DATA PREPARATION
# ==============================================================================

cox_data_tcga <- tcga_ov_final[, c(
  "bcr_patient_barcode",
  "Stage",
  "OS",
  "OS.time",
  "age_at_initial_pathologic_diagnosis",
  "ZNF429_survminer"
)]

cox_data_tcga <- cox_data_tcga %>%
  mutate(
    Sample = substr(
      as.character(bcr_patient_barcode),
      1,
      12
    ),
    
    OS = as.numeric(
      as.character(OS)
    ),
    
    OS.time = as.numeric(
      as.character(OS.time)
    ),
    
    Age = as.numeric(
      as.character(
        age_at_initial_pathologic_diagnosis
      )
    ),
    
    Stage_numeric = suppressWarnings(
      as.numeric(
        as.character(Stage)
      )
    ),
    
    Stage_group = case_when(
      Stage_numeric %in% c(1, 2) ~ "I/II",
      Stage_numeric %in% c(3, 4) ~ "III/IV",
      TRUE ~ NA_character_
    ),
    
    Stage_group = factor(
      Stage_group,
      levels = c(
        "I/II",
        "III/IV"
      )
    ),
    
    ZNF429_group = factor(
      as.character(
        ZNF429_survminer
      ),
      levels = c(
        "low",
        "high"
      )
    )
  ) %>%
  filter(
    complete.cases(
      OS,
      OS.time,
      Age,
      Stage_group,
      ZNF429_group
    )
  )

cat(
  "\n============================================================\n",
  "TCGA-OV DATA SUMMARY\n",
  "============================================================\n",
  sep = ""
)

cat(
  "Samples:",
  nrow(cox_data_tcga),
  "\n"
)

cat(
  "Events:",
  sum(cox_data_tcga$OS == 1),
  "\n"
)

print(
  table(
    cox_data_tcga$ZNF429_group,
    useNA = "ifany"
  )
)

print(
  table(
    cox_data_tcga$Stage_group,
    useNA = "ifany"
  )
)


# ==============================================================================
# 2. TCGA-OV MULTIVARIATE COX MODEL
# ==============================================================================

fit_cox_tcga <- coxph(
  Surv(OS.time, OS) ~
    ZNF429_group +
    Stage_group +
    Age,
  data = cox_data_tcga,
  na.action = na.omit,
  x = TRUE,
  y = TRUE
)

print(
  summary(
    fit_cox_tcga
  )
)

tcga_znf <- ggforest(
  fit_cox_tcga,
  data = model.frame(
    fit_cox_tcga
  ),
  main = "Multivariate Cox analysis — TCGA-OV",
  fontsize = 1
)

print(tcga_znf)


# ==============================================================================
# 3. INCA ADJUVANT DATA PREPARATION
# ==============================================================================

cox_data_adj <- PCR_INCA_ADj[, c(
  "Sample",
  "ZNF429",
  "OS_time5",
  "OS_status5",
  "Stage",
  "Idade_do_diagnostico"
)]

cox_data_adj <- cox_data_adj %>%
  transmute(
    Sample = as.character(Sample),
    
    OS = as.numeric(
      as.character(OS_status5)
    ),
    
    # Convert months to years.
    OS.time = as.numeric(
      as.character(OS_time5)
    ) / 12,
    
    Age = as.numeric(
      as.character(
        Idade_do_diagnostico
      )
    ),
    
    Stage_original = as.character(
      Stage
    ),
    
    Stage_group = case_when(
      Stage_original == "I/II" ~ "I/II",
      Stage_original %in% c(
        "III",
        "IV"
      ) ~ "III/IV",
      TRUE ~ NA_character_
    ),
    
    Stage_group = factor(
      Stage_group,
      levels = c(
        "I/II",
        "III/IV"
      )
    ),
    
    ZNF429_group = factor(
      as.character(ZNF429),
      levels = c(
        "low",
        "high"
      )
    )
  ) %>%
  filter(
    complete.cases(
      OS,
      OS.time,
      Age,
      Stage_group,
      ZNF429_group
    )
  )

cat(
  "\n============================================================\n",
  "INCA ADJUVANT DATA SUMMARY\n",
  "============================================================\n",
  sep = ""
)

cat(
  "Samples:",
  nrow(cox_data_adj),
  "\n"
)

cat(
  "Events:",
  sum(cox_data_adj$OS == 1),
  "\n"
)

print(
  table(
    cox_data_adj$ZNF429_group,
    useNA = "ifany"
  )
)

print(
  table(
    cox_data_adj$Stage_group,
    useNA = "ifany"
  )
)

print(
  table(
    cox_data_adj$OS,
    useNA = "ifany"
  )
)


# ==============================================================================
# 4. INCA ADJUVANT MULTIVARIATE COX MODEL
# ==============================================================================

fit_cox_adj <- coxph(
  Surv(OS.time, OS) ~
    ZNF429_group +
    Stage_group +
    Age,
  data = cox_data_adj,
  na.action = na.omit,
  x = TRUE,
  y = TRUE
)

print(
  summary(
    fit_cox_adj
  )
)

adj_znf <- ggforest(
  fit_cox_adj,
  data = model.frame(
    fit_cox_adj
  ),
  main = "Multivariate Cox analysis — INCA (Adjuvant)",
  fontsize = 1.2
)

print(adj_znf)


# ==============================================================================
# 5. INCA NEOADJUVANT DATA PREPARATION
# ==============================================================================

cox_data_neo <- PCR_INCA_NEO[, c(
  "Sample",
  "ZNF429",
  "OS_time5",
  "OS_status5",
  "Stage",
  "Idade_do_diagnostico"
)]

cox_data_neo <- cox_data_neo %>%
  transmute(
    Sample = as.character(Sample),
    
    OS = as.numeric(
      as.character(OS_status5)
    ),
    
    # Convert months to years.
    OS.time = as.numeric(
      as.character(OS_time5)
    ) / 12,
    
    Age = as.numeric(
      as.character(
        Idade_do_diagnostico
      )
    ),
    
    Stage_group = factor(
      as.character(Stage),
      levels = c(
        "III",
        "IV"
      )
    ),
    
    ZNF429_group = factor(
      as.character(ZNF429),
      levels = c(
        "low",
        "high"
      )
    )
  ) %>%
  filter(
    complete.cases(
      OS,
      OS.time,
      Age,
      Stage_group,
      ZNF429_group
    )
  )

cat(
  "\n============================================================\n",
  "INCA NEOADJUVANT DATA SUMMARY\n",
  "============================================================\n",
  sep = ""
)

cat(
  "Samples:",
  nrow(cox_data_neo),
  "\n"
)

cat(
  "Events:",
  sum(cox_data_neo$OS == 1),
  "\n"
)

print(
  table(
    cox_data_neo$ZNF429_group,
    useNA = "ifany"
  )
)

print(
  table(
    cox_data_neo$Stage_group,
    useNA = "ifany"
  )
)

print(
  table(
    cox_data_neo$OS,
    useNA = "ifany"
  )
)


# ==============================================================================
# 6. INCA NEOADJUVANT MULTIVARIATE COX MODEL
# ==============================================================================

fit_cox_neo <- coxph(
  Surv(OS.time, OS) ~
    ZNF429_group +
    Stage_group +
    Age,
  data = cox_data_neo,
  na.action = na.omit,
  x = TRUE,
  y = TRUE
)

print(
  summary(
    fit_cox_neo
  )
)

neo_znf <- ggforest(
  fit_cox_neo,
  data = model.frame(
    fit_cox_neo
  ),
  main = "Multivariate Cox analysis — INCA (Neoadjuvant)",
  fontsize = 1.2
)

print(neo_znf)



# ==============================================================================
# 7. EXTRACT MODEL RESULTS
# ==============================================================================

extract_cox_results <- function(model, cohort_name) {
  
  model_summary <- summary(model)
  
  data.frame(
    Cohort = cohort_name,
    Variable = rownames(
      model_summary$coefficients
    ),
    HR = model_summary$conf.int[
      ,
      "exp(coef)"
    ],
    Lower95 = model_summary$conf.int[
      ,
      "lower .95"
    ],
    Upper95 = model_summary$conf.int[
      ,
      "upper .95"
    ],
    Pvalue = model_summary$coefficients[
      ,
      "Pr(>|z|)"
    ],
    row.names = NULL
  ) %>%
    mutate(
      HR = round(HR, 3),
      Lower95 = round(Lower95, 3),
      Upper95 = round(Upper95, 3),
      Pvalue = signif(Pvalue, 3)
    )
}

cox_results_tcga <- extract_cox_results(
  fit_cox_tcga,
  "TCGA-OV"
)

cox_results_adj <- extract_cox_results(
  fit_cox_adj,
  "INCA Adjuvant"
)

cox_results_neo <- extract_cox_results(
  fit_cox_neo,
  "INCA Neoadjuvant"
)

cox_results_all <- bind_rows(
  cox_results_tcga,
  cox_results_adj,
  cox_results_neo
)

print(cox_results_all)


# ==============================================================================
# 9. OPTIONAL EXPORTS
# ==============================================================================

# Uncomment to save the combined forest plot.
#
# ggsave(
#   filename = "Forest_Cox_ZNF429_TCGA_INCA.png",
#   plot = forest_panel,
#   width = 12,
#   height = 18,
#   dpi = 300
# )
#
# ggsave(
#   filename = "Forest_Cox_ZNF429_TCGA_INCA.pdf",
#   plot = forest_panel,
#   width = 12,
#   height = 18
# )
#
# write.csv(
#   cox_results_all,
#   "Cox_ZNF429_TCGA_INCA_results.csv",
#   row.names = FALSE
# )
