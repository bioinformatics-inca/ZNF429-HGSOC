# ==============================================================================
# 03_TCGA_INCA_final_panel.R
# ==============================================================================
# Title:
#   Final figure panel combining TCGA internal performance and INCA external
#   validation
#
# Panel contents:
#
#   A–C — TCGA time-dependent ROC curves at 2, 3 and 4 years
#         AUC and Brier score are displayed for each model.
#
#   D–F — INCA external-validation ROC curves at 2, 3 and 4 years
#         Only AUC is displayed.
#
#   G   — TCGA Kaplan–Meier curve
#
#   H   — INCA Kaplan–Meier curve
#
# Important:
#   The Number at risk tables do not receive panel letters.
#
# Execution order:
#
#   source("01_TCGA_internal_validation.R")
#   source("02_INCA_external_validation.R")
#   source("03_TCGA_INCA_final_panel.R")
#
# Scripts 1, 2 and 3 must be run in the same R session.
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
# 1. TCGA ROC-DATA EXTRACTION FUNCTION
# ==============================================================================

extract_roc_df_panel <- function(
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


# ==============================================================================
# 2. TCGA ROC-PLOT FUNCTION
#
# TCGA panels display:
#   - time-dependent AUC
#   - time-dependent Brier score
#
# Models:
#   - Clinical
#   - ZNF429
#   - Combined
# ==============================================================================

make_tcga_roc_plot <- function(
    timepoint,
    auc_index,
    title_text
) {
  
  roc_df <- bind_rows(
    extract_roc_df_panel(
      roc_clinical,
      timepoint,
      "Clinical"
    ),
    
    extract_roc_df_panel(
      roc_molecular,
      timepoint,
      "ZNF429"
    ),
    
    extract_roc_df_panel(
      roc_combined,
      timepoint,
      "Combined"
    )
  )
  
  metric_label <- paste0(
    "Clinical: AUC = ",
    round(
      roc_clinical$AUC[auc_index],
      3
    ),
    " | Brier = ",
    round(
      brier_clinical[
        as.character(timepoint)
      ],
      3
    ),
    
    "\nZNF429: AUC = ",
    round(
      roc_molecular$AUC[auc_index],
      3
    ),
    " | Brier = ",
    round(
      brier_znf429[
        as.character(timepoint)
      ],
      3
    ),
    
    "\nCombined: AUC = ",
    round(
      roc_combined$AUC[auc_index],
      3
    ),
    " | Brier = ",
    round(
      brier_combined[
        as.character(timepoint)
      ],
      3
    )
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
      linewidth = 1.35
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
      y = 0.20,
      hjust = 0,
      size = 5,
      label = metric_label
    ) +
    labs(
      title = title_text,
      x = "False Positive Rate",
      y = "True Positive Rate",
      color = NULL
    ) +
    coord_cartesian(
      xlim = c(
        0,
        1
      ),
      ylim = c(
        0,
        1
      ),
      expand = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Clinical" = "#F8766D",
        "Combined" = "#00BA38",
        "ZNF429" = "#619CFF"
      )
    ) +
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
      ),
      legend.position = "bottom"
    )
}


# ==============================================================================
# 3. INCA ROC-DATA EXTRACTION FUNCTION
# ==============================================================================

extract_roc_inca_panel <- function(
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


# ==============================================================================
# 4. INCA ROC-PLOT FUNCTION
#
# INCA panels display only the external-validation AUC.
# ==============================================================================

make_inca_roc_plot <- function(
    timepoint,
    auc_index,
    title_text
) {
  
  roc_df <- extract_roc_inca_panel(
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
      linewidth = 1.5,
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
      size = 8,
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
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    coord_cartesian(
      xlim = c(
        0,
        1
      ),
      ylim = c(
        0,
        1
      ),
      expand = FALSE
    ) +
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


# ==============================================================================
# 5. BUILD TCGA ROC PANELS
#
# Letters are assigned manually so that only the intended plots receive tags.
# ==============================================================================

p_tcga_2 <- make_tcga_roc_plot(
  timepoint = 2,
  auc_index = 1,
  title_text = "TCGA: 2-Year ROC"
) +
  labs(
    tag = "A"
  )

p_tcga_3 <- make_tcga_roc_plot(
  timepoint = 3,
  auc_index = 2,
  title_text = "TCGA: 3-Year ROC"
) +
  labs(
    tag = "B"
  )

p_tcga_4 <- make_tcga_roc_plot(
  timepoint = 4,
  auc_index = 3,
  title_text = "TCGA: 4-Year ROC"
) +
  labs(
    tag = "C"
  )

tcga_roc_row <-
  p_tcga_2 |
  p_tcga_3 |
  p_tcga_4


# ==============================================================================
# 6. BUILD INCA ROC PANELS
# ==============================================================================

p_inca_2 <- make_inca_roc_plot(
  timepoint = 2,
  auc_index = 1,
  title_text = "INCA: 2-Year ROC"
) +
  labs(
    tag = "D"
  )

p_inca_3 <- make_inca_roc_plot(
  timepoint = 3,
  auc_index = 2,
  title_text = "INCA: 3-Year ROC"
) +
  labs(
    tag = "E"
  )

p_inca_4 <- make_inca_roc_plot(
  timepoint = 4,
  auc_index = 3,
  title_text = "INCA: 4-Year ROC"
) +
  labs(
    tag = "F"
  )

inca_roc_row <-
  p_inca_2 |
  p_inca_3 |
  p_inca_4


# ==============================================================================
# 7. BUILD TCGA KAPLAN–MEIER PANEL
#
# The letter G is assigned only to the survival curve.
# The Number at risk table explicitly receives no tag.
# ==============================================================================

tcga_km_plot <- p_km$plot +
  labs(
    title = "TCGA training cohort",
    tag = "G"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    plot.tag = element_text(
      face = "bold",
      size = 16
    ),
    legend.position = "bottom"
  )

tcga_risk_table <- p_km$table +
  labs(
    tag = NULL
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    plot.tag = element_blank()
  )

tcga_km_full <-
  tcga_km_plot /
  tcga_risk_table +
  plot_layout(
    heights = c(
      3,
      0.9
    )
  )


# ==============================================================================
# 8. BUILD INCA KAPLAN–MEIER PANEL
#
# The letter H is assigned only to the survival curve.
# The Number at risk table explicitly receives no tag.
# ==============================================================================

inca_km_plot <- p_km_inca$plot +
  labs(
    title = "INCA external validation cohort",
    tag = "H"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    plot.tag = element_text(
      face = "bold",
      size = 16
    ),
    legend.position = "bottom"
  )

inca_risk_table <- p_km_inca$table +
  labs(
    tag = NULL
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    plot.tag = element_blank()
  )

inca_km_full <-
  inca_km_plot /
  inca_risk_table +
  plot_layout(
    heights = c(
      3,
      0.9
    )
  )

km_row <-
  tcga_km_full |
  inca_km_full


# ==============================================================================
# 9. ASSEMBLE FINAL MULTIPANEL FIGURE
#
# plot_annotation(tag_levels = "A") is intentionally not used.
#
# This prevents Patchwork from assigning letters automatically to the
# Number at risk tables.
# ==============================================================================

final_panel_all <- (
  tcga_roc_row /
    inca_roc_row /
    km_row +
    plot_layout(
      heights = c(
        1.5,
        1.5,
        1.7
      )
    )
) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16
    )
  )

print(final_panel_all)


# ==============================================================================
# 10. EXPORT HIGH-RESOLUTION FIGURES
# ==============================================================================

ggsave(
  filename = "Figure_S11.svg",
  plot = final_panel_all,
  width = 20,
  height = 18,
  dpi = 300
)

ggsave(
  filename = "Figure_S11.pdf",
  plot = final_panel_all,
  width = 20,
  height = 18
)
