# ==============================================================================
# ZNF429 ONCOPRINT ANALYSIS IN TWO SAMPLE SETS
# ==============================================================================
#
# This script runs the same mutation/CNA workflow twice:
#
#   1. TCGA-selected patients:
#      samples_keep <- substr(tcga_ov_final$bcr_patient_barcode, 1, 12)
#
#   2. All patients available in annot_znf429:
#      samples_keep <- annot_znf429$Tumor_Sample_Barcode
#
# Separate output folders are created to prevent files from being overwritten:
#
#   outputs/TCGA_selected_patients/
#   outputs/All_annotated_patients/
#
#
# the TCGA-OV MAF is loaded
# with TCGAmutations::tcga_load(study = "OV").
# ==============================================================================


# ==============================================================================
# 0. LOAD DATA
# ==============================================================================

load("data/dados_oncoprint_ZNF429.RData")

library(TCGAmutations)
library(maftools)
library(dplyr)
library(tidyr)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplot2)

# ==============================================================================
# 1. INPUT PREPARATION
# ==============================================================================

normalize_tcga_barcode <- function(x) substr(as.character(x), 1, 12)

annot_znf429 <- annot_znf429 %>%
  mutate(Tumor_Sample_Barcode = normalize_tcga_barcode(Tumor_Sample_Barcode))

if (!exists("maf_raw")) {
  maf_ov_znf <- TCGAmutations::tcga_load(study = "OV")
  maf_raw <- maf_ov_znf@data
}

# ==============================================================================
# 2. FIXED PARAMETERS
# ==============================================================================

nonsyn_classes <- c(
  "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
  "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
  "Splice_Site", "Translation_Start_Site"
)

point_mut_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Splice_Site",
  "Nonstop_Mutation", "Translation_Start_Site"
)

blacklist_genes <- c("TTN", "MUC16", "OBSCN", "FLG", "RYR2", "RYR3")

genes_19p12_plot <- c(
  "LINC00664", "ZNF254", "ZNF429", "ZNF430", "ZNF486", "ZNF626",
  "ZNF682", "ZNF726", "ZNF737", "ZNF826P", "ZNF85", "ZNF90",
  "ZNF93", "ZNF431", "ZNF493", "ZNF708", "ZNF714", "ZNF738"
)

genes_plot <- c(
  "TP53", "BRCA1", "BRCA2", "CCNE1", "RB1", "NF1", "PTEN",
  "ARID1A", "CDK12", "PIK3CA", "NOTCH3", "MYC", "PAX8", "ZNF429"
)

genes_all_cna <- unique(c(genes_19p12_plot, genes_plot))

# ==============================================================================
# 3. COLORS AND ONCOPRINT DRAWING FUNCTIONS
# ==============================================================================

znf429_col <- c("high" = "#4a7c59", "low" = "#f79d65")

col <- c(
  "SNV" = "#00778A",
  "INDEL" = "#F4D2E5",
  "AMP" = "#4E475F",
  "DEL" = "#A30248"
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w * 0.96, h * 0.96,
              gp = gpar(fill = "#CFCFCF", col = "white"))
  },
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w * 0.90, h * 0.90,
              gp = gpar(fill = col["SNV"], col = NA))
  },
  INDEL = function(x, y, w, h) {
    grid.rect(x, y, w * 0.90, h * 0.90,
              gp = gpar(fill = col["INDEL"], col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y + h * 0.22, w * 0.90, h * 0.35,
              gp = gpar(fill = col["AMP"], col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x, y - h * 0.22, w * 0.90, h * 0.35,
              gp = gpar(fill = col["DEL"], col = NA))
  }
)

# ==============================================================================
# 4. MAIN ANALYSIS FUNCTION
# ==============================================================================

run_znf429_oncoprint_analysis <- function(samples_keep, analysis_id, run_19p12 = TRUE) {
  
  output_dir <- file.path("outputs", analysis_id)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_file <- function(filename) file.path(output_dir, filename)
  
  samples_keep <- unique(normalize_tcga_barcode(samples_keep))
  samples_keep <- samples_keep[!is.na(samples_keep) & samples_keep != ""]
  
  cat(
    "\n============================================================\n",
    "ANALYSIS: ", analysis_id, "\n",
    "Input samples: ", length(samples_keep), "\n",
    "============================================================\n",
    sep = ""
  )
  
  if (run_19p12) {
    writeLines(genes_19p12_plot, output_file("genes_19p12_plot.txt"))
  }
  
  writeLines(genes_plot, output_file("genes_oncodrivers_plot.txt"))
  
  # ---------------------------------------------------------------------------
  # 4.1 MUTATION DATA
  # ---------------------------------------------------------------------------
  
  maf_filtered <- maf_raw %>%
    mutate(
      Tumor_Sample_Barcode = normalize_tcga_barcode(Tumor_Sample_Barcode),
      VAF = t_alt_count / t_depth
    ) %>%
    filter(
      Tumor_Sample_Barcode %in% samples_keep,
      IMPACT %in% c("HIGH", "MODERATE"),
      Variant_Classification %in% nonsyn_classes,
      t_depth >= 30,
      VAF >= 0.05,
      is.na(ExAC_AF) | ExAC_AF < 0.01,
      is.na(GMAF) | GMAF < 0.01,
      !Hugo_Symbol %in% blacklist_genes
    )
  
  maf_final <- read.maf(maf = maf_filtered, rmFlags = TRUE)
  
  # ---------------------------------------------------------------------------
  # 4.2 CNA DATA
  # ---------------------------------------------------------------------------
  
  data_cna_long <- data_cna %>%
    as.data.frame() %>%
    filter(Hugo_Symbol %in% genes_all_cna) %>%
    pivot_longer(
      cols = starts_with("TCGA"),
      names_to = "Sample",
      values_to = "CNA_score"
    ) %>%
    mutate(
      Tumor_Sample_Barcode = normalize_tcga_barcode(Sample),
      sample_type_code = substr(Sample, 14, 15),
      CNA_score = as.integer(CNA_score),
      CN = case_when(
        CNA_score == 2 ~ "Amp",
        CNA_score == -2 ~ "Del",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      sample_type_code == "01",
      Tumor_Sample_Barcode %in% samples_keep,
      !is.na(CNA_score)
    )
  
  # ---------------------------------------------------------------------------
  # 4.3 COMMON SAMPLES
  # ---------------------------------------------------------------------------
  
  common_samples <- intersect(
    samples_keep,
    unique(c(
      maf_filtered$Tumor_Sample_Barcode,
      data_cna_long$Tumor_Sample_Barcode
    ))
  )
  
  common_samples <- intersect(
    common_samples,
    annot_znf429$Tumor_Sample_Barcode
  )
  
  maf_final@data <- maf_final@data %>%
    filter(Tumor_Sample_Barcode %in% common_samples)
  
  data_cna_long <- data_cna_long %>%
    filter(Tumor_Sample_Barcode %in% common_samples)
  
  cnTable <- data_cna_long %>%
    filter(CN %in% c("Amp", "Del")) %>%
    transmute(Hugo_Symbol, Entrez_Gene_Id, Tumor_Sample_Barcode, CN) %>%
    distinct()
  
  # ---------------------------------------------------------------------------
  # 4.4 FISHER TEST FOR CNA
  # ---------------------------------------------------------------------------
  
  run_fisher_cna_gene <- function(gene_name) {
    
    gene_cna_samples <- cnTable %>%
      filter(Hugo_Symbol == gene_name, CN %in% c("Amp", "Del")) %>%
      pull(Tumor_Sample_Barcode) %>%
      unique()
    
    fisher_gene_df <- annot_znf429 %>%
      filter(Tumor_Sample_Barcode %in% common_samples) %>%
      mutate(
        ZNF429_group = factor(ZNF429_group, levels = c("high", "low")),
        CNA_status = ifelse(
          Tumor_Sample_Barcode %in% gene_cna_samples,
          "CNA_altered",
          "CNA_neutral"
        ),
        CNA_status = factor(
          CNA_status,
          levels = c("CNA_altered", "CNA_neutral")
        )
      )
    
    tab_gene <- table(
      fisher_gene_df$ZNF429_group,
      fisher_gene_df$CNA_status
    )
    
    high_altered <- tab_gene["high", "CNA_altered"]
    high_neutral <- tab_gene["high", "CNA_neutral"]
    low_altered <- tab_gene["low", "CNA_altered"]
    low_neutral <- tab_gene["low", "CNA_neutral"]
    
    high_total <- sum(tab_gene["high", ])
    low_total <- sum(tab_gene["low", ])
    
    high_pct <- 100 * high_altered / high_total
    low_pct <- 100 * low_altered / low_total
    
    if (sum(tab_gene[, "CNA_altered"]) == 0) {
      fisher_p <- NA_real_
      fisher_or <- NA_real_
      ci_low <- NA_real_
      ci_high <- NA_real_
    } else {
      fisher_res <- fisher.test(tab_gene)
      fisher_p <- fisher_res$p.value
      fisher_or <- as.numeric(fisher_res$estimate)
      ci_low <- fisher_res$conf.int[1]
      ci_high <- fisher_res$conf.int[2]
    }
    
    data.frame(
      Hugo_Symbol = gene_name,
      high_total = high_total,
      high_CNA_altered = high_altered,
      high_CNA_neutral = high_neutral,
      high_pct_CNA_altered = high_pct,
      low_total = low_total,
      low_CNA_altered = low_altered,
      low_CNA_neutral = low_neutral,
      low_pct_CNA_altered = low_pct,
      odds_ratio = fisher_or,
      ci_low = ci_low,
      ci_high = ci_high,
      p_value = fisher_p,
      stringsAsFactors = FALSE
    )
  }
  
  fisher_cna_19p12_results <- NULL
  
  if (run_19p12) {
    fisher_cna_19p12_results <- bind_rows(
      lapply(genes_19p12_plot, run_fisher_cna_gene)
    ) %>%
      mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
      arrange(p_value)
    
    fwrite(
      fisher_cna_19p12_results,
      output_file("Fisher_CNA_19p12_genes_ZNF429_high_vs_low.tsv"),
      sep = "\t"
    )
  }
  
  # ---------------------------------------------------------------------------
  # 4.5 CNA FREQUENCIES
  # ---------------------------------------------------------------------------
  
  freq_gene_cna <- data_cna_long %>%
    group_by(Hugo_Symbol, Entrez_Gene_Id) %>%
    summarise(
      n_samples = n_distinct(Tumor_Sample_Barcode),
      AMP_n = sum(CNA_score == 2, na.rm = TRUE),
      DEL_n = sum(CNA_score == -2, na.rm = TRUE),
      GAIN_n = sum(CNA_score == 1, na.rm = TRUE),
      LOSS_n = sum(CNA_score == -1, na.rm = TRUE),
      AMP_pct = 100 * mean(CNA_score == 2, na.rm = TRUE),
      DEL_pct = 100 * mean(CNA_score == -2, na.rm = TRUE),
      GAIN_pct = 100 * mean(CNA_score == 1, na.rm = TRUE),
      LOSS_pct = 100 * mean(CNA_score == -1, na.rm = TRUE),
      Altered_strong_pct = 100 * mean(CNA_score %in% c(2, -2), na.rm = TRUE),
      Altered_any_pct = 100 * mean(CNA_score != 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(Altered_strong_pct))
  
  cna_recurrent_strong <- freq_gene_cna %>%
    filter(AMP_pct >= 10 | DEL_pct >= 5) %>%
    arrange(desc(Altered_strong_pct))
  
  # ---------------------------------------------------------------------------
  # 4.6 CREATE ONCOPRINT MATRIX
  # ---------------------------------------------------------------------------
  
  make_oncoprint_data <- function(genes_use, group_use) {
    
    samples_group <- annot_znf429 %>%
      filter(ZNF429_group == group_use) %>%
      pull(Tumor_Sample_Barcode)
    
    samples_group <- intersect(samples_group, common_samples)
    
    mut_long <- maf_final@data %>%
      filter(
        Hugo_Symbol %in% genes_use,
        Tumor_Sample_Barcode %in% samples_group
      ) %>%
      mutate(
        Alteration = case_when(
          Variant_Classification %in% c(
            "Missense_Mutation", "Nonsense_Mutation", "Splice_Site",
            "Nonstop_Mutation", "Translation_Start_Site"
          ) ~ "SNV",
          Variant_Classification %in% c(
            "Frame_Shift_Del", "Frame_Shift_Ins",
            "In_Frame_Del", "In_Frame_Ins"
          ) ~ "INDEL",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Alteration)) %>%
      select(Hugo_Symbol, Tumor_Sample_Barcode, Alteration) %>%
      distinct()
    
    cnv_plot <- cnTable %>%
      filter(
        Hugo_Symbol %in% genes_use,
        Tumor_Sample_Barcode %in% samples_group
      ) %>%
      mutate(
        Alteration = case_when(
          CN == "Amp" ~ "AMP",
          CN == "Del" ~ "DEL",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Alteration)) %>%
      select(Hugo_Symbol, Tumor_Sample_Barcode, Alteration) %>%
      distinct()
    
    alter_long <- bind_rows(mut_long, cnv_plot) %>%
      group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
      summarise(
        alteration = paste(unique(Alteration), collapse = ";"),
        .groups = "drop"
      )
    
    onco_mat <- alter_long %>%
      complete(
        Hugo_Symbol = genes_use,
        Tumor_Sample_Barcode = samples_group,
        fill = list(alteration = "")
      ) %>%
      pivot_wider(
        names_from = Tumor_Sample_Barcode,
        values_from = alteration
      ) %>%
      as.data.frame()
    
    rownames(onco_mat) <- onco_mat$Hugo_Symbol
    onco_mat$Hugo_Symbol <- NULL
    onco_mat <- as.matrix(onco_mat)
    onco_mat <- onco_mat[genes_use, , drop = FALSE]
    
    sample_order <- colnames(onco_mat)[
      order(-colSums(onco_mat != ""), na.last = TRUE)
    ]
    
    onco_mat <- onco_mat[, sample_order, drop = FALSE]
    
    cat(
      "Group: ", group_use,
      " | Samples: ", length(samples_group),
      " | Altered cells: ", sum(onco_mat != ""),
      "\n",
      sep = ""
    )
    
    list(
      alter_long = alter_long,
      onco_mat = onco_mat,
      sample_order = sample_order,
      samples_group = samples_group
    )
  }
  
  # ---------------------------------------------------------------------------
  # 4.7 PAIRED ONCOPRINT FUNCTIONS
  # ---------------------------------------------------------------------------
  
  force_gene_order <- function(genes_use) {
    if ("ZNF429" %in% genes_use) {
      genes_use <- c("ZNF429", setdiff(genes_use, "ZNF429"))
    }
    genes_use
  }
  
  make_total_bar_mat <- function(high_data, low_data, gene_order_fixed) {
    
    onco_mat_total <- cbind(high_data$onco_mat, low_data$onco_mat)
    onco_mat_total <- as.matrix(onco_mat_total)
    onco_mat_total <- onco_mat_total[gene_order_fixed, , drop = FALSE]
    
    alt_types <- c("SNV", "INDEL", "AMP", "DEL")
    
    total_bar_mat <- sapply(alt_types, function(alteration_type) {
      mat_logical <- matrix(
        grepl(alteration_type, onco_mat_total, fixed = TRUE),
        nrow = nrow(onco_mat_total),
        ncol = ncol(onco_mat_total),
        dimnames = dimnames(onco_mat_total)
      )
      rowSums(mat_logical, na.rm = TRUE)
    })
    
    total_bar_mat <- as.matrix(total_bar_mat)
    rownames(total_bar_mat) <- rownames(onco_mat_total)
    colnames(total_bar_mat) <- alt_types
    
    total_bar_mat[
      gene_order_fixed,
      alt_types,
      drop = FALSE
    ]
  }
  
  make_paired_oncoprint_panel <- function(
    onco_data,
    panel_title,
    gene_order_fixed,
    show_gene_names = TRUE,
    pct_side_use = "left",
    show_legend = TRUE,
    right_bar_mat = NULL
  ) {
    
    row_order_fixed <- gene_order_fixed[
      gene_order_fixed %in% rownames(onco_data$onco_mat)
    ]
    
    # Always create one logical matrix for each alteration type.
    # This prevents oncoPrint() from failing when one group has no alterations.
    # The drawing functions below are vectorized for efficient rendering.
    alteration_types <- c("SNV", "INDEL", "AMP", "DEL")
    
    onco_mat_list <- lapply(alteration_types, function(alteration_type) {
      matrix(
        grepl(alteration_type, onco_data$onco_mat, fixed = TRUE),
        nrow = nrow(onco_data$onco_mat),
        ncol = ncol(onco_data$onco_mat),
        dimnames = dimnames(onco_data$onco_mat)
      )
    })
    
    names(onco_mat_list) <- alteration_types
    
    right_annot <- NULL
    show_builtin_pct <- TRUE
    
    if (pct_side_use == "right") {
      
      pct_values <- round(
        100 * rowSums(onco_data$onco_mat != "") / ncol(onco_data$onco_mat)
      )
      
      pct_labels <- paste0(pct_values, "%")
      names(pct_labels) <- rownames(onco_data$onco_mat)
      pct_labels <- pct_labels[rownames(onco_data$onco_mat)]
      
      if (!is.null(right_bar_mat)) {
        
        right_bar_mat <- right_bar_mat[
          rownames(onco_data$onco_mat),
          c("SNV", "INDEL", "AMP", "DEL"),
          drop = FALSE
        ]
        
        right_annot <- rowAnnotation(
          pct = anno_text(
            pct_labels,
            just = "left",
            location = 0,
            gp = gpar(fontsize = 9),
            width = max_text_width(pct_labels) + unit(3, "mm")
          ),
          `All cohort` = anno_barplot(
            right_bar_mat,
            beside = FALSE,
            border = FALSE,
            width = unit(2.7, "cm"),
            gp = gpar(
              fill = col[colnames(right_bar_mat)],
              col = NA
            ),
            axis_param = list(
              side = "top",
              gp = gpar(fontsize = 8)
            )
          ),
          annotation_name_side = "top",
          annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
          show_annotation_name = c(
            pct = FALSE,
            `All cohort` = TRUE
          )
        )
        
      } else {
        
        right_annot <- rowAnnotation(
          pct = anno_text(
            pct_labels,
            just = "left",
            location = 0,
            gp = gpar(fontsize = 9),
            width = max_text_width(pct_labels) + unit(3, "mm")
          ),
          show_annotation_name = FALSE
        )
      }
      
      show_builtin_pct <- FALSE
    }
    
    oncoPrint(
      onco_mat_list,
      alter_fun = alter_fun,
      alter_fun_is_vectorized = TRUE,
      col = col,
      top_annotation = NULL,
      left_annotation = NULL,
      right_annotation = right_annot,
      column_order = onco_data$sample_order,
      row_order = row_order_fixed,
      column_title = panel_title,
      column_title_gp = gpar(fontsize = 13, fontface = "bold"),
      show_column_names = FALSE,
      show_row_names = show_gene_names,
      row_names_side = "right",
      row_names_gp = gpar(fontsize = 10, fontface = "italic"),
      row_names_max_width = max_text_width(gene_order_fixed) + unit(4, "mm"),
      show_pct = show_builtin_pct,
      pct_side = pct_side_use,
      pct_gp = gpar(fontsize = 9),
      remove_empty_columns = FALSE,
      remove_empty_rows = FALSE,
      show_heatmap_legend = FALSE,
      heatmap_legend_param = list(
        title = "Alterations",
        at = c("SNV", "INDEL", "AMP", "DEL"),
        labels = c("SNV", "INDEL", "Amplification", "Deep deletion")
      )
    )
  }
  
  make_paired_oncoprint <- function(
    genes_use,
    file_prefix,
    width_use = 14,
    height_use = 6
  ) {
    
    gene_order_fixed <- force_gene_order(genes_use)
    
    high_data <- make_oncoprint_data(
      genes_use = gene_order_fixed,
      group_use = "high"
    )
    
    low_data <- make_oncoprint_data(
      genes_use = gene_order_fixed,
      group_use = "low"
    )
    
    high_data$onco_mat <- high_data$onco_mat[
      gene_order_fixed, , drop = FALSE
    ]
    
    low_data$onco_mat <- low_data$onco_mat[
      gene_order_fixed, , drop = FALSE
    ]
    
    total_bar_mat <- make_total_bar_mat(
      high_data = high_data,
      low_data = low_data,
      gene_order_fixed = gene_order_fixed
    )
    
    p_high <- make_paired_oncoprint_panel(
      onco_data = high_data,
      panel_title = paste0(
        "ZNF429 high (N = ",
        ncol(high_data$onco_mat),
        ")"
      ),
      gene_order_fixed = gene_order_fixed,
      show_gene_names = TRUE,
      pct_side_use = "left",
      show_legend = TRUE
    )
    
    p_low <- make_paired_oncoprint_panel(
      onco_data = low_data,
      panel_title = paste0(
        "ZNF429 low (N = ",
        ncol(low_data$onco_mat),
        ")"
      ),
      gene_order_fixed = gene_order_fixed,
      show_gene_names = FALSE,
      pct_side_use = "right",
      show_legend = FALSE,
      right_bar_mat = total_bar_mat
    )
    
    p_combined <- p_high + p_low
    
    # Use one explicit shared legend instead of merging internal legends.
    # This avoids NULL legend grobs in some ComplexHeatmap/grid versions.
    alteration_legend <- Legend(
      title = "Alterations",
      at = c("SNV", "INDEL", "AMP", "DEL"),
      labels = c("SNV", "INDEL", "Amplification", "Deep deletion"),
      legend_gp = gpar(
        fill = unname(col[c("SNV", "INDEL", "AMP", "DEL")]),
        col = NA
      ),
      nrow = 1
    )
    
    pdf(
      output_file(
        paste0(
          file_prefix,
          "_ZNF429_high_vs_low_single_right_barplot.pdf"
        )
      ),
      width = width_use,
      height = height_use
    )
    
    draw(
      p_combined,
      ht_gap = unit(5, "mm"),
      merge_legends = FALSE,
      heatmap_legend_list = list(alteration_legend),
      heatmap_legend_side = "bottom"
    )
    
    dev.off()
    
    png(
      output_file(
        paste0(
          file_prefix,
          "_ZNF429_high_vs_low_single_right_barplot.png"
        )
      ),
      width = 6500,
      height = 2800,
      res = 500
    )
    
    draw(
      p_combined,
      ht_gap = unit(5, "mm"),
      merge_legends = FALSE,
      heatmap_legend_list = list(alteration_legend),
      heatmap_legend_side = "bottom"
    )
    
    dev.off()
    
    total_bar_df <- data.frame(
      Hugo_Symbol = rownames(total_bar_mat),
      total_bar_mat,
      check.names = FALSE
    )
    
    fwrite(
      total_bar_df,
      output_file(
        paste0(
          file_prefix,
          "_total_cohort_right_barplot_counts.tsv"
        )
      ),
      sep = "\t"
    )
    
    list(
      high_data = high_data,
      low_data = low_data,
      total_bar_mat = total_bar_mat,
      p_high = p_high,
      p_low = p_low,
      p_combined = p_combined,
      gene_order_fixed = gene_order_fixed
    )
  }
  
  # ---------------------------------------------------------------------------
  # 4.8 GENERATE ONCOPRINTS
  # ---------------------------------------------------------------------------
  
  drivers_pair <- make_paired_oncoprint(
    genes_use = genes_plot,
    file_prefix = "OncoPrint_HGSOC_oncodrivers_data_cna",
    width_use = 14,
    height_use = 5.5
  )
  
  p19_pair <- NULL
  data_19p12_high <- NULL
  data_19p12_low <- NULL
  
  if (run_19p12) {
    p19_pair <- make_paired_oncoprint(
      genes_use = genes_19p12_plot,
      file_prefix = "OncoPrint_19p12_data_cna",
      width_use = 14,
      height_use = 6
    )
    
    data_19p12_high <- p19_pair$high_data
    data_19p12_low <- p19_pair$low_data
  }
  
  drivers_high_data <- drivers_pair$high_data
  drivers_low_data <- drivers_pair$low_data
  
  # ---------------------------------------------------------------------------
  # 4.9 SAVE TABLES
  # ---------------------------------------------------------------------------
  
  fwrite(
    drivers_high_data$alter_long,
    output_file("OncoPrint_HGSOC_oncodrivers_ZNF429_high_alter_long.tsv"),
    sep = "\t"
  )
  
  fwrite(
    drivers_low_data$alter_long,
    output_file("OncoPrint_HGSOC_oncodrivers_ZNF429_low_alter_long.tsv"),
    sep = "\t"
  )
  
  if (run_19p12) {
    fwrite(
      data_19p12_high$alter_long,
      output_file("OncoPrint_19p12_ZNF429_high_alter_long.tsv"),
      sep = "\t"
    )
    
    fwrite(
      data_19p12_low$alter_long,
      output_file("OncoPrint_19p12_ZNF429_low_alter_long.tsv"),
      sep = "\t"
    )
  }
  
  fwrite(
    cnTable,
    output_file("TCGA_OV_data_cna_AMP_DEL_filtered.tsv"),
    sep = "\t"
  )
  
  fwrite(
    freq_gene_cna,
    output_file("TCGA_OV_data_cna_frequency_filtered.tsv"),
    sep = "\t"
  )
  
  fwrite(
    cna_recurrent_strong,
    output_file("TCGA_OV_data_cna_recurrent_strong.tsv"),
    sep = "\t"
  )
  
  # ---------------------------------------------------------------------------
  # 4.10 CNA FISHER TEST FOR ONCODRIVERS
  # ---------------------------------------------------------------------------
  
  fisher_cna_oncodrivers_results <- bind_rows(
    lapply(genes_plot, run_fisher_cna_gene)
  ) %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
  
  fwrite(
    fisher_cna_oncodrivers_results,
    output_file("Fisher_CNA_oncodrivers_genes_ZNF429_high_vs_low.tsv"),
    sep = "\t"
  )
  
  # ---------------------------------------------------------------------------
  # 4.11 POINT-MUTATION FISHER TEST FOR ONCODRIVERS
  # ---------------------------------------------------------------------------
  
  run_fisher_point_mut_gene <- function(gene_name) {
    
    gene_mut_samples <- maf_final@data %>%
      filter(
        Hugo_Symbol == gene_name,
        Tumor_Sample_Barcode %in% common_samples,
        Variant_Classification %in% point_mut_classes
      ) %>%
      pull(Tumor_Sample_Barcode) %>%
      unique()
    
    fisher_gene_df <- annot_znf429 %>%
      filter(Tumor_Sample_Barcode %in% common_samples) %>%
      mutate(
        ZNF429_group = factor(ZNF429_group, levels = c("high", "low")),
        POINT_MUT_status = ifelse(
          Tumor_Sample_Barcode %in% gene_mut_samples,
          "POINT_mutated",
          "POINT_wildtype"
        ),
        POINT_MUT_status = factor(
          POINT_MUT_status,
          levels = c("POINT_mutated", "POINT_wildtype")
        )
      )
    
    tab_gene <- table(
      fisher_gene_df$ZNF429_group,
      fisher_gene_df$POINT_MUT_status
    )
    
    high_mutated <- tab_gene["high", "POINT_mutated"]
    high_wildtype <- tab_gene["high", "POINT_wildtype"]
    low_mutated <- tab_gene["low", "POINT_mutated"]
    low_wildtype <- tab_gene["low", "POINT_wildtype"]
    
    high_total <- sum(tab_gene["high", ])
    low_total <- sum(tab_gene["low", ])
    
    high_pct <- 100 * high_mutated / high_total
    low_pct <- 100 * low_mutated / low_total
    
    if (sum(tab_gene[, "POINT_mutated"]) == 0) {
      fisher_p <- NA_real_
      fisher_or <- NA_real_
      ci_low <- NA_real_
      ci_high <- NA_real_
    } else {
      fisher_res <- fisher.test(tab_gene)
      fisher_p <- fisher_res$p.value
      fisher_or <- as.numeric(fisher_res$estimate)
      ci_low <- fisher_res$conf.int[1]
      ci_high <- fisher_res$conf.int[2]
    }
    
    data.frame(
      Hugo_Symbol = gene_name,
      high_total = high_total,
      high_POINT_mutated = high_mutated,
      high_POINT_wildtype = high_wildtype,
      high_pct_POINT_mutated = high_pct,
      low_total = low_total,
      low_POINT_mutated = low_mutated,
      low_POINT_wildtype = low_wildtype,
      low_pct_POINT_mutated = low_pct,
      odds_ratio = fisher_or,
      ci_low = ci_low,
      ci_high = ci_high,
      p_value = fisher_p,
      stringsAsFactors = FALSE
    )
  }
  
  fisher_point_mut_oncodrivers_results <- bind_rows(
    lapply(genes_plot, run_fisher_point_mut_gene)
  ) %>%
    mutate(p_adj_BH = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
  
  fwrite(
    fisher_point_mut_oncodrivers_results,
    output_file(
      "Fisher_POINT_MUT_oncodrivers_genes_ZNF429_high_vs_low.tsv"
    ),
    sep = "\t"
  )
  
  # ---------------------------------------------------------------------------
  # 4.12 SUMMARY
  # ---------------------------------------------------------------------------
  
  cat(
    "\nFinalizado: ", analysis_id, "\n",
    "Amostras de entrada: ", length(samples_keep), "\n",
    "Amostras em comum: ", length(common_samples), "\n",
    "High: ", length(drivers_high_data$sample_order), "\n",
    "Low: ", length(drivers_low_data$sample_order), "\n",
    "Genes oncodrivers: ", length(genes_plot), "\n",
    "Genes 19p12: ", ifelse(run_19p12, length(genes_19p12_plot), 0), "\n",
    "Mutações filtradas: ", nrow(maf_filtered), "\n",
    "CNA AMP/DEL filtradas: ", nrow(cnTable), "\n",
    "Diretório de saída: ", output_dir, "\n",
    sep = ""
  )
  
  list(
    analysis_id = analysis_id,
    run_19p12 = run_19p12,
    samples_keep = samples_keep,
    common_samples = common_samples,
    maf_filtered = maf_filtered,
    maf_final = maf_final,
    data_cna_long = data_cna_long,
    cnTable = cnTable,
    freq_gene_cna = freq_gene_cna,
    cna_recurrent_strong = cna_recurrent_strong,
    fisher_cna_19p12_results = fisher_cna_19p12_results,
    fisher_cna_oncodrivers_results = fisher_cna_oncodrivers_results,
    fisher_point_mut_oncodrivers_results =
      fisher_point_mut_oncodrivers_results,
    drivers_pair = drivers_pair,
    p19_pair = p19_pair
  )
}

# ==============================================================================
# 5. RUN 1 — PATIENTS IN tcga_ov_final
# ==============================================================================

samples_keep_tcga <- unique(
  normalize_tcga_barcode(tcga_ov_final$bcr_patient_barcode)
)

results_tcga_selected <- run_znf429_oncoprint_analysis(
  samples_keep = samples_keep_tcga,
  analysis_id = "TCGA_selected_patients",
  run_19p12 = FALSE
)

# ==============================================================================
# 6. RUN 2 — ALL PATIENTS IN annot_znf429
# ==============================================================================

samples_keep_all <- unique(
  normalize_tcga_barcode(annot_znf429$Tumor_Sample_Barcode)
)

results_all_annotated <- run_znf429_oncoprint_analysis(
  samples_keep = samples_keep_all,
  analysis_id = "All_annotated_patients",
  run_19p12 = TRUE
)

# ==============================================================================
# 7. COMPARE BOTH RUNS
# ==============================================================================

run_summary <- data.frame(
  Analysis = c(
    "TCGA_selected_patients",
    "All_annotated_patients"
  ),
  Input_samples = c(
    length(results_tcga_selected$samples_keep),
    length(results_all_annotated$samples_keep)
  ),
  Common_samples = c(
    length(results_tcga_selected$common_samples),
    length(results_all_annotated$common_samples)
  ),
  Filtered_mutations = c(
    nrow(results_tcga_selected$maf_filtered),
    nrow(results_all_annotated$maf_filtered)
  ),
  Strong_CNA_events = c(
    nrow(results_tcga_selected$cnTable),
    nrow(results_all_annotated$cnTable)
  )
)

print(run_summary)

fwrite(
  run_summary,
  file.path("outputs", "ZNF429_analysis_run_summary.tsv"),
  sep = "\t"
)