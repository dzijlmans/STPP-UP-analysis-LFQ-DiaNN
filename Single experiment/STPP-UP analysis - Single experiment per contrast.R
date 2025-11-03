# ------------------------------------------------------------
#  STPP-UP analysis (Single experiment, per contrast)
#
#  Description:
#  This script performs a full analysis of DIA-NN processed STPP-UP proteomics data:
#    - Import and filtering of protein data
#    - DEP normalization and differential analysis
#    - Carrier correction to adjust for baseline expression differences
#    - Plotting of QC and volcano plots
#    - Exports results to browsable Excel file
#
#  Requirements:
#    Input files (in working directory):
#      - pg_matrix.tsv                      (protein data; DIA-NN output)
#      - pr_matrix.tsv                      (peptide data; DIA-NN output)
#      - analysis_parameters_stppup.xlsx    (parameter file)
#  
#  Considerations:
#  This script was designed to be used for the analysis of a single STPPUP 
#  experiment, with the only experimental variables being heating type,
#  drug treatment and replicates. In this script, the data is separated 
#  by contrast (comparison) instead of being analyzed as a whole.
# ------------------------------------------------------------


# ------------------------------------------------------------
# Loading libraries and setting parameters
# ------------------------------------------------------------

# Required packages
library(tidyverse)
library(DEP)
library(openxlsx)
library(ggrepel)
library(tcltk)
library(PTXQC)

# Set working directory
tkmessageBox(message="Choose your working directory")
setwd(choose.dir())

# User instructions
tkmessageBox(message = paste(
  "Please ensure that abundance column names (and thus the .raw files)",
  "contain the following information, separated by underscores:\n\n",
  "- heating conditions (e.g. Carrier)\n",
  "- compounds used (e.g. C1)\n",
  "- replicate  (e.g. 1)\n\n",
  "Example: 20250205_Astral_AUR_DZ114_Carrier_C1_1.raw"
))

# Define parameters
parameters <- read.xlsx("analysis_parameters_stppup.xlsx")


# ------------------------------------------------------------
# Data import and preprocessing
# ------------------------------------------------------------

# Import data
proteinGroups <- read.delim(list.files(pattern = "pg_matrix.tsv"))
peptides <- read.delim(list.files(pattern = "pr_matrix.tsv"))

# Filter proteins with more than N number of peptides
filtered_proteinGroups <- data.frame(Protein.Group = names(table(unique(peptides[, c('Protein.Group', 'Stripped.Sequence')])$Protein.Group)),
                                     Peptide.Count = as.integer(table(unique(peptides[, c('Protein.Group', 'Stripped.Sequence')])$Protein.Group))) %>%
  filter(Peptide.Count > (parameters$unique_peptides - 1))

proteinGroups <- proteinGroups[proteinGroups$Protein.Group %in% filtered_proteinGroups$Protein.Group, ]

# Remove contaminants and keratins
proteinGroups <- proteinGroups %>%
  make_unique("Genes", "Protein.Group", delim = ";") %>% 
  filter(
    !grepl("Cont_", ID),
      !grepl("KRT", name),
      !grepl("cRAP", name),
      !grepl("cRAP", Protein.Group)
  )

# Tidy column names
value_columns <- which(str_detect(names(proteinGroups), parameters$common_name_value_columns))
prefix <- LCSn(colnames(proteinGroups[, value_columns]))

colnames(proteinGroups) <- gsub(
  ".raw", "",
    gsub(paste0(".*", prefix), "", colnames(proteinGroups)))

# Transform data to long format
proteinGroups <-
  proteinGroups %>% pivot_longer(
    cols = all_of(value_columns),
    names_to = c("heating", "compound", "replicate"),
    names_sep = "_",
    values_to = "values"
  ) %>% pivot_wider(names_from = c("compound", "replicate"),
                    values_from = "values")

# Split dataset into carrier & heating conditions
datalist <- split(proteinGroups, f = ~ proteinGroups$heating)                   
datalist_carrier <- datalist[names(datalist) %in% c("Carrier")]
datalist_other <- datalist[!names(datalist) %in% c("Carrier")]

results_list <- list()


# ------------------------------------------------------------
# Contrast extraction
# ------------------------------------------------------------

# Helper function
extract_conditions <- function(df) {
  sample_cols <- colnames(df)[grepl("_[0-9]+$", colnames(df))]
  # remove trailing _replicate digits
  conds <- unique(sub("_[0-9]+$", "", sample_cols))
  return(conds)
}

contrasts <- NULL
if (parameters$comparison == "manual") {
  raw <- parameters$contrasts
  contrasts <- str_split(raw, "[;]")[[1]] %>% trimws()
} else if (parameters$comparison == "control") {
  ctrl <- parameters$control
  conds <- extract_conditions(datalist_carrier[[1]])
  contrasts <- paste0(conds[!conds %in% ctrl], "_vs_", ctrl)
} else if (tolower(parameters$comparison) == "all") {
  conds <- extract_conditions(datalist_carrier[[1]])
  pairs <- combn(conds, 2, simplify = FALSE)
  contrasts <- vapply(pairs, function(x) paste0(x[1], "_vs_", x[2]), character(1))
}

message("Contrasts to run: ", paste(contrasts, collapse = ", "))


# ------------------------------------------------------------
# Per contrast analysis loop
# ------------------------------------------------------------

for (contrast in contrasts) {
  
  message("Analyzing contrast ", contrast, "...")
  
  base_condition  <- strsplit(contrast, "_vs_")[[1]][2]
  treat_condition <- strsplit(contrast, "_vs_")[[1]][1]

  
  # ------------------------------------------------------------
  # DEP Analysis (Carrier)
  # ------------------------------------------------------------
  
  data_carrier_full <- datalist_carrier[[1]]
  keep_cols <- grep(paste0("^(", treat_condition, "|", base_condition, ")_[0-9]+$"),
                            colnames(data_carrier_full), value = TRUE)
  
  data_carrier <- data_carrier_full[, c("name", "ID", keep_cols)]
  value_columns_carrier <- which(sapply(data_carrier, is.numeric))
  
  # Create experimental design
  experimental_design_carrier <- data.frame(label = colnames(data_carrier[, value_columns_carrier])) %>%
    separate_wider_delim(
      label,
      delim = "_",
      names = c("condition", "replicate"),
      cols_remove = FALSE
    ) %>% select(label, condition, replicate)
  
  # DEP filtering and normalization
  data_carrier_se <- make_se(data_carrier, value_columns_carrier, experimental_design_carrier)
  data_carrier_filt <- filter_proteins(data_carrier_se, type = parameters$filtering_type, thr = 0)
  data_carrier_norm <- suppressMessages(normalize_vsn(data_carrier_filt))
  if (parameters$filtering_type == "condition") {
    data_carrier_imp <- impute(data_carrier_norm, fun = "MinProb", q = 0.01)          #choose imputation method [NOT RECOMMENDED]
  } else {
    data_carrier_imp <- data_carrier_norm
  }
  
  # Differential enrichment analysis - Carrier
  data_carrier_diff <-suppressMessages(test_diff(data_carrier_imp, type = "manual", test = contrast))
  dep_carrier <- add_rejections(data_carrier_diff, alpha = parameters$p_value, lfc = parameters$log2_FC)
  results_list[[paste0("Carrier_", contrast)]] <- get_results(dep_carrier)
  
  # QC and volcano plot - Carrier
  pdf(paste0("Carrier_", contrast, ".pdf"))
  print(plot_frequency(data_carrier_se))
  print(plot_numbers(data_carrier_filt))
  if (parameters$filtering_type == "condition") {
    print(plot_detect(data_carrier_filt))
  }
  if (parameters$filtering_type == "condition") {
    print(plot_normalization(data_carrier_filt, data_carrier_norm, data_carrier_imp))
  } else {
    print(plot_normalization(data_carrier_filt, data_carrier_norm))
  }
  if (parameters$filtering_type == "condition") {
    print(plot_imputation(data_carrier_norm, data_carrier_imp))
  }
  print(plot_volcano(dep_carrier, contrast = contrast, label_size = 2, add_names = TRUE))
  
  dev.off()
  
  
  # ------------------------------------------------------------
  # DEP Analysis (Heated)
  # ------------------------------------------------------------
  
  for (j in seq_along(datalist_other)) {
    data_full_heat <- datalist_other[[j]]
    heat_name <- names(datalist_other[j])
    
    data_heat <- data_full_heat[, c("name", "ID", keep_cols)]
    value_columns_heat <- which(sapply(data_heat, is.numeric))
    
    # Create experimental design
    experimental_design_heat <- data.frame(label = colnames(data_heat[, value_columns_heat])) %>%
      separate_wider_delim(
        label,
        delim = "_",
        names = c("condition", "replicate"),
        cols_remove = FALSE
      ) %>% select(label, condition, replicate)
    
    # DEP filtering and normalization
    data_heat_se <- make_se(data_heat, value_columns_heat, experimental_design_heat)
    data_heat_filt <- filter_proteins(data_heat_se, type = parameters$filtering_type, thr = 0)
    data_heat_norm <- suppressMessages(normalize_vsn(data_heat_filt))
    if (parameters$filtering_type == "condition") {
      data_heat_imp <- impute(data_heat_norm, fun = "MinProb", q = 0.01)          #choose imputation method [NOT RECOMMENDED]
    } else {
      data_heat_imp <- data_heat_norm
    }
    
    # Differential enrichment analysis - Uncorrected
    data_heat_diff <-suppressMessages(test_diff(data_heat_imp, type = "manual", test = contrast))
    dep_not_cor <- add_rejections(data_heat_diff, alpha = parameters$p_value, lfc = parameters$log2_FC)
    results_list[[paste0(heat_name, "_", contrast, "_uncorrected")]] <- get_results(dep_not_cor)
    
    
    # ------------------------------------------------------------
    # Carrier correction
    # ------------------------------------------------------------
    
    data_heat_imp_cor <- data_heat_imp
    data_corr <- as.data.frame(data_heat_imp_cor@assays@data@listData[[1]])
    data_corr$name <- rownames(data_corr)
    
    treat_cols <- grep(treat_condition, colnames(data_corr), value = TRUE)
    
    carrier_df <- get_results(dep_carrier) %>%
      dplyr::select(name, ratio = paste0(contrast, "_ratio"))
    
    data_corr <- data_corr %>%
      left_join(carrier_df, by = "name") %>%
      mutate(ratio = ifelse(is.na(ratio), 0, ratio))
    
    data_corr[treat_cols] <- data_corr[treat_cols] - data_corr$ratio
    data_heat_imp_cor@assays@data@listData[[1]][, treat_cols] <- as.matrix(data_corr[treat_cols])
    
    # Differential enrichment analysis - Corrected
    data_heat_diff_cor <-suppressMessages(test_diff(data_heat_imp_cor, type = "manual", test = contrast))
    dep_cor <- add_rejections(data_heat_diff_cor, alpha = parameters$p_value, lfc = parameters$log2_FC)
    results_list[[paste0(heat_name, "_", contrast, "_corrected")]] <- get_results(dep_cor)
    
    
    # ------------------------------------------------------------
    # Volcano plots and QCs
    # ------------------------------------------------------------
    
    # Extract uncorrected results
    results_uncorrected <- get_results(dep_not_cor) %>%
      select(name, ID, contains(contrast) &
               (contains("p.val") | contains("ratio")))
    colnames(results_uncorrected) <- sub("ratio", "log2.FC", colnames(results_uncorrected))
    
    p_val_column   <- grep("_p.val$", names(results_uncorrected), value = TRUE)
    log2_fc_column <- grep("_log2.FC$", names(results_uncorrected), value = TRUE)
    
    # Select proteins to highlight
    if (parameters$volcano == "protein list") {
      proteins_to_highlight <- str_split(parameters$protein_list, ";")[[1]]
    }
    if (parameters$volcano == "TopN") {
      sorted_results <- results_uncorrected[order(abs(results_uncorrected[[log2_fc_column]]), decreasing = TRUE), ]
      proteins_to_highlight <- head(sorted_results[sorted_results[[p_val_column]] < parameters$p_value, ], parameters$TopN)$name
    }
    
    # Extract corrected results
    results_corrected <- get_results(dep_cor) %>%
      select(name, ID, contains(contrast) &
               (contains("p.val") | contains("ratio")))
    
    colnames(results_corrected) <- sub("ratio", "log2.FC", colnames(results_corrected))
  
    # Generate volcano plots and QC figures
    pdf(paste0(heat_name, "_", contrast, ".pdf"))
    print(plot_frequency(data_heat_se))
    print(plot_numbers(data_heat_filt))
    if (parameters$filtering_type == "condition") {
      print(plot_detect(data_heat_filt))
    }
    if (parameters$filtering_type == "condition") {
      print(plot_normalization(data_heat_filt, data_heat_norm, data_heat_imp))
    } else {
      print(plot_normalization(data_heat_filt, data_heat_norm))
    }
    if (parameters$filtering_type == "condition") {
      print(plot_imputation(data_heat_norm, data_heat_imp))
    }
    
    
    # Uncorrected volcano
    print(
      ggplot(results_uncorrected, aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]]))) +
        geom_point(color = "grey60", alpha = 0.5) +
        geom_point(
          data = results_uncorrected[results_uncorrected$name %in% proteins_to_highlight, ],
          color = "#E41A1C", size = 2
        ) +
        geom_text_repel(
          data = results_uncorrected[results_uncorrected$name %in% proteins_to_highlight, ],
          aes(label = name),
          max.overlaps = Inf, size = 3
        ) +
        geom_hline(yintercept = -log10(parameters$p_value), linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = c(-parameters$log2_FC, parameters$log2_FC), linetype = "dashed", color = "grey50") +
        theme_bw(base_size = 12) +
        xlab(paste0("log2 FC (", treat_condition, " / ", base_condition, ")")) +
        ylab("-log10(p-value)") +
        ggtitle(paste0(heat_name, ": ", treat_condition, " vs ", base_condition, " (uncorrected)"))
    )
    
    # Corrected volcano
    print(
      ggplot(results_corrected, aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]]))) +
        geom_point(color = "grey60", alpha = 0.5) +
        geom_point(
          data = results_corrected[results_corrected$name %in% proteins_to_highlight, ],
          color = "#377EB8", size = 2
        ) +
        geom_text_repel(
          data = results_corrected[results_corrected$name %in% proteins_to_highlight, ],
          aes(label = name),
          max.overlaps = Inf, size = 3
        ) +
        geom_hline(yintercept = -log10(parameters$p_value), linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = c(-parameters$log2_FC, parameters$log2_FC), linetype = "dashed", color = "grey50") +
        theme_bw(base_size = 12) +
        xlab(paste0("log2 FC (", treat_condition, " / ", base_condition, ")")) +
        ylab("-log10(p-value)") +
        ggtitle(paste0(heat_name, ": ", treat_condition, " vs ", base_condition, " (carrier corrected)"))
    )
    
    # Compare uncorrected vs carrier ratios
    merged_data <- merge(
      get_results(dep_carrier), get_results(dep_not_cor),
      by = c("name", "ID")
    ) %>%
      select(name, ID, contains(contrast) & contains("ratio"))
    colnames(merged_data) <- c("name", "ID", "carrier", "uncorrected")
    
    print(
      ggplot(merged_data, aes(x = carrier, y = uncorrected)) +
        geom_point(color = "grey60", alpha = 0.5) +
        geom_point(
          data = merged_data[merged_data$name %in% proteins_to_highlight, ],
          color = "#E41A1C", size = 2
        ) +
        geom_text_repel(
          data = merged_data[merged_data$name %in% proteins_to_highlight, ],
          aes(label = name), max.overlaps = Inf, size = 3
        ) +
        xlab("log2FC carrier") +
        ylab("log2FC STPP-UP (uncorrected)") +
        ggtitle(paste0(heat_name, ": ", contrast, " QC comparison")) +
        theme_bw(base_size = 12)
    )
    
    dev.off()

  }
  
  
  # ------------------------------------------------------------
  # Write Excel file per contrast
  # ------------------------------------------------------------
  
  wb <- createWorkbook()
  
  addWorksheet(wb, "Carrier")
  addWorksheet(wb, paste0(heat_name, "_uncorrected"))
  addWorksheet(wb, paste0(heat_name, "_corrected"))
  
  writeData(wb, "Carrier", results_list[[paste0("Carrier_", contrast)]])
  writeData(wb, paste0(heat_name, "_uncorrected"), results_list[[paste0(heat_name, "_", contrast, "_uncorrected")]])
  writeData(wb, paste0(heat_name, "_corrected"), results_list[[paste0(heat_name, "_", contrast, "_corrected")]])
  
  saveWorkbook(wb, paste0("DEP_results_", contrast, ".xlsx"), overwrite = TRUE)
  
}


# ------------------------------------------------------------
# Save workspace and session info
# ------------------------------------------------------------

save.image(file.path(getwd(), "/STPP-UP analysis workspace.RData"))
sessionInfo()