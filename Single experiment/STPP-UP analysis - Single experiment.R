# ------------------------------------------------------------
#  STPP-UP analysis (Single experiment)
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
#  drug treatment and replicates
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
# DEP Analysis (Carrier)
# ------------------------------------------------------------

data <- datalist_carrier[[1]]
value_columns <- which(sapply(data, is.numeric))

# Create experimental design
experimental_design <- data.frame(label = colnames(data[, value_columns])) %>%
  separate_wider_delim(
    label,
    delim = "_",
    names = c("condition", "replicate"),
    cols_remove = FALSE
  ) %>% select(label, condition, replicate)

# DEP filtering and normalization
data_se <- make_se(data, value_columns, experimental_design)
data_filt <- filter_proteins(data_se, type = parameters$filtering_type, thr = 0)
data_norm <- suppressMessages(normalize_vsn(data_filt))
if (parameters$filtering_type == "condition") {
  data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)          #choose imputation method [NOT RECOMMENDED]
} else {
  data_imp <- data_norm
}

# Differential enrichment analysis - Carrier
if (parameters$comparison == "control") {
  data_diff <- test_diff(data_imp, type = "control", control = parameters$control)
} else if (parameters$comparison == "manual") {
  manual_contrasts <- str_split(parameters$contrasts, ";")[[1]]
  data_diff <- test_diff(data_imp, type = "manual", test = manual_contrasts)
} else if (parameters$comparison == "all") {
  data_diff <- test_diff(data_imp, type = "all")
}

# Extract contrasts
contrasts <- sub("_diff", "", str_subset(names(data_diff@elementMetadata@listData), "_diff"))
dep_carrier <- add_rejections(data_diff, alpha = parameters$p_value, lfc = parameters$log2_FC)
results_list[["Carrier"]] <- get_results(dep_carrier)

# QC and volcano plots - Carrier
pdf("Carrier_QC.pdf")
print(plot_frequency(data_se))
print(plot_numbers(data_filt))
if (parameters$filtering_type == "condition") {
  print(plot_detect(data_filt))
}
if (parameters$filtering_type == "condition") {
  print(plot_normalization(data_filt, data_norm, data_imp))
} else {
  print(plot_normalization(data_filt, data_norm))
}
if (parameters$filtering_type == "condition") {
  print(plot_imputation(data_norm, data_imp))
}
for (contrast in contrasts) {
  print(plot_volcano(dep_carrier, contrast = contrast, label_size = 2, add_names = TRUE))
}
dev.off()


# ------------------------------------------------------------
# DEP Analysis (Heated)
# ------------------------------------------------------------

for (j in seq_along(datalist_other)) {
  data <- datalist_other[[j]]
  name <- names(datalist_other[j])
  value_columns <- which(sapply(data, is.numeric))
  
  # Create experimental design
  experimental_design <- data.frame(label = colnames(data[, value_columns])) %>%
    separate_wider_delim(
      label,
      delim = "_",
      names = c("condition", "replicate"),
      cols_remove = FALSE
    ) %>% select(label, condition, replicate)
  
  # DEP filtering and normalization
  data_se <- make_se(data, value_columns, experimental_design)
  data_filt <- filter_proteins(data_se, type = parameters$filtering_type, thr = 0)
  data_norm <- suppressMessages(normalize_vsn(data_filt))
  if (parameters$filtering_type == "condition") {
    data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)          #choose imputation method [NOT RECOMMENDED]
  } else {
    data_imp <- data_norm
  }
  
  # Differential enrichment analysis - Uncorrected
  if (parameters$comparison == "control") {
    data_diff <- test_diff(data_imp, type = "control", control = parameters$control)
  } else if (parameters$comparison == "manual") {
    manual_contrasts <- str_split(parameters$contrasts, ";")[[1]]
    data_diff <- test_diff(data_imp, type = "manual", test = manual_contrasts)
  } else if (parameters$comparison == "all") {
    data_diff <- test_diff(data_imp, type = "all")
  }
  
  contrasts <- sub("_diff", "", str_subset(names(data_diff@elementMetadata@listData), "_diff"))
  dep_not_cor <- add_rejections(data_diff, alpha = parameters$p_value, lfc = parameters$log2_FC)
  results_list[[paste0(name, "_uncorrected")]] <- get_results(dep_not_cor)
  
  # QC and volcano plots - Heated
  pdf(paste0(name, "_QC.pdf"))
  print(plot_frequency(data_se))
  print(plot_numbers(data_filt))
  if (parameters$filtering_type == "condition") {
    print(plot_detect(data_filt))
  }
  if (parameters$filtering_type == "condition") {
    print(plot_normalization(data_filt, data_norm, data_imp))
  } else {
    print(plot_normalization(data_filt, data_norm))
  }
  if (parameters$filtering_type == "condition") {
    print(plot_imputation(data_norm, data_imp))
  }
  dev.off()
  
  
  # ------------------------------------------------------------
  # Carrier correction (per contrast)
  # ------------------------------------------------------------
  
  for (contrast in contrasts) {
    data_imp_cor <- data_imp
    data_corr <- as.data.frame(data_imp_cor@assays@data@listData[[1]])
    data_corr$name <- rownames(data_corr)
    
    base_condition  <- strsplit(contrast, "_vs_")[[1]][2]
    treat_condition <- strsplit(contrast, "_vs_")[[1]][1]
    treat_cols <- grep(treat_condition, colnames(data_corr), value = TRUE)
    
    
    carrier_change <- get_results(dep_carrier) %>%
      dplyr::select(name, ratio = paste0(contrast, "_ratio"))
    
    data_corr <- data_corr %>%
      left_join(carrier_df, by = "name") %>%
      mutate(ratio = ifelse(is.na(ratio), 0, ratio))
    
    data_corr[treat_cols] <- data_corr[treat_cols] - data_corr$ratio
    data_imp_cor@assays@data@listData[[1]][, treat_cols] <- as.matrix(data_corr[treat_cols])
    
    # Differential enrichment analysis - Corrected
    if (parameters$comparison == "control") {
      data_diff_cor <- test_diff(data_imp_cor, type = "control", control = parameters$control)
    } else if (parameters$comparison == "manual") {
      manual_contrasts <- str_split(parameters$contrasts, ";")[[1]]
      data_diff_cor <- test_diff(data_imp_cor, type = "manual", test = manual_contrasts)
    } else if (parameters$comparison == "all") {
      data_diff_cor <- test_diff(data_imp_cor, type = "all")
    }
    
    dep_cor <- add_rejections(data_diff_cor, alpha = parameters$p_value, lfc = parameters$log2_FC)
    results_list[[paste0(name, "_", contrast, "_corrected")]] <- get_results(dep_cor)
    
    
    # ------------------------------------------------------------
    # Volcano plots and QC for this contrast
    # ------------------------------------------------------------
    
    # Extract uncorrected results
    results <- get_results(dep_not_cor) %>%
      select(name, ID, contains(contrast) &
               (contains("p.val") | contains("ratio")))
    colnames(results) <- sub("ratio", "log2.FC", colnames(results))

    p_val_column   <- grep("_p.val$", names(results), value = TRUE)
    log2_fc_column <- grep("_log2.FC$", names(results), value = TRUE)
    
    # Select proteins to highlight
    if (parameters$volcano == "protein list") {
      proteins_to_highlight <- str_split(parameters$protein_list, ";")[[1]]
    }
    if (parameters$volcano == "TopN") {
    sorted_results <- results[order(abs(results[[log2_fc_column]]), decreasing = TRUE), ]
    proteins_to_highlight <- head(sorted_results[sorted_results[[p_val_column]] < parameters$p_value, ], 30)$name
    }
    
    # Extract corrected results
    results_corrected <- get_results(dep_cor) %>%
      select(name, ID, contains(contrast) &
               (contains("p.val") | contains("ratio")))
    
    colnames(results_corrected) <- sub("ratio", "log2.FC", colnames(results_corrected))
    
    
    # Generate volcano plots and QC figures
    pdf(paste0(name, "_", contrast, "_volcano_and_QC.pdf"))
    
    # Uncorrected volcano
    print(
      ggplot(results, aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]]))) +
        geom_point(color = "grey60", alpha = 0.5) +
        geom_point(
          data = results[results$name %in% proteins_to_highlight, ],
          color = "#E41A1C", size = 2
        ) +
        geom_text_repel(
          data = results[results$name %in% proteins_to_highlight, ],
          aes(label = name),
          max.overlaps = Inf, size = 3
        ) +
        geom_hline(yintercept = -log10(parameters$p_value), linetype = "dashed", color = "grey50") +
        geom_vline(xintercept = c(-parameters$log2_FC, parameters$log2_FC), linetype = "dashed", color = "grey50") +
        theme_bw(base_size = 12) +
        xlab(paste0("log2 FC (", treat_condition, " / ", base_condition, ")")) +
        ylab("-log10(p-value)") +
        ggtitle(paste0(name, ": ", treat_condition, " vs ", base_condition, " (uncorrected)"))
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
        ggtitle(paste0(name, ": ", treat_condition, " vs ", base_condition, " (carrier corrected)"))
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
          aes(label = name), size = 3
        ) +
        xlab("log2FC carrier") +
        ylab("log2FC STPP-UP (uncorrected)") +
        ggtitle(paste0(name, ": ", contrast, " QC comparison")) +
        theme_bw(base_size = 12)
    )
    
    dev.off()
  }
}


# ------------------------------------------------------------
# Save workspace and session info
# ------------------------------------------------------------

write.xlsx(results_list, "DEP_results.xlsx")

save.image(file.path(getwd(), "/STPP-UP analysis workspace.RData"))
sessionInfo()