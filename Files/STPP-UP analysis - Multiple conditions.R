# ------------------------------------------------------------
#  STPP-UP analysis (Multiple conditions, per contrast)
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
#  This script was designed to be used for the analysis of a STPPUP 
#  experiment containing an additional grouping variable other than heating type,
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
library(SummarizedExperiment)

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

# Helper function
extract_conditions <- function(df) {
  sample_cols <- colnames(df)[grepl("_[0-9]+$", colnames(df))]
  # remove trailing _replicate digits
  conds <- unique(sub("_[0-9]+$", "", sample_cols))
  return(conds)
}

make_df <- function(se, suffix) {
  as.data.frame(assay(se)) %>%
    setNames(paste0(names(.), suffix)) %>%
    tibble::rownames_to_column("name")
}

impute_mixed <- function(SE,
                         MAR_method = "knn",
                         MNAR_method = "MinProb",
                         MAR_args = list(k = 10),
                         MNAR_args = list(q = 0.01)) {
  require(tidyverse)
  require(DEP)
  require(SummarizedExperiment)
  
  if (MAR_method == MNAR_method) {
    print("Imputing for MAR and MNAR together (non-mixed)")
    SE_final <- do.call(DEP::impute, c(list(SE, fun = MAR_method), MAR_args))
    return(SE_final)
  } else {
    df_long <- get_df_long(SE)
    
    proteins_MNAR <- df_long %>%
      group_by(name, condition) %>%
      summarize(all_missing = all(is.na(intensity)), .groups = "drop") %>%
      filter(all_missing) %>%
      pull(name) %>%
      unique()
    
    X <- assay(SE)
    MNAR_rows <- rownames(X) %in% proteins_MNAR
    mnar_matrix <- is.na(X) & MNAR_rows
    
    SE_MAR <- do.call(DEP::impute, c(list(SE, fun = MAR_method), MAR_args))
    Y <- assay(SE_MAR)
    Y[mnar_matrix] <- NA
    assay(SE_MAR) <- Y
    SE_final <- do.call(DEP::impute, c(list(SE_MAR, fun = MNAR_method), MNAR_args))
    return(SE_final)
    
  }
  
}


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
    names_to = c("heating", "condition_b", "compound", "replicate"),
    names_sep = "_",
    values_to = "values"
  ) %>% pivot_wider(names_from = c("compound", "replicate"),
                    values_from = "values")

# Split dataset by condition b
datalist_condition_b <- split(proteinGroups, f = proteinGroups$condition_b)

results_list <- list()

# ------------------------------------------------------------
# Analysis loop per secondary condition
# ------------------------------------------------------------

for (i in seq_along(datalist_condition_b)) {
  
  # Split dataset into carrier & heating conditions
  name_condition_b <- names(datalist_condition_b[i])
  dataset <- datalist_condition_b[[i]]
  datalist <- split(dataset, f = dataset$heating)                   #specify variable to be split on
  datalist_carrier <- datalist[names(datalist) %in% c("Carrier")]
  datalist_other <- datalist[!names(datalist) %in% c("Carrier")]
  
  # ------------------------------------------------------------
  # Contrast extraction
  # ------------------------------------------------------------
  
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
  
  message("Contrasts to run for ", name_condition_b, ": ", paste(contrasts, collapse = ", "))
  
  
  # ------------------------------------------------------------
  # Per contrast analysis loop
  # ------------------------------------------------------------
  
  for (contrast in contrasts) {
    
    set.seed(123)
    
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
      data_carrier_imp <- suppressMessages(impute_mixed(data_carrier_norm))          #choose imputation method [NOT RECOMMENDED]
    } else {
      data_carrier_imp <- data_carrier_norm
    }
    
    # Differential enrichment analysis - Carrier
    data_carrier_diff <-suppressMessages(test_diff(data_carrier_imp, type = "manual", test = contrast))
    dep_carrier <- add_rejections(data_carrier_diff, alpha = parameters$p_value, lfc = parameters$log2_FC)
    
    # Extract Carrier results
    results <- get_results(dep_carrier) %>%
      select(name, ID, contains(contrast) &
               (contains("p.val") | contains("ratio")))
    colnames(results) <- sub("ratio", "log2.FC", colnames(results))
    results <- results %>%
      mutate(!!paste0(contrast, "_p.adj") :=
               p.adjust(get(paste0(contrast, "_p.val")), method = "BH"))
    
    raw_df  <- make_df(data_carrier_filt, "_raw")
    norm_df <- make_df(data_carrier_norm, "_normalized")
    merged_data <- merge(raw_df, norm_df, by = "name")
    if (parameters$filtering_type == "condition") {
      imp_df <- make_df(data_carrier_imp, "_imputed")
      merged_data <- merge(merged_data, imp_df, by = "name")
    }
    results_list[[paste0(name_condition_b, "_Carrier_", contrast)]] <- merged_data %>%
      merge(results, by = "name")

    # QC and volcano plot - Carrier
    pdf(paste0(name_condition_b, "_Carrier_", contrast, ".pdf"))
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
    
    # Set up Excel file
    wb <- createWorkbook()
    addWorksheet(wb, "Carrier")
    writeData(wb, "Carrier", results_list[[paste0(name_condition_b, "_Carrier_", contrast)]])
    
    
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
        data_heat_imp <- impute_mixed(data_heat_norm)          #choose imputation method [NOT RECOMMENDED]
      } else {
        data_heat_imp <- data_heat_norm
      }
      
      # Differential enrichment analysis - Uncorrected
      data_heat_diff <-suppressMessages(test_diff(data_heat_imp, type = "manual", test = contrast))
      dep_not_cor <- add_rejections(data_heat_diff, alpha = parameters$p_value, lfc = parameters$log2_FC)

      
      # ------------------------------------------------------------
      # Carrier correction
      # ------------------------------------------------------------
      
      if (parameters$filtering_type == "complete") {
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
      }
      
      
      # ------------------------------------------------------------
      # Volcano plots and QCs
      # ------------------------------------------------------------
      
      # Extract uncorrected results
      results_uncorrected <- get_results(dep_not_cor) %>%
        select(name, ID, contains(contrast) &
                 (contains("p.val") | contains("ratio")))
      colnames(results_uncorrected) <- sub("ratio", "log2.FC", colnames(results_uncorrected))
      
      results_uncorrected <- results_uncorrected %>%
        mutate(!!paste0(contrast, "_p.adj") :=
                 p.adjust(get(paste0(contrast, "_p.val")), method = "BH"))
      
      raw_df  <- make_df(data_heat_filt, "_raw")
      norm_df <- make_df(data_heat_norm, "_normalized")
      merged_data <- merge(raw_df, norm_df, by = "name")
      if (parameters$filtering_type == "condition") {
        imp_df <- make_df(data_heat_imp, "_imputed")
        merged_data <- merge(merged_data, imp_df, by = "name")
      }
      results_list[[paste0(name_condition_b, "_", heat_name, "_", contrast, "_uncorrected")]] <- merged_data %>%
        merge(results_uncorrected, by = "name")
      
      # Imputation only - remove MNAR proteins for plotting
      if (parameters$filtering_type == "condition") {
        results_uncorrected_complete <- results_uncorrected
        proteins_MNAR <- get_df_long(data_heat_norm) %>%
          group_by(name, condition) %>%
          summarize(all_missing = all(is.na(intensity)), .groups = "drop") %>%
          filter(all_missing) %>%
          pull(name) %>%
          unique()
        results_uncorrected <- results_uncorrected[!results_uncorrected$name%in%proteins_MNAR, ]
      }
      
      # Select proteins to highlight
      p_val_column   <- grep("_p.val$", names(results_uncorrected), value = TRUE)
      log2_fc_column <- grep("_log2.FC$", names(results_uncorrected), value = TRUE)
      
      if (parameters$volcano == "protein list") {
        proteins_to_highlight <- str_split(parameters$protein_list, ";")[[1]]
      }
      if (parameters$volcano == "TopN") {
        sorted_results <- results_uncorrected[order(abs(results_uncorrected[[log2_fc_column]]), decreasing = TRUE), ]
        proteins_to_highlight <- head(sorted_results[sorted_results[[p_val_column]] < parameters$p_value, ], parameters$TopN)$name
      }
      
      # Extract corrected results
      if (parameters$filtering_type == "complete") {
      results_corrected <- get_results(dep_cor) %>%
        select(name, ID, contains(contrast) &
                 (contains("p.val") | contains("ratio")))
      
      colnames(results_corrected) <- sub("ratio", "log2.FC", colnames(results_corrected))
      
      results_corrected <- results_corrected %>%
        mutate(!!paste0(contrast, "_p.adj") :=
                 p.adjust(get(paste0(contrast, "_p.val")), method = "BH"))
      
      raw_df  <- make_df(data_heat_filt, "_raw")
      norm_df <- make_df(data_heat_norm, "_normalized")
      merged_data <- merge(raw_df, norm_df, by = "name")
      cor_df <- make_df(data_heat_imp_cor, "_corrected")
      merged_data <- merge(merged_data, cor_df, by = "name")
      results_list[[paste0(name_condition_b, "_", heat_name, "_", contrast, "_corrected")]] <- merged_data %>%
        merge(results_corrected, by = "name")
      }
      
      # Generate volcano plots and QC figures
      pdf(paste0(name_condition_b, "_", heat_name, "_", contrast, ".pdf"))
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
          ggtitle(paste0(name_condition_b, " - ", heat_name, ": ", treat_condition, " vs ", base_condition, " (uncorrected)"))
      )
      
      # Corrected volcano (Carrier correction)
      if (parameters$filtering_type == "complete") {
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
          ggtitle(paste0(name_condition_b, " - ", heat_name, ": ", treat_condition, " vs ", base_condition, " (carrier corrected)"))
      )
        }
      
      # Imputed highlighted volcano (Imputation)
      if (parameters$filtering_type == "condition") {
        proteins_imputed <- get_df_long(data_heat_norm) %>%
          group_by(name, condition) %>%
          summarize(NAs = any(is.na(intensity))) %>%
          filter(NAs) %>%
          pull(name) %>%
          unique()
        results_uncorrected$imputed <- ifelse(results_uncorrected$name %in% proteins_imputed, "Yes", "No")
        
        #MAR
        print(
          ggplot(results_uncorrected, aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]]), color = imputed, shape = imputed, fill = imputed)) +
            geom_point(alpha = 0.5) +
            scale_shape_manual(name = "Imputed", values = c("No" = 21, "Yes" = 22)) +
            scale_fill_manual(name = "Imputed", values = c("No" = "#72B173", "Yes" = "#492050")) +
            scale_color_manual(name = "Imputed", values = c("No" = "#72B173", "Yes" = "#492050")) +
            geom_text_repel(
              data = results_uncorrected[results_uncorrected$name %in% proteins_to_highlight, ],
              aes(label = name),
              max.overlaps = Inf, size = 3, color = "black"
            ) +
            geom_hline(yintercept = -log10(parameters$p_value), linetype = "dashed", color = "grey50") +
            geom_vline(xintercept = c(-parameters$log2_FC, parameters$log2_FC), linetype = "dashed", color = "grey50") +
            theme_bw(base_size = 12) +
            xlab(paste0("log2 FC (", treat_condition, " / ", base_condition, ")")) +
            ylab("-log10(p-value)") +
            ggtitle(paste0(name_condition_b, " - ", heat_name, ": ", treat_condition, " vs ", base_condition, " (uncorrected) - Imputed MAR only")) +
            theme(
              legend.position = c(0.90, 0.10),   # (x, y) coordinates in relative units (0–1)
              legend.background = element_rect(fill = alpha("white", 0.8), color = "black"),
              legend.key = element_blank()
            ) +
            guides(fill = guide_legend(override.aes = list(shape = c(21, 22))),
                   color = guide_legend(override.aes = list(color = c("#72B173", "#492050"))),
                   shape = guide_legend(override.aes = list(fill = c("#72B173", "#492050"))))
        )
      
      #MAR + MNAR
      if (parameters$volcano == "TopN") {
        sorted_results <- results_uncorrected_complete[order(abs(results_uncorrected_complete[[log2_fc_column]]), decreasing = TRUE), ]
        proteins_to_highlight_ <- head(sorted_results[sorted_results[[p_val_column]] < parameters$p_value, ], parameters$TopN)$name
      }
      results_uncorrected_complete$imputed <- ifelse(results_uncorrected_complete$name %in% proteins_imputed, "Yes", "No")
      
      print(
        ggplot(results_uncorrected_complete, aes(x = .data[[log2_fc_column]], y = -log10(.data[[p_val_column]]), color = imputed, shape = imputed, fill = imputed)) +
          geom_point(alpha = 0.5) +
          scale_shape_manual(name = "Imputed", values = c("No" = 21, "Yes" = 22)) +
          scale_fill_manual(name = "Imputed", values = c("No" = "#72B173", "Yes" = "#492050")) +
          scale_color_manual(name = "Imputed", values = c("No" = "#72B173", "Yes" = "#492050")) +
          geom_text_repel(
            data = results_uncorrected_complete[results_uncorrected_complete$name %in% proteins_to_highlight_, ],
            aes(label = name),
            max.overlaps = Inf, size = 3, color = "black"
          ) +
          geom_hline(yintercept = -log10(parameters$p_value), linetype = "dashed", color = "grey50") +
          geom_vline(xintercept = c(-parameters$log2_FC, parameters$log2_FC), linetype = "dashed", color = "grey50") +
          theme_bw(base_size = 12) +
          xlab(paste0("log2 FC (", treat_condition, " / ", base_condition, ")")) +
          ylab("-log10(p-value)") +
          ggtitle(paste0(name_condition_b, " - ", heat_name, ": ", treat_condition, " vs ", base_condition, " (uncorrected) - Imputed MAR & MNAR")) +
          theme(
            legend.position = c(0.90, 0.10),   # (x, y) coordinates in relative units (0–1)
            legend.background = element_rect(fill = alpha("white", 0.8), color = "black"),
            legend.key = element_blank()
          ) +
          guides(fill = guide_legend(override.aes = list(shape = c(21, 22))),
                 color = guide_legend(override.aes = list(color = c("#72B173", "#492050"))),
                 shape = guide_legend(override.aes = list(fill = c("#72B173", "#492050"))))
      )
    }
    
      # Compare uncorrected vs carrier ratios
      merged_data <- merge(
        get_results(dep_carrier), get_results(dep_not_cor),
        by = c("name", "ID")
      ) %>%
        select(name, ID, contains(contrast) & contains("ratio")) %>% filter(name %in% results_uncorrected$name)
      colnames(merged_data) <- c("name", "ID", "carrier", "uncorrected")
      
      if (parameters$filtering_type == "complete") {
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
          ggtitle(paste0(name_condition_b, " - ", heat_name, ": ", treat_condition, " vs ", base_condition, " QC comparison")) +
          theme_bw(base_size = 12)
      )
        }
      
      if (parameters$filtering_type == "condition") {
        #MAR
        merged_data$imputed <- ifelse(merged_data$name %in% proteins_imputed, "Yes", "No")
        print(
          ggplot(merged_data, aes(x = carrier, y = uncorrected, color = imputed, shape = imputed, fill  = imputed)) +
            geom_point(alpha = 0.5) +
            geom_point(
              data = merged_data[merged_data$name %in% proteins_to_highlight, ],
              size = 2
            ) +
            geom_text_repel(
              data = merged_data[merged_data$name %in% proteins_to_highlight, ],
              aes(label = name), max.overlaps = Inf, size = 3, color = "black"
            ) +
            scale_shape_manual(name = "Imputed", values = c("No" = 21, "Yes" = 22)) +
            scale_fill_manual(name = "Imputed", values = c("No" = "#72B173", "Yes" = "#492050")) +
            scale_color_manual(name = "Imputed", values = c("No" = "#72B173", "Yes" = "#492050")) +
            xlab("log2FC carrier") +
            ylab("log2FC STPP-UP (uncorrected)") +
            ggtitle(paste0(name_condition_b, " - ", heat_name, ": ", treat_condition, " vs ", base_condition, " QC comparison - Imputed MAR only")) +
            theme_bw(base_size = 12) +
            theme(
              legend.position = c(0.90, 0.10),   # (x, y) coordinates in relative units (0–1)
              legend.background = element_rect(fill = alpha("white", 0.8), color = "black"),
              legend.key = element_blank()
            ) +
            guides(fill = guide_legend(override.aes = list(shape = c(21, 22))),
                   color = guide_legend(override.aes = list(color = c("#72B173", "#492050"))),
                   shape = guide_legend(override.aes = list(fill = c("#72B173", "#492050"))))
        )
        
        #MAR + MNAR
        merged_data <- merge(
          get_results(dep_carrier), get_results(dep_not_cor),
          by = c("name", "ID")
        ) %>%
          select(name, ID, contains(contrast) & contains("ratio"))
        colnames(merged_data) <- c("name", "ID", "carrier", "uncorrected")
        merged_data$imputed <- ifelse(merged_data$name %in% proteins_imputed, "Yes", "No")
        print(
          ggplot(merged_data, aes(x = carrier, y = uncorrected, color = imputed, shape = imputed, fill  = imputed)) +
            geom_point(alpha = 0.5) +
            geom_text_repel(
              data = merged_data[merged_data$name %in% proteins_to_highlight, ],
              aes(label = name), max.overlaps = Inf, size = 3, color = "black"
            ) +
            scale_shape_manual(name = "Imputed", values = c("No" = 21, "Yes" = 22)) +
            scale_fill_manual(name = "Imputed", values = c("No" = "#72B173", "Yes" = "#492050")) +
            scale_color_manual(name = "Imputed", values = c("No" = "#72B173", "Yes" = "#492050")) +
            xlab("log2FC carrier") +
            ylab("log2FC STPP-UP (uncorrected)") +
            ggtitle(paste0(heat_name, ": ", treat_condition, " vs ", base_condition, " QC comparison - Imputed MAR & MNAR")) +
            theme_bw(base_size = 12) +
            theme(
              legend.position = c(0.90, 0.10),   # (x, y) coordinates in relative units (0–1)
              legend.background = element_rect(fill = alpha("white", 0.8), color = "black"),
              legend.key = element_blank()
            ) +
            guides(fill = guide_legend(override.aes = list(shape = c(21, 22))),
                   color = guide_legend(override.aes = list(color = c("#72B173", "#492050"))),
                   shape = guide_legend(override.aes = list(fill = c("#72B173", "#492050"))))
        )
        }
      
      dev.off()
      
      # Set up Excel file
      addWorksheet(wb, paste0(heat_name, "_uncorrected"))
      if (parameters$filtering_type == "complete") {
        addWorksheet(wb, paste0(heat_name, "_corrected"))
        }
      writeData(wb, paste0(heat_name, "_uncorrected"), results_list[[paste0(name_condition_b, "_", heat_name, "_", contrast, "_uncorrected")]])
      if (parameters$filtering_type == "complete") {
        writeData(wb, paste0(heat_name, "_corrected"), results_list[[paste0(name_condition_b, "_", heat_name, "_", contrast, "_corrected")]])
      }
      
    }
    
    
    # ------------------------------------------------------------
    # Write Excel file per contrast
    # ------------------------------------------------------------
    
    saveWorkbook(wb, paste0("DEP_results_", name_condition_b, "_", contrast, ".xlsx"), overwrite = TRUE)
    
  }
  
}


# ------------------------------------------------------------
# Save workspace and session info
# ------------------------------------------------------------

save.image(file.path(getwd(), "/STPP-UP analysis workspace.RData"))
sessionInfo()