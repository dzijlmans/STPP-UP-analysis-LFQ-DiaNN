#Starting up
library(tidyverse)
library(DEP)
library(openxlsx)
library(ggrepel)
library(tcltk)
library(PTXQC)
tkmessageBox(message="Choose your working directory")
setwd(tk_choose.dir())
#set parameters
filtering_type <- "complete"                                                     #choose from "complete" or "condition"
##load and clean data
tkmessageBox(message="Choose your DiaNN output file (*pg_matrix.tsv)")
proteinGroups <- read.delim(tk_choose.files())                                   #select pg_matrix DiaNN output file

tkmessageBox(message="Please ensure that the abundance column names (and thus the .raw files) contain the name of the mass spectrometer used (e.g. Astral, Exploris) and the following information, in order and separated by an underscore:

- culture conditions or cell lines used 
- heating conditions
- compounds used
- replicate
             
Example: 20250205_Astral_AUR_DZ114_Ad_Carrier_C1_1.raw")

proteinGroups <-
  make_unique(proteinGroups, "Genes", "Protein.Group", delim = ";")
proteinGroups <- proteinGroups %>% filter(!grepl("Cont_", ID) &
                                            !grepl("KRT", name))

###tidy column names
value_columns <- which(str_detect(names(proteinGroups), "Astral"))               #specify Mass spec (this should be in the value column names)
prefix <- LCSn(colnames(proteinGroups[, value_columns]))
colnames(proteinGroups) <- gsub(".raw", "", gsub(".*" %>% paste0(prefix), "", colnames(proteinGroups)))
proteinGroups <-
  proteinGroups %>% pivot_longer(
    cols = all_of(value_columns),
    names_to = c("culture", "heating", "compound", "replicate"),                 #change this depending on the amount of different variables used in your experiment
    names_sep = "_",
    values_to = "values"
  )

proteinGroups <-
  proteinGroups %>% pivot_wider(names_from = c("compound", "replicate"),         #change this if necessary
                                values_from = "values")

###split data according to desired variable
datalist_var1 <- split(proteinGroups, f = ~ proteinGroups$culture)               #specify variable to be split on
results_var1 <- list()
results_var2 <- list()

#run DEP analysis
for (i in 1:length(datalist_var1)) {
  #for each element of the list separate carrier and heated options first
  name_var1 <- names(datalist_var1[i])
  dataset <- datalist_var1[[i]]
  datalist_var2 <- split(dataset, f = ~ proteinGroups$heating)                   #specify variable to be split on
  datalist_carrier <- datalist_var2[names(datalist_var2) %in% c("Carrier")]
  datalist_other <- datalist_var2[!names(datalist_var2) %in% c("Carrier")]
  
  #run DEP analysis on carrier samples
  for (j in 1:length(datalist_carrier)) {
    data <- datalist_carrier[[j]]
    name_var2 <- names(datalist_carrier[j])
    value_columns <- which(sapply(data, is.numeric))
    
    experimental_design <- data.frame(label = colnames(data[, value_columns]))
    experimental_design <- experimental_design %>% separate_wider_delim(
      label,
      delim = "_",
      names = c("condition", "replicate"),
      cols_remove = FALSE
    )
    experimental_design <- experimental_design[, c("label", "condition", "replicate")]
    
    #Make object and do analysis
    data_se <- make_se(data, value_columns, experimental_design)
    data_filt <- filter_proteins(data_se, type = filtering_type, thr = 0)
    data_norm <- normalize_vsn(data_filt)
    if (filtering_type == "complete") {
      data_imp <- data_norm
    }
    if (filtering_type == "condition") {
      data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
    }
    data_diff_manual <- test_diff(data_imp, type = "control", control = "DMSO")   #specify control
    contrast <- sub("_diff", "", str_subset(
      names(data_diff_manual@elementMetadata@listData),
      "_diff"
    ))
    dep <- add_rejections(data_diff_manual, alpha = 0.05, lfc = log2(2))    
    results_var2[[paste0(name_var2, "_", contrast)]] <- get_results(dep)
    
    ####make pdf containing QCs
    pdf(paste0(name_var1, "_", name_var2, ".pdf"))
    
    print(plot_frequency(data_se))
    print(plot_numbers(data_filt))
    if (filtering_type == "condition") {
      print(plot_detect(data_filt))
    }
    print(plot_normalization(data_filt, data_norm, data_imp))
    print(plot_imputation(data_norm, data_imp))
    print(plot_volcano(
      dep,
      contrast = contrast,
      label_size = 2,
      add_names = TRUE
    ))
    dev.off()
    
  }
  
  #next, perform DEP analysis on heated samples, with and without normalization
  for (k in 1:length(datalist_other)) {
    data <- datalist_other[[k]]
    name_var2 <- names(datalist_other[k])
    value_columns <- which(sapply(data, is.numeric))
    
    experimental_design <- data.frame(label = colnames(data[, value_columns]))
    experimental_design <- experimental_design %>% separate_wider_delim(
      label,
      delim = "_",
      names = c("condition", "replicate"),
      cols_remove = FALSE
    )
    experimental_design <- experimental_design[, c("label", "condition", "replicate")]
    
    #Make object and do analysis
    data_se <- make_se(data, value_columns, experimental_design)
    data_filt <- filter_proteins(data_se, type = filtering_type, thr = 0)
    data_norm <- normalize_vsn(data_filt)
    if (filtering_type == "complete") {
      data_imp <- data_norm
    }
    if (filtering_type == "condition") {
      data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
    }

    #do carrier correction
    data_imp_not_corrected <- data_imp
    data_imp_corrected <- data_imp
    data_corr <- data.frame(data_imp_corrected@assays@data@listData[[1]])
    data_corr$name <- rownames(data_imp_corrected@assays@data@listData[[1]])
    data_corr <- merge(data_corr, dplyr::select(results_var2[[grep("Carrier", names(results_var2))]], c("name", "C1_vs_DMSO_ratio")), by = "name", all.x = TRUE)
    data_corr[is.na(data_corr)] <- 0
    data_corr <- data_corr[order(match(data_corr$name, rownames(data_imp_corrected@assays@data@listData[[1]]))),]
    data_corr[, c(9:11)] <- data_corr[, c(2:4)] + data_corr$C1_vs_DMSO_ratio
    data_imp_corrected@assays@data@listData[[1]][, c(1:3)] <- as.matrix(data_corr[, c(9:11)])
    
    #no carrier correction
    data_diff_manual_not_cor <- 
      test_diff(data_imp_not_corrected, type = "control", control = "DMSO")       #specify control
    contrast <- sub("_diff", "", str_subset(
      names(data_diff_manual_not_cor@elementMetadata@listData),
      "_diff"
    ))
    dep_not_cor <- add_rejections(data_diff_manual_not_cor, alpha = 0.05, lfc = log2(2))

    results_var2[[paste0(name_var2, "_", contrast, "_uncorrected")]] <- get_results(dep_not_cor)
    
    #carrier correction
    data_diff_manual_cor <- test_diff(data_imp_corrected, type = "control", control = "DMSO")
    contrast <- sub("_diff", "", str_subset(
      names(data_diff_manual_cor@elementMetadata@listData),
      "_diff"
    ))
    dep_cor <- add_rejections(data_diff_manual_cor, alpha = 0.05, lfc = log2(2))
    results_var2[[paste0(name_var2, "_", contrast, "_corrected")]] <- get_results(dep_cor)
    
    ####make pdf containing QCs
    pdf(paste0(name_var1, "_", name_var2, ".pdf"))
    
    print(plot_frequency(data_se))
    print(plot_numbers(data_filt))
    if (filtering_type == "condition") {
      print(plot_detect(data_filt))
    }
    print(plot_normalization(data_filt, data_norm, data_imp))
    print(plot_imputation(data_norm, data_imp))
    print(plot_volcano(
      dep_not_cor,
      contrast = contrast,
      label_size = 2,
      add_names = TRUE
    ))
    print(plot_volcano(
      dep_cor,
      contrast = contrast,
      label_size = 2,
      add_names = TRUE
    ))
    dev.off()
  }
  
  #QCs
  pdf(paste0(name_var1, "_QC.pdf"))
  merged_data <- merge(get_results(dep), get_results(dep_not_cor), by = c("name", "ID")) %>% select(name, ID, contains("ratio"))
  colnames(merged_data) <- c("name", "ID", "carrier", "corrected")
  print(ggplot(merged_data,
         aes(carrier, corrected)) +
    geom_point(color = "grey50") +
    geom_point(data = merged_data[merged_data$name%in%c("DHFR", "TYMS"), ], color = "#E41A1C") +
    geom_text_repel(data = merged_data[merged_data$name%in%c("DHFR", "TYMS"), ], aes(label = name)) +
    xlab("log2FC carrier") +
    ylab("log2FC STPP-UP") +
    theme_bw(base_size = 16))
  dev.off()
  write.xlsx(results_var2, paste0("DEP_results_", name_var1, ".xlsx"))
}

save.image(paste0(getwd(), "/DiaNN - DEP analysis.RData"))



## to do
#add option for multiple contrasts
