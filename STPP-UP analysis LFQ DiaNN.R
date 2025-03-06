#Starting up
library(tidyverse)
library(DEP)
library(openxlsx)
library(ggrepel)
library(tcltk)
library(PTXQC)
setwd(tk_choose.dir())
#set parameters
filtering_type <- "complete"                                                     #choose from "complete" or "condition"
##load and clean data
proteinGroups <- read.delim(tk_choose.files())                                   #select pg_matrix DiaNN output file
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
    print(plot_detect(data_filt))
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
    print(plot_detect(data_filt))
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





## to do
#validate that carrier correction is done correctly
#add option for multiple contrasts
#add ifelse statements for type of filtering (complete or condition)











proteinGroups %>% drop_na()








openxlsx::write.xlsx(proteinGroups_list, "raw_dataset_filtered.xlsx")

proteinGroups_list_carrier <- proteinGroups_list["Carrier"]
proteinGroups_list_STPPUP <- proteinGroups_list["STPPUP"]


#preparation for DEP analysis
experimental_design <- data.frame(label = unique(gsub("^[^_]*_", "", label_names)))
experimental_design <- experimental_design %>% separate_wider_delim(label, delim = "_", names = c("condition", "replicate"), cols_remove = FALSE)
experimental_design <- experimental_design[, c("label", "condition", "replicate")]

#DEP analysis - Carrier
DEP_results_carrier <- list()
for (i in 1:length(proteinGroups_list_carrier)) {
  data <- datalist[[1]]
    proteinGroups_list_carrier[[i]] %>% dplyr::select(-c(experiment)) #drop column to match input example DEP
  value_columns <- which(colnames(data)%in%experimental_design$label) # set value column numbers
  data_se <- make_se_parse(data, value_columns)
  data_se <- filter_proteins(data_se, "complete")
  #normalization
  data_norm <- normalize_vsn(data_se)
  plot_normalization(data_se, data_norm)
  ggsave(paste0("QC_normalization.", names(proteinGroups_list_carrier[i]), ".pdf"))
  
  data_diff <-
    test_diff(data_norm, type = "control", control = "DMSO")
  dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(0.5))
  plot_volcano(dep, contrast = "C1_vs_DMSO", label_size = 2, add_names = TRUE)
  ggsave(paste0("DEP.volcano.plot.", names(proteinGroups_list_carrier[i]), ".pdf" ))
  
  DEP_results_carrier[[i]] <- get_results(dep)
  names(DEP_results_carrier)[i] <- paste0(names(proteinGroups_list_carrier[i]), "_results")
}

write.xlsx(DEP_results_carrier, "DEP_results_carrier.xlsx")



#DEP analysis - STPPUP - uncorrected
DEP_results_STPPUP_uncorrected <- list()
for (i in 1:length(proteinGroups_list_STPPUP)) {
  data <-
    proteinGroups_list_STPPUP[[i]] %>% dplyr::select(-c(experiment)) #drop column to match input example DEP
  value_columns <- which(colnames(data)%in%experimental_design$label) # set value column numbers
  data_se <- make_se(data, value_columns, experimental_design)
  data_se <- filter_proteins(data_se, "complete")
    #normalization
  data_norm <- normalize_vsn(data_se)
  plot_normalization(data_se, data_norm)
  ggsave(paste0("QC_normalization.", names(proteinGroups_list_STPPUP[i]), ".pdf"))
  
  data_diff <-
    test_diff(data_norm, type = "control", control = "DMSO")
  dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(0.5))
  plot_volcano(dep, contrast = "C1_vs_DMSO", label_size = 2, add_names = TRUE)
  ggsave(paste0("DEP.volcano.plot.", names(proteinGroups_list_STPPUP[i]), ".uncorrected.pdf"))
  
  
  data_results <- get_results(dep)
  DEP_results_STPPUP_uncorrected[[i]] <- data_results
  names(DEP_results_STPPUP_uncorrected)[i] <- paste0(names(proteinGroups_list_STPPUP[i]), "_results")
}

write.xlsx(DEP_results_STPPUP_uncorrected, "DEP_results_STPPUP_uncorrected.xlsx")


#DEP analysis - STPPUP - corrected
#change data_corr column numbers according to replicate numbers
DEP_results_STPPUP_corrected <- list()
for (i in 1:length(proteinGroups_list_STPPUP)) {
  data <-
    proteinGroups_list_STPPUP[[i]] %>% dplyr::select(-c(experiment)) #drop column to match input example DEP
  value_columns <- which(colnames(data)%in%experimental_design$label) # set value column numbers
  data_se <- make_se(data, value_columns, experimental_design)
  data_se <- filter_proteins(data_se, "complete")
  #normalization
  data_norm <- normalize_vsn(data_se)
  plot_normalization(data_se, data_norm)
  ggsave(paste0("QC_normalization.", names(proteinGroups_list_STPPUP[i]), ".pdf"))
  
  #adjust for carrier
  data_corr <- data.frame(data_norm@assays@data@listData[[1]])
  data_corr$name <- rownames(data_norm@assays@data@listData[[1]])
  data_corr <- merge(data_corr, dplyr::select(DEP_results_carrier[[i]], c("name", "C1_vs_DMSO_ratio")), by = "name")
  data_corr <- data_corr[order(match(data_corr$name, rownames(data_norm@assays@data@listData[[1]]))),]
  data_corr[, c(9:11)] <- data_corr[, c(2:4)] + data_corr$C1_vs_DMSO_ratio
  data_norm@assays@data@listData[[1]][, c(1:3)] <- as.matrix(data_corr[, c(9:11)])
  
  data_diff <-
    test_diff(data_norm, type = "control", control = "DMSO")
  dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(0.5))
  plot_volcano(dep, contrast = "C1_vs_DMSO", label_size = 2, add_names = TRUE)
  ggsave(paste0("DEP.volcano.plot.", names(proteinGroups_list_STPPUP[i]), ".corrected.pdf"))
  
  
  data_results <- get_results(dep)
  DEP_results_STPPUP_corrected[[i]] <- data_results
  names(DEP_results_STPPUP_corrected)[i] <- paste0(names(proteinGroups_list_STPPUP[i]), "_results")
}

write.xlsx(DEP_results_STPPUP_corrected, "DEP_results_STPPUP_corrected.xlsx")


#QCs and Volcano plots
for (j in 1:length(proteinGroups_list_STPPUP)) {
  carrier <- DEP_results_carrier[[j]]
  uncorrected <- DEP_results_STPPUP_uncorrected[[j]]
  corrected <- DEP_results_STPPUP_corrected[[j]]
  
  #Set proteins to highlight, either by name or enrichment
  #proteins_to_highlight <- corrected_data$name[-log10(corrected_data$C1_vs_DMSO_p.val) > 3 & abs(corrected_data$C1_vs_DMSO_ratio) > 0.3]
  proteins_to_highlight <- c("DHFR", "TYMS")
  
  #QC - carrier vs uncorrected
  merged_data <- merge(carrier, uncorrected, by = c("name", "ID"))
    ggplot(merged_data,
         aes(C1_vs_DMSO_ratio.x, C1_vs_DMSO_ratio.y)) +
    geom_point(color = "grey50") +
    geom_point(data = merged_data[merged_data$name%in%c("DHFR", "TYMS"), ], color = "#E41A1C") +
   # geom_text_repel(data = merged_data[merged_data$name%in%proteins_to_highlight, ], aes(label = name)) +
    xlab("log2FC carrier") +
    ylab("log2FC STPP-UP") +
    theme_bw(base_size = 16)
  ggsave(paste0("QC_carrier_vs_uncorrected_", names(proteinGroups_list_STPPUP[j]), ".pdf"), height = 5.21, width = 6.92)
  
  #Volcano carrier
  ggplot(carrier, aes(C1_vs_DMSO_ratio, -log10(C1_vs_DMSO_p.val))) +
    geom_point(alpha = 0.1) +
    geom_point(data = carrier[carrier$name%in%proteins_to_highlight, ], color = "#E41A1C") +
    geom_text_repel(data = carrier[carrier$name%in%proteins_to_highlight, ], aes(label = name)) +
    theme_bw()
  ggsave("Volcano_carrier.pdf")
  
  #Volcano uncorrected
  ggplot(uncorrected, aes(C1_vs_DMSO_ratio, -log10(C1_vs_DMSO_p.val))) +
    geom_point(alpha = 0.1) +
    geom_point(data = uncorrected[uncorrected$name%in%proteins_to_highlight, ], color = "#E41A1C") +
    geom_text_repel(data = uncorrected[uncorrected$name%in%proteins_to_highlight, ], aes(label = name)) +
    theme_bw()
  ggsave("Volcano_uncorrected.pdf")
  
  #Volcano corrected
  ggplot(corrected, aes(C1_vs_DMSO_ratio, -log10(C1_vs_DMSO_p.val))) +
    geom_point(alpha = 0.1) +
    geom_point(data = corrected[corrected$name%in%proteins_to_highlight, ], color = "#E41A1C") +
    geom_text_repel(data = corrected[corrected$name%in%proteins_to_highlight, ], aes(label = name)) +
    theme_bw()
  ggsave("Volcano_corrected.pdf")
}


save.image(paste0(getwd(), "/DEP analysis.RData"))
