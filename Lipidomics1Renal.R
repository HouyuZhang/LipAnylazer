setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(janitor)

# ==============================================================================
# 0. Install necessary packages
# ==============================================================================
metanr_packages <- function(){
  metr_pkgs <- c("impute", "pcaMethods", "globaltest", "GlobalAncova", "Rgraphviz", 
                 "preprocessCore", "genefilter", "SSPA", "sva", "limma", 
                 "KEGGgraph", "siggenes","BiocParallel", "MSnbase", "multtest",
                 "RBGL","edgeR","fgsea","devtools","crmn")
  list_installed <- installed.packages()
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  if(length(new_pkgs)!=0){if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install(new_pkgs)
    print(c(new_pkgs, " packages added..."))
  }
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}
metanr_packages()
devtools::install_github("xia-lab/OptiLCMS", build = TRUE, build_vignettes = F)
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = F)

# ==============================================================================
# 1. Tidy data and run functions
# ==============================================================================
source("D:/OneDrive/Github_respiratories/Lipidomics_analysis/LipidomicsAnlysis.R")
MatchLipidMAPS(MeasurementsFile_standardized = "Lipidomics_5depots_rawdata.csv", 
               LipidMAPSDB = "lipid_database.csv",
               IdentificationFile = "Lipidomics_5depots_Identificaitons.csv", 
               OutUnmatchedLipids = T,
               MeasurementsFile_standardized_matched = "Lipidomics_5depots_rawdata_matched.csv"
)

Format_mummichog_input(MeasurementsFile_standardized = "Lipidomics_5depots_rawdata.csv", 
                       OutputPairWise = T, run_mummichog_Functional_analysis = T,
                       SelectedSample = "")

Refine_MZ_Identifications(MeasurementsFile_standardized_matched = "Lipidomics_5depots_rawdata_matched.csv",
                          RefinePlan = "PlanC")

Analyze_Lipids(MeasurementsFile_standardized_matched_refined = "Lipidomics_5depots_rawdata_matched_refined_PlanC.csv",
               OutputPairWise = T, KeywordSelected = "adi", run_RunMetaboAnalystR = T)



setwd("../From_Detain/")
source("../210923_Lipidomics/LipidomicsAnlysis.R")
MatchLipidMAPS(MeasurementsFile_standardized = "standardized_Measurements_RT_TN.csv", 
               LipidMAPSDB = "../210923_Lipidomics/lipid_database.csv",
               IdentificationFile = "CompoundIdentifications.csv", 
               OutUnmatchedLipids = T,
               MeasurementsFile_standardized_matched = "standardized_Measurements_RT_TN_matched.csv"
)

Format_mummichog_input(MeasurementsFile_standardized = "standardized_Measurements_RT_TN.csv", 
                       OutputPairWise = T, SelectedSample = "")

mummichog_Functional_analysis(pktablePath = "functional_hfdcd_ing_ma.csv", 
                              PvalueThreshold = 1, Species = "mmu", rowNormMet = "NULL", EnrichType = "MainClass",
                              SamplesDel = "hfd_ing_ma_3")

mummichog_Functional_analysis(pktablePath = "functional_hfdcd_ing_ma.csv",
                              PvalueThreshold = 0.15,
                              Species = "mmu",
                              rowNormMet = "NULL",
                              EnrichType = "KEGG",
                              SamplesDel = "hfd_ing_ma_3")

Refine_MZ_Identifications(MeasurementsFile_standardized_matched = "standardized_Measurements_RT_TN_matched.csv",
                          RefinePlan = "PlanC")

Analyze_Lipids(MeasurementsFile_standardized_matched_refined = "standardized_Measurements_RT_TN_matched_refined_PlanC.csv",
               OutputPairWise = T, KeywordSelected = "rt", run_RunMetaboAnalystR = T)

# ==============================================================================
# 1. Tidy data and run functions
# ==============================================================================
Level = c("main_class","sub_class","abbrev")[3]

  read_csv("Lipidomics_5depots_rawdata_matched_refined_PlanC.csv") %>%
  select(c(compound:adi_epi_3),starts_with("adi_peri")) %>%
  mutate_at(vars(starts_with("adi")), funs(.*100000/sum(.))) %>%
  filter(str_detect(main_class,"Ceramide")) %>%
  filter(str_detect(sub_class,"acylsphingosines")) %>%
  select_at(vars(Level, starts_with("adi"))) %>%
  group_by(!!as.name(Level)) %>% 
  mutate_if(is.character,as.numeric) %>% 
  summarise_all(c("sum")) %>%
  filter(!!as.name(Level) != "-") %>% 
  filter(rowSums(across(where(is.numeric))) != 0) %>%
  reshape2::melt(c(Level)) %>% 
  mutate(group = gsub("_[0-9]","",variable)) %>%
  ggplot(aes(x = group, y = value, fill = group)) + 
    geom_boxplot(notch=F, size=0.4) +
    geom_jitter(width = 0.1, alpha = 0.5,color = "black") + 
    ggpubr::stat_compare_means(method = "t.test", paired = F, 
                               comparisons = list(c("adi_epi", "adi_peri"))) +
    scale_fill_manual(values = c("#66A61E","#D95F02")) +
    facet_wrap(~abbrev, ncol =6, scales = "free") +
    labs(x = "", y = "Relative lipid signal") +
    theme_bw() +
    theme(
      axis.title.x = element_text(color="black", size=14, face="bold"),
      # axis.text.x = element_text(color="black", size=11, face="bold",angle = 40, hjust = 1),
      axis.text.x = element_blank(),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.text.y = element_text(color="black", size=11, face="bold"),
      legend.title = element_blank(),
      legend.text = element_text(color="black", size=14, face="bold"),
      
    )
ggsave("abbrev_22.pdf",width = 10, height = 6)


