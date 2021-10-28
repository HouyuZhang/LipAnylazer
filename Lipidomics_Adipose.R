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
source("../210923_Lipidomics/LipidomicsAnlysis.R")
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

mummichog_Functional_analysis(pktablePath = "standardized_Measurements_RT_TN_mummichogInput_tn_renal_ma+rt_renal_ma.csv", 
                              PvalueThreshold = 0.5, Species = "mmu", rowNormMet = "NULL", EnrichType = "KEGG",
                              SamplesDel = "rt_renal_ma_3")

Refine_MZ_Identifications(MeasurementsFile_standardized_matched = "standardized_Measurements_RT_TN_matched.csv",
                          RefinePlan = "PlanC")

Analyze_Lipids(MeasurementsFile_standardized_matched_refined = "standardized_Measurements_RT_TN_matched_refined_PlanC.csv",
               OutputPairWise = T, KeywordSelected = "rt", run_RunMetaboAnalystR = T)



