setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("D:/OneDrive - pku.edu.cn/0.5_Github_respiratories/Lipidomics_analysis/LipidomicsAnlysis.R")

# ==============================================================================
# 1. iWAT differentiation
# ==============================================================================
setwd("D:/OneDrive/4_PKU-CIBR/1_Rotation_CarlosLab/1_Lipidomics2Lab/2_iWAT_differentiation/")
Format_mummichog_input(MeasurementsFile_standardized = "iWATdiff_standardized_Measurements.csv", OutputPairWise = T)

mummichog_Functional_analysis(pktablePath = "iWATdiff_standardized_Measurements_mummichogInput_ing_cell+ing_ma.csv", 
                              Species = "mmu", rowNormMet = "NULL",
                              PvalueThreshold = 0.05, EnrichType = c("KEGG","MainClass")
)
mummichog_Functional_analysis(pktablePath = "iWATdiff_standardized_Measurements_mummichogInput_ing_svf+ing_cell.csv", 
                               Species = "mmu", rowNormMet = "NULL",
                               PvalueThreshold = 0.05, EnrichType = c("KEGG","MainClass")
)
mummichog_Functional_analysis(pktablePath = "iWATdiff_standardized_Measurements_mummichogInput_ing_svf+ing_ma.csv", 
                               Species = "mmu", rowNormMet = "NULL",
                               PvalueThreshold = 0.05, EnrichType = c("KEGG","MainClass")
)
mummichog_Functional_analysis(pktablePath = "iWATdiff_standardized_Measurements_mummichogInput_ing.csv", 
                              Species = "mmu", rowNormMet = "NULL",
                              PvalueThreshold = 0.05, EnrichType = c("KEGG","MainClass")
                              )
Replot_mummichog(mummichog_PEpath = "./")

# ==============================================================================
# 2_iBAT_differentiation_noIS
# ==============================================================================
setwd("D:/OneDrive/4_PKU-CIBR/1_Rotation_CarlosLab/1_Lipidomics2Lab/2_iBAT_differentiation_noIS/")
Format_mummichog_input(MeasurementsFile_standardized = "iBATdiff_Measurements.csv", OutputPairWise = T)

mummichog_Functional_analysis(pktablePath = "iBATdiff_Measurements_mummichogInput_ibat_cell+ibat_ma.csv", 
                              Species = "mmu", rowNormMet = "SumNorm",
                              PvalueThreshold = 0.05, EnrichType = c("KEGG","MainClass")
)
mummichog_Functional_analysis(pktablePath = "iBATdiff_Measurements_mummichogInput_ibat_svf+ibat_cell.csv", 
                              Species = "mmu", rowNormMet = "SumNorm",
                              PvalueThreshold = 0.05, EnrichType = c("KEGG","MainClass")
)
mummichog_Functional_analysis(pktablePath = "iBATdiff_Measurements_mummichogInput_ibat_svf+ibat_ma.csv", 
                              Species = "mmu", rowNormMet = "SumNorm",
                              PvalueThreshold = 0.05, EnrichType = c("KEGG","MainClass")
)
mummichog_Functional_analysis(pktablePath = "iBATdiff_Measurements_mummichogInput.csv", 
                              Species = "mmu", rowNormMet = "SumNorm",
                              PvalueThreshold = 0.05, EnrichType = c("KEGG","MainClass")
)
setwd("D:/OneDrive/4_PKU-CIBR/1_Rotation_CarlosLab/1_Lipidomics2Lab/")
Replot_mummichog(mummichog_PEpath = "./")

# ==============================================================================
# 3. OB/Lean moneky
# ==============================================================================
setwd("D:/OneDrive/4_PKU-CIBR/1_Rotation_CarlosLab/1_Lipidomics2Lab/3_OBLean_Monkey/")
Format_mummichog_input(MeasurementsFile_standardized = "OBLean_Monkey_standardized_Measurements.csv", OutputPairWise = F)

mummichog_Functional_analysis(pktablePath = "OBLean_Monkey_standardized_Measurements_mummichogInput.csv", 
                              Species = "hsa", rowNormMet = "NULL",
                              PvalueThreshold = 0.4, EnrichType = c("KEGG","MainClass")
)
Replot_mummichog(mummichog_PEpath = "./")

# ==============================================================================
# 4. CE/RT moneky
# ==============================================================================
setwd("D:/OneDrive/4_PKU-CIBR/1_Rotation_CarlosLab/1_Lipidomics2Lab/1_CE_RT_noIS/")
Format_mummichog_input(MeasurementsFile_standardized = "CE_RT_Measurements2.csv", OutputPairWise = T)

mummichog_Functional_analysis(pktablePath = "CE_RT_Measurements2_mummichogInput_ma_ce_ibat_exp1+ma_rt_ibat_exp1.csv", 
                              Species = "mmu", rowNormMet = "SumNorm", PvalueThreshold = 0.05, EnrichType = c("KEGG","MainClass")
)
mummichog_Functional_analysis(pktablePath = "CE_RT_Measurements2_mummichogInput_svf_ce_ibat_exp1+svf_rt_ibat_exp1.csv", 
                               Species = "mmu", rowNormMet = "SumNorm", PvalueThreshold = 1, EnrichType = c("KEGG","MainClass")
)

mummichog_Functional_analysis(pktablePath = "CE_RT_Measurements2_mummichogInput_ma_ce_ibat_exp2+ma_rt_ibat_exp2.csv", 
                              Species = "mmu", rowNormMet = "SumNorm", PvalueThreshold = 0.4, EnrichType = c("KEGG","MainClass")
)
mummichog_Functional_analysis(pktablePath = "CE_RT_Measurements2_mummichogInput_svf_ce_ibat_exp2+svf_rt_ibat_exp2.csv", 
                              Species = "mmu", rowNormMet = "SumNorm", PvalueThreshold = 0.1, EnrichType = c("KEGG","MainClass")
)


Replot_mummichog(mummichog_PEpath = "../")

