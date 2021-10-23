#=========================================================================================
# This script is for standardizing the identified Lipids concentration against 
# a list of Internal Standards, with specially taken the replicates consistency into consideration

#Version 1.2 created by Houyu Zhang on 2021/10/23
#Issue report on Hughiez047@gmail.com
#Copyright (c) 2021 __CarlosLab@PKU__. All rights reserved.
#=========================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

##' Standardize samples against internal controls within defined groups,
##' and map compound value (rt_m/z) against LipidMAPS database
##' 
##' @param MeasurementsFile LPMS output file containing compound signal across sample
##' @param ISMeasurementsFile LPMS output file containing counts for internal standards
##' @param IdentificationFile LPMS output file containing possible identification for this LPMS experiment 
##' @param LipidMAPSDB A lipid information database file downloaded from LIPIDMAPS website
##' @param Output_MeasurementsFile Specify a output file name for standardized-converted lipid signals
##' @param GroupLabel Defined groups for the excel sheet which will be standardized together
##' @examples
##'
##' standardize_by_standards(MeasurementsFile = "test1.csv", 
##'                         ISMeasurementsFile = "IS-Measurements-1.xls",
##'                         IdentificationFile = "Identifications-1.csv",
##'                         LipidMAPSDB = "lipid_database.csv",
##'                         Output_MeasurementsFile = "standardized_test1_Converted_zhy.csv",
##'                         GroupLabel = list(c(1),c(2:4),c(5:7),c(8:10),c(11:13),c(14:16),
##'                                                c(17:19),c(20:22),c(23:25),c(26:28),c(29:31),c(32:34),
##'                                                c(35:37),c(38:42),c(43:47),c(48:50),c(51:53)))
##'

standardize_by_standards <- function(MeasurementsFile = "", 
                                     ISMeasurementsFile = "", 
                                     IdentificationFile = "",
                                     LipidMAPSDB = "", 
                                     Output_MeasurementsFile = "", 
                                     GroupLabel = list()){
  
  MF_original <- read_csv(MeasurementsFile, skip = 2, show_col_types = F) %>% clean_names()
  #Count total sample numbers through GroupLabel
  SampleNumbers <- max(unlist(GroupLabel))
  #First 13 cols are MS meta-information
  SampleColStartFrom <- 13 + SampleNumbers 
  MF_standardized <- MF_original
  
  #Since Standards might differ across defined group, consistent Standards will be used
  ReplicatesStandards <- list()
  for (item in 1:length(GroupLabel)){
    tmpVC <- c()
    for (sampleNum in GroupLabel[[item]]){
      IS <- read_xls(ISMeasurementsFile, sheet = sampleNum) %>% 
        clean_names() %>% arrange(observed_rt_min)
      tmpVC <- c(tmpVC, IS$component_name)
    }
    tmpVC <- table(tmpVC) %>% as.data.frame()
    #Only use Standards shared in all samples
    ReplicatesStandards[[item]] <- as.vector(tmpVC[tmpVC$Freq == length(GroupLabel[[item]]),]$tmpVC)
  }
  
  for (sampleNum in 1:SampleNumbers){
    StandardsUsedHere <- ReplicatesStandards[[which(sapply(GroupLabel , `%in%`, x = sampleNum))]]
    cat("--->>> Start process sample",sampleNum,"...\n")
    IS <- read_xls(ISMeasurementsFile, sheet = sampleNum) %>% 
      clean_names() %>% arrange(observed_rt_min)
    IS <- IS[IS$component_name %in% StandardsUsedHere,]
    
    cat("Refine RT intervals...\n")
    IS$start <- NA; IS$end <- NA
    for (i in 1:(nrow(IS)-1)){
      RT1 = as.numeric(IS[i,"observed_rt_min"])
      RT2 = as.numeric(IS[i+1,"observed_rt_min"])
      HalfIntervalLen <- (RT2 - RT1)/2
      IS[i,"end"] = RT1 + HalfIntervalLen
      IS[i+1,"start"] = RT2 - HalfIntervalLen
    }
    IS[1,"start"] = 0; IS[nrow(IS),"end"] = 100
    
    cat("Standardize lipid concentration using corresponding internal control...\n")
    for (LipidNum in 1:nrow(MF_original)){
      rt = MF_original[LipidNum,]$retention_time_min
      index <- which(with(IS, start <= rt & end >= rt))
      MF_standardized[LipidNum,sampleNum + SampleColStartFrom] <- 
        as.numeric(MF_original[LipidNum,sampleNum + SampleColStartFrom])/as.numeric(IS[index,"detector_counts"])
    }
  }
  
  MF_standardized_Slected <- MF_standardized[,-c(14:(13+SampleNumbers))]
  identification <- read_csv(IdentificationFile, show_col_types = F) %>% clean_names()
  lipid_database <- read_csv(LipidMAPSDB, show_col_types = F) %>% clean_names()
  
  Converted_standardize_MeasurementsFile <- lipid_database %>% 
    select(compound_id, category, main_class, sub_class, abbrev) %>% 
    right_join(identification, by="compound_id") %>% 
    select(compound, category, main_class, sub_class, abbrev) %>% 
    right_join(MF_standardized_Slected, by="compound") %>% 
    filter(category != "NA") %>%
    select(-c(neutral_mass_da:minimum_cv_percent)) %>% 
    arrange(compound)
  
  write.csv(x = Converted_standardize_MeasurementsFile, file = Output_MeasurementsFile, row.names = F)

  cat("<<",MeasurementsFile,">> has been standardized using information in <<",ISMeasurementsFile, 
      ">> and wrote into <<",Output_MeasurementsFile,">>\n")
}

##' Convert Standard (STD) Lipid signal file to mummichog accepted format (m/z__rt)
##' 
##' @param STD_MeasurementsFile same as the Output_MeasurementsFile from standardize_by_standards()
##' @param OutputPairWise Besides output all samples, indicate whether output pairwise samples
##' @examples
##'
##' Format_mummichog_input(STD_MeasurementsFile = "standardized_Measurements1_Converted.csv", 
##'                      OutputPairWise = T)
##'                      
                     
Format_mummichog_input <- function(STD_MeasurementsFile = "", 
                                   OutputPairWise = TRUE){
  
  prefix <- tools::file_path_sans_ext(STD_MeasurementsFile)
  mummichog_MeasurementsFile <- read_csv(STD_MeasurementsFile, show_col_types = F) %>% clean_names() %>% 
    mutate(compound = gsub("[a-z]|/","",compound)) %>% 
    separate(compound, c("RT", "MZ"), "_") %>% 
    mutate(Label = paste0(MZ,"__",RT)) %>% 
    relocate(Label, .before = RT) %>% 
    select(-c(RT:abbrev)) %>%
    unique()
  
  labels <- gsub("_[0-9]+","",colnames(mummichog_MeasurementsFile))
  mummichog_MeasurementsFile <- rbind(labels, mummichog_MeasurementsFile)
  write.csv(file = paste0(prefix,"_mummichogInput.csv"), mummichog_MeasurementsFile, row.names = F)
  
 if(OutputPairWise){
   groups <- unique(gsub("_[0-9]+","",colnames(mummichog_MeasurementsFile))[-1])
   if (length(groups) > 1){
     pairwise <- combn(groups,2)
     for (i in 1:ncol(pairwise)){
       x = pairwise[1,i]
       y = pairwise[2,i]
       tmp <- mummichog_MeasurementsFile %>% select(Label,starts_with(x),starts_with(y))
       fileName <- paste0(prefix,"_mummichogInput_",x,"+",y,".csv")
       write.csv(file = fileName, tmp, row.names = F, quote = F)
     }
   }
  }
}

##' Do mummichog analysis on KEGG/MainClass/SubClass using m/z__rt information
##' 
##' @param ktablePath A pair-wise sample file output by Format_mummichog_input()
##' @param PvalueThreshold Threshold to selected lipids for enrichment analysis
##' @param EnrichType Define enrichment methods, this is multi-optional
##' @examples
##' 
##' mummichog_Functional_analysis(pktablePath = "standardized_Measurements1_Converted_mummichogInput_G1+G3.csv",
##'                              PvalueThreshold = 0.5, 
##'                              EnrichType = c("KEGG","MainClass","SubClass")[c(1,2)])

mummichog_Functional_analysis <- function(pktablePath = "", 
                                          PvalueThreshold = 0.2, 
                                          EnrichType = c("KEGG","MainClass","SubClass")[c(1,2)]){
  suppressPackageStartupMessages(library(MetaboAnalystR))
  prefix <- tools::file_path_sans_ext(pktablePath)
  
  FunAnalysis <- InitDataObjects("pktable", data.type = "stat", paired = FALSE);
  FunAnalysis <- Read.TextData(FunAnalysis, filePath = pktablePath, format = "colu", lbl.type = "disc")
  FunAnalysis <- SanityCheckData(FunAnalysis)
  FunAnalysis <- ReplaceMin(FunAnalysis)
  FunAnalysis <- FilterVariable(FunAnalysis, filter = "iqr", qcFilter = "F", rsd = 25)
  FunAnalysis <- PreparePrenormData(FunAnalysis)
  FunAnalysis <- Normalization(FunAnalysis, rowNorm = "SumNorm", transNorm = "NULL", scaleNorm = "NULL",
                               ratio = FALSE, ratioNum = 20)
  FunAnalysis <- Ttests.Anal(FunAnalysis, nonpar = F, threshp = 0.05, paired = FALSE, 
                             equal.var = T, pvalType = "raw", all_results = T)
  FunAnalysis <- Convert2Mummichog(FunAnalysis, rt = T) # mummichog_input: m.z, r.tm p.value, t.score
    
  FunAnalysis <- InitDataObjects("mass_table", anal.type = "mummichog", paired = FALSE)
  FunAnalysis <- SetPeakFormat(FunAnalysis, type = "mprt")
  # Ion Mode: Negative Mode; Mass Tolerance (ppm): 5.0
  FunAnalysis <- UpdateInstrumentParameters(FunAnalysis, instrumentOpt = 5, msModeOpt = "positive")
  dataValue <- as.Date(Sys.Date(), "%Y-%m-%d")
  FunAnalysis <- Read.PeakListData(FunAnalysis, filename = paste0("mummichog_input_",dataValue,".txt")) # change date
  FunAnalysis <- SanityCheckMummichogData(FunAnalysis) # Retention time tolerance
  FunAnalysis <- SetPeakEnrichMethod(FunAnalysis, algOpt = "mum", version = "v2")
  FunAnalysis <- SetMummichogPval(FunAnalysis, cutoff = PvalueThreshold)

  if ("KEGG" %in% EnrichType){
    #Do KEGG analysis (This need more significant lipid for enrichment)
    FunAnalysis <- PerformPSEA(FunAnalysis, lib = "mmu_kegg", libVersion = "current", minLib = 3, permNum = 50)
    FunAnalysis <- PlotPeaks2Paths(FunAnalysis, paste0(prefix,"_peaks_to_KEGG_"), "pdf", dpi = 100, width=NA)
  }
  if ("MainClass" %in% EnrichType){
    FunAnalysis <- PerformPSEA(FunAnalysis, lib = "main_lipid_class_mset", libVersion = "current", minLib = 3, permNum = 50)
    FunAnalysis <- PlotPeaks2Paths(FunAnalysis, paste0(prefix,"_peaks_to_MLC_"), format = "pdf", dpi = 100, width=NA)
  }
  if ("SubClass" %in% EnrichType){
    FunAnalysis <- PerformPSEA(FunAnalysis, lib = "sub_lipid_class_mset", libVersion = "current", minLib = 3, permNum = 50)
    FunAnalysis <- PlotPeaks2Paths(FunAnalysis, paste0(prefix,"_peaks_to_SLC_"), "pdf", dpi = 100, width=NA)
  }
  #Rename intermediate files
  for (file in c("mum_raw.qs","complete_norm.qs","row_norm.qs","prenorm.qs",
                 "preproc.qs","data_orig.qs",paste0("mummichog_input_",dataValue,".txt"),
                 "t_test.csv","t_test_all.csv","scattermum.json",
                 "mummichog_query.json","mummichog_pathway_enrichment.csv",
                 "mum_res.qs","mummichog_matched_compound_all.csv","initial_ecs.qs")){
    file.rename(file, paste0(prefix,"_",file))
  }
}

##' Replot mummichog dotplot across pair-wise comparisons on a figure
##' 
##'  @param mummichog_PEpath Path of directory contain mummichog enrichment results
##'  @examples
##'  Replot_mummichog(mummichog_PEpath = "./")
Replot_mummichog <- function(mummichog_PEpath = ""){

  PathwayList <- c()
  for (file in list.files(mummichog_PEpath, pattern = "mummichog_pathway_enrichment.csv")){
    PathwayList <- c(PathwayList, read.csv(file)[,1])
  }
  PathwayList <- unique(PathwayList)
  
  names_ <- c("MainClass","Condition","mLogadjP","EnrFactor")
  PathwayInte = data.frame(matrix(nrow=0, ncol = length(names_)))
  colnames(PathwayInte) <- names_
  
  for (file in list.files(mummichog_PEpath, pattern = "mummichog_pathway_enrichment.csv")){
    tmp <- read_csv(file) %>% rename_(MainClass = names(.)[1]) %>%
      mutate(Condition = gsub("_mummichog_pathway_enrichment.csv","",file), 
             mLogadjP = -log10(Gamma), EnrFactor = Hits.sig/Expected) %>%
      select(MainClass,Condition, mLogadjP, EnrFactor) %>%
      filter(mLogadjP > 1.3)
    added <- PathwayList[!PathwayList %in% tmp$MainClass]
    names_ <- c("MainClass","Condition","mLogadjP","EnrFactor")
    tmp2 = data.frame(matrix(nrow=length(added), ncol = length(names_)))
    colnames(tmp2) <- names_
    tmp2$MainClass <- added
    tmp2$Condition <- gsub("_mummichog_pathway_enrichment.csv","",file)
    tmp <- tmp %>% add_row(tmp2)
    PathwayInte <- rbind(PathwayInte, tmp)
  }
  PathwayInte$Condition <- gsub(".*_","",PathwayInte$Condition)
  pdf(paste0("PathwayEnrichement_compre_",Sys.Date(),".pdf"), width = 8, height = 8)
  p <- ggplot(PathwayInte,aes(x = Condition, y = MainClass, size = EnrFactor, color = mLogadjP)) +
    geom_point() + 
    scale_color_gradient(low="white", high="red") +
    theme_bw() +
    theme(
      axis.text.x = element_text(color="black", size=12,face="bold", angle=70, hjust=1),
      axis.text.y = element_text(color="black", size=12, face="bold"),
      legend.title = element_text(color="black", size=14, face="bold"), 
      legend.text = element_text(color="black", size=12, face="bold")
    )
  plot(p)
  dev.off()
}

##' Refine lipid identification based on given rules
##' @param STD_MeasurementsFile same as the Output_MeasurementsFile from standardize_by_standards()
##' @param RefinePlan Chose schemes to refine lipid identification (Recommend PlanC)
##' @examples
##' 
##'Refine_MZ_Identifications(STD_MeasurementsFile = "standardized_Measurements1_Converted.csv",
##'                          RefinePlan = c("PlanA","PlanB","PlanC")[3])

Refine_MZ_Identifications <- function(STD_MeasurementsFile = "",
                                      RefinePlan = c("PlanA","PlanB","PlanC")[3]){
  
  prefix <- tools::file_path_sans_ext(STD_MeasurementsFile)
  MF_std <- read_csv(STD_MeasurementsFile, show_col_types = F) %>% clean_names()
  
  if ("PlanA" %in% RefinePlan){
    # Refine based on category information
    Framebone1 <- MF_std %>% select(-main_class,-sub_class,-abbrev) %>% 
      group_by(compound, category) %>% summarise(categoryCounts = n())
    compoundList <- unique(Framebone1$compound)
    
    cns = c(colnames(MF_std)[-c(3:5)])
    Framebone2 = data.frame(matrix(nrow=0, ncol = length(cns)))
    colnames(Framebone2) <- cns
    
    for (i in compoundList){
      tmp <- filter(MF_std, compound == i)
      for (cate in unique(tmp$category)){
        tmp1 <- tmp %>% filter(category == cate) %>% select(-c(3:5))
        Framebone2 <- rbind(Framebone2,tmp1[1,])
      }
    }
    if (identical(Framebone1$category,Framebone2$category)){
      Framebone <- cbind(Framebone1,Framebone2[,-c(1,2)])
    }
    write.csv(x = Framebone, file = paste0(prefix,"_refined_PlanA_category.csv"), row.names = F)
    
    # Refine based on Main_class information
    Framebone1 <- MF_std %>% select(-category,-sub_class,-abbrev) %>% 
      group_by(compound, main_class) %>% summarise(main_classCounts = n())
    compoundList <- unique(Framebone1$compound)
    
    cns = c(colnames(MF_std)[-c(2,4,5)])
    Framebone2 = data.frame(matrix(nrow=0, ncol = length(cns)))
    colnames(Framebone2) <- cns
    
    for (i in compoundList){
      # i = compoundList[1]
      tmp <- filter(MF_std, compound == i)
      for (cate in unique(tmp$main_class)){
        tmp1 <- tmp %>% filter(main_class == cate) %>% select(-c(2,4,5))
        Framebone2 <- rbind(Framebone2,tmp1[1,])
      }
    }
    Framebone <- cbind(Framebone1,Framebone2[,-c(1,2)])
    write.csv(x = Framebone, file = paste0(prefix,"_refined_PlanA_main_class.csv"), row.names = F)
    
    # Refine based on sub_class information
    Framebone1 <- MF_std %>% select(-category,-main_class) %>% 
      group_by(compound, sub_class,abbrev) %>% summarise(Counts = n())
    compoundList <- unique(Framebone1$compound)
    
    cns = c(colnames(MF_std)[-c(2,3)])
    Framebone2 = data.frame(matrix(nrow=0, ncol = length(cns)))
    colnames(Framebone2) <- cns
    
    for (i in compoundList){
      # i = compoundList[1]
      tmp <- filter(MF_std, compound == i)
      item <- filter(Framebone1, compound == i)
      for (n in 1:nrow(item)){
        tmp1 <- tmp %>% filter(sub_class == item[n,]$sub_class & abbrev == item[n,]$abbrev) %>% select(-c(2,3))
        Framebone2 <- rbind(Framebone2, tmp1[1,])
      }
    }
    
    if (identical(Framebone1$sub_class,Framebone2$sub_class)){
      Framebone <- cbind(Framebone1, Framebone2[,-c(1:3)])
    }
    write.csv(x = Framebone, file = paste0(prefix,"_refined_PlanA_sub_class_Abbrv.csv"), row.names = F)
  }
  
  if ("PlanB" %in% RefinePlan){
    compoundList <- unique(MF_std$compound)
    #The most frequent label will be used to represent
    cns = c(colnames(MF_std), c("category_ratio","main_class_ratio","sub_class_ratio","abbrev_ratio"))
    refinedDF = data.frame(matrix(nrow=0, ncol = length(cns))) %>% as_tibble()
    colnames(refinedDF) <- cns
    
    for (i in compoundList){
      tmp <- filter(MF_std, compound == i)
      
      for (level in c("category","main_class","sub_class","abbrev")){
        
        if (! is.na(tmp[,level][[1]][1])){
          ranks <- table(tmp[,level]) %>% as.data.frame() %>% arrange(Freq)
          tmp[1,level] <- as.vector(ranks$Var1[1])
          ratio <- ranks[1,"Freq"]/unname(colSums(ranks[,2,drop = F]))
          tmp[,paste0(level,"_ratio")] <- ratio
        } else { 
          tmp[1,level] <- NA
          tmp[,paste0(level,"_ratio")] <- NA
        }
      }
      refinedDF <- rbind(refinedDF, tmp[1,])
    }
    
    refinedDF <- refinedDF %>% relocate(category_ratio, .after = category) %>% 
      relocate(main_class_ratio, .after = main_class) %>% 
      relocate(sub_class_ratio, .after = sub_class) %>%
      relocate(abbrev_ratio, .after = abbrev)
    
    write.csv(x = refinedDF, file = paste0(prefix,"_refined_PlanB.csv"), row.names = F)
  }
  
  if ("PlanC" %in% RefinePlan){
    compoundList <- unique(MF_std$compound)
    #The most frequent label will be used to represent
    refinedDF <- data.frame(matrix(nrow=0, ncol = length(colnames(MF_std)))) %>% as_tibble()
    colnames(refinedDF) <- colnames(MF_std)
    
    for (i in compoundList){
      tmp <- filter(MF_std, compound == i)
      if(nrow(tmp) == 1){
        refinedDF <- rbind(refinedDF, tmp)
      } else {
        # Only keep items with consistent category and main_class, merge abbrev when sub_class are the same
        if (length(unique(tmp$category)) == 1 & 
            length(unique(tmp$main_class)) == 1 & 
            length(unique(tmp$sub_class)) == 1){
            tmp$abbrev <- paste0(tmp$abbrev,collapse=";")
            refinedDF <- rbind(refinedDF, tmp[1,])
        }
      }
    }
    write.csv(x = refinedDF, file = paste0(prefix,"_refined_PlanC.csv"), row.names = F)
  }
}

##' Draw lipid signal on MainClass, SubClass, Abbrev level within each Category
##' @param DataFile File contain pairwise lipid signals, this should be a subset of the Output_MeasurementsFile from standardize_by_standards()
##' @param GroupAccordingTo Choose the level to plot
##' @param GroupScheme Define the replicates information
##' @param PlotSpecific Plot specific term in specific level
##' @param KeywordSelected All samples should start with this char
##' @examples
##' 
##' HeatMapDraw(DataFile = "Epi_Peri.csv", KeywordSelected = "G",
##'             PlotSpecific = c("main_class","Ceramides [SP02]"),
##'             GroupScheme = c(rep(c("Epi", "Peri"),c(3,3))),
##'             GroupAccordingTo = c("abbrev","sub_class","main_class")[2]
             
HeatMapDraw <- function(DataFile = "",
                        GroupAccordingTo = c("abbrev","sub_class")[2],
                        GroupScheme = c(rep(c("Epi", "Peri"),c(3,3))),
                        PlotSpecific = c("main_class","Ceramides [SP02]"),
                        KeywordSelected = "G"){
  
  cat("Processing",DataFile,"...\n")
  prefix <- tools::file_path_sans_ext(DataFile)
  MF_STD <- read_csv(DataFile) %>% arrange(category) %>% 
    filter(rowSums(across(where(is.numeric))) != 0)
  
  #Heatmaps of all features
  MF_STD_count <- MF_STD %>% group_by(category) %>% dplyr::summarise(Count = n())
  annotation_col = data.frame(Condition = factor(GroupScheme))
  annotation_row = data.frame(LipidCategory = factor(rep(c("FA", "GL", "GP","PK","PR","SL","SP", "ST"), 
                                                         as.vector(MF_STD_count$Count))))
  ann_colors = list(LipidCategory = c(FA = "#66C2A5", GL = "#FC8D62", GP = "#8DA0CB", PK = "#E78AC3",
                                      PR = "#A6D854",SL ="#E5C494", SP = "#FFD92F", ST = "#B3B3B3")
  ) 
  pdf(paste0(prefix,"_All_features_heatmap.pdf"), width = 10, height = 10)
  p1 <- ComplexHeatmap::pheatmap(MF_STD[,6:11], 
                                 scale = "row", cluster_rows = T, cluster_cols = T,
                                 show_rownames = FALSE, show_colnames = T,
                                 annotation_col = annotation_col, annotation_row = annotation_row,
                                 annotation_colors = ann_colors,
                                 annotation_legend = FALSE, annotation_names_col = FALSE,
                                 annotation_names_row = F,
                                 row_split = annotation_row$LipidCategory,
                                 column_split = annotation_col$Condition,
                                 fontsize = 12)
  dev.off()
  
  #Plot designed pathway
  pdf(paste0(prefix,"_",paste0(PlotSpecific,collapse = "_"),"_features_heatmap.pdf"), width = 10, height = 10)
  Forplot <- MF_STD %>% filter(!!as.name(PlotSpecific[1]) == PlotSpecific[2]) %>% as.data.frame()
  rownames(Forplot) <- make.names(Forplot$abbrev, unique = TRUE)
  p <- ComplexHeatmap::pheatmap(Forplot[,6:11], scale = "row", cluster_rows = T, cluster_cols = T,
                                show_rownames = T, show_colnames = T,
                                annotation_col = annotation_col, annotation_colors = ann_colors, 
                                annotation_legend = FALSE, annotation_names_col = FALSE,
                                annotation_names_row = F, column_split = annotation_col$Condition,
                                fontsize = 12)
  dev.off()
  
  #Heatmaps of category/main_class/sub_class/abbrev levels
  MF_STD_grouped <- MF_STD %>% group_split(category, .keep = TRUE)
  names(MF_STD_grouped) <- unique(MF_STD$category)
  
  for (i in names(MF_STD_grouped)){
    cat("Plotting on",GroupAccordingTo,"level for(",i,")...\n")
    heatmap_data <- MF_STD_grouped[[i]] %>% 
      na.omit() %>% 
      select_at(vars(GroupAccordingTo, starts_with(KeywordSelected))) %>%
      group_by(!!as.name(GroupAccordingTo)) %>% 
      filter(GroupAccordingTo != "-") %>% 
      mutate_if(is.character,as.numeric) %>% 
      summarise_all(c("sum")) %>%
      filter(!!as.name(GroupAccordingTo) != "-") %>% 
      as.data.frame()
    
    if (nrow(heatmap_data) > 0){
      rownames(heatmap_data) <- gsub("\\[.*]/","/",heatmap_data[,GroupAccordingTo])
      rownames(heatmap_data) <- gsub("\\[.*","",rownames(heatmap_data))
      
      heatmapdraw_matrix <- heatmap_data %>% select_at(vars(starts_with(KeywordSelected)))
      write.csv(heatmapdraw_matrix, paste0(prefix,"_",GroupAccordingTo,"_",i,"_value.csv"), row.names=TRUE)
      
      annotation_col = data.frame(Condition = factor(GroupScheme))
      pdf(paste0(prefix,"_",GroupAccordingTo,"_",i,"_heatmap.pdf"), width = 10, height = 10)
      p2 <- ComplexHeatmap::pheatmap(heatmapdraw_matrix, 
                                     cellwidth = 15, cellheight = 15,
                                     scale = "row", cluster_rows = T, cluster_cols = T,
                                     show_rownames = TRUE, show_colnames = T,
                                     annotation_col = annotation_col,
                                     annotation_legend = FALSE, annotation_names_col = FALSE,
                                     column_split = annotation_col$Condition,
                                     legend = TRUE, fontsize = 12)
      dev.off()
    }
  }
}

#=========================================================================================

volcano_plotting <- function(PlotFile = "volcano_fattyacyls.csv"){
  
  volcano_pk <- read_csv(PlotFile)
  volcano_pk$mob <- rowMeans(volcano_pk[,7:11])
  volcano_pk$ml <- rowMeans(volcano_pk[,12:15])
  
  vpk <- volcano_pk %>% group_by(abbrev) %>% 
    summarise(mob1 = sum(mob1), mob2 = sum(mob2), mob3 = sum(mob3),mob4 = sum(mob4),
              mob5 = sum(mob5), ml1 = sum(ml1), ml2 = sum(ml2),ml4= sum(ml4),ml5= sum(ml5))
  
  vpk <- transform(vpk, Group="")
  vpk <- transform(vpk, Group=case_when(vpk$fdr<0.05&vpk$log2_fc>1.2~"Up-regulated",
                                        vpk$fdr<0.05&vpk$log2_fc<(-1.2)~"Down-regulated",
                                        is.na(Group) ~ "Nonsignificant"))
  
  up <- subset(vpk, Group=="Up-regulated")
  # up <- up[order(up$log2_fc), ][1:10, ]
  down <- subset(vpk, Group=="Down-regulated")
  # down <- down[order(down$log2_fc), ][1:10, ]
  
  ggscatter(vpk, x="log2_fc", y="log10_fdr", color="Group", label=vpk$abbrev, 
            palette=c("#2f5688", "#BBBBBB", "#CC0000"), size=1, font.label = 8, repel=T) +
    theme_base() +
    xlab("log2FoldChange") +
    ylab("-log10(Adjusted P-value)") +
    theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12)) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10)) +
    theme(legend.title=element_text(size=12)) +
    theme(legend.text=element_text(size=11)) +
    geom_hline(yintercept=1, linetype="dashed") +
    geom_vline(xintercept=c((-1.2),1.2), linetype="dashed")
  
}


#=========================================================================================
# Run example
#=========================================================================================
standardize_by_standards(MeasurementsFile = "test1.csv", 
                         ISMeasurementsFile = "IS-Measurements-1.xls",
                         IdentificationFile = "Identifications-1.csv",
                         LipidMAPSDB = "lipid_database.csv",
                         Output_MeasurementsFile = "standardized_test1_Converted_zhy.csv",
                         GroupLabel = list(c(1),c(2:4),c(5:7),c(8:10),c(11:13),c(14:16),
                                                c(17:19),c(20:22),c(23:25),c(26:28),c(29:31),c(32:34),
                                                c(35:37),c(38:42),c(43:47),c(48:50),c(51:53)))

Format_mummichog_input(STD_MeasurementsFile = "standardized_Measurements1_Converted.csv", OutputPairWise = T)

mummichog_Functional_analysis(pktablePath = "standardized_Measurements1_Converted_mummichogInput_G1+G3.csv",
                              PvalueThreshold = 0.5, 
                              EnrichType = c("KEGG","MainClass","SubClass")[c(1,2)])

Replot_mummichog(mummichog_PEpath = "./")

Refine_MZ_Identifications(STD_MeasurementsFile = "standardized_Measurements1_Converted.csv",
                          RefinePlan = c("PlanA","PlanB","PlanC")[3])

HeatMapDraw(DataFile = "Epi_Peri.csv", KeywordSelected = "G",
            PlotSpecific = c("main_class","Ceramides [SP02]"),
            GroupScheme = c(rep(c("Epi", "Peri"),c(3,3))),
            GroupAccordingTo = c("abbrev","sub_class","main_class")[2]
)




