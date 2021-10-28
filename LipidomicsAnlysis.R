#=========================================================================================
# This script contain main functions used for lipidomics analysis

# Version 1.4 created by Houyu Zhang on 2021/10/27
# Issue report on Hughiez047@gmail.com
# Copyright (c) 2021 __CarlosLab@PKU__. All rights reserved.
#=========================================================================================

# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(janitor))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(MetaboAnalystR))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(pls))
suppressPackageStartupMessages(library(tools))

#=========================================================================================
# Released functions
#=========================================================================================

##' Standardize samples against internal controls within defined groups
##' 
##' @param MeasurementsFile LPMS output file containing compound signal across sample
##' @param ISMeasurementsFile LPMS output file containing counts for internal standards
##' @param MeasurementsFile_standardized Specify a output csv file name for standardized-converted lipid signals
##' @param GroupLabel Defined groups for the excel sheet which will be standardized together
##' @examples
##'
##' standardize_by_standards(MeasurementsFile = "test1.csv",
##'                          ISMeasurementsFile = "IS-Measurements-1.xls",
##'                          MeasurementsFile_standardized = "test1_standardized.csv",
##'                          GroupLabel = list(c(1),c(2:4),c(5:7),c(8:10),c(11:13),c(14:16),
##'                                                c(17:19),c(20:22),c(23:25),c(26:28),c(29:31),c(32:34),
##'                                                c(35:37),c(38:42),c(43:47),c(48:50),c(51:53)))

standardize_by_standards <- function(MeasurementsFile = "", 
                                     ISMeasurementsFile = "", 
                                     MeasurementsFile_standardized = "", 
                                     GroupLabel = list()){
  
  MF_original <- read_csv(MeasurementsFile, skip = 2, show_col_types = F) %>% clean_names() %>% 
    filter(retention_time_min > 1.5)
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
  write.csv(MF_standardized[,-c(14:(13+SampleNumbers))], MeasurementsFile_standardized, row.names = F)

  cat("<<",MeasurementsFile,">> has been standardized using information in <<",ISMeasurementsFile, 
      ">> and wrote into <<",MeasurementsFile_standardized,">>\n")
}

##' Map standardized compounds(rt_m/z) against LipidMAPS database
##' 
##' @param MeasurementsFile_standardized A standardized lipid file output by standardize_by_standards() 
##' @param IdentificationFile LPMS output file containing possible identification for this LPMS experiment 
##' @param LipidMAPSDB A lipid information database file downloaded from LIPIDMAPS website
##' @param OutUnmatchedLipids Whether output unmatched lipids through the LIPIDMAPS database (identification = 0)
##' @param MeasurementsFile_standardized_matched Output matched file
##' @examples
##' 
##' MatchLipidMAPS(MeasurementsFile_standardized = "test1_standardized.csv",
##'                IdentificationFile = "Identifications-1.csv",
##'                LipidMAPSDB = "lipid_database.csv",
##'                OutUnmatchedLipids = T,
##'                MeasurementsFile_standardized_matched = "test1_standardized_matched.csv"
##' )

MatchLipidMAPS <- function(MeasurementsFile_standardized = "",
                           IdentificationFile = "",
                           LipidMAPSDB = "",
                           OutUnmatchedLipids = T,
                           MeasurementsFile_standardized_matched = ""){
  
  MF_standardized <- read_csv(MeasurementsFile_standardized, show_col_types = F) %>% clean_names()
  identification <- read_csv(IdentificationFile, show_col_types = F) %>% clean_names()
  lipid_database <- read_csv(LipidMAPSDB, show_col_types = F) %>% clean_names()
  
  if(OutUnmatchedLipids){
    prefix <- file_path_sans_ext(MeasurementsFile_standardized)
    UM <- MF_standardized %>% filter(identifications == 0) %>% select(-c(neutral_mass_da:minimum_cv))
    write.csv(x = UM, file = paste0(prefix,"_unmatched.csv"), row.names = F)
  }

  MF_standardized_matched <- lipid_database %>% 
    select(compound_id, mass, category, main_class, sub_class, abbrev) %>% 
    right_join(identification, by="compound_id") %>% 
    select(compound, category, main_class, sub_class, abbrev) %>% 
    right_join(MF_standardized, by="compound") %>% 
    filter(category != "NA") %>%
    select(-c(neutral_mass_da:minimum_cv)) %>% 
    arrange(compound) 
  
  write.csv(MF_standardized_matched, MeasurementsFile_standardized_matched, row.names = F)
  cat("Finished matching LipidMAPS database!\n")
}

##' Convert Standard (STD) Lipid signal file to mummichog accepted format (m/z__rt)
##' 
##' @param STD_MeasurementsFile same as the Output_MeasurementsFile from standardize_by_standards()
##' @param OutputPairWise Besides output all samples, indicate whether output pairwise samples
##' @param SelectedSample Specify samples head with these strings for output in a file
##' @examples
##'
##' Format_mummichog_input(MeasurementsFile_standardized = "standardized_Measurements1_Converted.csv",
##'                        OutputPairWise = T, SelectedSample = c("epi","Peri","ibat"))
##' 
                     
Format_mummichog_input <- function(MeasurementsFile_standardized = "", 
                                   OutputPairWise = TRUE,
                                   SelectedSample = ""){
  
  prefix <- file_path_sans_ext(MeasurementsFile_standardized)
  mummichog_MeasurementsFile <- read_csv(MeasurementsFile_standardized, show_col_types = F) %>% clean_names() %>% 
    mutate(compound = gsub("[a-z]|/","",compound)) %>% 
    separate(compound, c("RT", "MZ"), "_") %>% 
    mutate(Label = paste0(MZ,"__",RT)) %>% 
    relocate(Label, .before = RT) %>% 
    select(-c(RT:minimum_cv)) %>%
    unique()
  
  # colnames(mummichog_MeasurementsFile) <- c("Label",c("G1_1","G1_2","G1_3","G2_1","G2_2","G2_3",
  #                                                     "G3_1","G3_2","G3_3","G4_1","G4_2","G4_3"))
  
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
        tmp <- mummichog_MeasurementsFile %>% select(Label, starts_with(x), starts_with(y))
        fileName <- paste0(prefix,"_mummichogInput_",x,"+",y,".csv")
        write.csv(file = fileName, tmp, row.names = F, quote = F)
      }
    }
  }
  
  if(SelectedSample != ""){
    tmp <- mummichog_MeasurementsFile %>% select(Label,starts_with(SelectedSample))
    fileName <- paste0(prefix,"_mummichogInput_",paste0(SelectedSample, collapse = "+"),".csv")
    write.csv(file = fileName, tmp, row.names = F, quote = F)
  }
  cat("Finished converting compound to mummichog supported format!\n")
}

##' Do mummichog analysis on KEGG/MainClass/SubClass using m/z__rt information
##' 
##' @param ktablePath A pair-wise sample file output by Format_mummichog_input()
##' @param PvalueThreshold Threshold to selected lipids for enrichment analysis
##' @param Species Specify species for KEGG analysis, support "mmu" and "hsa" now
##' @param rowNormMet method for column normalization, recommend using NULL for data has internal controls
##' @param EnrichType Define enrichment methods, this is multi-optional
##' @param SamplesDel samples to delete, these samples might be sudo-samples for MetaboAnalyst checking
##' @examples
##' 
##' mummichog_Functional_analysis(pktablePath = "standardized_Measurements1_Converted_mummichogInput_G1+G3.csv",
##'                               PvalueThreshold = 0.5, Species = c("mmu","hsa")[2], rowNormMet = c("SumNorm","NULL")[2],
##'                               EnrichType = c("KEGG","MainClass","SubClass")[c(1,2)])

mummichog_Functional_analysis <- function(pktablePath = "", 
                                          PvalueThreshold = 0.2, 
                                          Species = c("mmu","hsa")[2],
                                          rowNormMet = c("SumNorm","NULL")[2],
                                          EnrichType = c("KEGG","MainClass","SubClass")[c(1,2)],
                                          SamplesDel = ""){

  prefix <- file_path_sans_ext(pktablePath)
  
  FunAnalysis <- InitDataObjects("pktable", data.type = "stat", paired = FALSE);
  FunAnalysis <- Read.TextData(FunAnalysis, filePath = pktablePath, format = "colu", lbl.type = "disc")
  FunAnalysis <- SanityCheckData(FunAnalysis)
  FunAnalysis <- ReplaceMin(FunAnalysis)
  FunAnalysis <- FilterVariable(FunAnalysis, filter = "iqr", qcFilter = "F", rsd = 25)
  
  if(SamplesDel != ""){
    FunAnalysis <- GetGroupNames(FunAnalysis, "")
    feature.nm.vec <- c("")
    smpl.nm.vec <- SamplesDel
    grp.nm.vec <- c("")
    FunAnalysis <- UpdateData(FunAnalysis)
  }
  
  FunAnalysis <- PreparePrenormData(FunAnalysis)
  FunAnalysis <- Normalization(FunAnalysis, rowNorm = rowNormMet, transNorm = "NULL", scaleNorm = "NULL",
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
    FunAnalysis <- PerformPSEA(FunAnalysis, lib = paste0(Species,"_kegg"), libVersion = "current", minLib = 3, permNum = 50)
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
##'  
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
##' @param MeasurementsFile_standardized_matched same as the Output_MeasurementsFile from standardize_by_standards()
##' @param RefinePlan Chose schemes to refine lipid identification (Recommend PlanC)
##' @examples
##' 
##'Refine_MZ_Identifications(MeasurementsFile_standardized_matched = "standardized_Measurements1_Converted.csv",
##'                          RefinePlan = c("PlanA","PlanB","PlanC")[3])

Refine_MZ_Identifications <- function(MeasurementsFile_standardized_matched = "",
                                      RefinePlan = c("PlanA","PlanB","PlanC")[c(3)]){
  
  prefix <- file_path_sans_ext(MeasurementsFile_standardized_matched)
  MF_std <- read_csv(MeasurementsFile_standardized_matched, show_col_types = F) %>% clean_names()
  
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
    refinedDF <- data.frame(matrix(nrow=0, ncol = length(colnames(MF_std)))) %>% as_tibble()
    colnames(refinedDF) <- colnames(MF_std)
    
    for (i in compoundList){
      CurrentCompound <- filter(MF_std, compound == i)
      if(nrow(CurrentCompound) == 1){
        refinedDF <- rbind(refinedDF, CurrentCompound)
      } else {
        # Only keep items with consistent category and main_class, merge abbrev when sub_class are the same
        if (length(unique(CurrentCompound$category)) == 1 &  length(unique(CurrentCompound$main_class)) == 1){
            if(length(unique(CurrentCompound$sub_class)) == 1){
              CurrentCompound$abbrev <- paste0(unique(CurrentCompound$abbrev),collapse="/")
              refinedDF <- rbind(refinedDF, CurrentCompound[1,])
            } else {
              tbl <- CurrentCompound %>% group_by(abbrev) %>% plyr::count() %>% arrange(desc(freq)) %>%
                mutate_at(vars(freq), funs(./sum(.)))
              if(tbl$freq[1] >= 0.8){
                KeyAbbrev <- tbl[1,]$abbrev
                CurrentCompound <- filter(CurrentCompound, abbrev == KeyAbbrev) %>% 
                  mutate(sub_class = paste0(unique(sub_class),collapse = "/"))
                refinedDF <- rbind(refinedDF, CurrentCompound[1,])
              }
            }
        }
      }
    }
    write.csv(x = refinedDF, file = paste0(prefix,"_refined_PlanC.csv"), row.names = F)
  }
  cat("Finished refining lipids using",RefinePlan,"\n")
}

##' Draw costume Volcano Plot using ggplot
##' @param PlotFile The file returned by Volcano.Anal()
##' @param ThresholdFC Threshold for fold-change
##' @param ThresholdSig Threshold for significance
##' @param TopUpDownShown Slected top up,down significant features for shown
##' @examples
##' 
##' volcano_plotting(PlotFile = "refinedDF_5depots_abbrev_adi_epi+adi_ibat_volcano.csv",
##'                  ThresholdFC = 1.5, ThresholdSig = 0.05, TopUpDownShown = c(10,10))
##'                  
volcano_plotting <- function(PlotFile = "",
                             ThresholdFC = 1.5,
                             ThresholdSig = 0.05,
                             TopUpDownShown = c(10,10)){
  volcano_pk <- read_csv(PlotFile, show_col_types = FALSE) %>% 
    clean_names() %>% 
    transform(Group=NA) %>%
    transform(Group=case_when(raw_pval < ThresholdSig & log2_fc > ThresholdFC ~ "Sig.Up",
                              raw_pval < ThresholdSig & log2_fc < -ThresholdFC ~ "Sig.Down",
                              is.na(Group) ~ "Nonsig."))
  
  tbl <- table(volcano_pk$Group) 
  volcano_pk <- volcano_pk %>% 
    transform(Group=case_when(Group == "Sig.Up" ~ paste0("Sig.Up [",ifelse(is.na(tbl["Sig.Up"]),0,tbl[["Sig.Up"]]),"]"),
                              Group == "Sig.Down" ~ paste0("Sig.Down [",ifelse(is.na(tbl["Sig.Down"]),0,tbl[["Sig.Down"]]),"]"),
                              Group == "Nonsig." ~ paste0("Nonsig. [",ifelse(is.na(tbl["Nonsig."]),0,tbl[["Nonsig."]]),"]")))
  
  volcano_pk$Label <- ""
  
  JudgeNumUp <- ifelse(is.na(tbl["Sig.Up"]),0,tbl[["Sig.Up"]])
  if(JudgeNumUp !=0 ){
    UpShownNum <- ifelse(JudgeNumUp < TopUpDownShown[1], JudgeNumUp, TopUpDownShown[1])
    UpShownItem <- volcano_pk %>% arrange(desc(log10_p)) %>% filter(str_detect(Group, "Sig.Up")) %>% slice(1:UpShownNum)
    volcano_pk$Label[volcano_pk$x1 %in% UpShownItem$x1] <- UpShownItem$x1
  }

  JudgeNumDown <- ifelse(is.na(tbl["Sig.Down"]),0,tbl[["Sig.Down"]])
  if(JudgeNumDown != 0){
    DownShownNum <- ifelse(JudgeNumDown < TopUpDownShown[2], JudgeNumDown, TopUpDownShown[2])
    DownShownItem <- volcano_pk %>% arrange(desc(log10_p)) %>% filter(str_detect(Group, "Sig.Down")) %>% slice(1:DownShownNum)
    volcano_pk$Label[volcano_pk$x1 %in% DownShownItem$x1] <- DownShownItem$x1
  }
  
  prefix <- file_path_sans_ext(PlotFile)
  pdf(paste0(prefix,"_volcanoPlot_Custome_FC",ThresholdFC,"_Sig",ThresholdSig,".pdf"), height = 8, width = 12)
  p <- ggplot(volcano_pk, aes(x=log2_fc, y=log10_p, color=Group, label=Label)) +
    geom_point(shape=21) + 
    theme_bw() + geom_text_repel() + 
    scale_color_manual(values=c("grey", "#1F78B4", "#E31A1C")) +
    geom_vline(xintercept=c(-ThresholdFC, ThresholdFC), col="black", linetype="dashed") +
    geom_hline(yintercept=-log10(ThresholdSig), col="black", linetype="dashed") + 
    labs(x = "log2(FoldChange)", y = "-log10(P-value)") + 
    theme(
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.text.x = element_text(color="black", size=12, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.text.y = element_text(color="black", size=12, face="bold"),
      legend.title = element_blank(), 
      legend.text = element_text(color="black", size=12, face="bold")
    )
  plot(p)
  dev.off()
}

##' Draw costume scatter plot for pair-wise T.test or multi-group anova result
##' @param PlotFile The file returned by Volcano.Anal()
##' @param ThresholdSig Threshold for significance
##' @param SigIndex select using raw p-value (p_value) or adjusted p-value (fdr)
##' @param TopShown Top significant items for labeling
##' @examples
##' 
##' T_Anavo_plotting(PlotFile = "refinedDF_5depots_abbrev_anova_posthoc.csv",
##'                  SigIndex = c("p_value","fdr")[1],ThresholdSig = 0.005,TopShown = 10)
##'                  
T_Anavo_plotting <- function(PlotFile = "",
                             SigIndex = c("p_value","fdr")[1],
                             ThresholdSig = 0.005,
                             TopShown = 10){
  prefix <- file_path_sans_ext(PlotFile)
  
  T_Anavo_pk <- read_csv(PlotFile, show_col_types = FALSE) %>% clean_names() %>%
    transform(Group=case_when(eval(parse(text = SigIndex)) < ThresholdSig ~ "Significant",
                              eval(parse(text = SigIndex)) >= ThresholdSig  ~ "Nonsignificant"))
  
  tbl <- table(T_Anavo_pk$Group)
  T_Anavo_pk <- T_Anavo_pk %>% 
    transform(Group=case_when(Group == "Significant" ~ paste0("Significant [",ifelse(is.na(tbl["Significant"]),0,tbl[["Significant"]]),"]"),
                              Group == "Nonsignificant" ~ paste0("Nonsignificant [",ifelse(is.na(tbl["Nonsignificant"]),0,tbl[["Nonsignificant"]]),"]"))) %>%
    arrange(eval(parse(text = SigIndex)))
  
  T_Anavo_pk$Label <- ""
  JudgeNum <- ifelse(is.na(tbl["Significant"]),0,tbl[["Significant"]])[[1]]
  
  if(JudgeNum != 0){
    TopShownNum <- ifelse(JudgeNum < TopShown, JudgeNum, TopShown)
    TopShownItem <- T_Anavo_pk %>% filter(str_detect(Group, "Significant")) %>% slice(1:TopShownNum)
    T_Anavo_pk$Label[T_Anavo_pk$x1 %in% TopShownItem$x1] <- TopShownItem$x1
  }
  
  # T_Anavo_pk$x1 <- factor(T_Anavo_pk$x1, levels = T_Anavo_pk$x1)
  
  pdf(paste0(prefix,"_TAnavoPlot_Custome_Sig",ThresholdSig,".pdf"), height = 7, width = 10)
  p <- ggplot(T_Anavo_pk, aes(x=x1, y=log10_p, color=Group, label=Label)) +
    geom_point(shape=21) + 
    theme_bw() + geom_text_repel() + 
    scale_color_manual(values=c("grey", "#1F78B4", "#E31A1C")) +
    geom_hline(yintercept=-log10(ThresholdSig), col="black", linetype="dashed") + 
    labs(x = "Lipid features", y = "-log10(P-value)") + 
    theme(
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top",
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.text.y = element_text(color="black", size=12, face="bold"),
      legend.title = element_blank(), 
      legend.text = element_text(color="black", size=12, face="bold")
    )
  plot(p)
  dev.off()
}

##' Run MetaboAnalystR codes
##' @param pktablePath A standard MetaboAnalystR input file
##' @param rowNormMet method for column normalization, recommend using NULL for data has internal controls
##' @examples
##' 
##' RunMetaboAnalystR(pktablePath = "refinedDF_5depots_abbrev_adi_peri+adi_ibat.csv",
##' rowNormMet = c("SumNorm","NULL")[1])
##' 
RunMetaboAnalystR <- function(pktablePath = "",
                              rowNormMet = c("SumNorm","NULL")[1]){

  prefix <- file_path_sans_ext(pktablePath)
  # Step1. create the mSet Object, specifying that the data to be uploaded
  # is a peak table ("pktable") and that statistical analysis will be performed ("stat").
  mSet <- InitDataObjects(data.type = "pktable", anal.type = "stat", paired = FALSE)
  # read in the filtered peak list
  mSet <- Read.TextData(mSetObj = mSet, filePath = pktablePath, format = "colu", lbl.type = "disc")
  
  # Step2. perform data processing (filtering/normalization)
  mSet <- SanityCheckData(mSetObj = mSet)
  # Perform data processing - Minimum Value Replacing
  mSet <- ReplaceMin(mSetObj = mSet)
  mSet <- SanityCheckData(mSetObj = mSet)
  mSet <- FilterVariable(mSetObj = mSet, filter = "none", qcFilter = "F", rsd = 20)
  mSet <- PreparePrenormData(mSetObj = mSet)
  mSet <- Normalization(mSetObj = mSet, rowNorm = rowNormMet, transNorm = "NULL", scaleNorm = "NULL", ratio=FALSE, ratioNum=20)
  mSet <- PlotNormSummary(mSetObj = mSet, imgName = paste0(prefix,"_NormFeature_"), format="pdf", dpi = 100, width=NA)
  mSet <- PlotSampleNormSummary(mSetObj = mSet, imgName=paste0(prefix,"_NormSample_"), format="pdf", dpi = 100, width=NA)
  
  #Step3. Do differential lipid detection
  if(length(unique(mSet$dataSet$prenorm.cls)) == 2){
    cat("Only 2 groups detected, will do t.test...\n")
    mSet <- Ttests.Anal(mSet, nonpar = F, threshp = 0.05, paired = F, equal.var = F, pvalType = "raw", all_results = T)
    # mSet <- PlotTT(mSet, imgName = paste0(prefix,"_Ttest_"), format = "pdf", dpi = 100, width=NA)
    file.rename("fold_change.csv", paste0(prefix,"_fold_change.csv"))
    file.rename("t_test.csv", paste0(prefix,"_t_test.csv"))
    file.rename("t_test_all.csv", paste0(prefix,"_t_test_all.csv"))
    T_Anavo_plotting(PlotFile = paste0(prefix,"_t_test_all.csv"), ThresholdSig = 0.05)
    
    mSet <- Volcano.Anal(mSet, paired = FALSE, fcthresh = 1, cmpType = 0, nonpar = F, threshp = 1, 
                         equal.var = FALSE, pval.type = "raw")
    file.rename("volcano.csv", paste0(prefix,"_volcano.csv"))
    volcano_plotting(PlotFile = paste0(prefix,"_volcano.csv"), ThresholdFC = 1.5, ThresholdSig = 0.05)
    # mSet <- PlotVolcano(mSet, paste0(prefix,"_volcanoPlot_"), plotLbl = 1, format = "pdf", dpi = 100, width=NA)
    
  } else if (length(unique(mSet$dataSet$prenorm.cls)) > 2){
    cat("More than 2 groups detected, will do Anova test...\n")
    mSet <- ANOVA.Anal(mSetObj = mSet, nonpar = F, thresh = 0.05, post.hoc = "fisher", all_results = FALSE)
    # mSet <- PlotANOVA(mSetObj = mSet, imgName = paste0(prefix,"_anova_"), format = "pdf", dpi = 100, width=NA)
    file.rename("anova_posthoc.csv", paste0(prefix,"_anova_posthoc.csv"))
    T_Anavo_plotting(PlotFile = paste0(prefix,"_anova_posthoc.csv"), ThresholdSig = 0.005)
  }
  # Step4. Plot overall heatmap view
  mSet <- PlotSubHeatMap(mSet, imgName = paste0(prefix,"_FeatureSignalHeatmap_"), format = "pdf", dpi = 100, width=NA, 
                         dataOpt = "norm", scaleOpt = "row", smplDist = "euclidean",clstDist = "ward.D",
                         palette = "bwm", method.nm = "tanova", top.num = 100, viewOpt = "overview", 
                         rowV = T, colV = T, border = T, grp.ave = F)
  mSet <- PlotCorrHeatMap(mSet, imgName = paste0(prefix,"_CorrSample_"), format = "pdf", dpi = 100, width=NA, 
                          target = "row", cor.method = "pearson", colors = "bwm", viewOpt = "overview", fix.col = T, 
                          no.clst = F, corrCutoff = "0")
  file.rename("correlation_table.csv", paste0(prefix,"_CorrSample_table.csv"))
  mSet <- PlotCorrHeatMap(mSet, imgName = paste0(prefix,"_CorrFeature_"), format = "pdf", dpi = 100, width=NA, 
                          target = "col", cor.method = "pearson", colors = "bwm", viewOpt = "overview", fix.col = T, 
                          no.clst = F, corrCutoff = "0")
  file.rename("correlation_table.csv", paste0(prefix,"_CorrFeature_table.csv"))
  
  # Step5. perform PCA
  mSet <- PCA.Anal(mSetObj = mSet)
  mSet <- PlotPCAPairSummary(mSetObj = mSet, imgName=paste0(prefix,"_pca_pair_"), format="pdf", dpi = 100, width=NA, pc.num = 5)
  mSet <- PlotPCAScree(mSetObj = mSet, imgName=paste0(prefix,"_pca_scree_"), format="pdf", dpi = 100, width=NA, scree.num = 5)
  mSet <- PlotPCA2DScore(mSetObj=mSet, imgName=paste0(prefix,"_pca_score2d_"), format="pdf", 
                         72, width=NA, pcx=1, pcy=2, reg=0.95, show=1, grey.scale=0)
  mSet <- PlotPCALoading(mSetObj = mSet, imgName=paste0(prefix,"_pca_loading_"), format="pdf", dpi = 100, width=NA, inx1 = 1,inx2 = 2)
  mSet <- PlotPCABiplot(mSetObj = mSet, imgName=paste0(prefix,"_pca_biplot_"), format="pdf", dpi = 100, width=NA, inx1 = 1,inx2 = 2)
  mSet <- PlotPCA3DScoreImg(mSet, imgName = paste0(prefix,"_pca_score3d_"), "pdf", dpi =100, width=NA, 1,2,3, angl = 40)
  
  # Step6. perform PLS-DA
  mSet <- PLSR.Anal(mSet, reg=TRUE)
  mSet <- PlotPLSPairSummary(mSet, paste0(prefix,"_pls_pair_"), format = "pdf", dpi = 100, width=NA, pc.num = 5)
  mSet <- PlotPLS2DScore(mSet, paste0(prefix,"_pls_score2d_"), format = "pdf", dpi = 100, width=NA,
                         inx1 = 1,inx2 = 2,reg = 0.95,show = 1,grey.scale = 0)
  library(pls)
  mSet <- PlotPLS3DScoreImg(mSet, paste0(prefix,"_pls_score3d_"), format = "pdf", dpi = 100, width=NA,1,2,3,40)
  mSet <- PlotPLSLoading(mSet, paste0(prefix,"_pls_loading_"), format = "pdf", dpi = 100, width=NA, 1, 2)
  mSet <- PLSDA.CV(mSet, methodName = "L", compNum = 3, choice = "Q2")
  mSet <- PlotPLS.Classification(mSet, paste0(prefix,"_pls_cv_"), format = "pdf", dpi = 100, width=NA)
  mSet <- PlotPLS.Imp(mSet, paste0(prefix,"_pls_imp_"), format = "pdf", dpi = 100, width=NA,
                      type = "vip", feat.nm = "Comp. 1", feat.num = 30, color.BW = FALSE)
  
  for (fileName in c("data_orig.qs","preproc.qs","prenorm.qs","row_norm.qs","complete_norm.qs",
                     "pca_loadings.csv","pca_score.csv","plsda_coef.csv","plsda_loadings.csv",
                     "plsda_score.csv","plsda_vip.csv"
  )){
    file.rename(fileName, paste0(prefix,"_",fileName))
  }
}

##' Unify the Lipid signal at each level and return all/pairwise samples
##' @param Refined_MeasurementsFile Refined lipid signal file output from Refine_MZ_Identifications()
##' @param OutputPairWise Besides output all samples, indicate whether output pairwise samples
##' @param KeywordSelected All samples should start with this char
##' @param run_RunMetaboAnalystR Whether run MetaboAnalystR analysis in the meanwhile 
##' @examples
##' 
##' Analyze_Lipids(MeasurementsFile_standardized_matched_refined = "refinedDF_5depots.csv",
##'                OutputPairWise = T,
##'                KeywordSelected = "adi",
##'                run_RunMetaboAnalystR = T)
            
Analyze_Lipids <- function(MeasurementsFile_standardized_matched_refined = "",
                           OutputPairWise = T,
                           KeywordSelected = "adi",
                           run_RunMetaboAnalystR = T){
  
  LipidLevels <- c("compound","category","main_class","sub_class","abbrev")
  prefix <- file_path_sans_ext(MeasurementsFile_standardized_matched_refined)
  
  for (level in LipidLevels){
    slectedrefinedDF <- read_csv(MeasurementsFile_standardized_matched_refined, show_col_types = FALSE) %>% 
      select_at(vars(level, starts_with(KeywordSelected))) %>%
      group_by(!!as.name(level)) %>% 
      summarise_all(c("sum")) %>% 
      filter(!!as.name(level) != "-")
    
    labels <- sub(level,"Label",gsub("_[1-3]","",colnames(slectedrefinedDF)))
    slectedrefinedDF <- rbind(labels, slectedrefinedDF)
    
    fileName <- paste0(prefix,"_",level,".csv")
    cat("Processing",fileName,"...\n")
    write.csv(file = fileName, slectedrefinedDF, row.names = F)
    if(run_RunMetaboAnalystR){RunMetaboAnalystR(pktablePath = fileName, rowNormMet = "SumNorm")}
    
    if(OutputPairWise){
      groups <- unique(gsub("_[0-9]+","",colnames(slectedrefinedDF))[-1])
      if (length(groups) > 1){
        pairwise <- combn(groups,2)
        for (i in 1:ncol(pairwise)){
          x = pairwise[1,i]
          y = pairwise[2,i]
          tmp <- slectedrefinedDF %>% select(level,starts_with(x),starts_with(y))
          
          fileName <- paste0(prefix,"_",level,"_",x,"+",y,".csv")
          cat("Processing",fileName,"...\n")
          write.csv(file = fileName, tmp, row.names = F)
          if(run_RunMetaboAnalystR){RunMetaboAnalystR(pktablePath = fileName, rowNormMet = "SumNorm")}
        }
      }
    }
  }
}

##' Draw lipid signal on MainClass, SubClass, Abbrev level within each Category
##' @param MeasurementsFile_standardized_matched_refined File contain pairwise lipid signals, this can be a subset of the Refine_MZ_Identifications()
##' @param Level Choose the level to plot
##' @param GroupScheme Define the replicates information
##' @param PlotSpecific Plot specific term in specific level
##' @param KeywordSelected All samples should start with this char
##' @examples
##' 
##' HeatMapDraw(MeasurementsFile_standardized_matched_refined = "Epi_Peri.csv", KeywordSelected = "G",
##'             PlotSpecific = c("main_class","Ceramides [SP02]"),
##'             GroupScheme = c(rep(c("Epi", "Peri"),c(3,3))),
##'             Level = c("abbrev","sub_class","main_class")[2]

HeatMapDraw <- function(MeasurementsFile_standardized_matched_refined = "",
                        Level = c("main_class","sub_class","abbrev")[2],
                        GroupScheme = c(rep(c("Epi", "Peri"),c(3,3))),
                        PlotSpecific = c("main_class","Ceramides [SP02]"),
                        KeywordSelected = "adi"){
  
  cat("Processing",Refined_MeasurementsFile,"...\n")
  prefix <- file_path_sans_ext(Refined_MeasurementsFile)
  MF_STD <- read_csv(Refined_MeasurementsFile) %>% arrange(category) %>% 
    filter(rowSums(across(where(is.numeric))) != 0)
  
  #Heatmaps of all features
  MF_STD_count <- MF_STD %>% group_by(category) %>% dplyr::summarise(Count = n())
  annotation_col = data.frame(Condition = factor(GroupScheme))
  ExistCategory <- gsub(".*\\[|\\]","",unique(MF_STD$category))
  annotation_row = data.frame(LipidCategory = factor(rep(ExistCategory, as.vector(MF_STD_count$Count))))
  ann_colors = list(Condition = structure(c("#6A3D9A","#B15928"), names = as.vector(unique(annotation_col$Condition))),
                    LipidCategory = structure(brewer.pal(length(ExistCategory), "Paired"), names = ExistCategory))
  
  Forplot <- MF_STD %>% select_at(vars(starts_with(KeywordSelected)))
  pdf(paste0(prefix,"_All_features_heatmap.pdf"), width = 10, height = 10)
  p1 <- ComplexHeatmap::pheatmap(Forplot, 
                                 scale = "row", cluster_rows = T, cluster_cols = T,
                                 show_rownames = F, show_colnames = T,
                                 annotation_col = annotation_col, annotation_row = annotation_row,
                                 annotation_colors = ann_colors,
                                 annotation_legend = F, annotation_names_col = F,
                                 annotation_names_row = F,
                                 row_split = annotation_row$LipidCategory,
                                 column_split = annotation_col$Condition,
                                 fontsize = 12)
  dev.off()
  
  #Plot designed pathway
  Forplot <- MF_STD %>% filter(!!as.name(PlotSpecific[1]) == PlotSpecific[2]) %>% 
    select_at(vars(starts_with(KeywordSelected))) %>% as.data.frame()
  rownames(Forplot) <- make.names(Forplot$abbrev, unique = TRUE)
  
  pdf(paste0(prefix,"_",paste0(PlotSpecific,collapse = "_"),"_features_heatmap.pdf"), width = 10, height = 10)
  p <- ComplexHeatmap::pheatmap(Forplot, scale = "row", cluster_rows = T, cluster_cols = T,
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
    cat("Plotting on",Level,"level for(",i,")...\n")
    heatmap_data <- MF_STD_grouped[[i]] %>% 
      na.omit() %>% 
      select_at(vars(Level, starts_with(KeywordSelected))) %>%
      group_by(!!as.name(Level)) %>% 
      mutate_if(is.character,as.numeric) %>% 
      summarise_all(c("sum")) %>%
      filter(!!as.name(Level) != "-") %>% 
      as.data.frame()
    
    if (nrow(heatmap_data) > 0){
      rownames(heatmap_data) <- gsub("\\[.*]/","/",heatmap_data[,Level])
      rownames(heatmap_data) <- gsub("\\[.*","",rownames(heatmap_data))
      
      heatmapdraw_matrix <- heatmap_data %>% select_at(vars(starts_with(KeywordSelected)))
      write.csv(heatmapdraw_matrix, paste0(prefix,"_",Level,"_",i,"_value.csv"), row.names=TRUE)
      
      annotation_col = data.frame(Condition = factor(GroupScheme))
      pdf(paste0(prefix,"_",Level,"_",i,"_heatmap.pdf"), width = 10, height = 10)
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
# Developing functions
#=========================================================================================
.RelativeDiffPlot <- function(MeasurementsFile_standardized_matched_refined = "refinedDF_5depots.csv",
                             OutputPairWise = T,
                             SelectedSample = c("adi_epi","adi_peri")){
  
  LipidLevels <- c("main_class","sub_class","abbrev")
  prefix <- file_path_sans_ext(MeasurementsFile_standardized_matched_refined)
  
  for (level in LipidLevels){
    # slectedrefinedDF <- 
      read_csv(MeasurementsFile_standardized_matched_refined, show_col_types = FALSE) %>% 
      select_at(vars(category, level, starts_with(SelectedSample))) %>%
      mutate(MergedLevel = paste0(category,"__",!!as.name(level))) %>%
      relocate(MergedLevel, .before = category) %>%
      group_by(MergedLevel) %>% 
      select(-c(category,level)) %>%
      summarise_all(c("sum")) %>% 
      separate(MergedLevel, c("category", level), "__") %>%
      mutate_at(vars(starts_with(SelectedSample)), funs(./sum(.))) %>%
      rowwise() %>%
      mutate(TreatMean = mean())
    
    # Relative Difference = (experimental - control)/control *100%

  }
  
  if(SelectedSample != ""){
    tmp <- mummichog_MeasurementsFile %>% select(Label,starts_with(SelectedSample))
    fileName <- paste0(prefix,"_mummichogInput_",paste0(SelectedSample, collapse = "+"),".csv")
    write.csv(file = fileName, tmp, row.names = F, quote = F)
  }
  

  # lipid_cert_exp1_ma <- read_csv("lipid_cert_exp1_ma.csv")
  # certmaexp1_pvalue_distinct <- read_csv("certmaexp1_pvalue_distinct.csv")
  # 
  # # data_bubble <- lipid_cert_exp1_ma %>% 
  # #   select(compound, ma_ce_ibat_1, ma_ce_ibat_2, ma_ce_ibat_3, ma_rt_ibat_1, ma_rt_ibat_2, ma_rt_ibat_3) %>% 
  # #   right_join(certmaexp1_pvalue_distinct, by = "compound")
  # # 
  # # write.csv(data_bubble, file = "data_bubble.csv")
  
  # Y_axis stands for relative percentage difference in abundance of all lipid species between cold-treated and control mice. 
  # Each dot represents a lipid species.
  # The dot size indicates significance.
  # The different lipid categories are color-coded.
  bubble <- read_csv("cert_ma_exp1_bubble.csv")
  ggplot(bubble, aes(x=relative_difference, y=category)) +
    # labs(title="The relative percentage difference in detected lipid species between cold-treated and control mice") +
    # theme(axis.title.x=element_text(vjust=1, size=20)) +
    xlab("Relative Difference (%)") +
    ylab("Lipid Category") +
    # theme(axis.title.x=element_text(margin=margin(t=10, r=10, b=10, l=10))) +
    # scale_size_area(breaks=c(0.25, 0.50, 0.75)) +
    scale_fill_manual(breaks=c(0.75,0.50,0.25)) +
    geom_point(aes(size=p_value, color=category), alpha=0.6) +
    scale_color_brewer(palette = "Accent") +
    theme_test() +
    scale_size(name="P_value")
}

.LipidAbundancePlot <- function(MeasurementsFile_standardized_matched_refined = ""){
  
  mergedDf1 <- mergedDf[,c(70,16:30)] %>% as.data.frame()
  for (i in 2:16){
    for (j in 1:nrow(mergedDf1)){
      if (mergedDf1[j,i] != 0){mergedDf1[j,i] <- mergedDf1[j,1]}
      else {mergedDf1[j,i] <- "ND"}
    }
  }
  
  LipidClass <- mergedDf1[,-1] %>% map(table) %>% as.data.frame() 
  rownames(LipidClass) <- LipidClass$adi_epi_1.Var1
  LipidClass <- LipidClass[,seq(2,30,2)] %>% t() %>% as.data.frame()
  LipidClass$stage <- sub(".Freq","",rownames(LipidClass))
  LipidClass <- reshape2::melt(LipidClass, id = "stage")
  LipidClass <- LipidClass[-grep("ND",LipidClass$variable),]
  LipidClass$stage <- factor(LipidClass$stage, levels = unique(LipidClass$stage))
  LipidClass$variable <- factor(LipidClass$variable,levels = unique(LipidClass$variable)[c(4:8,1,3,2)])
  
  pdf("Detected lipid species in each stage.pdf", height = 7, width = 9)
  ggplot(LipidClass, aes(x=stage, y=value, fill=variable)) +
    geom_bar(color="black", stat="identity") +
    theme_bw() +
    labs(x="", y = "Identified number lipid species") +
    scale_fill_manual(values = brewer.pal(n = 8, name = "Dark2")) +
    theme(
      axis.text.x = element_text(color="black", size=14, face="bold",angle = 45, hjust=1),
      axis.title.y = element_text(color="black", size=20, face="bold"),
      axis.text.y = element_text(color="black", size=14, face="bold"),
      legend.title = element_blank(), 
      legend.text = element_text(color="black", size=12, face="bold")
    )
  dev.off()
}

#=========================================================================================
# Run example functions
#=========================================================================================
# standardize_by_standards(MeasurementsFile = "test1.csv",
#                          ISMeasurementsFile = "IS-Measurements-1.xls",
#                          MeasurementsFile_standardized = "test1_standardized.csv",
#                          GroupLabel = list(c(1),c(2:4),c(5:7),c(8:10),c(11:13),c(14:16),
#                                                c(17:19),c(20:22),c(23:25),c(26:28),c(29:31),c(32:34),
#                                                c(35:37),c(38:42),c(43:47),c(48:50),c(51:53)))
# 
# MatchLipidMAPS(MeasurementsFile_standardized = "test1_standardized.csv",
#                IdentificationFile = "Identifications-1.csv",
#                LipidMAPSDB = "lipid_database.csv",
#                OutUnmatchedLipids = T,
#                MeasurementsFile_standardized_matched = "test1_standardized_matched.csv"
# )
# 
# Format_mummichog_input(MeasurementsFile_standardized = "standardized_Measurements1_Converted.csv",
#                        OutputPairWise = T, SelectedSample = c("epi","Peri","ibat"))
# 
# mummichog_Functional_analysis(pktablePath = "standardized_Measurements1_Converted_mummichogInput_G1+G3.csv",
#                               PvalueThreshold = 0.5, Species = c("mmu","hsa")[2], rowNormMet = c("SumNorm","NULL")[2],
#                               EnrichType = c("KEGG","MainClass","SubClass")[c(1,2)])
# 
# Replot_mummichog(mummichog_PEpath = "./")
# 
# Refine_MZ_Identifications(STD_MeasurementsFile = "standardized_Measurements1_Converted.csv",
#                           RefinePlan = c("PlanA","PlanB","PlanC")[3])
# 
# volcano_plotting(PlotFile = "refinedDF_5depots_abbrev_adi_epi+adi_ibat_volcano.csv",
#                  ThresholdFC = 1.5, ThresholdSig = 0.05, TopUpDownShown = c(10,10))
# 
# T_Anavo_plotting(PlotFile = "refinedDF_5depots_abbrev_anova_posthoc.csv",
#                  SigIndex = c("p_value","fdr")[1],ThresholdSig = 0.005,TopShown = 10)
# 
# HeatMapDraw(DataFile = "Epi_Peri.csv",
#             KeywordSelected = "G",
#             PlotSpecific = c("main_class","Ceramides [SP02]"),
#             GroupScheme = c(rep(c("Epi", "Peri"),c(3,3))),
#             GroupAccordingTo = c("abbrev","sub_class","main_class")[2]
# )
# 
# Analyze_Lipids(STD_MeasurementsFile = "refinedDF_5depots.csv",
#               OutputPairWise = T,
#               KeywordSelected="adi",
#               run_RunMetaboAnalystR = T)
