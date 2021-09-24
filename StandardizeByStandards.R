#=========================================================================================
# This script is for standardizing the identified Lipids concentration against 
# a list of Internal Standards, with specially taken the replicates consistency into consideration

#Version 1.1 created by Houyu Zhang on 2021/09/23
#Issue report on Hughiez047@gmail.com
#Copyright (c) 2021 __CarlosLab@PKU__. All rights reserved.
#=========================================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(readxl)

standardize_by_standards <- function(MeasurementsFile = "", ISMeasurementsFile = "",
                                     standardize_MeasurementsFile = "", ReplicatesLabel = list()){
  
  MF <- read_csv(MeasurementsFile, skip = 2)
  #Count total sample numbers through replicatesLabel
  SampleNumbers <- max(unlist(ReplicatesLabel))
  SampleColStartFrom <- 13 + SampleNumbers
  standardized_MF <- MF
  
  #Since Standards might differ across replicates, we need to find consistent Standards for replicates
  ReplicatesStandards <- list()
  for (item in 1:length(ReplicatesLabel)){
    tmpVC <- c()
    for (sampleNum in ReplicatesLabel[[item]]){
      IS <- read_xls(ISMeasurementsFile, sheet = sampleNum) %>% arrange(`Observed RT (min)`)
      tmpVC <- c(tmpVC, IS$`Component name`)
    }
    tmpVC <- table(tmpVC) %>% as.data.frame()
    #Only store Standards shared in all replicates
    ReplicatesStandards[[item]] <- as.vector(tmpVC[tmpVC$Freq == length(ReplicatesLabel[[item]]),]$tmpVC)
  }

  for (sampleNum in 1:SampleNumbers){
    StandardsUsedHere <- ReplicatesStandards[[which(sapply(ReplicatesLabel , `%in%`, x = sampleNum))]]
    cat("--->>> Start process sample",sampleNum,"...\n")
    IS <- read_xls(ISMeasurementsFile, sheet = sampleNum) %>% arrange(`Observed RT (min)`)
    IS <- IS[IS$`Component name` %in% StandardsUsedHere,]
    
    cat("Refine RT intervals...\n")
    IS$start <- NA; IS$end <- NA
    for (i in 1:(nrow(IS)-1)){
      RT1 = as.numeric(IS[i,"Observed RT (min)"])
      RT2 = as.numeric(IS[i+1,"Observed RT (min)"])
      HalfIntervalLen <- (RT2 - RT1)/2
      IS[i,"end"] = RT1 + HalfIntervalLen
      IS[i+1,"start"] = RT2 - HalfIntervalLen
    }
    IS[1,"start"] = 0; IS[nrow(IS),"end"] = 100
    
    cat("Normalize lipid concentration using corresponding Standard...\n")
    for (LipidNum in 1:nrow(MF)){
      rt = MF[LipidNum,]$`Retention time (min)`
      index <- which(with(IS, start <= rt & end >= rt))
      standardized_MF[LipidNum,sampleNum + SampleColStartFrom] <- as.numeric(MF[LipidNum,sampleNum + SampleColStartFrom])/as.numeric(IS[index,"Detector counts"])
    }
  }
  
  write.csv(file = standardize_MeasurementsFile, x = standardized_MF[,-c(14:(13+SampleNumbers))], row.names = F)
  cat("<<",MeasurementsFile,">> has been standardized using information in <<",ISMeasurementsFile, 
    ">> and wrote into <<",standardize_MeasurementsFile,">>\n")
}

MatchDatabase <- function(LipidMAPSDB = "",standardize_MeasurementsFile = "",
                          IdentificationFile = "", Converted_standardize_MeasurementsFile = ""
                          ){
  standardized_MF <- read_csv(standardize_MeasurementsFile)
  identification <- read_csv(IdentificationFile)
  lipid_database <- read_csv(LipidMAPSDB)
  
  standardized_MF <- janitor::clean_names(standardized_MF)
  identification <- janitor::clean_names(identification)
  
  colnames(identification)[2] <- "compound_id"
  
  lipidMAPS_identification <- lipid_database %>% 
    select(compound_id, category, main_class, sub_class, abbrev) %>% 
    right_join(identification, by="compound_id") 
  
  colnames(standardized_MF)[1] <- "compound"
  Converted_standardized_MF <- lipidMAPS_identification %>% 
    select(compound, category, main_class, sub_class, abbrev) %>% 
    right_join(standardized_MF ,by="compound") 
    # unique() %>% 
    # filter(category != "NA")
  write.csv(x = Converted_standardized_MF, file = Converted_standardize_MeasurementsFile, row.names = F)
}

# standardize_by_standards(MeasurementsFile = "test1.csv",ISMeasurementsFile = "IS-Measurements-1.xls",
#                          standardize_MeasurementsFile = "standardized_test1.csv",
#                          ReplicatesLabel = list(c(1),c(2:4),c(5:7),c(8:10),c(11:13),c(14:16),
#                                                 c(17:19),c(20:22),c(23:25),c(26:28),c(29:31),c(32:34),
#                                                 c(35:37),c(38:42),c(43:47),c(48:50),c(51:53)))
# MatchDatabase(LipidMAPSDB = "lipid_database.csv",
#               standardize_MeasurementsFile = "standardized_test1.csv",
#               IdentificationFile = "Identifications-1.csv",
#               Converted_standardize_MeasurementsFile = "standardized_test1_Converted.csv")

standardize_by_standards(MeasurementsFile = "Measurements-1.csv",ISMeasurementsFile = "IS-Measurements-1.xls",
                         standardize_MeasurementsFile = "standardized_Measurements1.csv",
                         ReplicatesLabel = list(c(1),c(2:4),c(5:7),c(8:10),c(11:13),c(14:16),
                                                c(17:19),c(20:22),c(23:25),c(26:28),c(29:31),c(32:34),
                                                c(35:37),c(38:42),c(43:47),c(48:50),c(51:53)))
MatchDatabase(LipidMAPSDB = "lipid_database.csv",
                          standardize_MeasurementsFile = "standardized_Measurements1.csv",
                          IdentificationFile = "Identifications-1.csv",
                          Converted_standardize_MeasurementsFile = "standardized_Measurements1_Converted.csv")

standardize_by_standards(MeasurementsFile = "Measurements-2.csv",ISMeasurementsFile = "IS-Measurements-2.xls",
                         standardize_MeasurementsFile = "standardized_Measurements2.csv",
                         ReplicatesLabel = list(c(1:3),c(4:6),c(7:9),c(10:12),c(13),c(14),
                                                c(15),c(16),c(17),c(18),c(19),c(20)))
MatchDatabase(LipidMAPSDB = "lipid_database.csv",
              standardize_MeasurementsFile = "standardized_Measurements2.csv",
              IdentificationFile = "Identifications-2.csv",
              Converted_standardize_MeasurementsFile = "standardized_Measurements2_Converted.csv")

# standardize_by_standards(MeasurementsFile = "test2.csv",ISMeasurementsFile = "IS-Measurements-2.xls",
#                          standardize_MeasurementsFile = "standardized_test2.csv",
#                          ReplicatesLabel = list(c(1:3),c(4:6),c(7:9),c(10:12),c(13),c(14),
#                                                 c(15),c(16),c(17),c(18),c(19),c(20)))

