#=========================================================================================
# This script is for standardizing the identified Lipids concentration against 
# a list of Internal Standards, with specially taken the replicates consistency into consideration

#Version 1.1 created by Houyu Zhang on 2021/09/28
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
    right_join(standardized_MF ,by="compound") %>% 
    filter(category != "NA") %>%
    arrange(compound)
    # unique()

  write.csv(x = Converted_standardized_MF, file = Converted_standardize_MeasurementsFile, row.names = F)
}

RefineMZIdentifications <- function(Converted_standardize_MeasurementsFile = "standardized_Measurements1_Converted.csv",
                                    RefinePlan = c("Freq","seperate")[1]){
  
  prefix <- tools::file_path_sans_ext(Converted_standardize_MeasurementsFile)
  Converted_standardized_MF <- read_csv(Converted_standardize_MeasurementsFile)
  
  if (RefinePlan == "seperate"){
    
    # Refine based on category information
    Framebone1 <- Converted_standardized_MF %>% select(-main_class,-sub_class,-abbrev) %>% 
      group_by(compound, category) %>% summarise(categoryCounts = n())
    compoundList <- unique(Framebone1$compound)
    
    cns = c(colnames(Converted_standardized_MF)[-c(3:5)])
    Framebone2 = data.frame(matrix(nrow=0, ncol = length(cns)))
    colnames(Framebone2) <- cns
    
    for (i in compoundList){
      # i = compoundList[1]
      tmp <- filter(Converted_standardized_MF, compound == i)
      for (cate in unique(tmp$category)){
        tmp1 <- tmp %>% filter(category == cate) %>% select(-c(3:5))
        Framebone2 <- rbind(Framebone2,tmp1[1,])
      }
    }
    if (identical(Framebone1$category,Framebone2$category)){
      Framebone <- cbind(Framebone1,Framebone2[,-c(1,2)])
    }
    write.csv(x = Framebone, file = paste0(prefix,"_category.csv"), row.names = F)
    
    # Refine based on Main_class information
    Framebone1 <- Converted_standardized_MF %>% select(-category,-sub_class,-abbrev) %>% 
      group_by(compound, main_class) %>% summarise(main_classCounts = n())
    compoundList <- unique(Framebone1$compound)
    
    cns = c(colnames(Converted_standardized_MF)[-c(2,4,5)])
    Framebone2 = data.frame(matrix(nrow=0, ncol = length(cns)))
    colnames(Framebone2) <- cns
    
    for (i in compoundList){
      # i = compoundList[1]
      tmp <- filter(Converted_standardized_MF, compound == i)
      for (cate in unique(tmp$main_class)){
        tmp1 <- tmp %>% filter(main_class == cate) %>% select(-c(2,4,5))
        Framebone2 <- rbind(Framebone2,tmp1[1,])
      }
    }
    Framebone <- cbind(Framebone1,Framebone2[,-c(1,2)])
    write.csv(x = Framebone, file = paste0(prefix,"_main_class.csv"), row.names = F)
    
    # Refine based on sub_class information
    Framebone1 <- Converted_standardized_MF %>% select(-category,-main_class) %>% 
      group_by(compound, sub_class,abbrev) %>% summarise(Counts = n())
    compoundList <- unique(Framebone1$compound)
    
    cns = c(colnames(Converted_standardized_MF)[-c(2,3)])
    Framebone2 = data.frame(matrix(nrow=0, ncol = length(cns)))
    colnames(Framebone2) <- cns
    
    for (i in compoundList){
      # i = compoundList[1]
      tmp <- filter(Converted_standardized_MF, compound == i)
      item <- filter(Framebone1, compound == i)
      for (n in 1:nrow(item)){
          tmp1 <- tmp %>% filter(sub_class == item[n,]$sub_class & abbrev == item[n,]$abbrev) %>% select(-c(2,3))
          Framebone2 <- rbind(Framebone2, tmp1[1,])
      }
    }
    
    if (identical(Framebone1$sub_class,Framebone2$sub_class)){
      Framebone <- cbind(Framebone1,Framebone2[,-c(1:3)]) %>% filter(m_z != "NA")
    }
    write.csv(x = Framebone, file = paste0(prefix,"_sub_class_Abbrv.csv"), row.names = F)
  }
  
  if (RefinePlan == "Freq"){
    compoundList <- unique(Converted_standardized_MF$compound)
    #The most frequent label will be used to represent
    cns = c(colnames(Converted_standardized_MF), 
            c("category_ratio","main_class_ratio","sub_class_ratio","abbrev_ratio"))
    refinedDF = data.frame(matrix(nrow=0, ncol = length(cns))) %>% as_tibble()
    colnames(refinedDF) <- cns
    
    for (i in compoundList){
      tmp <- filter(Converted_standardized_MF, compound == i)
      
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

      #Merge result
      refinedDF <- rbind(refinedDF, tmp[1,])
    }
    
    refinedDF <- refinedDF %>% relocate(category_ratio, .after = category) %>% 
      relocate(main_class_ratio, .after = main_class) %>% 
      relocate(sub_class_ratio, .after = sub_class) %>%
      relocate(abbrev_ratio, .after = abbrev)
    
    write.csv(x = refinedDF, file = paste0(prefix,"_refinedAll.csv"), row.names = F)
  }
}

standardize_by_standards(MeasurementsFile = "test1.csv",ISMeasurementsFile = "IS-Measurements-1.xls",
                         standardize_MeasurementsFile = "standardized_test1.csv",
                         ReplicatesLabel = list(c(1),c(2:4),c(5:7),c(8:10),c(11:13),c(14:16),
                                                c(17:19),c(20:22),c(23:25),c(26:28),c(29:31),c(32:34),
                                                c(35:37),c(38:42),c(43:47),c(48:50),c(51:53)))
MatchDatabase(LipidMAPSDB = "lipid_database.csv",
              standardize_MeasurementsFile = "standardized_test1.csv",
              IdentificationFile = "Identifications-1.csv",
              Converted_standardize_MeasurementsFile = "standardized_test1_Converted.csv")

RefineMZIdentifications(Converted_standardize_MeasurementsFile = "standardized_test1_Converted.csv",
                        RefinePlan = c("Freq","seperate")[1])


#=========================================================================================

standardize_by_standards(MeasurementsFile = "Measurements-1.csv",ISMeasurementsFile = "IS-Measurements-1.xls",
                         standardize_MeasurementsFile = "standardized_Measurements1.csv",
                         ReplicatesLabel = list(c(1),c(2:4),c(5:7),c(8:10),c(11:13),c(14:16),
                                                c(17:19),c(20:22),c(23:25),c(26:28),c(29:31),c(32:34),
                                                c(35:37),c(38:42),c(43:47),c(48:50),c(51:53)))
MatchDatabase(LipidMAPSDB = "lipid_database.csv",
              standardize_MeasurementsFile = "standardized_Measurements1.csv",
              IdentificationFile = "Identifications-1.csv",
              Converted_standardize_MeasurementsFile = "standardized_Measurements1_Converted.csv")

RefineMZIdentifications(Converted_standardize_MeasurementsFile = "standardized_Measurements1_Converted.csv",
                        RefinePlan = c("Freq","seperate")[2])

standardize_by_standards(MeasurementsFile = "Measurements-2.csv",ISMeasurementsFile = "IS-Measurements-2.xls",
                         standardize_MeasurementsFile = "standardized_Measurements2.csv",
                         ReplicatesLabel = list(c(1:3),c(4:6),c(7:9),c(10:12),c(13),c(14),
                                                c(15),c(16),c(17),c(18),c(19),c(20)))
MatchDatabase(LipidMAPSDB = "lipid_database.csv",
              standardize_MeasurementsFile = "standardized_Measurements2.csv",
              IdentificationFile = "Identifications-2.csv",
              Converted_standardize_MeasurementsFile = "standardized_Measurements2_Converted.csv")

RefineMZIdentifications(Converted_standardize_MeasurementsFile = "standardized_Measurements2_Converted.csv",
                        RefinePlan = c("Freq","seperate")[2])



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
  
  
  




