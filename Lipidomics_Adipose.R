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

devtools::install_github("xia-lab/OptiLCMS", build = TRUE, build_vignettes = F)
devtools::install_github("xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = F)

# ==============================================================================
# 1. Define costume MetaboAnalystR functions
# ==============================================================================
RunLipidomics_mutiGroups <- function(pktablePath = ""){
  
  prefix <- tools::file_path_sans_ext(pktablePath)
  library(MetaboAnalystR)
  remove(mSet)
  # First step is to create the mSet Object, specifying that the data to be uploaded
  # is a peak table ("pktable") and that statistical analysis will be performed ("stat").
  mSet <- InitDataObjects(data.type = "pktable", anal.type = "stat", paired = FALSE)
  # Second step is to read in the filtered peak list
  mSet <- Read.TextData(mSetObj = mSet, filePath = pktablePath, format = "colu", lbl.type = "disc")
  # The third step is to perform data processing using MetaboAnalystR (filtering/normalization)
  mSet <- SanityCheckData(mSetObj = mSet)
  # Perform data processing - Minimum Value Replacing
  mSet <- ReplaceMin(mSetObj = mSet)
  mSet <- SanityCheckData(mSetObj = mSet)
  mSet <- FilterVariable(mSetObj = mSet, filter = "none", qcFilter = "F", rsd = 20)
  mSet <- PreparePrenormData(mSetObj = mSet)
  mSet <- Normalization(mSetObj = mSet, rowNorm = "SumNorm", transNorm = "NULL", scaleNorm = "NULL", ratio=FALSE, ratioNum=20)
  mSet <- PlotNormSummary(mSetObj = mSet, imgName = paste0(prefix,"_lnorm_"), format="pdf", dpi = 100, width=NA)
  mSet <- PlotSampleNormSummary(mSetObj = mSet, imgName=paste0(prefix,"_snorm_"), format="pdf", dpi = 100, width=NA)
  
  # The fourth step is to perform fold-change analysis
  mSet <- ANOVA.Anal(mSetObj = mSet, nonpar = F, thresh = 0.05, post.hoc = "fisher", all_results = FALSE)
  mSet <- PlotANOVA(mSetObj = mSet, imgName = paste0(prefix,"_anova_"), format = "pdf", dpi = 100, width=NA)
  
  # The sixth step is to perform PCA
  mSet <- PCA.Anal(mSetObj = mSet)
  mSet <- PlotPCAPairSummary(mSetObj = mSet, imgName=paste0(prefix,"_pca_pair_"), format="pdf", dpi = 100, width=NA, pc.num = 5)
  mSet <- PlotPCAScree(mSetObj = mSet, imgName=paste0(prefix,"_pca_scree_"), format="pdf", dpi = 100, width=NA, scree.num = 5)
  mSet <- PlotPCA2DScore(mSetObj=mSet, imgName=paste0(prefix,"_pca_score2d_"), format="pdf", 
                         72, width=NA, pcx=1, pcy=2, reg=0.95, show=1, grey.scale=0)
  mSet <- PlotPCALoading(mSetObj = mSet, imgName=paste0(prefix,"_pca_loading_"), format="pdf", dpi = 100, width=NA, inx1 = 1,inx2 = 2)
  mSet <- PlotPCABiplot(mSetObj = mSet, imgName=paste0(prefix,"_pca_biplot_"), format="pdf", dpi = 100, width=NA, inx1 = 1,inx2 = 2)
  mSet <- PlotPCA3DScoreImg(mSet, imgName = paste0(prefix,"_pca_score3d_"), "pdf", dpi =100, width=NA, 1,2,3, angl = 40)
  
  mSet <- PlotSubHeatMap(mSet, paste0(prefix,"_heatmap_"), "pdf", 100, width=NA, dataOpt = "norm", 
                         scaleOpt = "row", smplDist = "euclidean",clstDist = "ward.D",palette = "bwm", 
                         method.nm = "tanova", top.num = 100, viewOpt = "overview", T, T, border = T, F)
  
  # The seventh step is to perform PLS-DA
  mSet <- PLSR.Anal(mSet, reg=TRUE)
  mSet <- PlotPLSPairSummary(mSet, paste0(prefix,"_pls_pair_"), format = "pdf", dpi = 100, width=NA, pc.num = 5)
  mSet <- PlotPLS2DScore(mSet, paste0(prefix,"_pls_score2d_"), format = "pdf", dpi = 100, width=NA, 
                         inx1 = 1,inx2 = 2,reg = 0.95,show = 1,grey.scale = 0)
  library(pls)
  mSet <- PlotPLS3DScoreImg(mSet, paste0(prefix,"_pls_score3d_"), format = "pdf", dpi = 100, width=NA,1,2,3,40)
  mSet <- PlotPLSLoading(mSet, paste0(prefix,"_pls_loading_"), format = "pdf", dpi = 100, width=NA, 1, 2)
  mSet <- PLSDA.CV(mSet, methodName = "L",compNum = 5, choice = "Q2")
  mSet <- PlotPLS.Classification(mSet, paste0(prefix,"_pls_cv_"), format = "pdf", dpi = 100, width=NA)
  mSet <- PlotPLS.Imp(mSet, paste0(prefix,"_pls_imp_"), format = "pdf", dpi = 100, width=NA, 
                      type = "vip", feat.nm = "Comp. 1", feat.num = 30, color.BW = FALSE)
  
  for (file in c("row_norm.qs","prenorm.qs","preproc.qs","complete_norm.qs","data_orig.qs",
                 "anova_posthoc.csv","pca_loadings.csv","pca_score.csv","plsda_coef.csv",
                 "plsda_loadings.csv","plsda_score.csv","plsda_vip.csv")){
    file.rename(file, paste0(prefix,"_",file))
  }
}

RunLipidomics_pairwise_Functional_analysis <- function(pktablePath = ""){
  library(MetaboAnalystR)
  prefix <- tools::file_path_sans_ext(pktablePath)

  FunAnalysis<-InitDataObjects("pktable", "stat", FALSE);
  FunAnalysis<-Read.TextData(FunAnalysis, pktablePath, "colu", "disc")
  FunAnalysis<-SanityCheckData(FunAnalysis)
  FunAnalysis<-ReplaceMin(FunAnalysis);
  FunAnalysis<-FilterVariable(FunAnalysis, "iqr", "F", 25)
  FunAnalysis<-PreparePrenormData(FunAnalysis)
  FunAnalysis<-Normalization(FunAnalysis, "SumNorm", "NULL", "NULL", ratio=FALSE, ratioNum=20)
  FunAnalysis<-Ttests.Anal(FunAnalysis, nonpar = F, threshp = 0.05, paired = FALSE, 
                           equal.var = T, pvalType = "raw", all_results = T)
  FunAnalysis<-Convert2Mummichog(FunAnalysis, rt=TRUE) # mummichog_input: m.z, r.tm p.value, t.score
  
  FunAnalysis<-InitDataObjects("mass_table", "mummichog", FALSE)
  FunAnalysis<-SetPeakFormat(FunAnalysis, "mprt");
  FunAnalysis<-UpdateInstrumentParameters(FunAnalysis, 5, "positive") # Ion Mode: Negative Mode; Mass Tolerance (ppm): 5.0
  FunAnalysis<-Read.PeakListData(FunAnalysis, "mummichog_input_2021-10-11.txt"); # change date
  FunAnalysis<-SanityCheckMummichogData(FunAnalysis) # Retention time tolerance
  FunAnalysis<-SetPeakEnrichMethod(FunAnalysis, "mum", "v2")
  
  FunAnalysis<-SetMummichogPval(FunAnalysis, 0.25)
  # FunAnalysis<-PerformPSEA(FunAnalysis, "mmu_kegg", libVersion = "current", permNum = 100)
  # FunAnalysis<-PlotPeaks2Paths(FunAnalysis, paste0(prefix,"_peaks_to_KEGG_"), "pdf", 72, width=NA)
  
  FunAnalysis<-PerformPSEA(FunAnalysis, "main_lipid_class_mset", libVersion = "current", minLib = 3, permNum = 50)
  FunAnalysis<-PlotPeaks2Paths(FunAnalysis, paste0(prefix,"_peaks_to_MLC_"), "pdf", 72, width=NA)

  # FunAnalysis<-PerformPSEA(FunAnalysis, "sub_lipid_class_mset", libVersion = "current", minLib = 3, permNum = 100)
  # FunAnalysis<-PlotPeaks2Paths(FunAnalysis, paste0(prefix,"_peaks_to_SLC_"), "pdf", 72, width=NA)
  
  for (file in c("mum_raw.qs","complete_norm.qs","row_norm.qs","prenorm.qs",
                 "preproc.qs","data_orig.qs","mummichog_input_2021-10-11.txt",
                 "t_test.csv","t_test_all.csv","scattermum.json",
                 "mummichog_query.json","mummichog_pathway_enrichment.csv",
                 "mum_res.qs","mummichog_matched_compound_all.csv","initial_ecs.qs")){
    file.rename(file, paste0(prefix,"_",file))
  }
  # mSet<-InitDataObjects("mass_table", "mummichog", FALSE)
  # mSet<-SetPeakFormat(mSet, "pvalue")
  # mSet<-UpdateInstrumentParameters(mSet, 5.0, "positive", "yes", 0.02);
  # mSet<-SetRTincluded(mSet, "minutes")
  # mSet<-Read.TextData(mSet, pktablePath, "colu", "disc");
  # mSet<-SanityCheckMummichogData(mSet)
  # mSet<-ReplaceMin(mSet)
  # mSet<-SanityCheckMummichogData(mSet)
  # mSet<-FilterVariable(mSet, "none", "F", 25)
  # mSet<-PreparePrenormData(mSet)
  # mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
  # mSet<-PlotNormSummary(mSet, paste0(prefix,"_norm_"), "pdf", 72, width=NA)
  # mSet<-PlotSampleNormSummary(mSet, paste0(prefix,"_snorm_"), "pdf", 72, width=NA)
  # mSet<-SetPeakEnrichMethod(mSet, "mum", "v2")
  # mSet<-DoPeakConversion(mSet)
  # mSet<-SetMummichogPval(mSet, 0.05)
  # mSet<-PerformPSEA(mSet, "mmu_kegg", "current", 3 , 100)
  # mSet<-PlotPeaks2Paths(mSet, paste0(prefix,"_peaks_to_paths_"), "pdf", 72, width=NA)
  # 
  # mSet<-PerformPSEA(mSet, "main_lipid_class_mset", "current", 3 , 100)
  # mSet<-PlotPeaks2Paths(mSet, paste0(prefix,"_peaks_to_MLC_"), "pdf", 72, width=NA)
}

RunLipidomics_pairwise_Functional_analysis(pktablePath = "Compound_adi_epi+adi_ibat_depots.csv")

# ==============================================================================
# 2. Tidy data and run functions
# ==============================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)

lipid_df <- read_csv("lipid_df.csv")
identifications_df <- read_csv("Identificaitons.csv")
lipid_database <- read_csv("lipid_database.csv")

### 2.1 Pair-wise analysis ####
#Store matrix based on m/s_RT information mainly for pathway analysis
lipid_df2 <- lipid_df %>% mutate(Compound = gsub("[a-z]|/","",compound)) %>% 
  separate(Compound, c("RT", "MZ"), "_") %>% 
  mutate(Compound = paste0(MZ,"__",RT)) %>% 
  select(Compound,starts_with("adi"))

labels <- sub("Compound","Label", gsub("_[1-3]","",colnames(lipid_df2)))
lipid_df2 <- rbind(labels, lipid_df2)
#Output all depots
write.csv(file = "Compound_5depots.csv", lipid_df2, row.names = F, quote = F)
RunLipidomics_mutiGroups(pktablePath = "Compound_5depots.csv")

#Output pair-wise depot
pairwise <- combn(unique(gsub("_[1-3]","",colnames(lipid_df2))[-1]),2)
for (i in (1:ncol(pairwise))[-8]){
  x = pairwise[1,i]
  y = pairwise[2,i]
  tmp <- lipid_df2 %>% select(Compound,starts_with(x),starts_with(y))
  fileName <- paste0("Compound_",x,"+",y,"_depots.csv")
  write.csv(file = fileName, tmp, row.names = F, quote = F)
  if (!file.exists(paste0("Compound_",x,"+",y,"_depots_peaks_to_MLC_dpi72.pdf"))){
    RunLipidomics_pairwise_Functional_analysis(pktablePath = fileName)
  }
}

### Pathway enrichment analysis
PathwayList <- c()
for (file in list.files("./", pattern = "mummichog_pathway_enrichment.csv")){
  PathwayList <- c(PathwayList, read.csv(file)[,1])
}
PathwayList <- unique(PathwayList)

names_ <- c("MainClass","Condition","mLogadjP","EnrFactor")
PathwayInte = data.frame(matrix(nrow=0, ncol = length(names_)))
colnames(PathwayInte) <- names_

for (file in list.files("./", pattern = "mummichog_pathway_enrichment.csv")){
  tmp <- read_csv(file) %>% rename_(MainClass = names(.)[1]) %>%
    mutate(Condition = gsub("_depots_mummichog_pathway_enrichment.csv|Compound_|adi_","",file), 
           mLogadjP = -log10(Gamma), EnrFactor = Hits.sig/Expected) %>%
    select(MainClass,Condition, mLogadjP, EnrFactor) %>%
    filter(mLogadjP > 1.3)
  added <- PathwayList[!PathwayList %in% tmp$MainClass]
  names_ <- c("MainClass","Condition","mLogadjP","EnrFactor")
  tmp2 = data.frame(matrix(nrow=length(added), ncol = length(names_)))
  colnames(tmp2) <- names_
  tmp2$MainClass <- added
  tmp2$Condition <- gsub("_depots_mummichog_pathway_enrichment.csv|Compound_|adi_","",file)
  tmp <- tmp %>% add_row(tmp2)
  PathwayInte <- rbind(PathwayInte, tmp)
}
PathwayInte$Condition <- factor(PathwayInte$Condition, 
                                   levels = unique(PathwayInte$Condition)[c(3,6,8,1,5,9,2,7,4)])

pdf("MainClass_enrichemnt_P0.05.pdf", width = 8, height = 8)
ggplot(PathwayInte,aes(x = Condition, y = MainClass, size = EnrFactor, color = mLogadjP)) +
  geom_point() + 
  scale_color_gradient(low="white", high="red") +
  theme_bw() +
  theme(
    axis.text.x = element_text(color="black", size=12,face="bold", angle=70, hjust=1),
    axis.text.y = element_text(color="black", size=12, face="bold"),
    legend.title = element_text(color="black", size=14, face="bold"), 
    legend.text = element_text(color="black", size=12, face="bold")
  )
dev.off()

### 2.2 Multi-groups analysis ####
#NOTE: An ID might correspond to many rt_m/z values and vice verse
colnames(identifications_df)[2] <- "compound_id"
colnames(lipid_df)[1] <- "Compound"

CompreLipids <- lipid_database %>% 
  select(compound_id, category, main_class, sub_class, abbrev) %>% 
  right_join(identifications_df, by="compound_id") %>%
  select(Compound, category, main_class, sub_class, abbrev) %>% 
  right_join(lipid_df, by="Compound") %>% 
  select(Compound:abbrev,starts_with("adi")) %>% 
  filter(category != "NA") %>%
  arrange(Compound)

CompreLipids %>% filter(Compound == "13.60_898.7808m/z") %>% data.table::as.data.table()

compoundList <- unique(CompreLipids$Compound)
#The most frequent label will be used to represent
cns = c(colnames(CompreLipids), c("category_ratio","main_class_ratio","sub_class_ratio","abbrev_ratio"))
refinedDF = data.frame(matrix(nrow=0, ncol = length(cns))) %>% as_tibble()
colnames(refinedDF) <- cns

for (i in compoundList){
  tmp <- filter(CompreLipids, Compound == i)
  
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

for (level in c("category","main_class","sub_class","abbrev")){
  slectedrefinedDF <- refinedDF %>% 
    filter(paste0(level,"_ratio") > 0.3) %>%
    select_at(vars(level, starts_with("adi"))) %>%
    group_by(!!as.name(level)) %>% 
    summarise_all(c("sum")) %>% 
    filter(level != "-") %>% 
    mutate_at(vars(adi_epi_1:adi_ibat_3), funs(./sum(.)))
  
  labels <- sub(level,"Label",gsub("_[1-3]","",colnames(slectedrefinedDF)))
  slectedrefinedDF <- rbind(labels, slectedrefinedDF)
  
  write.csv(file = paste0(level,"_5depots.csv"), slectedrefinedDF, row.names = F)
  write.csv(file = paste0(level,"_4depots.csv"), slectedrefinedDF[,1:13], row.names = F)
  RunLipidomics_mutiGroups(pktablePath = paste0(level,"_5depots.csv"))
  RunLipidomics_mutiGroups(pktablePath = paste0(level,"_4depots.csv"))
}

# ==============================================================================
# 2. Lipid abundance analysis
# ==============================================================================
library(RColorBrewer)
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




