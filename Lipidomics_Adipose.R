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


# ==============================================================================
# 1. Tidy data and run functions
# ==============================================================================
lipid_df <- read_csv("lipid_df.csv") %>% clean_names()
identifications_df <- read_csv("Identificaitons.csv") %>% clean_names()
lipid_database <- read_csv("lipid_database.csv") %>% clean_names()

# ==============================================================================
# 2. Multi-groups analysis
# ==============================================================================
#Make annotations
CompreLipids <- lipid_database %>% 
  select(compound_id, category, main_class, sub_class, abbrev) %>% 
  right_join(identifications_df, by="compound_id") %>%
  select(compound, category, main_class, sub_class, abbrev) %>% 
  right_join(lipid_df, by="compound") %>% 
  select(compound:abbrev,starts_with("adi")) %>% 
  filter(category != "NA") %>%
  arrange(compound)

#A small function for search 
CompreLipids %>% filter(compound == "13.70_856.7474n") %>% data.table::as.data.table()

#Refine lipids based on frequency (The most frequent label will be used to represent)
compoundList <- unique(CompreLipids$compound)
cns = c(colnames(CompreLipids), c("category_ratio","main_class_ratio","sub_class_ratio","abbrev_ratio"))
refinedDF = data.frame(matrix(nrow=0, ncol = length(cns))) %>% as_tibble()
colnames(refinedDF) <- cns
LipidLevels <- c("category","main_class","sub_class","abbrev")

for (i in compoundList){
  tmp <- filter(CompreLipids, compound == i)
  for (level in LipidLevels){
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
  refinedDF <- rbind(refinedDF, tmp[1,])  #Merge result
}
write.csv(file = "refinedDF_5depots.csv", refinedDF, row.names = F)

#Run MetabolanalystR
for (level in LipidLevels){
  slectedrefinedDF <- refinedDF %>% 
    filter(paste0(level,"_ratio") > 0.3) %>%
    select_at(vars(level, starts_with("adi"))) %>%
    group_by(!!as.name(level)) %>% 
    summarise_all(c("sum")) %>% 
    filter(level != "-")
  
  labels <- sub(level,"Label",gsub("_[1-3]","",colnames(slectedrefinedDF)))
  slectedrefinedDF <- rbind(labels, slectedrefinedDF)
  
  write.csv(file = paste0(level,"_5depots.csv"), slectedrefinedDF, row.names = F)
  write.csv(file = paste0(level,"_4depots.csv"), slectedrefinedDF[,1:13], row.names = F)
  RunLipidomics_mutiGroups(pktablePath = paste0(level,"_5depots.csv"))
  RunLipidomics_mutiGroups(pktablePath = paste0(level,"_4depots.csv"))
}

# ==============================================================================
# 2. Pair-wise analysis
# ==============================================================================
#Store matrix based on m/s_RT information mainly for pathway analysis
lipid_df2 <- lipid_df %>% mutate(Compound = gsub("[a-z]|/","",compound)) %>% 
  separate(compound, c("RT", "MZ"), "_") %>% 
  mutate(compound = paste0(MZ,"__",RT)) %>% 
  select(compound,starts_with("adi"))

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

### Pathway enrichment plot
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

# ==============================================================================
# 3. Specific Pair analysis
# ==============================================================================
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
Epi_Peri <- read_csv("refinedDF_5depots.csv") %>%
  select_at(vars(category,main_class,sub_class,abbrev,starts_with("adi_epi"),starts_with("adi_peri"))) %>%
  na.omit() %>%
  filter(rowSums(across(where(is.numeric))) != 0) %>%
  mutate_at(vars(adi_epi_1:adi_peri_3), funs(.*1000/sum(.))) 

write.csv(file = "Epi_Peri.csv", Epi_Peri, row.names = F)

HeatMapDraw <- function(DataFile = "Epi_Peri.csv",
                        GroupAccordingTo = c("abbrev","sub_class")[2],
                        KeywordSelected = "adi"){
  
  cat("Processing",DataFile,"...\n")
  sig_monkey <- read.csv(DataFile) %>% arrange(category) # remove mouse_obese6
  
  # 2.1 Monkey: Heatmaps of significant features
  count_monkey <- sig_monkey %>% group_by(category) %>% count() 
  annotation_col = data.frame(Condition = factor(rep(c("Epi", "Peri"),c(3,3)))) # variable 
  annotation_row = data.frame(LipidCategory = factor(rep(c("FA", "GL", "GP","PK","PR","SL","SP", "ST"), 
                                                         as.vector(count_monkey$n)))) # variable
  ann_colors = list(
    Condition = c(Epi = "#1B9E77", Peri = "#D95F02"),
    LipidCategory = c(FA = "#66C2A5", GL = "#FC8D62", GP = "#8DA0CB", PK = "#E78AC3", 
                      PR = "#A6D854",SL ="#E5C494", SP = "#FFD92F", ST = "#B3B3B3")
  ) # variable
  prefix <- tools::file_path_sans_ext(DataFile)
  pdf(paste0(prefix,"_All_features_heatmap.pdf"), width = 10, height = 10)
  p1 <- ComplexHeatmap::pheatmap(sig_monkey[,5:10], # variable, 2-11 columns are numeric data
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
  pdf(paste0(prefix,"_Cerimides_features_heatmap.pdf"), width = 10, height = 10)
  Forplot <- sig_monkey %>% filter(main_class == "Ceramides [SP02]") %>% select(starts_with("adi"))
  p <- ComplexHeatmap::pheatmap(Forplot, scale = "row", cluster_rows = T, cluster_cols = T,
                                 show_rownames = T, show_colnames = T,
                                 annotation_col = annotation_col, annotation_colors = ann_colors, 
                                 annotation_legend = FALSE, annotation_names_col = FALSE,
                                 annotation_names_row = F, column_split = annotation_col$Condition,
                                 fontsize = 12)
  dev.off()
  # 3. Heatmaps of category/main_class/sub_class/abbrev levels
  # 3.1 Split group into category-level
  sig_monkey_grouped <- sig_monkey %>% group_split(category, .keep = TRUE)
  names(sig_monkey_grouped) <- unique(sig_monkey$category)
  
  for (i in names(sig_monkey_grouped)){
    
    cat("Plotting on",GroupAccordingTo,"level for(",i,")...\n")
    heatmap_data <- sig_monkey_grouped[[i]] %>% 
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
      
      annotation_col = data.frame(Condition = factor(rep(c("Epi", "Peri"),c(3,3))))
      ann_colors = list(Condition = c(Epi = "#1B9E77", Peri = "#D95F02")) # ann_colors of monkeys
      
      pdf(paste0(prefix,"_",GroupAccordingTo,"_",i,"_heatmap.pdf"), width = 10, height = 10)
      p2 <- ComplexHeatmap::pheatmap(heatmapdraw_matrix, 
                                     cellwidth = 15, cellheight = 15,
                                     scale = "row", cluster_rows = T, cluster_cols = T,
                                     show_rownames = TRUE, show_colnames = T,
                                     annotation_col = annotation_col, annotation_colors = ann_colors, 
                                     annotation_legend = FALSE, annotation_names_col = FALSE,
                                     column_split = annotation_col$Condition,
                                     legend = TRUE, fontsize = 12)
      dev.off()
    }
  }
}

HeatMapDraw(DataFile = "Epi_Peri.csv", KeywordSelected = "adi",
            GroupAccordingTo = c("abbrev","sub_class","main_class")[2]
            )

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




