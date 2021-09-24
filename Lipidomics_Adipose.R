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
# 1. Tidy data
# ==============================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)

#Merge all lipid MAPS data
LipidMAPS_df <- list.files(path="Lipid_MAPS/", full.names = TRUE) %>% 
  map_df(~read_csv(., col_types = cols(.default = "c")))

#LM_ID: 46292 unique
#COMMON_NAME: 41190 unique, 4947 "-", 155 others
#SYSTEMATIC_NAME: 40549 unique, 5676 "-", 67 others
#ABBREV: 6447 unique
length(unique(LipidMAPS_df$LM_ID)) == 46292
length(unique(LipidMAPS_df$COMMON_NAME))
LipidMAPS_df$COMMON_NAME <- ifelse(LipidMAPS_df$COMMON_NAME == "-", LipidMAPS_df$LM_ID, LipidMAPS_df$COMMON_NAME)

lipid_df <- read_csv("dataSet/lipid_df.csv")
identifications_df <- read_csv("dataSet/Identificaitons.csv")

summary(lipid_df)
summary(identifications_df)
summary(LipidMAPS_df)

#NOTE: An ID might correspond to many rt_m/z values and vice verse
slectedIdentifications <- identifications_df[match(lipid_df$compound, identifications_df$Compound),]
slectedLipidMAPS <- LipidMAPS_df[match(slectedIdentifications$`Compound ID`, LipidMAPS_df$LM_ID),]

length(unique(slectedIdentifications$`Compound ID`))

if (identical(lipid_df$compound, slectedIdentifications$Compound)){
  mergedDf <- cbind(lipid_df, slectedIdentifications, slectedLipidMAPS) %>% as_tibble()
}
mergedDf <- mergedDf[!is.na(mergedDf$LM_ID),]

colnames(mergedDf)
table(mergedDf$`Compound ID`) %>% as.data.frame() %>% arrange(desc(Freq))
zz <- mergedDf[grep("LMGL03010350",mergedDf$`Compound ID`),]

searchMZ <- function(mz=NA){
  return(mergedDf[grep(mz,mergedDf$compound),c(1,66:73)] %>% data.table::as.data.table())
  }
searchMZ(mz = "13.70_856.7474")
searchMZ(mz = "13.60_898.7808")
searchMZ(mz = "13.67_830.7328")
searchMZ(mz = "11.31_759.5750")
searchMZ(mz = "13.11_888.7627")
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

# ==============================================================================
# 3. Custome analysis 
# ==============================================================================
mergedDf_1 <- mergedDf[,c(1,16:30)]
mergedDf_1 <- cbind(mergedDf_1[,1],sweep(mergedDf_1[,2:16],2,colSums(mergedDf_1[,2:16]),`/`))
labels <- sub("compound","Label",gsub("_[1-3]","",colnames(mergedDf_1)))
mergedDf_1 <- rbind(labels, mergedDf_1)

write.csv(file = "dataSet/mergedDf_5depots.csv", mergedDf_1, row.names = F)
write.csv(file = "dataSet/mergedDf_4depots.csv", mergedDf_1[,1:13], row.names = F)

RunLipidomics <- function(pktablePath = ""){
  
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
  mSet <- PlotNormSummary(mSetObj = mSet, imgName = paste0(prefix,"_lnorm_0_"), format="pdf", dpi = 100, width=NA)
  mSet <- PlotSampleNormSummary(mSetObj = mSet, imgName=paste0(prefix,"_snorm_0_"), format="pdf", dpi = 100, width=NA)
  
  # The fourth step is to perform fold-change analysis
  mSet <- ANOVA.Anal(mSetObj = mSet, nonpar = F, thresh = 0.05, post.hoc = "fisher", all_results = FALSE)
  mSet <- PlotANOVA(mSetObj = mSet, imgName = paste0(prefix,"_anova_0_"), format = "pdf", dpi = 100, width=NA)
  
  # The sixth step is to perform PCA
  mSet <- PCA.Anal(mSetObj = mSet)
  mSet <- PlotPCAPairSummary(mSetObj = mSet, imgName=paste0(prefix,"_pca_pair_0_"), format="pdf", dpi = 100, width=NA, pc.num = 5)
  mSet <- PlotPCAScree(mSetObj = mSet, imgName=paste0(prefix,"_pca_scree_0_"), format="pdf", dpi = 100, width=NA, scree.num = 5)
  mSet <- PlotPCA2DScore(mSetObj=mSet, imgName=paste0(prefix,"_pca_score2d_0_"), format="pdf", 
                         72, width=NA, pcx=1, pcy=2, reg=0.95, show=1, grey.scale=0)
  mSet <- PlotPCALoading(mSetObj = mSet, imgName=paste0(prefix,"_pca_loading_0_"), format="pdf", dpi = 100, width=NA, inx1 = 1,inx2 = 2)
  mSet <- PlotPCABiplot(mSetObj = mSet, imgName=paste0(prefix,"_pca_biplot_0_"), format="pdf", dpi = 100, width=NA, inx1 = 1,inx2 = 2)
  mSet <- PlotPCA3DScoreImg(mSet, imgName = paste0(prefix,"_pca_score3d_0_"), "pdf", dpi =100, width=NA, 1,2,3, angl = 40)
  
  # The seventh step is to perform PLS-DA
  mSet <- PLSR.Anal(mSet, reg=TRUE)
  mSet <- PlotPLSPairSummary(mSet, paste0(prefix,"_pls_pair_0_"), format = "pdf", dpi = 100, width=NA, pc.num = 5)
  mSet <- PlotPLS2DScore(mSet, paste0(prefix,"_pls_score2d_0_"), format = "pdf", dpi = 100, width=NA, 
                         inx1 = 1,inx2 = 2,reg = 0.95,show = 1,grey.scale = 0)
  library(pls)
  mSet <- PlotPLS3DScoreImg(mSet, paste0(prefix,"_pls_score3d_0_"), format = "pdf", dpi = 100, width=NA,1,2,3,40)
  mSet <- PlotPLSLoading(mSet, paste0(prefix,"_pls_loading_0_"), format = "pdf", dpi = 100, width=NA, 1, 2)
  mSet <- PLSDA.CV(mSet, methodName = "L",compNum = 5, choice = "Q2")
  mSet <- PlotPLS.Classification(mSet, paste0(prefix,"_pls_cv_0_"), format = "pdf", dpi = 100, width=NA)
  mSet <- PlotPLS.Imp(mSet, paste0(prefix,"_pls_imp_0_"), format = "pdf", dpi = 100, width=NA, 
                      type = "vip", feat.nm = "Comp. 1", feat.num = 30, color.BW = FALSE)
  for (file in c("row_norm.qs","prenorm.qs","preproc.qs","complete_norm.qs","data_orig.qs",
                 "anova_posthoc.csv","pca_loadings.csv","pca_score.csv","plsda_coef.csv",
                 "plsda_loadings.csv","plsda_score.csv","plsda_vip.csv")){
    file.rename(file, paste0(prefix,"_",file))
  }
}

RunLipidomics(pktablePath = "dataSet/mergedDf_5depots.csv")
RunLipidomics(pktablePath = "dataSet/mergedDf_4depots.csv")













