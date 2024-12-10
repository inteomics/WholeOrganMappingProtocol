SampleTableExtSel = readRDS('SampleTableExtSel_v1.RDS')
source('MDA_BLD_plotting_functions.R')
#install.packages("matrixStats")
library("matrixStats")
### General parameters:

FileVer = '10'

PLOT_CHR=F
SIG_LEV=1

SCALE_CUT = 5
SCALE_CUT_CHR = 5
SCALE_CUT_CHR_3D = 5

### Map-specific parameters:

#FilterA
FRAC_SAMP_A = 0.8 
LFC_CRIT_A = 1.5

#FilterB
FRAC_SAMP_B = 1
FRAC_SAMP_NEG_B = 0.2
LFC_CRIT_B = 1.5

#FilterC
FRAC_SAMP_C = 1
FRAC_SAMP_NEG_C = 0.2
LFC_CRIT_C = 1.2

rowZScores <- function(X, na.rm=FALSE){
  (X - rowMeans(X, na.rm=na.rm)) / matrixStats::rowSds(X, na.rm=na.rm)
}

#######################


mapID = "Map19"

  
  ########## data preparation ########### 
  DEGFile = paste0('MDA_BLD1_RNAseq_DESeq2_',mapID,'_v10.txt')
  AnnotColumns = c('SampleType')

  tSampleTableExtSel=SampleTableExtSel[SampleTableExtSel$MapID==mapID,]
  tSampleTableExtSel = rbind(tSampleTableExtSel[tSampleTableExtSel$SampleType=="N",],tSampleTableExtSel[tSampleTableExtSel$SampleType!="N",])
  
  NewID = paste(tSampleTableExtSel$SampleID_exp,tSampleTableExtSel$SampleType,sep="_")
  NewID = gsub('_T_T','_T',NewID)
  tSampleTableExtSel = cbind(FileNameTemp=NewID,tSampleTableExtSel)
  rownames(tSampleTableExtSel) = NewID
   
  ResTab = read.table(DEGFile, header = T, sep="\t", dec=",", stringsAsFactors = F,check.names=FALSE,quote="#")
  tID = ResTab$gene_id
  IDTab = table(ResTab$gene_name)
  DupIDs = IDTab[IDTab>1]
  for(z in 1:length(DupIDs)) {
    ID = DupIDs[z]
    idx = ResTab$gene_name==names(ID)
    ResTab$gene_name[idx] = paste0(ResTab$gene_name[idx],' (',1:ID,')')
  }
  ResTab[,1]= ResTab$gene_name
  ResTab$gene_name = tID
  colnames(ResTab)[7] = 'ensembl_id'
  ResTab$IsMonotonic = NULL
  

  ########## p-value filtering ########### 
  fp1 = ResTab[,paste0('Class1.',mapID,'.NU_LGIN_HGIN_UC_vs_NU.padj')]<SIG_LEV
  fp2 = ResTab[,paste0('Class2.',mapID,'.HGIN_UC_vs_NU.padj')]<SIG_LEV
  fp3 = ResTab[,paste0('Class3.',mapID,'.UC_vs_NU.padj')]<SIG_LEV
  Filt_padj = data.frame(padj.Combined=sum(fp1|fp2|fp3), padj.NU_LGIN_HGIN_UC=sum(fp1), padj.HGIN_UC=sum(fp2), padj.UC=sum(fp3))
  

  ########## LFC filtering ###########
  Samp_NU_LGIN = as.character(tSampleTableExtSel$FileNameTemp[tSampleTableExtSel$Class3=="NU_LGIN"])
  Samp_HGIN = as.character(tSampleTableExtSel$FileNameTemp[tSampleTableExtSel$Class3=="HGIN"])
  Samp_UC = as.character(tSampleTableExtSel$FileNameTemp[tSampleTableExtSel$Class3=="UC"])
  Samp_ALL = as.character(tSampleTableExtSel$FileNameTemp)
  Samp_HGIN_UC = as.character(tSampleTableExtSel$FileNameTemp[tSampleTableExtSel$Class3 %in% c("HGIN","UC") ])
  Samp_NU_LGIN_HGIN = as.character(tSampleTableExtSel$FileNameTemp[tSampleTableExtSel$Class3 %in% c("NU_LGIN","HGIN") ])
  
  # Detect field effect genes, i.e., LG & HG & UC: check the genes with their log2 ratios 
  # of all samples from HG and UC that pass threshold LFC_CRIT_A, and more than FRAC_SAMP_A of the samples from the LG group pass the same threshold;
  Num_NU_LGIN = length(Samp_NU_LGIN)*FRAC_SAMP_A
  Num_HGIN = length(Samp_HGIN)
  Num_UC = length(Samp_UC)
  ffA1 = rowSums(ResTab[,Samp_NU_LGIN]>LFC_CRIT_A)>=Num_NU_LGIN | rowSums(ResTab[,Samp_NU_LGIN]< -LFC_CRIT_A)>=Num_NU_LGIN
  ffA2 = rowSums(ResTab[,Samp_HGIN]>LFC_CRIT_A)>=Num_HGIN | rowSums(ResTab[,Samp_HGIN]< -LFC_CRIT_A)>=Num_HGIN
  ffA3 = rowSums(ResTab[,Samp_UC]>LFC_CRIT_A)>=Num_UC | rowSums(ResTab[,Samp_UC]< -LFC_CRIT_A)>=Num_UC
  
  # Detect genes in HG & UC: check the genes with FRAC_SAMP_B samples from HG and UC that pass threshold LFC_CRIT_B, 
  # and less than FRAC_SAMP_NEG_B of the samples from the LG group pass the same threshold
  Num_NU_LGIN = length(Samp_NU_LGIN)*FRAC_SAMP_NEG_B
  Num_HGIN_UC = length(Samp_HGIN_UC)*FRAC_SAMP_B
  ffB1 = rowSums(ResTab[,Samp_HGIN_UC]>LFC_CRIT_B)>=Num_HGIN_UC | rowSums(ResTab[,Samp_HGIN_UC]< -LFC_CRIT_B)>=Num_HGIN_UC
  ffB2 = rowSums(abs(ResTab[,Samp_NU_LGIN])>LFC_CRIT_B)<Num_NU_LGIN 
  
  #	Detect genes in UC only: check the genes with FRAC_SAMP_C samples from UC that pass threshold LFC_CRIT_C, and less than 
  # FRAC_SAMP_NEG_C of the samples from the LG and HG groups pass the same threshold
  Num_NU_LGIN_HGIN = length(Samp_NU_LGIN_HGIN)*FRAC_SAMP_NEG_C
  Num_UC = length(Samp_UC)*FRAC_SAMP_C
  ffC1 = rowSums(ResTab[,Samp_UC]>LFC_CRIT_C)>=Num_UC | rowSums(ResTab[,Samp_UC]< -LFC_CRIT_C)>=Num_UC
  ffC2 = rowSums(abs(ResTab[,Samp_NU_LGIN_HGIN])>1.5)<Num_NU_LGIN_HGIN 
        

  # Actual data filtration 
  Filter_A = ffA1 & ffA2 & ffA3
  Filter_B = ffB1 & ffB2
  Filter_C = ffC1 & ffC2
  Filter_cmb = Filter_A|Filter_B|Filter_C
  ResTab_monotonic = ResTab[Filter_cmb,]  
        
  
  ########## Table sorting ###########
  # Sort the data based on average expression level
  tAnnotRows=data.frame(Filter_C=as.character(Filter_C),Filter_B=as.character(Filter_B),Filter_A=as.character(Filter_A),stringsAsFactors = F)[Filter_cmb,]
  dataOrderAll = order(rowSums(ResTab_monotonic[,Samp_ALL]))
  ResTab_monotonic = ResTab_monotonic[dataOrderAll,]  
  tAnnotRows = tAnnotRows[dataOrderAll,]
  
  # Additionally sort part of the data that doesn't pass the A filter
  tmp = ResTab_monotonic[tAnnotRows$Filter_A==F,]  
  dataOrderField = order(rowSums(tmp[,Samp_HGIN_UC]))
  tmp2 = tmp[dataOrderField,]
  ResTab_monotonic[tAnnotRows$Filter_A==F,] = tmp2 
  
  
  ########## Create heatmap ###########
  FilePrefix = paste0('Heatmaps/MDA_BLD1_RNAseq_DESeq2_',mapID,'_Heatmap_DEG',
                      '_LFC_',LFC_CRIT_A,';',LFC_CRIT_B,';',LFC_CRIT_C,
                      '_FRC_',FRAC_SAMP_A,';',FRAC_SAMP_B,';',FRAC_SAMP_C,
                      '_FRCN_',FRAC_SAMP_NEG_B,';',FRAC_SAMP_NEG_C,
                      '_ScaleCut_',SCALE_CUT)
  
  plotDEGheatmap(ResTab_monotonic, tSampleTableExtSel, AnnotColumns, FilePrefix, FileVer=FileVer, featureSelection='none', sigLev=1.01, FactorCol="Class3", 
                 skipDEGprefix='control', cutree_cols=NA, cutree_rows=NA, plotWidth=12, plotHeight=10, xlabel="", ylabel="", scaleCut=SCALE_CUT, isLOG=T,
                 noNovel=F,cluster_col=T,cluster_row=F, colorscheme="RdBu",showLFC = F,showAVG = T,LFCcut=0,repelRowLabels=F,zscoreTransform=F,
                 reorderColumns = F, reverseReorderedColumns=F, clusterInGroup=F,Format="PDF")
  
  
  ########## Selection of Top10 genes ###########
  Ngenes = dim(ResTab_monotonic)[1]
  
  idxField = tAnnotRows$Filter_A==T
  idxPlus = rowSums(ResTab_monotonic[,Samp_HGIN_UC]) >0
  excluded = c('HSPD1P6','SCN11A')
  idxExcluded = !ResTab_monotonic$gene_id %in% excluded
  
  fieldMinus = head((1:Ngenes)[idxField & !idxPlus & idxExcluded],10)
  fieldPlus = tail((1:Ngenes)[idxField & idxPlus & idxExcluded],10)
  MidMinus = head((1:Ngenes)[!idxField & !idxPlus & idxExcluded],10)
  MidPlus = tail((1:Ngenes)[!idxField & idxPlus & idxExcluded],10)
  HiMinus = tail((1:Ngenes)[!idxField & !idxPlus & idxExcluded],10)
  HiPlus = head((1:Ngenes)[!idxField & idxPlus & idxExcluded],10)
  
  ResTab_monotonic$isField_Top10 = 1:Ngenes %in% c(fieldMinus,fieldPlus)
  ResTab_monotonic$isHGIN_UC_Top10 = 1:Ngenes %in% c(MidMinus,MidPlus)
  ResTab_monotonic$isUC_Top10 = 1:Ngenes %in% c(HiMinus,HiPlus) 
  
  top10_idx = unique(c(fieldMinus,MidMinus,HiMinus,HiPlus,MidPlus,fieldPlus))
  ResTab_sel_Top10 = ResTab_monotonic[top10_idx,]
  tAnnotRows_Top10 = tAnnotRows[top10_idx,]
  
  plotDEGheatmap(ResTab_sel_Top10, tSampleTableExtSel, AnnotColumns, paste0(FilePrefix,'_Top10'), 
                 FileVer=FileVer, featureSelection='none', sigLev=1.01, FactorCol="Class3", 
                 skipDEGprefix='control', cutree_cols=2, cutree_rows=NA, plotWidth=12, plotHeight=10, xlabel="", ylabel="", scaleCut=SCALE_CUT,isLOG=T,
                 noNovel=F,cluster_col=T,cluster_row=F, colorscheme="RdBu",showLFC = F,showAVG = T,LFCcut=0,repelRowLabels=F,zscoreTransform=F,
                 reorderColumns = F, reverseReorderedColumns=F, clusterInGroup=F,Format="PDF") 
  
  
  #save each of the table to an excel file
  fileExcel = paste0('Heatmaps/MDA_BLD1_RNAseq_DESeq2_',mapID,'_HeatmapsData_DEG',
                              '_LFC_',LFC_CRIT_A,';',LFC_CRIT_B,';',LFC_CRIT_C,
                              '_FRC_',FRAC_SAMP_A,';',FRAC_SAMP_B,';',FRAC_SAMP_C,
                              '_FRCN_',FRAC_SAMP_NEG_B,';',FRAC_SAMP_NEG_C,'_',FileVer,'.xlsx')
  
  library('openxlsx')    
  wb <- createWorkbook()
  addWorksheet(wb = wb, sheetName = 'AllGenes', gridLines = T);  writeDataTable(wb = wb, sheet = 1, x = ResTab, withFilter=F)
  addWorksheet(wb = wb, sheetName = 'MonotonicGenes', gridLines = T);  writeDataTable(wb = wb, sheet = 2, x = ResTab_monotonic, withFilter=F)
  addWorksheet(wb = wb, sheetName = 'MonotonicGenesTop10', gridLines = T);  writeDataTable(wb = wb, sheet = 3, x = ResTab_sel_Top10, withFilter=F)
  saveWorkbook(wb, fileExcel, overwrite = TRUE)


