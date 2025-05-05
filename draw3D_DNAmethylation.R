#################################################################################
#Draw 3D Heatmaps by geographical locations and chromosome-genomic locations 
#for DNA Methylation, RNAseq, and CNV Data
#Created by Ziqiao Wang, PhD (wzqjanet@gmail.com)
##############################
###### load utility functions for 3D plots ###### 
source("drawMap_smooth_color.R")
source("drawBaseMap.R")
source("buildUnits.R")
source("drawUnit.R")
################################
###### load mapInfo files ###### 
# mapInfo
mapInfo <- read.table("Map24Info.csv",header=TRUE,sep=",",row.names="Sample")

####################################################
###### define output file names and locations ###### 

### PB ###
#directory of CNV data, 
cnv.PB.dir <- "~map3D"
cnv.PB.out <- "~map3D"

###### length of each chromosome ######
# positive in red, negative in blue, zero and NA in white
col.scale = function(x, max=0.25) {
  cols = NULL
  for(tmp in x) {
    if(is.na(tmp)) cols = c(cols,rgb(1,1,1)) else {
      #scaling
      tmp = tmp/max
      if(tmp>1) tmp  =  1
      if(-tmp>1) tmp = -1
      #stopifnot(abs(tmp)<=1,'tmp not in scale')
      if(tmp>0) cols=c(cols, rgb(1,1-tmp,1-tmp)) else
        cols = c(cols,rgb(1+tmp, 1+tmp,1)) 
    }
  }
  cols
}

chrLengths <- c(247.3, # chr 1: 247,249,719
                243.0, # chr 2: 242,951,149
                199.6, # chr 3: 199,501,827
                191.3, # chr 4: 191,273,063
                180.9, # chr 5: 180,857,866
                170.9, # chr 6: 170,899,992
                158.9, # chr 7: 158,821,424
                146.3, # chr 8: 146,274,826
                140.3, # chr 9: 140,273,252
                135.4, # chr10: 135,374,737
                134.5, # chr11: 134,452,384
                132.4, # chr12: 132,349,534
                114.2, # chr13: 114,142,980
                106.4, # chr14: 106,368,585
                100.4, # chr15: 100,338,915
                88.9,  # chr16:  88,827,254
                78.8,  # chr17:  78,774,742
                76.2,  # chr18:  76,117,153
                63.9,  # chr19:  63,811,651
                62.5,  # chr20:  62,435,964
                47.0,  # chr21:  46,944,323
                49.7,  # chr22:  49,691,432
                155.0, # chr X: 154,913,754
                57.8)  # chr Y:  57,772,954


names(chrLengths) <- c(as.character(1:22),"X","Y")

#########################################
# FUNCTION: drawCNV
#    Input: chr    - the number of chromosome
#           Start  - the starting base position (Mb) of the plot region (Default = 0.0 Mb).
#           Stop   - the ending base position (Mb) of the plot region (Default = the length of the chromsome).
#           Height - the height of the plot (Default = 8 inches).
#           label  - whether to print the label (Default = FALSE)
#   Output: postscript file
#########################################


drawCNV <- function(Data, chr, info, Start=0.0, Stop=chrLengths[chr], Height=8, label=FALSE,output.name,cnv.PB.out) {
  
  
  #### output file name
  
  colnames(Data)[1:3] <- c("FeatureNum", "Chr", "BPstart")
  
  pdf(file=file.path(cnv.PB.out, output.name),width=6,height=11)
  drawMap(info, Data, chr=chr, mbStart=Start, mbStop=Stop, plotHeight=Height)
  if(label == TRUE) {
    text(1,(Height+0.2), output.name, cex=0.45)
  }
  dev.off()
}

load("heatmapbychromosome_log2ratio.RData")

annEPIC=input_gviz_map24[,c(1,2,3,4)]
anno_bychr=list()
for(i in 1:22){
  anno_bychr[[i]]=annEPIC[which(annEPIC$chromosome==paste0("chr",i)),]
  anno_bychr[[i]]=anno_bychr[[i]][order(anno_bychr[[i]]$start),]
}

load("twomap_3Dplot.RData")

#Test for one example
chr="chr11"
start_chr=(min(map24_3d$start[which(map24_3d$chromosome=="chr11")])-5000)/1000000
end_chr=(max(map24_3d$start[which(map24_3d$chromosome=="chr11")])+5000)/1000000
drawCNV(info=info_map24_3d,map24_3d, chr, Start=start_chr, Stop=end_chr,output.name= paste0("map24_methylation_", chr,".pdf"))


