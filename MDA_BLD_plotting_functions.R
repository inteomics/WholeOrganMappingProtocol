plotDEGheatmap = function(DEGFile, SampleTable, FactorCol, AnnotColumns, FilePrefix, FileVer, featureSelection='DEG', sigLev=0.05, skipDEGprefix="#",
                          cutree_cols=1, cutree_rows=1, plotWidth=12, plotHeight=8, xlabel="Patients", ylabel="Genes",isLOG=F,show_rownames='auto',show_colnames=T, zscoreTransform=T,
                          scaleCut=3,LFCcut=0,noNovel=T,cluster_col=T, cluster_row=T,colorscheme="RdYlBu",showLFC = T,showAVG = T,maxLFC = 4, repelRowLabels=F,
                          reorderColumns = F, reverseReorderedColumns=F, clusterInGroup=F, AnnotationRows=NULL,AnnotationRows_col=NULL,scaleSymetry=T,Format='PNG') {
  #featureSelection='DEG'; sigLev=0.05; skipDEGprefix="#"; cutree_cols=1; cutree_rows=1; plotWidth=12; plotHeight=8; xlabel="Patients"; ylabel="Genes";isLOG=F;show_rownames='auto';show_colnames=T; zscoreTransform=T; scaleCut=3;LFCcut=0;noNovel=T;cluster_col=T; cluster_row=T;colorscheme="RdYlBu";showLFC = T;showAVG = T;maxLFC = 4; repelRowLabels=F; reorderColumns = F; reverseReorderedColumns=F; clusterInGroup=F; AnnotationRows=NULL;AnnotationRows_col=NULL;scaleSymetry=T;Format='PNG'; FilePrefix='TEST'; FileVer=1
  require('grid')
  
  #read the DEG table
  if (is.character(DEGFile)) {
    ResTab = read.table(DEGFile, header = T, sep="\t", dec=",", stringsAsFactors = F,check.names=FALSE,quote="#")
  } else {
    ResTab = DEGFile
  }
  
  
  ###### remove novel miRNAs
  if (noNovel) {
    ResTab = ResTab[!grepl('novel_miRNA',ResTab$ID),]
  }
  
  expressionCols = as.character(SampleTable[!is.na(SampleTable[,FactorCol]),1])
  expressionCols = as.character(expressionCols[expressionCols %in% colnames(ResTab)])
  
  #get the log2 expression levels for each sample and DEG miRNA
  if (isLOG) {
    ResTabExp = ResTab[,expressionCols]
  } else {
    ResTabExp = log2(ResTab[,expressionCols] + 1)
  }
  rownames(ResTabExp) = rownames(ResTab) = ResTab[,1]
  
  ######### remove rows with NA only (caused by Combat) ##########
  #ResTab = ResTab[!is.na(rowSums(ResTabExp)),]
  #ResTabExp = ResTabExp[!is.na(rowSums(ResTabExp)),]
  
  
  #row annotations
  AVGreads = rowMeans(ResTab[,expressionCols])
  if (!isLOG) {
    AVGreads = log2(AVGreads+1)
  }
  annotation_row = NULL
  
  if (showAVG) {
    avgColsID=grepl('AverageExprs(norm)',colnames(ResTab),fixed=T) | grepl('AverageExprs (norm counts)',colnames(ResTab),fixed=T)
    avgCols = avgColsID & !grepl(skipDEGprefix,colnames(ResTab))
    rawdata = ResTab[,avgCols,drop=F]
    #transdata = adply(rawdata,1,function(x) { y=rescale(as.numeric(x), to = c(-1, 1)) y=})
    if (ncol(rawdata)>1) {
      avgData = data.frame(rowZScores(rawdata),stringsAsFactors = F,check.names=F)
    } else {
      avgData = rawdata
    }
    colnames(avgData) = gsub('.AverageExprs(norm)','',colnames(avgData),fixed=T)
    rownames(avgData) = rownames(ResTab)
    annotation_row = avgData
  }
  
  
  lfcColsID=grepl('.logFC',colnames(ResTab))
  lfcCols = lfcColsID & !grepl(skipDEGprefix,colnames(ResTab))
  if (showLFC) {
    logFCdata = data.frame(ResTab[,lfcCols],stringsAsFactors = F,check.names=F)
    if (dim(logFCdata)[2]==1) {
      colnames(logFCdata) = 'logFC'
    } else {
      colnames(logFCdata) = paste0('logFC: ',gsub('.logFC','',colnames(ResTab)[lfcCols]))
    }
    
    rownames(logFCdata) = rownames(ResTab)
    logFCdata[logFCdata>maxLFC] = maxLFC
    logFCdata[logFCdata< -maxLFC] = -maxLFC
    logFCdata_conv = round(logFCdata*10)/10
    for (i in 1:dim(logFCdata_conv)[2]) {
      logFCdata_conv[,i] = (logFCdata_conv[,i])
    }
    
    if (showAVG) {
      annotation_row = data.frame(avgData,logFCdata_conv,'log2 avg reads'=AVGreads, stringsAsFactors = F,check.names=FALSE)
    } else {
      annotation_row = data.frame(logFCdata_conv,'log2 avg reads'=AVGreads, stringsAsFactors = F,check.names=FALSE)
    }
  }
  
  #add custom annotation rows
  if(!is.null(AnnotationRows)) {
    if(!is.null(annotation_row)) {
      annotation_row = cbind(annotation_row,AnnotationRows)
    } else {
      annotation_row = AnnotationRows
    }
  }
  #annotation_row$chr=ResTab$chrom
  
  
  #get factor labels from the samples table
  annotation_col = data.frame(SampleTable[colnames(ResTabExp),unique(c(FactorCol,AnnotColumns))], stringsAsFactors = F,check.names=FALSE)
  rownames(annotation_col) = colnames(ResTabExp)
  colnames(annotation_col) = unique(c(FactorCol,AnnotColumns))
  
  
  library("RColorBrewer")
  #adjust the color hue
  gg_color_hue <- function(n) {
    #hues = seq(15, 375, length = n + 1)
    #hcl(h = hues, l = 65, c = 100)[1:n]
    if (n>9) {
      myPal <- colorRampPalette(brewer.pal(12, "Set3"))
      cols <- myPal(n)
    }
    else if (n>2) {
      cols = brewer.pal(n, "Greys")
    } else if (n==2) {
      cols = c("#F0F0F0","#636363")
    } else {
      cols = 'black'
    }
    return(cols)
  }
  
  
  
  
  ##set the order of factors (control first if it exists)
  ControlFactor = "N" #"Control"
  ancolors = list()
  for (i in 1:dim(annotation_col)[2]) {
    # categorial variables
    if(!is.numeric(annotation_col[,i])) {
      factorLev = as.character(unique(annotation_col[,i]))
      if (sum(factorLev==ControlFactor,na.rm = T)>0) {
        factorLev = sort(factorLev[factorLev!=ControlFactor])
        annotation_col[,i]=factor(annotation_col[,i],levels=c(ControlFactor,factorLev))
        #assign color to each factor
        ancolors[i] = list(gg_color_hue(length(factorLev)+1))
        names(ancolors[[i]])=c(ControlFactor,factorLev)
      } else {
        #assign color to each factor
        ancolors[i] = list(gg_color_hue(length(factorLev)))
        names(ancolors[[i]])=factorLev
      }
      #numeric variables
    } else {
      factorLev = annotation_col[,i]
      labs = as.character(c(0,0.05,max(factorLev)))
      ancolors[[i]] = c("red","white","black")
      names(ancolors[[i]])=labs
    }
  }
  names(ancolors) = unique(c(FactorCol,AnnotColumns))
  
  
  library("RColorBrewer")
  if (showLFC) {
    FClab = as.character(round(seq(-maxLFC,maxLFC,0.1)*10)/10)
    FCcolors = colorRampPalette(c("#006837","#FFFFFF","#A50026"))(length(FClab))
    names(FCcolors) = FClab
    ancolors[colnames(logFCdata)] = list(FCcolors)
  }
  if (showAVG) {
    AVGlab = (round(seq(-2,2,0.01)*100)/100)
    AVGcolors = colorRampPalette(c("#006837","#FFFFFF","#A50026"))(length(AVGlab))
    names(AVGcolors) = AVGlab
    ancolors[colnames(avgData)] = list(c('white','#267029'))# list(AVGcolors)
  }
  ancolors$'log2 avg reads' = c('0'='white','15'='#4070B0')
  
  #add colors for custom row annotation
  if (!is.null(AnnotationRows_col)) {
    ancolors = c(ancolors,AnnotationRows_col)
  }
  
  
  ## filter features based on DEG or ANOVA
  if (featureSelection=='DEG') {
    
    #select DEG columns skip columns only if more than one exists
    degColsID = grepl('.padj',colnames(ResTab),fixed=T)
    if (sum(degColsID)>1) {
      degCols = degColsID & !grepl(skipDEGprefix,colnames(ResTab)) & !grepl('ANOVA',colnames(ResTab))
      print(paste0('Skipping ',sum(degColsID & grepl(skipDEGprefix,colnames(ResTab))),' columns with ',skipDEGprefix,' prefix'))
    } else {
      degCols = degColsID
    }
    
    #select features with DEG status
    if (sum(degCols)>1) {
      ResTabDEG = rowSums(!is.na(ResTab[,degCols]) & ResTab[,degCols]<sigLev) > 0  & rowSums(abs(ResTab[,lfcCols])>LFCcut) > 0
    } else {
      ResTabDEG = !is.na(ResTab[,degCols]) & ResTab[,degCols]<sigLev  & abs(ResTab[,lfcCols])>LFCcut
    }
    #PlotData = t(rowZScores(t(ResTabExp[ResTabDEG,])))
    #z-score transformation
    if(zscoreTransform) {
      PlotData = rowZScores(ResTabExp[ResTabDEG,])
    } else {
      PlotData = ResTabExp[ResTabDEG,]
    }
    PlotFileName = paste0(FilePrefix,'_DEG_',sigLev,'_v',FileVer)
    
  } else if (featureSelection=='ANOVA') {
    
    ResTab_data = ResTabExp
    #z-score transformation
    if(zscoreTransform) {
      ResTab_data_zs = t(rowZScores(t(ResTab_data)))
      ResTab_data_zs = rowZScores(ResTab_data_zs)
    } else {
      ResTab_data_zs = ResTab_data
    }
    rownames(ResTab_data_zs)=rownames(ResTabExp)
    PlotData = ResTab_data_zs[ResTab$ANOVA.padj<sigLev & !is.na(ResTab$ANOVA.padj),]
    
    #ResTabAnova = makeANOVAexpr(ResTabExp,SampleTable,FactorCol)
    #PlotData = ResTab_data_zs[ResTabAnova$ANOVA_padj<sigLev,]
    
    PlotFileName = paste0(FilePrefix,'_ANOVA_',sigLev,'_v',FileVer)
  } else {
    
    #z-score transformation
    if(zscoreTransform) {
      PlotData = rowZScores(ResTabExp)
    } else {
      PlotData = ResTabExp
    }
    
    PlotFileName = paste0(FilePrefix,'_v',FileVer)
  }
  
  
  #automatic rownames on/off
  if (show_rownames=='auto') {
    if (dim(PlotData)[1]>300 & repelRowLabels==F) {
      show_rownames=F
    } else {
      show_rownames=T
    }
  }
  
  if (show_rownames==T & dim(PlotData)[1]>80 & repelRowLabels==F) {
    fontsize_row = 4
  } else if (repelRowLabels) {
    fontsize_row = 12
  } else {
    fontsize_row = 8
  }
  
  PlotData[is.na(PlotData)] = 0

  if (reorderColumns) {
    hc <- hclust(dist(t(PlotData)))
    dd <- as.dendrogram(hc)
    FactorVal = annotation_col[colnames(PlotData),FactorCol]
    labels = as.numeric(factor(FactorVal,levels = unique(FactorVal)))
    dd.reorder <- reorder(dd, labels)
    if(reverseReorderedColumns){
      dd.reorder <- rev(dd.reorder)
    }
    cluster_col = as.hclust(dd.reorder)
  }
  
  
  if(clusterInGroup) {
    FactorVal = annotation_col[colnames(PlotData),FactorCol]
    UniqFactors = unique(FactorVal)
    OrderLabels = 1:length(FactorVal)
    for (i in 1:length(UniqFactors)) {
      idx = FactorVal==UniqFactors[i]
      if (sum(idx)>1) {
        hc <- hclust(dist(t(PlotData[,idx])))
        dd <- as.dendrogram(hc)
        
        PlotFileName_dend = paste0(PlotFileName,'_DEND_',UniqFactors[i],'.pdf')
        pdf(PlotFileName_dend)
        plot(dd)
        dev.off()
        
        tOrder = order.dendrogram(dd)
      } else {
        tOrder = 1
      }
      OrderLabels[idx] = OrderLabels[idx][tOrder]
    }
    PlotData = PlotData[,OrderLabels]
  }
  
  #create and save the plot
  if (Format=="PNG") {
    PlotFileName = paste0(PlotFileName,'.png')
    png(PlotFileName,width = plotWidth, height = plotHeight, units = "in", res=300)
  } else {
    PlotFileName = paste0(PlotFileName,'.pdf')
    pdf(PlotFileName,width = plotWidth, height = plotHeight)
  }
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
  symetricalHeatmap(PlotData, scaleCut=scaleCut, colorscheme=colorscheme, cluster_col=cluster_col, cluster_row=cluster_row, show_rownames=show_rownames, show_colnames = show_colnames,
                    cutree_cols=cutree_cols, cutree_rows=cutree_rows, repelRowLabels=repelRowLabels,annotation_col=annotation_col, annotation_row=annotation_row,
                    annotation_colors=ancolors,fontsize_row=fontsize_row,scaleSymetry=scaleSymetry)
  setHook("grid.newpage", NULL, "replace")
  grid.text(ylabel, x=1-0.15, rot=270, gp=gpar(fontsize=16))
  grid.text(xlabel, x=0.35, y=0.03,  gp=gpar(fontsize=16))
  dev.off()
  
}




symetricalHeatmap = function(Data,scaleCut=NULL,colorscheme="RdYlBu",cluster_col=F,cluster_row=F,show_rownames=T,show_colnames=T,cutree_rows=NA,cutree_cols=NA,annotation_col=NA,
                             annotation_row=NA,annotation_colors=NA,display_numbers=F,fontsize_row=10,fontsize_col=10,repelRowLabels=F,scaleSymetry=T) {
  
  #Plots heatmap with colorscale centered arround zero
  #  repelRowLabels = T - plot all row labels repeled
  #                 = F - do not use row repel
  #                 = character_vector - plot only selected labels - repeled
  
  
  library("RColorBrewer")
  library("pheatmap")
  library("dendsort")
  
  maxval = max(abs(Data),na.rm=T)
  
  if(scaleSymetry) {
    range = c(-scaleCut,scaleCut)
    maxrange = c(-maxval,maxval)
  } else {
    range = maxrange = c(min(Data),max(Data))
  }
  
  
  bk = seq(range[1],range[2],scaleCut/20)
  if (colorscheme=="RdBu2") {
    colPalete = c('#535CA8','white','#E44E4D')
  } else if (colorscheme=="Jet") {
    colPalete = c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
    #jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    #colPalete (volcano, color = jet.colors, asp = 1, nlevels=100)
  } else if (colorscheme=="RdBu") {
    colPalete = rev(brewer.pal(n = 8, name ="RdBu"))
  } else {
    colPalete = rev(brewer.pal(n = brewer.pal.info[colorscheme,'maxcolors'], name =colorscheme))
  }
  col = colorRampPalette(colPalete)(length(bk)) #rev
  
  if (scaleCut<maxval) {
    bk_min = seq(maxrange[1]-0.01,range[1]-0.01,0.01)
    bk_max = seq(range[2]+0.01,maxrange[2]+0.01,0.01)
    bk_mod=c(bk_min,bk,bk_max)
    
    col_min = rep(col[1],length(bk_min))
    col_max = rep(col[length(col)],length(bk_max))
    col_mod = c(col_min,col,col_max)
  } else {
    bk_mod = bk
    col_mod = col
  }

  
  if (cluster_col=="sort") {
    mat_cluster_cols <- hclust(dist(t(Data)))
    sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
    mat_cluster_cols <- sort_hclust(mat_cluster_cols)
  } else {
    mat_cluster_cols=cluster_col
  }
  
  heat = pheatmap(Data, cluster_col=mat_cluster_cols, cluster_row=cluster_row, show_rownames=show_rownames, show_colnames=show_colnames, breaks=bk_mod, color=col_mod,
                  cutree_rows=cutree_rows, annotation_col=annotation_col, annotation_row=annotation_row, cutree_cols=cutree_cols,
                  annotation_colors=annotation_colors, display_numbers=display_numbers,fontsize_row=fontsize_row,fontsize_col=fontsize_col,border_color = "grey60") #
  
  #add repeled row labels
  if (is.logical(repelRowLabels)) {
    if (repelRowLabels) {
      #keep all labels but add separation
      add.flag(heat,kept.labels=rownames(Data),repel.degree = 1.5)
    }
  } else {
    #use only selected row labels
    add.flag(heat,kept.labels=repelRowLabels,repel.degree = 1.5)
  }
  
  return(heat)
}
