plot_3Dmap <- function(DEGFile, SampleTable, chr, mapID, output.name, Width=6, Height=8, label=FALSE, FillMissingData=F, scaleCut=2, FileType='PNG') {

  #### sample information
  tSamples=SampleTable[SampleTable$Map==mapID,]
  mapInfo = tSamples[,c('SampleType','MapPosX','MapPosY')]
  letters = c('A','B','C','D','E','F','G','H','I','J')
  mapInfo$MapPosX=match(mapInfo$MapPosX,letters)
  rownames(mapInfo) = tSamples$NewID
  
  #### read counts 
  #read the CNV data in log2(CN/2) scale
  #required columns: gene_id, chrom, start, sample1, sample2, ...
  if (is.character(DEGFile)) {
    tData = read.table(DEGFile, header = T, sep="\t", dec=",", stringsAsFactors = F,check.names=FALSE,quote="#")  
  } else {
    tData = DEGFile
  }
  
  
  SampleList = rownames(mapInfo)
  Data = tData[,c('gene_id','chrom','start',SampleList)]
  # change the colnames in order to fit the graphics function
  colnames(Data)[1:3] <- c("FeatureNum", "Chr", "BPstart")
	
  library(BSgenome.Hsapiens.UCSC.hg38)
  mbStart=0
  mbStop=seqlengths(Hsapiens)[[chr]]/1000000
		
	
  if(FileType=="PNG") {
    png(file=output.name,width = Width, height = Height+2, units = "in", res=600)   
  } else if(FileType=="PDF") {
    pdf(file=output.name,width = Width, height = Height+2)
  }
  
  drawMap(mapInfo, rawData=Data, chr=chr, mbStart=mbStart, mbStop=mbStop, plotHeight=Height,scaleCut = scaleCut,labelRows=F)
  if(label == TRUE) {
 	text(1,(Height+0.2), output.name, cex=0.45)
  }
  dev.off()
}




drawUnit <- function(xOffset, yOffset, unitPanels, unitColors){

  polygon(xOffset + unitPanels[,"Panel1X"],
          yOffset + unitPanels[,"Panel1Y"],border=NA,
          col=unitColors[1])
  polygon(xOffset + unitPanels[,"Panel2X"],
          yOffset + unitPanels[,"Panel2Y"],border=NA,
          col=unitColors[2])
  polygon(xOffset + unitPanels[,"Panel3X"],
          yOffset + unitPanels[,"Panel3Y"],border=NA,
          col=unitColors[3])
  
}


buildUnits <-
  function(theta=pi/6, phi=-pi/6,
           xyzDivs=c(0.1,0.2,0.02), xyzGaps=c(0.01,0.02,0.0)){

  xDiv <- xyzDivs[1]
  yDiv <- xyzDivs[2]
  zDiv <- xyzDivs[3]

  xGap <- xyzGaps[1]
  yGap <- xyzGaps[2]
  zGap <- xyzGaps[3]
  
  unitDims <- matrix(0,nrow=3,ncol=2)
  rownames(unitDims) <- c("X","Y","Z")
  colnames(unitDims) <- c("X","Y")
  
  unitDims["X","X"] <- (xDiv - xGap) * cos(theta) * cos(phi)
  unitDims["X","Y"] <- (xDiv - xGap) * sin(theta) * sin(phi)
  #unitDims["Y","X"] <- (yDiv - yGap) * cos(theta) 
  #unitDims["Y","Y"] <- (yDiv - yGap) * sin(theta)
  unitDims["Y","Y"] <- (yDiv - yGap) * cos(theta) 
  unitDims["Y","X"] <- (yDiv - yGap) * sin(theta)
  unitDims["Z","X"] <- 0
  unitDims["Z","Y"] <- -zDiv
  
  unitPanels <- matrix(0,nrow=4,ncol=6)
  colnames(unitPanels) <-
    c("Panel1X", "Panel1Y", "Panel2X", "Panel2Y", "Panel3X", "Panel3Y")
  
  unitPanels[,"Panel1X"] <-
    c(0,
      unitDims["X","X"],
      unitDims["X","X"] + unitDims["Y","X"],
      unitDims["Y","X"])
  unitPanels[,"Panel1Y"] <-
    c(0,
      unitDims["X","Y"],
      unitDims["X","Y"] + unitDims["Y","Y"],
      unitDims["Y","Y"])
  
  unitPanels[,"Panel2X"] <-
    c(0,
      unitDims["X","X"],
      unitDims["X","X"] + unitDims["Z","X"],
      unitDims["Z","X"])
  unitPanels[,"Panel2Y"] <-
    c(0,
      unitDims["X","Y"],
      unitDims["X","Y"] + unitDims["Z","Y"],
      unitDims["Z","Y"])
  
  unitPanels[,"Panel3X"] <-
    c(unitDims["X","X"],
      unitDims["X","X"] + unitDims["Y","X"],
      unitDims["X","X"] + unitDims["Y","X"] + unitDims["Z","X"],
      unitDims["X","X"] + unitDims["Z","X"])
  unitPanels[,"Panel3Y"] <-
    c(unitDims["X","Y"],
      unitDims["X","Y"] + unitDims["Y","Y"],
      unitDims["X","Y"] + unitDims["Y","Y"] + unitDims["Z","Y"],
      unitDims["X","Y"] + unitDims["Z","Y"])
  
  unitDims["X","X"] <- xDiv * cos(theta) * cos(phi)
  unitDims["X","Y"] <- xDiv * sin(theta) * sin(phi)
  unitDims["Y","Y"] <- yDiv * cos(theta) 
  unitDims["Y","X"] <- yDiv * sin(theta)

  list(mDivs = c(xDiv, yDiv, zDiv), unitDims=unitDims, unitPanels=unitPanels)
  
}




drawBaseMap <-
  function(mapInfo, theta=pi/6, phi=-pi/6, xyzDivs=c(0.1,0.2,0.02),
           zoffset=0, sampleIsUsed=NULL, drawBoundingBox=FALSE){
  
  ## DrawBaseMap
  ##
  ## Draws a 2d WOHM with shading indicating severity.

  ######
  ## Define the gray shades to be used
  ######

  gradeGrays <- c(0.9, 0.8, 0.7, 0.6, 0.5, 0.3, 0)
  gradeGrays <- gray(gradeGrays)
  names(gradeGrays) <- c(paste("D", 0:5, sep=""),'T')

  mapInfo[,"Grade"] <- as.character(mapInfo[,"Grade"])
  mapInfo[mapInfo[,"Grade"]=="N","Grade"] <- "D0"
  
  
  ######
  ## Define the size of a base map panel in inches
  ######

  xDiv <- xyzDivs[1]
  yDiv <- xyzDivs[2]

  unitDims <- matrix(0,nrow=2,ncol=2)
  rownames(unitDims) <- c("X","Y")
  colnames(unitDims) <- c("X","Y")
  unitDims["X","X"] <- xDiv * cos(theta) * cos(phi)
  unitDims["X","Y"] <- xDiv * sin(theta) * sin(phi)
  unitDims["Y","Y"] <- yDiv * cos(theta) 
  unitDims["Y","X"] <- yDiv * sin(theta)

  k <- 1/2
  #k <- 1
  unitDims["X","Y"] <- unitDims["X","Y"]*k
  unitDims["Y","Y"] <- unitDims["Y","Y"]*k


  
  mapPanel <- matrix(0,nrow=4,ncol=2)
  colnames(mapPanel) <- c("X","Y")
  mapPanel[,"X"] <-
    c(0,
      unitDims["X","X"],
      unitDims["X","X"] + unitDims["Y","X"],
      unitDims["Y","X"])
  mapPanel[,"Y"] <-
    c(0,
      unitDims["X","Y"],
      unitDims["X","Y"] + unitDims["Y","Y"],
      unitDims["Y","Y"])

  ######
  ## Get the x and y limits of the plotting area
  ######

  mapInfo[,"Row"] <- mapInfo[,"Row"] - min(mapInfo[,"Row"]) + 1
  mapInfo[,"Column"] <- mapInfo[,"Column"] - min(mapInfo[,"Column"]) + 1
  
  maxRow <- max(mapInfo[,"Row"])    
  maxCol <- max(mapInfo[,"Column"]) 

  boundingBox <-
    matrix(c(0,0,maxRow,0,maxRow,maxCol,0,maxCol),
           nrow=4,ncol=2,byrow=TRUE) %*% unitDims

  #print(boundingBox)
  
  ######
  ## Plot the bounding box and add the map units
  ######

  #plot(boundingBox[c(1:4,1),1],
  #     boundingBox[c(1:4,1),2],
  #     type="l",col="red",asp=2/3)

  if(drawBoundingBox){
  
    lines(boundingBox[c(1:4,1),1],
          boundingBox[c(1:4,1),2]+zoffset,
          type="l",col="black")

  }
    
  offsetMat <- (as.matrix(mapInfo[,c("Row","Column")]) - 1) %*% unitDims

  if(is.null(sampleIsUsed)){
    sampleIsUsed = rep(TRUE,dim(mapInfo)[1])
  }
  
  for(i1 in 1:dim(mapInfo)[1]){

    polygon(offsetMat[i1,1] + mapPanel[,"X"],
            offsetMat[i1,2] + mapPanel[,"Y"] + zoffset,
            col=gradeGrays[mapInfo[i1,"Grade"]],
            border=NA
            )
    if(!sampleIsUsed[i1]){
      polygon(offsetMat[i1,1] + mapPanel[,"X"],
              offsetMat[i1,2] + mapPanel[,"Y"] + zoffset,
              col="red", density=30, angle=105,
              border=NA
              )
    }
    #else{
    #  polygon(offsetMat[i1,1] + mapPanel[,"X"],
    #          offsetMat[i1,2] + mapPanel[,"Y"] + zoffset,
    #          col=gradeGrays[mapInfo[i1,"Grade"]],
    #          border=NA, density=30
    #          )      
    #}
  }
  
}



drawMap <- function(mapInfo, rawData, chr, mbStart, mbStop,scaleCut = 3,labelRows=F,
                    theta=pi/6, phi=-pi/6, plotHeight=6, drawBase=TRUE,
                    xyzDivs=c(0.1,0.2,0.02), xyzGaps=c(0.01,0.02,0.0),
                    xOrigin=0, yOrigin=0,
                    positiveRGBColor = matrix(c(0.4,0.4,1),1,3),
                    negativeRGBColor = matrix(c(1,0.4,0.4),1,3)){


  bStart <- mbStart * 1000000
  bStop  <- mbStop  * 1000000
  
  rawData <- rawData[rawData[,"BPstart"] > bStart & rawData[,"BPstart"] < bStop,]
  rawData <- rawData[as.character(rawData[,"Chr"]) == as.character(chr),]
  
  
  ## sort the data rows by BP, with the largest value at the bottom
  rawData <- rawData[order(-rawData[,"BPstart"]),]

  rawInfo    <- rawData[,c("FeatureNum","BPstart","Chr")]
  rawNumbers <- rawData[,-match(c("FeatureNum","BPstart","Chr"),
                                colnames(rawData))]

  #mapInfo <- mapInfo[match(colnames(rawNumbers),rownames(mapInfo)),]
  
  mapInfo[,"Row"] <- mapInfo[,"Row"] - min(mapInfo[,"Row"]) + 1
  mapInfo[,"Column"] <- mapInfo[,"Column"] - min(mapInfo[,"Column"]) + 1
  
  maxRow <- max(mapInfo[,"Row"])    
  maxCol <- max(mapInfo[,"Column"]) 
  

  ## sort the numerical data columns
  rawNumbers <-
    rawNumbers[,order(mapInfo[match(colnames(rawNumbers),
                                    rownames(mapInfo)),"Row"])]
  rawNumbers <-
    rawNumbers[,order(-mapInfo[match(colnames(rawNumbers),
                                     rownames(mapInfo)),"Column"])]
  sampleRow <- mapInfo[match(colnames(rawNumbers),
                             rownames(mapInfo)),"Row"]
  sampleCol <- mapInfo[match(colnames(rawNumbers),
                             rownames(mapInfo)),"Column"]

  
  mbDistance <- mbStop - mbStart
  scalingFactor <- plotHeight / mbDistance

  names(xyzDivs) <- c("x","y","z")
  names(xyzGaps) <- c("x","y","z")

  temp <- buildUnits(theta, phi, xyzDivs, xyzGaps)
  unitDims   <- temp$unitDims
  unitPanels <- temp$unitPanels

  k <- 1/2
  unitDims["X","Y"] <- unitDims["X","Y"]*k
  unitDims["Y","Y"] <- unitDims["Y","Y"]*k
  unitPanels[,"Panel1Y"] <- unitPanels[,"Panel1Y"]*k 
  unitPanels[,"Panel2Y"] <- unitPanels[,"Panel2Y"]*k 
  unitPanels[,"Panel3Y"] <- unitPanels[,"Panel3Y"]*k 
  
  ## ok, now that the data has been ordered, we want to
  ## actually plot the data. Part of the challenge here
  ## is simply to define the y-range over which the plot
  ## extends. To do this, we need to know the size of the
  ## bounding box encompassing the base.
  
  xOffset <-
    (sampleCol - 1) * unitDims["Y","X"] +
      (sampleRow - 1) * unitDims["X","X"]
  yBaseOffset <-
    (sampleCol - 1) * unitDims["Y","Y"] +
      (sampleRow - 1) * unitDims["X","Y"]

    
  boundingBox <-
    matrix(c(0,0,maxRow,0,maxRow,maxCol,0,maxCol),
           nrow=4,ncol=2,byrow=TRUE) %*% unitDims[1:2,]

  bbXmax <- max(boundingBox[,"X"])
  bbXmin <- min(boundingBox[,"X"])
  bbYmax <- max(boundingBox[,"Y"])
  bbYmin <- min(boundingBox[,"Y"])
  minX <- bbXmin
  maxX <- bbXmax
  minY <- bbYmin - (bbYmax - bbYmin)
  maxY <- bbYmax + plotHeight

  
  
  oldPar <- par(mai=c(0,0,0,0),omi=c(0,0,0,0),fin=c(maxX-minX,maxY-minY))
  
  plot(c(minX,maxX),c(minY,maxY),type="n",xaxt="n",yaxt="n",
       xlab="",ylab="",bty="n")
  
  
  
  maxval = round(max(abs(rawNumbers),na.rm=T),2)+0.01
  range = c(-scaleCut,scaleCut)
  maxrange = c(-maxval,maxval)
  
  bk = seq(range[1],range[2],0.01)
  library(RColorBrewer)
  colPalete = c("#3B53A3","white", "#EC2225") 
  #colPalete = rev(brewer.pal(n = 3, name ="RdBu"))
  
  col = colorRampPalette(colPalete,alpha = FALSE)(length(bk)) #rev  
  
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
  
  
  
  
  
  for(i1 in 1:dim(rawNumbers)[1]){
    
    probePosition <- rawInfo[i1,"BPstart"]/1000000 - mbStop
    yOffset <- yBaseOffset - probePosition * scalingFactor 

    #map19 1.1
    #map24 1.25
    if(labelRows) {
      text(1.1,max(yOffset),labels=rawInfo[i1,"FeatureNum"],cex=0.3, adj = 0)
    }
    
    for(i2 in 1:dim(rawNumbers)[2]){
      if (rawNumbers[i1,i2]!=0) {
        tColor = col_mod[as.character(round(bk_mod,2)) == as.character(round(rawNumbers[i1,i2],2))] 
        if (length(tColor)==0) {
          print(c(as.character(bk_mod), round(rawNumbers[i1,i2],2)))
        } else if (tColor!="#FFFFFF") {
          drawUnit(xOffset[i2], yOffset[i2], unitPanels, tColor)
        }
      }
    }
  }

  sampleIsUsed <- rep(FALSE, dim(mapInfo)[1])
  sampleIsUsed[match(colnames(rawNumbers),rownames(mapInfo))] <- TRUE

  drawBaseMap(mapInfo,zoffset=-(bbYmax-bbYmin),sampleIsUsed=sampleIsUsed,drawBoundingBox=T)

  ## add the genomic axis with ticks

  mbTicks <- pretty(c(mbStart,mbStop),n=10)
  mbTickLocs <- (mbTicks - mbStart) * scalingFactor
  tickSpacing <- abs(mbTickLocs[1]-mbTickLocs[2]) 
  nTicks <- length(mbTicks)
  
  lines(c(0,0), c(min(mbTickLocs),max(mbTickLocs)))
  for(i1 in 1:nTicks){
    segments(0, mbTickLocs[i1], 0.05, mbTickLocs[i1])
  }
  #text(0,mbTickLocs[1],labels=paste(mbTicks[nTicks],"Mb"),pos=1,adj=0)
  #text(0,mbTickLocs[nTicks],labels=paste(mbTicks[1],"Mb"),pos=3,adj=0)
  text(0.01,mbTickLocs[1]-0.25*tickSpacing,
       labels=paste(mbTicks[nTicks],"Mb"),adj=c(0,1),offset=1,cex=0.5)
  text(0.01,mbTickLocs[nTicks]+0.25*tickSpacing,
       labels=paste(mbTicks[1],"Mb"),adj=c(0,0),offset=,cex=0.5)
  #print(c(minX, minY, maxX, maxY))
  
}
