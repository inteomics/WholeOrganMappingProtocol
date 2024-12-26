drawMap <- function(mapInfo, rawData, chr, mbStart, mbStop,
                    theta=pi/6, phi=-pi/6, plotHeight=6, drawBase=TRUE,
                    xyzDivs=c(0.1,0.2,0.02), xyzGaps=c(0.01,0.02,0.0),
                    xOrigin=0, yOrigin=0,
                    positiveRGBColor = matrix(c(0.4,0.4,1),1,3),
                    negativeRGBColor = matrix(c(1,0.4,0.4),1,3)){

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
                57.8)  # chr Y:  57,772,954, 
  names(chrLengths) <- c(as.character(1:22),"X","Y")

  ##rawData <- read.table(rawDataFile,
  ##                      header=TRUE,sep="\t")

  bStart <- mbStart * 1000000
  bStop  <- mbStop  * 1000000
  
  rawData <- rawData[rawData[,"BPstart"] > bStart,]
  ##print(dim(rawData))
  rawData <- rawData[rawData[,"BPstart"] < bStop,]
  ##print(dim(rawData))
  ##print(rawData[1,])
  rawData <- rawData[as.character(rawData[,"Chr"]) == as.character(chr),]
  ##print(dim(rawData))

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
  
#  print(dim(rawInfo))
#  print(dim(rawNumbers))
  
  ## sort the numerical data columns

  rawNumbers <-
    rawNumbers[,order(mapInfo[match(colnames(rawNumbers),
                                    rownames(mapInfo)),"Row"])]
#  print(rawNumbers[1,])
  rawNumbers <-
    rawNumbers[,order(-mapInfo[match(colnames(rawNumbers),
                                     rownames(mapInfo)),"Column"])]
#  print(rawNumbers[1,])

  sampleRow <- mapInfo[match(colnames(rawNumbers),
                             rownames(mapInfo)),"Row"]
  sampleCol <- mapInfo[match(colnames(rawNumbers),
                             rownames(mapInfo)),"Column"]

#  print(sampleRow)
#  print(sampleCol)

  ##print(dim(rawData))
  ##print(rawData[1,])
  
  mbDistance <- mbStop - mbStart
  scalingFactor <- plotHeight / mbDistance

  names(xyzDivs) <- c("x","y","z")
  names(xyzGaps) <- c("x","y","z")
  
  positiveColors <- c(rgb(positiveRGBColor),
                      rgb(positiveRGBColor - 0.1),
                      rgb(positiveRGBColor - 0.1))
  negativeColors <- c(rgb(negativeRGBColor),
                      rgb(negativeRGBColor - 0.1),
                      rgb(negativeRGBColor - 0.1))

					  
  positiveColors2 <- c(rgb(positiveRGBColor*0.4),
                      rgb((positiveRGBColor - 0.1)*0.4),
                      rgb((positiveRGBColor - 0.1)*0.4))
  negativeColors2 <- c(rgb(negativeRGBColor*0.4),
                      rgb((negativeRGBColor - 0.1)*0.4),
                      rgb((negativeRGBColor - 0.1)*0.4))

  ##print(positiveColors)
  ##print(negativeColors)
  
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
#  +
#        plotHeight - mbStart * scalingFactor
#        yMax - mbStart * scalingFactor

#  print(yBaseOffset)
  
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

  #plot(c(minX,maxX),c(minY,maxY),type="n",xaxt="n",yaxt="n",
  #     xlab="",ylab="")

  # starting off, we have mai = c(1.02, 0.82, 0.82, 0.42)
  # fin=c(7,7), omi= c(0,0,0,0) -- the last one seems to have
  # worked. 
  
  oldPar <- par(mai=c(0,0,0,0),omi=c(0,0,0,0),fin=c(maxX-minX,maxY-minY))
  
  plot(c(minX,maxX),c(minY,maxY),type="n",xaxt="n",yaxt="n",
       xlab="",ylab="",bty="n")
  
  for(i1 in 1:dim(rawNumbers)[1]){
    
    probePosition <- rawInfo[i1,"BPstart"]/1000000 - mbStop
    yOffset <- yBaseOffset - probePosition * scalingFactor 

    #print(yOffset)
    
    for(i2 in 1:dim(rawNumbers)[2]){
      
      drawUnit(xOffset[i2], yOffset[i2], unitPanels, col.scale(rawNumbers[i1,i2]))
 #     if(rawNumbers[i1,i2] == 2){
 #      drawUnit(xOffset[i2], yOffset[i2], unitPanels, positiveColors2)
 #     } else if (rawNumbers[i1,i2] == 1) {
 #	    drawUnit(xOffset[i2], yOffset[i2], unitPanels, positiveColors)
#	  }
#      if(rawNumbers[i1,i2]  == -2){
#        drawUnit(xOffset[i2], yOffset[i2], unitPanels, negativeColors2)
#      } else if (rawNumbers[i1,i2] == -1) {
#	    drawUnit(xOffset[i2], yOffset[i2], unitPanels, negativeColors)
#	  }
      
    }
  }

  sampleIsUsed <- rep(FALSE, dim(mapInfo)[1])
  sampleIsUsed[match(colnames(rawNumbers),rownames(mapInfo))] <- TRUE

  drawBaseMap(mapInfo,zoffset=-(bbYmax-bbYmin),sampleIsUsed=sampleIsUsed)

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
