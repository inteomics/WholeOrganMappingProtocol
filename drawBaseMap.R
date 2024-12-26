drawBaseMap <-
  function(mapInfo, theta=pi/6, phi=-pi/6, xyzDivs=c(0.1,0.2,0.02),
           zoffset=0, sampleIsUsed=NULL, drawBoundingBox=FALSE){
  
  ## DrawBaseMap
  ##
  ## Draws a 2d WOHM with shading indicating severity.

  ######
  ## Define the gray shades to be used
  ######

  gradeGrays <- c(0.9, 0.8, 0.7, 0.6, 0.5, 0.3)
  gradeGrays <- gray(gradeGrays)
  names(gradeGrays) <- paste("D", 0:5, sep="")

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
          type="l",col="red")

  }
    
  offsetMat <- (as.matrix(mapInfo[,c("Row","Column")]) - 1) %*% unitDims

  if(is.null(sampleIsUsed)){
    sampleIsUsed = rep(TRUE,dim(mapInfo)[1])
  }
  
  for(i1 in 1:dim(mapInfo)[1]){

    polygon(offsetMat[i1,1] + mapPanel[,"X"],
            offsetMat[i1,2] + mapPanel[,"Y"] + zoffset,
            col=gradeGrays[mapInfo[i1,"Grade"]],
            border="black"
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
