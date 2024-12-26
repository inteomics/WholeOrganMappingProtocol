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
