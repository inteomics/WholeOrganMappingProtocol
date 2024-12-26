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
