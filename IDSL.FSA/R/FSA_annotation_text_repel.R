FSA_annotation_text_repel <- function(FSAspectra, nGridX, nGridY) {
  ##
  plotLabels <- FSAspectra
  nplotLabels <- nrow(plotLabels)
  ##
  gridLabelX <- seq((min(plotLabels[, 1]) - 10), (max(plotLabels[, 1]) + 10), length.out = nGridX)
  gridLabelY <- seq(min(plotLabels[, 2]), max(plotLabels[, 2]), length.out = nGridY)
  ##
  nGridX1 <- nGridX - 1
  xLabel <- do.call(c, lapply(1:(nGridY - 1), function(j) {
    xj <- which((plotLabels[, 2] >= gridLabelY[j]) & (plotLabels[, 2] <= gridLabelY[j + 1]))
    ##
    if (length(xj) > 0) {
      do.call(c, lapply(1:nGridX1, function(i) {
        xi <- which((plotLabels[xj, 1] >= gridLabelX[i]) & (plotLabels[xj, 1] <= gridLabelX[i + 1]))
        ##
        if (length(xi) > 0) {
          x <- xj[xi]
          xMax <- which.max(abs(plotLabels[x, 2]))
          x[xMax[1]]
        }
      }))
    }
  }))
  ##
  xLabel <- unique(xLabel)
  ##
  seqLabel <- seq(1, nplotLabels, 1)
  emptyLabels <- seqLabel[!(seqLabel %in% xLabel)]
  ##
  plotLabelsMZ <- round(FSAspectra[, 1], digits = 5)
  plotLabelsMZ[emptyLabels] <- ""
  ##
  return(plotLabelsMZ)
}