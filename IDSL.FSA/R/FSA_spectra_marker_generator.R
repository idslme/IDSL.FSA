FSA_spectra_marker_generator <- function(FSdb, ratio2basePeak4nSpectraMarkers = 0, aggregationLevel = NA) {
  ##
  if (ratio2basePeak4nSpectraMarkers > 0) {
    ratio2basePeak4nSpectraMarkersCheck <- TRUE
  } else {
    ratio2basePeak4nSpectraMarkersCheck <- FALSE
  }
  ##
  iSeqFSdb = seq(1, length(FSdb[["Num Peaks"]]), 1)
  ##
  ##############################################################################
  ##
  markersList <- lapply(iSeqFSdb, function(i) {
    ##
    if (FSdb[["Num Peaks"]][i] > 0) {
      FragList <- FSdb[["FragmentList"]][[i]]
      ##
      if (ratio2basePeak4nSpectraMarkersCheck) {
        xSpectraMarkers <- which(FragList[, 2] >= max(FragList[, 2])*ratio2basePeak4nSpectraMarkers)
        nSpectraMarkers <- length(xSpectraMarkers)
      } else {
        xSpectraMarkers <- seq(1, FSdb[["Num Peaks"]][i], 1)
        nSpectraMarkers <- FSdb[["Num Peaks"]][i]
      }
      ##
      list(spectraMarkersList = cbind(FragList[xSpectraMarkers, 1], rep(i, nSpectraMarkers)), nSM = nSpectraMarkers)
    } else {
      ##
      list(spectraMarkersList = c(0, i), nSM = 0)
    }
  })
  ##
  spectraMarkerMass <- do.call(rbind, lapply(iSeqFSdb, function(i) {
    markersList[[i]][[1]]
  }))
  ##
  nSpectraMarkers <- do.call(c, lapply(iSeqFSdb, function(i) {
    markersList[[i]][[2]]
  }))
  nSpectraMarkers <- as.numeric(nSpectraMarkers)
  ##
  markersList <- NULL
  x_non0 <- which(spectraMarkerMass[, 1] > 0)
  spectraMarkerMass <- matrix(spectraMarkerMass[x_non0, ], ncol = 2)
  ##
  ##############################################################################
  ##
  if (!is.na(aggregationLevel)) {
    spectraMarkerMass[, 1] <- round(spectraMarkerMass[, 1], digits = aggregationLevel)
    spectraMarkerMass <- matrix(spectraMarkerMass[order(spectraMarkerMass[, 1], decreasing = FALSE), ], ncol = 2)
    xDiffSpectraMarkers <- c(0, which(diff(spectraMarkerMass[, 1]) > 0), dim(spectraMarkerMass)[1])
    LDiffSpectraMarkers <- length(xDiffSpectraMarkers) - 1
    namesSpectraMarkerMass <- spectraMarkerMass[(xDiffSpectraMarkers[1:LDiffSpectraMarkers] + 1), 1]
    ##
    spectraMarkerMass <- lapply(1:LDiffSpectraMarkers, function(i) {
      unique(spectraMarkerMass[seq((xDiffSpectraMarkers[i] + 1), xDiffSpectraMarkers[i + 1], 1), 2])
    })
    ##
    names(spectraMarkerMass) <- namesSpectraMarkerMass
    ##
  } else {
    spectraMarkerMass <- data.frame(spectraMarkerMass)
    rownames(spectraMarkerMass) <- NULL
    colnames(spectraMarkerMass) <- c("spectraMarkerMass", "ID")
  }
  ##
  ##############################################################################
  ##
  IDlist <- list(spectraMarkerMass, nSpectraMarkers)
  ##
  return(IDlist)
}