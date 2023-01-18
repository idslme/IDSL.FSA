spectra_ion_filter <- function(spectraList, indexSpectraList = length(spectraList), massError, minPercentageDetectedScans = 10, rsdCutoff = 0, pearsonRHOthreshold = NA) {
  ##
  minRatioDetectedScans <- minPercentageDetectedScans/100
  rsdCutoffRatio <- rsdCutoff/100
  ##
  nScans <- length(indexSpectraList)
  nScans1 <- nScans + 1
  if (nScans < 3) {
    stop("The `spectra_ion_filter` function requires at least 3 indices for the `indexSpectraList` vector!")
  }
  ##
  imzRTint <- do.call(rbind, lapply(1:nScans, function(i) {
    specScan <- spectraList[[indexSpectraList[i]]]
    LspecScan <- dim(specScan)[1]
    cbind(rep((i + 1), LspecScan), specScan)
  }))
  ##
  imzRTint <- imzRTint[order(imzRTint[, 2], decreasing = FALSE), ]
  xDiffMZ <- c(0, which(diff(imzRTint[, 2]) > massError), nrow(imzRTint))
  LxDiffMZ <- length(xDiffMZ) - 1
  ##
  ##############################################################################
  ##
  corPeakTable <- do.call(rbind, lapply(1:LxDiffMZ, function(q) {
    ##
    counterMZ <- 0
    nSb <- xDiffMZ[q + 1] - xDiffMZ[q]
    if (nSb >= 3) { # At least 3 data points are necessary to calculate Pearson's correlation!
      if (nSb/nScans >= minRatioDetectedScans) {
        sbImzRTint <- imzRTint[seq((xDiffMZ[q] + 1), xDiffMZ[q + 1], 1), ]
        ##
        if (nSb == 1) {
          sbImzRTint <- matrix(sbImzRTint, nrow = 1)
          orderINT <- 1
        } else {
          orderINT <- order(sbImzRTint[, 3], decreasing = TRUE)
        }
        sbCorPeakTable <- matrix(rep(0, nScans1*nSb), ncol = nScans1)
        ##
        for (i in orderINT) {
          if (sbImzRTint[i, 1] != 0) {
            ##
            x <- c(i, which((abs(sbImzRTint[i, 2] - sbImzRTint[, 2]) <= massError) &
                              (sbImzRTint[i, 1] != sbImzRTint[, 1])))
            Lx <- length(x)
            if (Lx >= 3) { # At least 3 data points are necessary to calculate Pearson's correlation!
              if (Lx/nScans >= minRatioDetectedScans) {
                x <- FSA_aggregate(idVec = sbImzRTint[x, 1], variableVec = sbImzRTint[x, 2],
                                   indexVec = x, targetVar = sbImzRTint[i, 2])
                Lx <- length(x)
                if (Lx >= 3) { # At least 3 data points are necessary to calculate Pearson's correlation!
                  if (Lx/nScans >= minRatioDetectedScans) {
                    xSbCorPeakTable <- sbImzRTint[x, 3]
                    rsd <- sd(xSbCorPeakTable)/mean(xSbCorPeakTable)
                    if (rsd >= rsdCutoffRatio) {
                      ##
                      counterMZ <- counterMZ + 1
                      sbCorPeakTable[counterMZ, 1] <- sbImzRTint[i, 2]
                      ##
                      for (j in 1:Lx) {
                        sbCorPeakTable[counterMZ, sbImzRTint[x[j], 1]] <- sbImzRTint[x[j], 3]
                      }
                      sbImzRTint[x, ] <- 0
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    ##
    if (counterMZ > 0) {
      sbCorPeakTable <- sbCorPeakTable[1:counterMZ, ]
    } else {
      sbCorPeakTable <- NULL
    }
    ##
    return(sbCorPeakTable)
  }))
  ##
  imzRTint <- NULL
  ##
  ##############################################################################
  ##
  if (!is.null(corPeakTable)) {
    corPeakTable <- matrix(corPeakTable, ncol = nScans1)
    ##
    nAlignedPeaks <- nrow(corPeakTable)
    if (nAlignedPeaks == 1) {
      xmzID <- 1
    } else {
      if (!is.na(pearsonRHOthreshold)) {
        correlationMatrix <- suppressWarnings(cor(t(corPeakTable[, 2:nScans1]), method = "pearson"))
        ##
        xmzID <- do.call(c, lapply(1:(nAlignedPeaks - 1), function(i) {
          do.call(c, lapply((i + 1):nAlignedPeaks, function(j) {
            if (!is.na(correlationMatrix[i, j])) {
              if (correlationMatrix[i, j] >= pearsonRHOthreshold) {
                c(i, j)
              }
            }
          }))
        }))
        ##
        if (!is.null(xmzID)) {
          xmzID <- unique(xmzID)
        }
        ##
      } else {
        xmzID <- seq(1, nAlignedPeaks, 1)
      }
    }
    ##
    if (!is.null(xmzID)) {
      ##
      MS2_spectra <- do.call(rbind, lapply(xmzID, function(i) {
        c(corPeakTable[i, 1], sum(corPeakTable[i, 2:nScans1]))
      }))
    } else {
      MS2_spectra <- matrix(nrow = 0, ncol = 2)
    }
  } else {
    MS2_spectra <- matrix(nrow = 0, ncol = 2)
  }
  ##
  return(MS2_spectra)
}