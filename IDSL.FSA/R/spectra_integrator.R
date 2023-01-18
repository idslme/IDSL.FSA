spectra_integrator <- function(stackedSpectra, massError = 0, allowedNominalMass = FALSE) {
  ##############################################################################
  if (allowedNominalMass) {
    stackedSpectra[, 1] <- round(stackedSpectra[, 1], digits = 0)
    stackedSpectra <- matrix(stackedSpectra[order(stackedSpectra[, 1], decreasing = FALSE), ], ncol = 2)
    xDiff <- c(0, which(diff(stackedSpectra[, 1]) > 0), dim(stackedSpectra)[1])
    LxDiff <- (length(xDiff) - 1)
    integratedSpectra <- do.call(rbind, lapply(1:LxDiff, function(j) {
      MZ <- stackedSpectra[(xDiff[j] + 1), 1]
      INT <- sum(stackedSpectra[seq((xDiff[j] + 1), xDiff[j + 1], 1), 2])
      c(MZ, INT)
    }))
    ##
    if (LxDiff == 1) {
      integratedSpectra <- matrix(integratedSpectra, ncol = 2)
    }
    ##
  } else {
    ############################################################################
    stackedSpectra <- matrix(stackedSpectra[order(stackedSpectra[, 2], decreasing = TRUE), ], ncol = 2)
    LstackedSpectra <- dim(stackedSpectra)[1]
    integratedSpectra <- matrix(rep(0, 2*LstackedSpectra), ncol = 2)
    counter <- 0
    ##
    for (j in 1:LstackedSpectra) {
      if (stackedSpectra[j, 1] > 0) {
        counter <- counter + 1
        x <- which(abs(stackedSpectra[j, 1] - stackedSpectra[, 1]) <= massError)
        if (length(x) == 1) {
          integratedSpectra[counter, 1] <- stackedSpectra[x, 1]
          integratedSpectra[counter, 2] <- stackedSpectra[x, 2]
          ##
        } else {
          sumINT <- sum(stackedSpectra[x, 2])
          meanMZ <- sum(stackedSpectra[x, 1]*stackedSpectra[x, 2])/sumINT
          ##
          integratedSpectra[counter, 1] <- meanMZ
          integratedSpectra[counter, 2] <- sumINT
        }
        ##
        stackedSpectra[x, ] <- 0
      }
    }
    ##
    integratedSpectra <- integratedSpectra[1:counter, ]
    if (counter == 1) {
      integratedSpectra <- matrix(integratedSpectra, ncol = 2)
    }
  }
  ##
  return(integratedSpectra)
}