spectra_1A1B_mixer <- function(PEAK_A, PEAK_B, massError = 0, allowedNominalMass = FALSE) {
  ##
  nA <- nrow(PEAK_A)
  nB <- nrow(PEAK_B)
  stackedAB <- rbind(cbind(rep(1, nA), PEAK_A), cbind(rep(2, nB), PEAK_B))
  stackedAB <- stackedAB[order(stackedAB[, 3], decreasing = TRUE), ]
  nStackedAB <- nA + nB
  ##
  spectra1A1B <- matrix(rep(0, 2*nStackedAB), ncol = 2)
  counter1A1B <- 0
  ##############################################################################
  if (allowedNominalMass) {
    ##
    for (j in 1:nStackedAB) {
      if (stackedAB[j, 1] > 0) {
        counter1A1B <- counter1A1B + 1
        x <- which((stackedAB[j, 2] == stackedAB[, 2]) & (stackedAB[j, 1] != stackedAB[, 1]))
        Lx <- length(x)
        if (Lx == 0) {
          spectra1A1B[counter1A1B, 1] <- stackedAB[j, 2]
          spectra1A1B[counter1A1B, 2] <- stackedAB[j, 3]
          ##
          x <- j
          ##
        } else {
          ##
          if (Lx > 1) {
            xMax <- which.max(stackedAB[x, 3])
            x <- x[xMax[1]]
          }
          ##
          x <- c(j, x)
          ##
          spectra1A1B[counter1A1B, 1] <- stackedAB[j, 2]
          spectra1A1B[counter1A1B, 2] <- sum(stackedAB[x, 3])
        }
        ##
        stackedAB[x, ] <- 0
      }
    }
    ##
  } else {
    ############################################################################
    for (j in 1:nStackedAB) {
      if (stackedAB[j, 1] > 0) {
        counter1A1B <- counter1A1B + 1
        x <- which((abs(stackedAB[j, 2] - stackedAB[, 2]) <= massError) & (stackedAB[j, 1] != stackedAB[, 1]))
        Lx <- length(x)
        if (Lx == 0) {
          spectra1A1B[counter1A1B, 1] <- stackedAB[j, 2]
          spectra1A1B[counter1A1B, 2] <- stackedAB[j, 3]
          ##
          x <- j
          ##
        } else {
          ##
          if (Lx > 1) {
            xMin <- which.min(abs(stackedAB[j, 2] - stackedAB[x, 2]))
            x <- x[xMin[1]]
          }
          ##
          x <- c(j, x)
          ##
          sumINT <- sum(stackedAB[x, 3])
          meanMZ <- sum(stackedAB[x, 2]*stackedAB[x, 3])/sumINT
          ##
          spectra1A1B[counter1A1B, 1] <- meanMZ
          spectra1A1B[counter1A1B, 2] <- sumINT
        }
        ##
        stackedAB[x, ] <- 0
      }
    }
  }
  ##############################################################################
  spectra1A1B <- spectra1A1B[1:counter1A1B, ]
  if (counter1A1B == 1) {
    spectra1A1B <- matrix(spectra1A1B, ncol = 2)
  }
  ##
  return(spectra1A1B)
}