spectral_entropy_calculator <- function(FragmentList, allowedWeightedSpectralEntropy = TRUE, noiseRemovalRatio = 0.01) {
  ##
  xNon0 <- which(FragmentList[, 2]/max(FragmentList[, 2]) >= noiseRemovalRatio)
  LxNon0 <- length(xNon0)
  if (LxNon0 > 0) {
    FragmentList <- FragmentList[xNon0, ]
    if (LxNon0 == 1) {
      FragmentList <- matrix(FragmentList, nrow = 1)
    }
    ##
    FragmentList[, 2] <- FragmentList[, 2]/sum(FragmentList[, 2])
    ##
    spectralEntropy <- -do.call(sum, lapply(FragmentList[, 2], function(j) {
      j*log(j)
    }))
    if (allowedWeightedSpectralEntropy) {
      if (spectralEntropy < 3) {
        w <- 0.25 + spectralEntropy * 0.25
        FragmentList[, 2] <- FragmentList[, 2]^w
        ##
        FragmentList[, 2] <- FragmentList[, 2]/sum(FragmentList[, 2])
        ##
        spectralEntropy <- -do.call(sum, lapply(FragmentList[, 2], function(j) {
          j*log(j)
        }))
      }
    }
    ##
    listSpectralEntropy <- list(spectralEntropy, LxNon0, FragmentList)
  } else {
    listSpectralEntropy <- list(0, 0, matrix(nrow = 0, ncol = 2))
  }
  return(listSpectralEntropy)
}