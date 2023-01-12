spectral_entropy_similarity_score <- function(PEAK_A, S_PEAK_A, PEAK_B, S_PEAK_B, massError = 0, allowedNominalMass = FALSE) {
  ## To create 1:1 mixing AB spectral entropy
  spectra1A1B <- spectra_1A1B_mixer(PEAK_A, PEAK_B, massError, allowedNominalMass)
  ##
  S_spectra1A1B <- spectral_entropy_calculator(spectra1A1B, allowedWeightedSpectralEntropy = FALSE, noiseRemovalRatio = 0)[[1]]
  ##
  entropySimilarity <- 1 - (2*S_spectra1A1B - (S_PEAK_A + S_PEAK_B))/1.38629436111989  ## log(4) = 1.38629436111989
  ##
  if (is.numeric(entropySimilarity)) {
    if (!is.infinite(entropySimilarity)) {
      if (!is.nan(entropySimilarity)) {
        if (entropySimilarity > 0) {
          if (entropySimilarity > 1) {
            entropySimilarity <- 1 - entropySimilarity %% 1
          }
        } else {
          entropySimilarity <- 0
        }
      } else {
        entropySimilarity <- 0
      }
    } else {
      entropySimilarity <- 0
    }
  } else {
    entropySimilarity <- 0
  }
  ##
  return(entropySimilarity)
}