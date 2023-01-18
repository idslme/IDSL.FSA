\name{spectral_entropy_similarity_score}
\alias{spectral_entropy_similarity_score}
\title{
Spectral Entropy Calculator
}
\description{
This module measures similarity of spectral entropies between `PEAK_A` and `PEAK_B` fragment spectra using a method described by the reference paper.
}
\usage{
spectral_entropy_similarity_score(PEAK_A, S_PEAK_A, PEAK_B, S_PEAK_B, massError,
allowedNominalMass = FALSE)
}
\arguments{
  \item{PEAK_A}{
A matrix (m/z, int) of fragmentation spectra
}
  \item{S_PEAK_A}{
Spectral entropy of PEAK_A
}
  \item{PEAK_B}{
A matrix (m/z, int) of fragmentation spectra
}
  \item{S_PEAK_B}{
Spectral entropy of PEAK_B
}
  \item{massError}{
Mass accuracy in Da
}
  \item{allowedNominalMass}{
c(TRUE, FALSE). Select `TRUE` only for nominal mass analysis.
}
}
\value{
spectral entropy similarity between 0 - 1
}
\references{
Li, Y., Kind, T., Folz, J., Vaniya, A., Mehta, S.S. and Fiehn, O. (2021). Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. \emph{Nature methods}, 18(12), 1524-1531, \doi{10.1038/s41592-021-01331-z}
}
\examples{
allowedWeightedSpectralEntropy <- TRUE
##
A <- cbind(seq(50, 160, length.out = 10), seq(10, 90, length.out = 10))
sA <- spectral_entropy_calculator(A, allowedWeightedSpectralEntropy)
S_PEAK_A <- sA[[1]]
PEAK_A <- sA[[3]]
##
B <- cbind(seq(50, 160, length.out = 10), seq(50, 60, length.out = 10))
sB <- spectral_entropy_calculator(A, allowedWeightedSpectralEntropy)
S_PEAK_B <- sB[[1]]
PEAK_B <- sB[[3]]
##
allowedNominalMass = TRUE
entropyScore <- spectral_entropy_similarity_score(PEAK_A, S_PEAK_A, PEAK_B,
S_PEAK_B, allowedNominalMass)
}