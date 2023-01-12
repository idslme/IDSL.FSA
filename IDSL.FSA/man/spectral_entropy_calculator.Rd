\name{spectral_entropy_calculator}
\alias{spectral_entropy_calculator}
\title{
Spectral Entropy Calculator
}
\description{
This module calculates spectral entropy for a fragmentation pattern using a method described by the reference paper.
}
\usage{
spectral_entropy_calculator(FragmentList, allowedWeightedSpectralEntropy = TRUE,
noiseRemovalRatio = 0.01)
}
\arguments{
  \item{FragmentList}{
A matrix (m/z, int) of fragmentation pattern after intensity adjustment
}
  \item{allowedWeightedSpectralEntropy}{
c(TRUE, FALSE). Weighted entropy to transform low abundant signals prior to calculating entropy similarity score. Please see the reference for details on weight transformation.
}
  \item{noiseRemovalRatio}{
noise removal ratio ([0 - 1])relative to the basepeak to measure entropy similarity score.
}
}
\value{
  \item{spectralEntropy}{
spectral entropy
}
  \item{NumPeaks}{
NumPeaks
}
  \item{FragmentList}{
A matrix of two-columns after intensity normalization relative to summation of intensities AND entropy weight transformation when is selected.
}
}
\note{
noise removal on intensities should be performed prior to feeding to this function
}
\references{
Li, Y., Kind, T., Folz, J., Vaniya, A., Mehta, S.S. and Fiehn, O. (2021). Spectral entropy outperforms MS/MS dot product similarity for small-molecule compound identification. \emph{Nature methods}, 18(12), 1524-1531, \doi{10.1038/s41592-021-01331-z}
}
\examples{
FragmentList <- cbind(seq(50, 600, length.out = 10), seq(10, 90, length.out = 10))
SE <- spectral_entropy_calculator(FragmentList)
print(SE[[1]])
}