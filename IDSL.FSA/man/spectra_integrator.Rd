\name{spectra_integrator}
\alias{spectra_integrator}
\title{
Spectra Integrator
}
\description{
This function integrates individual m/z peaks from multiple chromatogram scans (spectra) into summed m/z peaks using a mass accuracy or nominal masses.
}
\usage{
spectra_integrator(stackedSpectra, massError = 0, allowedNominalMass = FALSE)
}
\arguments{
  \item{stackedSpectra}{
A matrix of two columns of the stacked spectra. First and second columns should represent m/z and intensity, respectively.
}
  \item{massError}{
Mass accuracy in Da
}
  \item{allowedNominalMass}{
c(TRUE, FALSE). Select `TRUE` only for nominal mass analysis.
}
}
\value{
A matrix of integrated spectra. First and second columns represent intensity-weighted average m/z and cumulated intensity, respectively.
}
\examples{
data(stackedSpectra)
massError <- 0.005 # Da
Integrated_spectra <- spectra_integrator(stackedSpectra[, 1:2], massError)
}