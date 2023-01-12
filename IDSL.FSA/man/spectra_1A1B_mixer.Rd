\name{spectra_1A1B_mixer}
\alias{spectra_1A1B_mixer}
\title{
Mixer 1:1 spectra A and B
}
\description{
This function creates 1:1 mixed AB spectra for spectral entropy calculation
}
\usage{
spectra_1A1B_mixer(PEAK_A, PEAK_B, massError = 0, allowedNominalMass = FALSE)
}
\arguments{
  \item{PEAK_A}{
A matrix (m/z, int) of fragmentation spectra
}
  \item{PEAK_B}{
A matrix (m/z, int) of fragmentation spectra
}
  \item{massError}{
Mass accuracy in Da
}
  \item{allowedNominalMass}{
c(TRUE, FALSE). Select `TRUE` only for nominal mass analysis.
}
}
\value{
A matrix of 1:1 mixing spectra. First and second columns represent intensity-weighted average m/z and cumulated intensity, respectively.
}