\name{FSA_spectra_marker_generator}
\alias{FSA_spectra_marker_generator}
\title{
FSA Spectra Marker Generator
}
\description{
This function generates spectra markers
}
\usage{
FSA_spectra_marker_generator(FSdb, ratio2basePeak4nSpectraMarkers = 0,
aggregationLevel = NA)
}
\arguments{
  \item{FSdb}{
FSdb object from the 'msp2FSdb' module
}
  \item{ratio2basePeak4nSpectraMarkers}{
Ratio of peaks in fragmentation spectra to the basepeak to calculate minimum qualified number of matched abundant peaks
}
  \item{aggregationLevel}{
c(NA, 0, 1, 2, 3). When `NA`, this function returns a matrix for the spectra markers. When integer numbers are used, the ion marker masses are grouped by a rounding digit equal to this number.
}
}
\value{
  \item{spectraMarkerMass}{
a grouped or a matrix of ion marker masses corresponding to FSdb ids
}
  \item{nSpectraMarkers}{
number of spectra markers for each FSdb id
}
}