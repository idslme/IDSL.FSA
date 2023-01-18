\name{spectra_ion_filter}
\alias{spectra_ion_filter}
\title{
Spectra Ion Filter
}
\description{
This function can detect m/z peaks that are related to each other across selected spectra lists.
}
\usage{
spectra_ion_filter(spectraList, indexSpectraList = length(spectraList), massError,
minPercentageDetectedScans = 10, rsdCutoff = 0, pearsonRHOthreshold = NA)
}
\arguments{
  \item{spectraList}{
a list of matrices of m/z and intensity values for each chromatogram scan
}
  \item{indexSpectraList}{
a vector of spectra indices for the analysis. This vector should have at least 3 elements to run this function.
}
  \item{massError}{
required mass error for m/z values
}
  \item{rsdCutoff}{
Relative standard deviations (in percent) to remove constant peaks (usually noisy peaks)
}
  \item{minPercentageDetectedScans}{
Minimum percentage of detected scans for an m/z peak
}
  \item{pearsonRHOthreshold}{
A threshold for pairwise Pearson's correlation coefficient across the selected spectra lists. This feature is recommended to find co-occurring peaks within a chromatographic peak. This feature may be used to eliminate instrument noises from MS2 data channels within an MS1 chromatographic peak for DDA analysis.
}
}
\value{
A matrix of m/z and cumulated intensities across the `indexSpectraList` spectra
}