\name{msp2TrainingMatrix}
\alias{msp2TrainingMatrix}
\title{
msp to Fragmentation Spectra DataBase (FSDB)
}
\description{
This function creates an aligned table from the spectra in the .msp file
}
\usage{
msp2TrainingMatrix(path, MSPfile = "", minDetectionFreq = 100,
selectedFSdbIDs = NULL, dimension = "wide", massAccuracy = 0.01,
allowedNominalMass = FALSE, allowedWeightedSpectralEntropy = TRUE,
noiseRemovalRatio = 0.01, number_processing_threads = 1)
}
\arguments{
  \item{path}{
Address of .msp file or an FSDB
}
  \item{MSPfile}{
A .msp file name or FSDB in .Rdata format
}
  \item{minDetectionFreq}{
A minimum detection frequency for an ion across the entire spectra
}
  \item{selectedFSdbIDs}{
selected MSP block/FSDB IDs to limit the screening to specific ion blocks
}
  \item{dimension}{
c("wide", "long"). *wide* or *long* alignment matrix output
}
  \item{massAccuracy}{
A mass accuracy (Da)
}
  \item{allowedNominalMass}{
c(TRUE, FALSE). Select `TRUE` only for nominal mass analysis.
}
  \item{allowedWeightedSpectralEntropy}{
c(TRUE, FALSE). Weighted entropy to transform low abundant signals prior to calculating entropy similarity score. Please see the reference for details on wight transformation.
}
  \item{noiseRemovalRatio}{
noise removal ratio ([0 - 1])relative to the basepeak to measure entropy similarity score.
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A FSDB file (.Rdata) and aligned spectra table (.csv) are stored in the same directory.
}
\examples{
\donttest{
temp_wd <- tempdir() # just a temporary folder
path_extdata <- system.file("extdata", package = "IDSL.FSA")
MSPfile <- "Kynurenine_Kynurenic_acid.msp"
file.copy(from = paste0(path_extdata, "/", MSPfile), to = temp_wd)
msp2TrainingMatrix(path = temp_wd, MSPfile, minDetectionFreq = 1)
}
}