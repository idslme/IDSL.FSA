\name{plotFSdb2SpectraCore}
\alias{plotFSdb2SpectraCore}
\title{
plot spectra from FSdb core
}
\description{
This function plots spectra figures from FSdb objects generated using the `msp2FSdb` function.
}
\usage{
plotFSdb2SpectraCore(FSdb, index)
}
\arguments{
  \item{FSdb}{
FSdb
}
  \item{index}{
index
}
}
\value{
spectra_figure object
}
\examples{
\donttest{
## To create the FSdb object
temp_wd <- tempdir() # just a temporary folder
path_extdata <- system.file("extdata", package = "IDSL.FSA")
MSPfile <- c("Kynurenine_Kynurenic_acid.msp")
file.copy(from = paste0(path_extdata, "/", MSPfile), to = temp_wd)
FSdb <- msp2FSdb(path = temp_wd, MSPfile)
## To plot spectra
index <- 1
plotFSdb2SpectraCore(FSdb, index)
}
}