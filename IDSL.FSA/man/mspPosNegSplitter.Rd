\name{mspPosNegSplitter}
\alias{mspPosNegSplitter}
\title{
MSP Pos/Neg Splitter
}
\description{
This function separates the positive and negative MSP blocks.
}
\usage{
mspPosNegSplitter(path, MSPfile = "", number_processing_threads = 1)
}
\arguments{
  \item{path}{
address of the .msp file.
}
  \item{MSPfile}{
name of the file with the .msp extension.
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
The .msp files are saved in the same location with `_Neg.msp` and `_Pos.msp` extensions.
}
\examples{
\donttest{
temp_wd <- tempdir() # just a temporary folder
path_extdata <- system.file("extdata", package = "IDSL.FSA")
MSPfile <- "Kynurenine_Kynurenic_acid.msp"
file.copy(from = paste0(path_extdata, "/", MSPfile), to = temp_wd)
mspPosNegSplitter(temp_wd, MSPfile)
}
}