\name{mgf2msp}
\alias{mgf2msp}
\title{
MGF to MSP
}
\description{
This function converts .mgf (Mascot generic format) files into the .msp (mass spectra) format.
}
\usage{
mgf2msp(path, MGFfile = "")
}
\arguments{
  \item{path}{
address of the .mgf file.
}
  \item{MGFfile}{
name of the file with the .mgf extension.
}
}
\value{
The .msp files are saved in the same location.
}
\examples{
\donttest{
temp_wd <- tempdir() # just a temporary folder
path_extdata <- system.file("extdata", package = "IDSL.FSA")
MGFfile <- "Training_000.mgf"
file.copy(from = paste0(path_extdata, "/", MGFfile), to = temp_wd)
mgf2msp(path = temp_wd, MGFfile)
}
}