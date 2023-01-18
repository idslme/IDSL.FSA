\name{xlsx2msp}
\alias{xlsx2msp}
\title{
xlsx to MSP
}
\description{
This function creates .msp files from an organized spreadsheet of fragmentation data.
}
\usage{
xlsx2msp(path, xlsxFileName = "", number_processing_threads = 1)
}
\arguments{
  \item{path}{
address of the spreadsheet
}
  \item{xlsxFileName}{
name of the file with the .xlsx extension.
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\note{
The spreadsheet should have only one column for the following headers (case-sensitive): c('ID', 'mz_fragment', 'int_fragment', 'Name')
}
\value{
The .msp files are saved in the same location.
}
\examples{
\donttest{
temp_wd <- tempdir() # just a temporary folder
path_extdata <- system.file("extdata", package = "IDSL.FSA")
xlsxFileName <- "PFAS_MSe.xlsx"
file.copy(from = paste0(path_extdata, "/", xlsxFileName), to = temp_wd)
xlsx2msp(temp_wd, xlsxFileName)
}
}