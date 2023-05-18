\name{FSdb2msp}
\alias{FSdb2msp}
\title{
Fragmentation Spectra DataBase (FSDB) to MSP
}
\description{
This function converts FSDB R objects into .msp standard files.
}
\usage{
FSdb2msp(path, FSdbFileName = "", UnweightMSP = FALSE,
number_processing_threads = 1)
}
\arguments{
  \item{path}{
address of .msp file(s)
}
  \item{FSdbFileName}{
name of the FSDB library name including `.Rdata` extension
}
  \item{UnweightMSP}{
to unweight fragmentation patterns
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
The .msp file is stored in the same folder
}