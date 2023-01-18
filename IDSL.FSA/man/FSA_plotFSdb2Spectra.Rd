\name{FSA_plotFSdb2Spectra}
\alias{FSA_plotFSdb2Spectra}
\title{
plot FSdb to Spectra
}
\description{
plot FSdb to Spectra
}
\usage{
FSA_plotFSdb2Spectra(path, allowedUnlink = TRUE, annexName = "", FSdb,
selectedFSdbIDs = NULL, number_processing_threads = 1, allowedVerbose = TRUE)
}
\arguments{
  \item{path}{
Address of .msp file(s)
}
  \item{allowedUnlink}{
allowedUnlink
}
  \item{annexName}{
annexName
}
  \item{FSdb}{
FSdb
}
  \item{selectedFSdbIDs}{
selected FSdb IDs. When `NULL`, the entire FSDB blocks are plotted. 
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
  \item{allowedVerbose}{
c(TRUE, FALSE). A `TRUE` allowedVerbose provides messages about the flow of the function.
}
}
\value{
spectra_figure object
}