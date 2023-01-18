\name{FSdb2precursorType}
\alias{FSdb2precursorType}
\title{
Precursor Types from Fragmentation Spectra DataBase (FSDB)
}
\description{
This function finds potential ionization pathways for molecular formulas using a vector of InChIKey values from an FSDB. This function only searches for the first 14 InChIKey letters; and therefore, may result with multiple potential precursor types.
}
\usage{
FSdb2precursorType(InChIKeyVector, libFSdb, tableIndicator = "Frequency",
number_processing_threads = 1)
}
\arguments{
  \item{InChIKeyVector}{
A vector of InChIKey values. This value may contain whole InChIKey strings or first 14 InChIKey letters.
}
  \item{libFSdb}{
A converted MSP library reference file using the `msp2FSdb` module which is an FSDB produced by the IDSL.FSA package.
}
  \item{tableIndicator}{
c("Frequency", "PrecursorMZ"). To show frequency or a median of 'PrecursorMZ' values in the output dataframe for each precursor type.
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A matrix of frequency for each InChIKey in the FSDB. The matrix column headers represent precursor types.
}
\examples{
address_input_msp <- system.file("extdata", package = "IDSL.FSA")
MSPfile_vector <- c("Kynurenine_Kynurenic_acid.msp")
libFSdb <- msp2FSdb(path = address_input_msp, MSPfile_vector)
##
InChIKeyVector <- c("HCZHHEIFKROPDY-UHFFFAOYSA-N", "YGPSJZOEDVAXAB-QMMMGPOBSA-N")
precursor_type_table <- FSdb2precursorType(InChIKeyVector, libFSdb,
tableIndicator = "Frequency", number_processing_threads = 1)
}