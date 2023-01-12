\name{FSA_msp_annotator}
\alias{FSA_msp_annotator}
\title{
FSA msp annotator
}
\description{
This function arranges the parameters for the annotation process
}
\usage{
FSA_msp_annotator(PARAM_SPEC, libFSdb, address_input_msp, output_path,
allowedVerbose = TRUE)
}
\arguments{
  \item{PARAM_SPEC}{
a parameter driven from the `FSA_SpectraSimilarity_xlsxAnalyzer` module.
}
  \item{libFSdb}{
a converted .msp library reference files (FSDB) using the `msp2FSdb` module
}
  \item{address_input_msp}{
address of the .msp files
}
  \item{output_path}{
output path
}
  \item{allowedVerbose}{
c(TRUE, FALSE). A `TRUE` allowedVerbose provides messages about the flow of the function.
}
}
\value{
A dataframe of matched annotated spectra stored in the output directory.
}