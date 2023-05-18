\name{FSA_SpectraSimilarity_xlsxAnalyzer}
\alias{FSA_SpectraSimilarity_xlsxAnalyzer}
\title{
FSA SpectraSimilarity xlsx Analyzer
}
\description{
This function processes the spreadsheet of the `SpectraSimilarity` tab to ensure the parameter inputs are consistent with the requirements of the IDSL.FSA pipeline.
}
\usage{
FSA_SpectraSimilarity_xlsxAnalyzer(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
FSA spreadsheet
}
}
\value{
This function returns the FSA SpectraSimilarity parameters to feed the `FSA_msp_annotator` module.
}