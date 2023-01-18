\name{FSA_FSdb_xlsxAnalyzer}
\alias{FSA_FSdb_xlsxAnalyzer}
\title{
FSA FSdb xlsx Analyzer
}
\description{
This function processes the spreadsheet of the `FSDB` tab to ensure the parameter inputs are consistent with the requirements of the IDSL.FSA pipeline.
}
\usage{
FSA_FSdb_xlsxAnalyzer(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
FSA spreadsheet
}
}
\value{
This function returns the FSDB parameters to feed the 'FSdb_file_generator' function.
}
