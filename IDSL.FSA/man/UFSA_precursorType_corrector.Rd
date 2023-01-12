\name{UFSA_precursorType_corrector}
\alias{UFSA_precursorType_corrector}
\title{
UFA Precursor Type Corrector
}
\description{
Precursor type corrector from MSP files. This function initially attempts to standardize the precursor types to be consistent with the `ionization_pathway_deconvoluter` module of the IDSL.SUFA package.
}
\usage{
UFSA_precursorType_corrector(precursorType, ionMode = NULL)
}
\arguments{
  \item{precursorType}{
precursorType
}
  \item{ionMode}{
ionMode
}
}
\value{
correctedPrecursorType
}
\examples{
uncorrectedPrecursorType <- c("[M]+", "[M+H]+", "[2M-Cl]-", "[3M+COO-H2O+Na-KO2+HCl-NH4]-")
precursorType <- UFSA_precursorType_corrector(uncorrectedPrecursorType, ionMode = NULL)
}