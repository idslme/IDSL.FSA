\name{stackedSpectra}
\alias{stackedSpectra}
\docType{data}
\title{
Example for a stacked spectra
}
\description{
A data to test the `spectra_integrator` function.
}
\usage{data("stackedSpectra")}
\format{
  \describe{
    \item{\code{mz}}{a numeric vector of m/z values}
    \item{\code{int}}{a numeric vector of intensities}
    \item{\code{scan_number}}{a numeric vector of chromatogram scan numbers}
  }
}
\details{
The `scan_number` column is not necessary to test the `spectra_integrator` function.
}
\examples{
data(stackedSpectra)
}
\keyword{datasets}
