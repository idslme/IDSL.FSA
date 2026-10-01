\name{FSA_uniqueTagPeakXcol}
\alias{FSA_uniqueTagPeakXcol}
\title{Integrate Untargeted MSP Tags with Aligned Peak Table}
\description{
\code{FSA_uniqueTagPeakXcol} integrates untargeted mass spectrometry (MSP) tags with an aligned peak table (\code{peakXcol}). It processes clustered peak IDs across multiple mzML files, filters them based on detection frequency, and performs frequent item-set mining. Finally, it annotates the aligned peaks with metadata for the top candidate compounds and saves the annotated matrix to disk.
}
\usage{
FSA_uniqueTagPeakXcol(uniqueMSPtagsUntargeted_folder,
                      peak_alignment_folder,
                      sortingMetavariable = "matchedRank",
                      metaVariables = "inchikey",
                      nCandidateCompounds = 5,
                      number_processing_threads = 1)
}
\arguments{
  \item{uniqueMSPtagsUntargeted_folder}{A character string specifying the directory path containing the untargeted MSP tags.}
  \item{peak_alignment_folder}{A character string specifying the directory path containing the aligned peak table.}
  \item{sortingMetavariable}{A character string indicating the metadata variable used to sort and rank the candidate compounds. Defaults to \code{"matchedRank"}.}
  \item{metaVariables}{A character vector of metadata variables to include in the final annotation. Defaults to \code{c("inchikey", "name", "precursor_type")}.}
  \item{nCandidateCompounds}{An integer specifying the maximum number of top candidate compounds to retain per peak. Defaults to \code{5}.}
  \item{number_processing_threads}{An integer specifying the number of threads to use for parallel processing. Defaults to \code{1}.}
}
\details{
This function parses data across multiple mzML files to find consistent clustered peak IDs. By filtering detection frequency and leveraging frequent item-set mining, it is able to confidently assign metadata variables to aligned MS features. The resulting annotated matrix is automatically saved to the target directory.
}
\value{
Returns \code{NULL} invisibly. The function is primarily called for its side effect of writing the annotated aligned peak matrix to disk within the specified working directories.
}