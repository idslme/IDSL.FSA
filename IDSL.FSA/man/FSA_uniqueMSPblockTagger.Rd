\name{FSA_uniqueMSPblockTagger}
\alias{FSA_uniqueMSPblockTagger}
\title{
FSA Unique MSP Block Tagger
}
\description{
This function removes similar MSP blocks. This function aggregates MSP blocks based on the `Name` values.
}
\usage{
FSA_uniqueMSPblockTagger(path, MSPfile = "", aggregateBy = "Name",
massError = 0.01, RTtolerance = NA, minEntropySimilarity = 0.75,
noiseRemovalRatio = 0.01, allowedNominalMass = FALSE,
allowedWeightedSpectralEntropy = TRUE, plotSpectra = FALSE,
number_processing_threads = 1)
}
\arguments{
  \item{path}{
Address of .msp file or an FSDB
}
  \item{MSPfile}{
name of .msp file
}
  \item{aggregateBy}{
a variable to aggregate the MSP blocks based on
}
  \item{massError}{
Mass accuracy in Da
}
  \item{RTtolerance}{
Retention time tolerance (min) to match msp blocks. Select \emph{NA} to ignore retention time match.
}
  \item{minEntropySimilarity}{
Minimum entropy similarity score
}
  \item{noiseRemovalRatio}{
noise removal ratio relative to the basepeak to measure entropy similarity score (in percent)
}
  \item{allowedNominalMass}{
c(TRUE, FALSE). Select `TRUE` only for nominal mass analysis.
}
  \item{allowedWeightedSpectralEntropy}{
c(TRUE, FALSE). Weighted entropy to measure entropy similarity score.
}
  \item{plotSpectra}{
c(TRUE, FALSE)
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
a list of similar MSP blocks is returned at the end and a subsetted .msp and FSDB files are saved in the `path` directory.
}