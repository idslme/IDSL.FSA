\name{FSA_uniqueMSPblockTaggerUntargeted}
\alias{FSA_uniqueMSPblockTaggerUntargeted}
\title{
FSA_uniqueMSPblockTaggerUntargeted
}
\description{
FSA_uniqueMSPblockTaggerUntargeted
}
\usage{
FSA_uniqueMSPblockTaggerUntargeted(path, MSPfile_vector = "",
peak_alignment_folder = NA, minCSAdetectionFrequency = 20, massError = 0.01,
massErrorPrecursor = 0.01, RTtolerance = 0.05, minEntropySimilarity = 0.75,
noiseRemovalRatio = 0.01, minCosineSimilarity = 0.75,
allowedNominalMass = FALSE, allowedWeightedSpectralEntropy = TRUE,
plotSpectra = FALSE, number_processing_threads = 1)
}
\arguments{
  \item{path}{
Address of .msp file(s)
}
  \item{MSPfile_vector}{
A vector of names of .msp files or one .msp file name. Default is "".
}
  \item{peak_alignment_folder}{
peak alignment folder from the IDSL.IPA analysis
}
  \item{minCSAdetectionFrequency}{
minimum CSA detection frequency across the processed batches of samples provided at \code{MSPfile_vector}
}
  \item{massError}{
Mass accuracy in Da
}
  \item{massErrorPrecursor}{
Mass accuracy (Da) to find precursor m/z in .msp files
}
  \item{RTtolerance}{
Retention time tolerance (min)
}
  \item{minEntropySimilarity}{
Minimum entropy similarity score
}
  \item{noiseRemovalRatio}{
noise removal ratio ([0 - 1]) relative to the basepeak to measure entropy similarity score.
}
  \item{minCosineSimilarity}{
minimum cosine similarity score for spectral comparison.
}
  \item{allowedNominalMass}{
c(TRUE, FALSE). Select `TRUE` only for nominal mass analysis.
}
  \item{allowedWeightedSpectralEntropy}{
c(TRUE, FALSE). Weighted entropy to transform low abundant signals prior to calculating entropy similarity score. Please see the reference for details on weight transformation.
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