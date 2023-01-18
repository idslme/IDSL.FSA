\name{fragmentation_spectra_annotator}
\alias{fragmentation_spectra_annotator}
\title{
Fragmentation Spectra Annotator
}
\description{
This module annotates fragmentation spectra from .MSP files.
}
\usage{
fragmentation_spectra_annotator(path, MSPfile = "", libFSdb,
libFSdbIDlist, targetedPrecursorType = NA, ratio2basePeak4nSpectraMarkers = 0,
allowedNominalMass = FALSE, allowedWeightedSpectralEntropy = TRUE,
noiseRemovalRatio = 0.01, roundingDigitPrefiltering = 1, minMatchedNumPeaks = 1,
massError = 0, maxNEME = 0, minIonRangeDifference = 0, minCosineSimilarity,
minEntropySimilarity, minRatioMatchedNspectraMarkers,
spectralEntropyDeviationPrefiltering, massErrorPrecursor = NA, RTtolerance = NA,
exportSpectraParameters = NULL, number_processing_threads = 1)
}
\arguments{
  \item{path}{
Address of .msp file(s)
}
  \item{MSPfile}{
name of the .msp file
}
  \item{libFSdb}{
A converted .msp library reference file using the `msp2FSdb` module which is an FSDB produced by the IDSL.FSA package.
}
  \item{libFSdbIDlist}{
Ion markers object from the FSDB reference
}
  \item{targetedPrecursorType}{
A vector of targeted precursor types
}
  \item{ratio2basePeak4nSpectraMarkers}{
Ratio of peaks in fragmentation spectra to the basepeak to calculate minimum qualified number of matched abundant peaks
}
  \item{allowedNominalMass}{
c(TRUE, FALSE). Select `TRUE` only for nominal mass analysis.
}
  \item{allowedWeightedSpectralEntropy}{
c(TRUE, FALSE). Weighted entropy to transform low abundant signals prior to calculating entropy similarity score. Please see the reference for details on wight transformation.
}
  \item{noiseRemovalRatio}{
noise removal ratio ([0 - 1])relative to the basepeak to measure entropy similarity score.
}
  \item{roundingDigitPrefiltering}{
Level of pre-filtering
}
  \item{minMatchedNumPeaks}{
Minimum matched number of peaks
}
  \item{massError}{
Mass accuracy in Da
}
  \item{maxNEME}{
Maximum value for Normalized Euclidean Mass Error (NEME) in mDa
}
  \item{minIonRangeDifference}{
Minimum distance (Da) between lowest and highest matched m/z to prevent matching only isotopic envelopes
}
  \item{minCosineSimilarity}{
Minimum cosine similarity score
}
  \item{minEntropySimilarity}{
Minimum entropy similarity score
}
  \item{minRatioMatchedNspectraMarkers}{
Minimum percentage of detection of abundant library peaks in percentage
}
  \item{spectralEntropyDeviationPrefiltering}{
Spectral entropy deviation for pre-filtering
}
  \item{massErrorPrecursor}{
Mass accuracy (Da) to find precursor m/z in .msp files
}
  \item{RTtolerance}{
Retention time tolerance (min)
}
  \item{exportSpectraParameters}{
Parameters for export MS/MS match figures
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
A dataframe of matched spectra
}
