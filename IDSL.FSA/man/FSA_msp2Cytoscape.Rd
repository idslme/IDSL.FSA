\name{FSA_msp2Cytoscape}
\alias{FSA_msp2Cytoscape}
\title{
FSA Cytoscape Files Generator
}
\description{
This function generates necessary files from pairwise MSP blocks analysis to create Cytoscape networks.
}
\usage{
FSA_msp2Cytoscape(path, MSPfile = "", mspVariableVector = NULL,
mspNodeID = NULL, massError = 0.01, RTtolerance = NA, minEntropySimilarity = 0.75,
allowedNominalMass = FALSE, allowedWeightedSpectralEntropy = TRUE,
noiseRemovalRatio = 0.01, number_processing_threads = 1)
}
\arguments{
  \item{path}{
address of .msp file or an FSDB
}
  \item{MSPfile}{
name of .msp file
}
  \item{mspVariableVector}{
a vector of msp variables
}
  \item{mspNodeID}{
msp Node ID which is the ID that is required for the `specsim' ID generation
}
  \item{massError}{
Mass accuracy in Da
}
  \item{RTtolerance}{
Retention time tolerance (min) to match msp blocks. Select \emph{NA} to ignore retention time match. This option is so helpful to find co-occurring compounds.
}
  \item{minEntropySimilarity}{
Minimum entropy similarity score
}
  \item{allowedNominalMass}{
c(TRUE, FALSE). Select `TRUE` only for nominal mass analysis.
}
  \item{allowedWeightedSpectralEntropy}{
c(TRUE, FALSE). Weighted entropy to measure entropy similarity score.
}
  \item{noiseRemovalRatio}{
noise removal ratio relative to the basepeak to measure entropy similarity score (in percent)
}
  \item{number_processing_threads}{
Number of processing threads for multi-threaded processing
}
}
\value{
  \item{node_attributes_dataFrame}{
node_attributes dataframe. A string to store using `writeTable` function of R after a tab separation.
}
  \item{edge_dataFrame}{
edge dataframe. A string to store using the `writeTable` function of R after a tab separation.
}
  \item{correlation_network}{
correlation_network dataframe. A string to store using the `writeTable` function of R after a tab separation.
}
  \item{FSDB}{
Fragmentation spectra database (FSDB) object
}
  \item{exclusionMSPnoideid}{
A vector of MSP node ids which can be excluded to create a library of unique MSP blocks. 
}
 \item{filteredNetworkSIF}{
A filtered network in the cytoscape SIF format that does not have redundant MSP blocks within a RT window. 
}
}
\references{
Shannon, P., Markiel, A., Ozier, O., Baliga, N.S., Wang, J.T., Ramage, D., Amin, N., Schwikowski, B. and Ideker, T., (2003). Cytoscape: a software environment for integrated models of biomolecular interaction networks. \emph{Genome research}, 13(11), 2498-2504, \doi{10.1101/gr.1239303}
}
\examples{
\donttest{
path_extdata <- system.file("extdata", package = "IDSL.FSA")
mspFileName <- "Kynurenine_Kynurenic_acid.msp"
##
listCytoscape <- FSA_msp2Cytoscape(path = path_extdata,
MSPfile = mspFileName, mspVariableVector = c("Name", "Collision_energy"),
mspNodeID = NULL, massError = 0.01, RTtolerance = NA, minEntropySimilarity = 0,
noiseRemovalRatio = 0, allowedNominalMass = FALSE,
allowedWeightedSpectralEntropy = TRUE, number_processing_threads = 1)
##
FSDB <- listCytoscape[["FSDB"]]
##
temp_wd <- tempdir() # just a temporary folder to save results
##
write.table(listCytoscape[["node_attributes_dataFrame"]], paste0(temp_wd,
"/node_attributes_dataFrame.txt"), quote = FALSE, sep = "\t", row.names = FALSE,
col.names = FALSE)
##
write.table(listCytoscape[["correlation_network"]], paste0(temp_wd,
"/correlation_network.sif"), quote = FALSE, sep = "\t", row.names = FALSE,
col.names = FALSE)
##
write.table(listCytoscape[["edge_dataFrame"]], paste0(temp_wd,
"/edge_dataFrame.txt"), quote = FALSE, sep = "\t", row.names = FALSE,
col.names = FALSE)
##
}
}