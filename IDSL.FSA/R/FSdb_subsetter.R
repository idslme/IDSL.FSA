FSdb_subsetter <- function(FSdb, inclusionIDs = NULL, exclusionIDs = NULL) {
  ##
  if (is.null(inclusionIDs) & is.null(exclusionIDs)) {
    stop("Select `inclusionIDs` or `exclusionIDs`!")
  } else {
    nFSdb <- length(FSdb[["Num Peaks"]])
    if (is.null(inclusionIDs)) {
      IDs <- seq(1, nFSdb, 1)
      inclusionIDs <- IDs[!(IDs %in% exclusionIDs)]
    }
    inclusionIDs <- sort(unique(inclusionIDs))
    if (max(inclusionIDs) > nFSdb) {
      stop("`inclusionIDs` are out of FSdb dimensions!")
    }
  }
  ##
  if (length(inclusionIDs) > 0) {
    FSdb[["PrecursorMZ"]] <- FSdb[["PrecursorMZ"]][inclusionIDs]
    FSdb[["Precursor Type"]] <- FSdb[["Precursor Type"]][inclusionIDs]
    FSdb[["Retention Time"]] <- FSdb[["Retention Time"]][inclusionIDs]
    FSdb[["Num Peaks"]] <- FSdb[["Num Peaks"]][inclusionIDs]
    FSdb[["Spectral Entropy"]] <- FSdb[["Spectral Entropy"]][inclusionIDs]
    FSdb[["FragmentList"]] <- FSdb[["FragmentList"]][inclusionIDs]
    ##
    ############################################################################
    ##
    if (length(inclusionIDs) > 1) {
      FSdb[["MSPLibraryParameters"]] <- FSdb[["MSPLibraryParameters"]][inclusionIDs, ]
    } else {
      FSdb[["MSPLibraryParameters"]] <- data.frame(FSdb[["MSPLibraryParameters"]][inclusionIDs, ])
    }
    ##
    dimMSP <- dim(FSdb[["MSPLibraryParameters"]])
    xNULL <- do.call(c, lapply(1:dimMSP[2], function(i) {
      x0 <- which(FSdb[["MSPLibraryParameters"]][, i] == "")
      if (length(x0) == dimMSP[1]) {
        i
      }
    }))
    #
    if (length(xNULL) > 0) {
      FSdb[["MSPLibraryParameters"]] <- FSdb[["MSPLibraryParameters"]][, -xNULL]
    }
    ##
    ############################################################################
    ##
  } else {
    FSdb <- NULL
  }
  ##
  return(FSdb)
}