FSA_uniqueTagPeakXcol <- function(uniqueMSPtagsUntargeted_folder, peak_alignment_folder, sortingMetavariable = "matchedRank",
                                  metaVariables = "inchikey", nCandidateCompounds = 5, number_processing_threads = 1) {
  ##
  FSA_logRecorder("Started mapping the annotated compounds from the unique spectra tags onto the aligned table.")
  ##
  metaVariables <- tolower(metaVariables)
  if (any(metaVariables == "inchikey14")) {
    metaVariables[metaVariables == "inchikey14"] = "inchikey"
  }
  ##
  FSDBmetaVariables <- do.call(c, lapply(metaVariables, function(f) {paste0("FSDB_", f)}))
  nMetaVariables <- length(metaVariables)
  ##
  ##############################################################################
  ## 
  peakXcol <- FSA_loadRdata(paste0(peak_alignment_folder, "/peakXcol.Rdata"))
  nAlignedPeaks <- nrow(peakXcol)
  nSamples <- ncol(peakXcol) - 3
  ##
  ##############################################################################
  ## 
  uniqueMSPtagsUntargeted <- FSA_loadRdata(paste0(uniqueMSPtagsUntargeted_folder, "/uniqueMSPtagsUntargeted.Rdata"))
  rownames(uniqueMSPtagsUntargeted[["MSPLibraryParameters"]]) <- NULL
  uniqueClusterInfo = uniqueMSPtagsUntargeted[['uniqueClusterInfo']]
  lengthUniqueClusterInfo <- length(uniqueClusterInfo[["ClusterIndices"]])
  minCSAdetectionFrequency <- uniqueClusterInfo[["minCSAdetectionFrequency"]]
  ##
  CSAmode <- FALSE
  if (tolower(uniqueMSPtagsUntargeted[["MSPLibraryParameters"]][["msp_mode"]][1]) == "csa") {
    CSAmode <- TRUE
  }
  ##
  ##############################################################################
  ##
  if (minCSAdetectionFrequency == 0) {
    minCSAdetectionFrequency <- 1e-16
  }
  ##
  call_uniqueCluster <- function(u) {
    ##
    uniqueCluster <- uniqueClusterInfo[["ClusterIndices"]][[u]]
    ##
    if (uniqueClusterInfo[["lengthClusterIndices"]][u] < minCSAdetectionFrequency) {
      return(NULL)
    }
    ##
    return(table(do.call(c, lapply(uniqueCluster, function(i) {
      ##
      mzMLfile <- uniqueClusterInfo[["mzMLfilename"]][i]
      IPA_collective_peakids <- uniqueClusterInfo[["CollectivePeakIDs"]][i]
      IPA_collective_peakids <- strsplit(IPA_collective_peakids, ",")[[1]]
      IPA_collective_peakids <- as.numeric(IPA_collective_peakids[IPA_collective_peakids != "0"])
      ##
      which(peakXcol[, mzMLfile]%in%IPA_collective_peakids)
    }))))
  }
  ##
  ##############################################################################
  ##
  ## Processing OS
  osType <- Sys.info()[['sysname']]
  ##
  if ((number_processing_threads == 1) || (osType == "Windows")) {
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = lengthUniqueClusterInfo, initial = 0, style = 3)
    ##
    uniqueAlignedTable <- lapply(seq(1, lengthUniqueClusterInfo, 1), function(u) {
      setTxtProgressBar(progressBARboundaries, u)
      ##
      call_uniqueCluster(u)
    })
    ##
    close(progressBARboundaries)
    ##
  } else {
    ##
    ############################################################################
    ##
    uniqueAlignedTable <- mclapply(seq(1, lengthUniqueClusterInfo, 1), function(u) {
      call_uniqueCluster(u)
    }, mc.cores = number_processing_threads)
    ##
    closeAllConnections()
    ##
    ############################################################################
    ##
  }
  ##
  ##############################################################################
  ##
  uniqueXRow <- matrix(0, nrow = nAlignedPeaks, ncol = 2)
  increasingUniqueOrder <- order(uniqueClusterInfo[["lengthClusterIndices"]], decreasing = FALSE)
  ##
  for (u in increasingUniqueOrder) {
    alignedTable <- uniqueAlignedTable[[u]]
    if (length(alignedTable) > 0) {
      peakIDstr <- names(alignedTable)
      peakIDnum <- as.numeric(peakIDstr)
      for (i in 1:length(peakIDstr)) {
        freqU <- alignedTable[[peakIDstr[i]]]
        if (uniqueXRow[peakIDnum[i], 2] < freqU) {
          uniqueXRow[peakIDnum[i], 1] <- u
          uniqueXRow[peakIDnum[i], 2] <- freqU
        }
      }
    }
  }
  ##
  ##############################################################################
  ## Frequent Item-set Mining
  tXRow <- base::tapply(seq(1, nAlignedPeaks, 1), uniqueXRow[, 1], FUN = 'c', simplify = FALSE)
  tXRow["0"] <- NULL
  ##
  freqRatio <- minCSAdetectionFrequency/nSamples
  ##
  uniqueXRowFiltered <- matrix(0, nrow = nAlignedPeaks, ncol = 2)
  ##
  for (u in names(tXRow)) {
    ##
    alignedPeakIDs <- tXRow[[u]]
    freqs <- uniqueXRow[alignedPeakIDs, 2]
    ##
    ## Calculate the dynamic cutoff for this cluster
    freqCutoff <- max(freqs) * freqRatio
    ##
    xFreq <- which(freqs >= freqCutoff)
    if (!CSAmode) {
      ##
      alignedPeakIDsTopkmean <- alignedPeakIDs[xFreq]
      uniqueXRowFiltered[alignedPeakIDsTopkmean, 1] <- as.numeric(u)
      uniqueXRowFiltered[alignedPeakIDsTopkmean, 2] <- freqs[xFreq]
      ##
    } else if (length(xFreq) > 1) {
      ##
      alignedPeakIDsTopkmean <- alignedPeakIDs[xFreq]
      uniqueXRowFiltered[alignedPeakIDsTopkmean, 1] <- as.numeric(u)
      uniqueXRowFiltered[alignedPeakIDsTopkmean, 2] <- freqs[xFreq]
    }
  }
  tXRowFiltered <- base::tapply(seq(1, nAlignedPeaks, 1), uniqueXRowFiltered[, 1], FUN = 'c', simplify = FALSE)
  tXRowFiltered["0"] <- NULL
  ##
  ##############################################################################
  ##
  annotatedUniqueTag <- FSA_loadRdata(paste0(uniqueMSPtagsUntargeted_folder, "/annotated_spectra_tables/SpectraAnnotationTable_uniqueMSPtagsUntargeted.msp.Rdata"))
  ##
  annType = NULL
  if (!any(colnames(annotatedUniqueTag) == "analyte_fsa_unique_id")) {
    FSA_message("Unable to aggregate the annotated unique MSP file. This function is primarily designed to analyze data from IDSL suite results!", failedMessage = TRUE)
    return()
  }
  ##
  nAnnotatedTag <- nrow(annotatedUniqueTag)
  ##
  tUniqueTag <- base::tapply(seq(1, nAnnotatedTag, 1), annotatedUniqueTag[, "analyte_fsa_unique_id"], FUN = 'c', simplify = FALSE)
  ##
  ##############################################################################
  ##
  FSA_Unique_ID <- as.character(uniqueMSPtagsUntargeted[["MSPLibraryParameters"]][["fsa_unique_id"]])
  ##
  ##############################################################################
  ##
  uniqueXRowAnnotated <- matrix("", nrow = nAlignedPeaks, ncol = nMetaVariables*nCandidateCompounds)
  ##
  progressBARboundaries <- txtProgressBar(min = 0, max = length(FSA_Unique_ID), initial = 0, style = 3)
  uCounter <- 0
  ##
  for (u in FSA_Unique_ID) {
    if (!is.null(tUniqueTag[[u]])) {
      a <- annotatedUniqueTag[tUniqueTag[[u]], ]
      a <- a[1:max(a[, sortingMetavariable]), ] ## This step is needed for the CSA .msp annotations that have redundant annotations for co-occurring IDS.IPA peak IDs
      ##
      xRowID <- tXRowFiltered[[u]]
      if (!is.null(xRowID)) {
        ##
        rowProperties <- rep("", nMetaVariables*nCandidateCompounds)
        iCounter <- 0
        for (i in 1:min(c(nCandidateCompounds, nrow(a)))) {
          for (j in FSDBmetaVariables) {
            iCounter <- iCounter + 1
            if (!is.null(a[i, j])) {
              rowProperties[iCounter] <- a[i, j]
            }
          }
        }
        ##
        for (i in xRowID) {
          uniqueXRowAnnotated[i, ] <- rowProperties
        }
        ##
      }
    }
    ##
    uCounter <- uCounter + 1
    setTxtProgressBar(progressBARboundaries, uCounter)
    ##
  }
  ##
  close(progressBARboundaries)
  ##
  ##############################################################################
  ##
  metaVariables <- do.call(c, lapply(metaVariables, function(f) {sub("^FSDB_", "", f)}))
  ##
  peakXRowAnnotated <- cbind(peakXcol[, c("mz", "RT","frequencyPeakXcol")], uniqueXRowFiltered, uniqueXRowAnnotated)
  ##
  colnames(peakXRowAnnotated) <- c(c("mz", "RT","frequencyPeakXcol"), c("uniqueID", "freqUniqueID"),
                                   do.call(c, lapply(1:nCandidateCompounds, function(i) {
                                     do.call(c, lapply(metaVariables, function(j) {
                                       paste0(j, "_", i)
                                     }))
                                   })))
  ##
  alignedUniqueDir <- paste0(uniqueMSPtagsUntargeted_folder, "/aligned_unique_peaks/")
  FSA_dir.create(alignedUniqueDir, allowedUnlink = FALSE)
  ##
  save(peakXRowAnnotated, file = paste0(alignedUniqueDir, "/peakXRowAnnotated.Rdata"))
  write.csv(peakXRowAnnotated, file = paste0(alignedUniqueDir, "/peakXRowAnnotated.csv"), row.names = TRUE)
  ##
  FSA_logRecorder(paste0("Completed mapping the annotated compounds from the unique spectra tags onto the aligned table. Results are stored in the `", alignedUniqueDir,"` folder!"))
  ##
  ##############################################################################
  ##
}