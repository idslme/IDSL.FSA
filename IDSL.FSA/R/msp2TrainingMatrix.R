msp2TrainingMatrix <- function(path, MSPfile = "", minDetectionFreq = 1, selectedFSdbIDs = NULL, dimension = "wide",
                               massAccuracy = 0.01, allowedNominalMass = FALSE, allowedWeightedSpectralEntropy = TRUE,
                               noiseRemovalRatio = 0.01, number_processing_threads = 1) {
  ##
  if (allowedNominalMass) {
    massAccuracy <- 0
  }
  ##
  dimension <- tolower(dimension)
  if (!((dimension == "wide") | (dimension == "long"))) {
    stop(FSA_logRecorder("The `dimension` variable must be either `wide` or `long`!"))
  }
  ##
  mspFileLocation <- paste0(path, "/", MSPfile)
  mspFileLocation <- gsub("\\", "/", mspFileLocation, fixed = TRUE)
  strmspFileLocation <- strsplit(mspFileLocation, "/")[[1]]
  LstrmspFileLocation <- length(strmspFileLocation)
  mspFileName <- strmspFileLocation[LstrmspFileLocation]
  mspFilePath <- paste0(strmspFileLocation[1:(LstrmspFileLocation - 1)], collapse = "/")
  ##
  ##############################################################################
  ##
  alignedMSPcsvFile <- paste0("/", mspFilePath, "/", dimension, "Aligned_", mspFileName, ".csv")
  if (file.exists(alignedMSPcsvFile)) {
    tryCatch(unlink(alignedMSPcsvFile), error = function(e){stop(paste0("Can't delete the `", alignedMSPcsvFile, "`!"))})
  }
  ##
  FSA_logRecorder("The alignment procedure may take a few minutes to several hours depending on the complexity of the .msp file!")
  ##
  ##############################################################################
  ##############################################################################
  ##
  strmspFileName <- strsplit(mspFileName, "[.]")[[1]]
  formatFileName <- tolower(strmspFileName[length(strmspFileName)])
  if (formatFileName == "msp") {
    FSA_logRecorder("Initiated deconvoluting the .msp file!")
    libFSdb <- msp2FSdb(path = mspFilePath, MSPfile_vector = mspFileName, massIntegrationWindow = massAccuracy, allowedNominalMass,
                        allowedWeightedSpectralEntropy, noiseRemovalRatio, number_processing_threads)
    FSA_logRecorder("Completed deconvoluting the .msp file!")
    ##
    save(libFSdb, file = paste0(mspFilePath, "/FSdb_", mspFileName, ".Rdata"))
    FSA_logRecorder("Stored the FSDB in the same address!")
  } else if (formatFileName == "rdata") {
    FSA_logRecorder("Loading the FSDB file in the `.Rdata` format!")
    libFSdb <- FSA_loadRdata(paste0(mspFilePath, "/", mspFileName))
  } else {
    stop(FSA_logRecorder("Allowed inputs : FSDB in `.Rdata` or MSP files in `.msp` formats!"))
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  FSA_logRecorder("Initiated aligning the fragmentation spectra!")
  if (is.null(selectedFSdbIDs)) {
    selectedFSdbIDs <- seq(1, length(libFSdb[["Num Peaks"]]), 1)
  } else if (length(which(!is.numeric(selectedFSdbIDs) | is.na(selectedFSdbIDs) | is.nan(selectedFSdbIDs) | (selectedFSdbIDs <= 0))) > 0) {
    stop(FSA_logRecorder("The `selectedFSdbIDs` variable contains non-numeric or non-standard values!"))
  } else if (max(selectedFSdbIDs) > length(libFSdb[["Num Peaks"]])) {
    stop(FSA_logRecorder("The `selectedFSdbIDs` variable contains indices greater than the size of the msp/FSdb file!"))
  }
  ##
  mainImzInt <- do.call(rbind, lapply(selectedFSdbIDs, function(i) {
    cbind(rep(i, libFSdb[["Num Peaks"]][i]), libFSdb[["FragmentList"]][[i]])
  }))
  libFSdb <- NULL
  ##
  nPeaks <- nrow(mainImzInt)
  mainImzInt <- mainImzInt[order(mainImzInt[, 2], decreasing = FALSE), ]
  xDiffMZ <- c(0, which(diff(mainImzInt[, 2]) > massAccuracy), nPeaks)
  LxDiffMZ <- length(xDiffMZ) - 1
  ##
  ##############################################################################
  ##############################################################################
  ##
  call_msp2alignedSpectraMatrix <- function(i) {
    ##
    nSb <- xDiffMZ[i + 1] - xDiffMZ[i]
    idSeedMZfragments <- rep(0, nSb)
    seedMZfragments <- idSeedMZfragments
    Counter <- 0
    ##
    if (nSb >= minDetectionFreq) {         # Minimum frequency of detection
      mainImzInt.sb <- mainImzInt[seq((xDiffMZ[i] + 1), xDiffMZ[i + 1], 1), ]
      if (nSb == 1) {
        mainImzInt.sb <- matrix(mainImzInt.sb, nrow = 1)
        orderINT <- 1
      } else {
        orderINT <- order(mainImzInt.sb[, 3], decreasing = TRUE)
      }
      ##
      for (j in orderINT) {
        if (mainImzInt.sb[j, 1] != 0) {
          if (allowedNominalMass) {
            x1 <- which(mainImzInt.sb[, 2] == mainImzInt.sb[j, 2])
          } else {
            x1 <- which(abs(mainImzInt.sb[, 2] - mainImzInt.sb[j, 2]) <= massAccuracy)
          }
          if (length(x1) >= minDetectionFreq) {         # Minimum frequency of detection
            ## To remove repeated scans
            x1 <- FSA_aggregate(idVec = mainImzInt.sb[x1, 1], variableVec = mainImzInt.sb[x1, 2],
                                indexVec = x1, targetVar = mainImzInt.sb[j, 2])
            ##
            if (length(x1) >= minDetectionFreq) {         # Minimum frequency of detection
              Counter <- Counter + 1
              idSeedMZfragments[x1] <- Counter
              seedMZfragments[Counter] <- mainImzInt.sb[j, 2]
              mainImzInt.sb[x1, ] <- 0
            }
          }
        }
      }
    }
    ##
    if (Counter > 0) {
      seedMZfragments <- seedMZfragments[1:Counter]
    } else {
      seedMZfragments <- NULL
    }
    ##
    list(seedMZfragments, idSeedMZfragments)
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = LxDiffMZ, initial = 0, style = 3)
    ##
    listSeedMZfragments <- lapply(1:LxDiffMZ, function(i) {
      setTxtProgressBar(progressBARboundaries, i)
      ##
      call_msp2alignedSpectraMatrix(i)
    })
    ##
    close(progressBARboundaries)
    ##
  } else {
    ## Processing OS
    osType <- Sys.info()[['sysname']]
    ##
    ############################################################################
    ##
    if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust", "LxDiffMZ")), envir = environment())
      ##
      listSeedMZfragments <- parLapply(clust, 1:LxDiffMZ, function(i) {
        call_msp2alignedSpectraMatrix(i)
      })
      ##
      stopCluster(clust)
      ##
      ##########################################################################
      ##
    } else {
      ##
      listSeedMZfragments <- mclapply(1:LxDiffMZ, function(i) {
        call_msp2alignedSpectraMatrix(i)
      }, mc.cores = number_processing_threads)
      ##
      closeAllConnections()
      ##
      ##########################################################################
      ##
    }
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  seedMZfragments <- do.call(c, lapply(listSeedMZfragments, function(i) {i[[1]]}))
  nSeedMZfragments <- length(seedMZfragments)
  ##
  if (nSeedMZfragments > 0) {
    ##
    idSeedMZfragments <- rep(0, nPeaks)
    counterSeedMZ <- 0
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = LxDiffMZ, initial = 0, style = 3)
    ##
    for (i in 1:LxDiffMZ) {
      if (!is.null(listSeedMZfragments[[i]][[1]])) {
        ##
        sbIdSeedMZfragments <- listSeedMZfragments[[i]][[2]] + counterSeedMZ
        ##
        idSeedMZfragments[seq((xDiffMZ[i] + 1), xDiffMZ[i + 1], 1)] <- sbIdSeedMZfragments
        ##
        counterSeedMZ <- max(sbIdSeedMZfragments)
      }
      setTxtProgressBar(progressBARboundaries, i)
    }
    ##
    close(progressBARboundaries)
    ##
    listSeedMZfragments <- NULL
    sbIdSeedMZfragments <- NULL
    ##
    mainImzInt <- cbind(mainImzInt, idSeedMZfragments)
    mainImzInt <- mainImzInt[(mainImzInt[, 4] != 0), ]
    ##
    FSA_logRecorder("Completed aligning the fragmentation spectra!")
    ##
    ############################################################################
    ############################################################################
    ##
    FSA_logRecorder(paste0("Initiated saving the `", dimension, "` aligned fragmentation spectra table!"))
    ##
    if (dimension == "long") {
      mainImzInt <- mainImzInt[order(mainImzInt[, 4], decreasing = FALSE), ]
      colnames(mainImzInt) <- c("spectraID", "mz", "intensity", "ionID")
      ##
      write.csv(mainImzInt, file = alignedMSPcsvFile, row.names = FALSE)
      ##
    } else if (dimension == "wide") {
      ##
      mainImzInt <- mainImzInt[order(mainImzInt[, 1], decreasing = FALSE), ]
      xDiff <- c(0, which(diff(mainImzInt[, 1]) > 0), nrow(mainImzInt))
      LxDiff <- length(xDiff) - 1
      ##
      orderSeedMZfragments <- order(seedMZfragments, decreasing = FALSE)
      headerRowSeedMZfragments <- paste0("spectraID,", paste0(seedMZfragments[orderSeedMZfragments], collapse = ","))
      ##
      write(headerRowSeedMZfragments, file = alignedMSPcsvFile, append = FALSE, sep = "\n")
      ##
      progressBARboundaries <- txtProgressBar(min = 0, max = LxDiff, initial = 0, style = 3)
      for (i in 1:LxDiff) {
        xID <- seq((xDiff[i] + 1), xDiff[i + 1], 1)  
        fragIDi <- rep("", nSeedMZfragments)
        fragIDi[mainImzInt[xID, 4]] <- as.character(mainImzInt[xID, 3])
        ##
        fragIDi <- paste0(mainImzInt[xID[1], 1], ",", paste0(fragIDi[orderSeedMZfragments], collapse = ","))
        ##
        write(fragIDi, file = alignedMSPcsvFile, append = TRUE, sep = "\n")
        ##
        setTxtProgressBar(progressBARboundaries, i)
      }
      close(progressBARboundaries)
    }
    ##
    FSA_logRecorder(paste0("Completed saving the `", dimension, "` aligned fragmentation spectra table!"))
    ##
  } else {
    stop(FSA_logRecorder("No common ion was detected!"))
  }
  return()
}