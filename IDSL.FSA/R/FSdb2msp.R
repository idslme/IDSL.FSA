FSdb2msp <- function(path, FSdbFileName = "", UnweightMSP = FALSE, number_processing_threads = 1) {
  ##
  FSdbFileLocation <- paste0(path, "/", FSdbFileName)
  strFSdbFileLocation <- strsplit(FSdbFileLocation, "/")[[1]]
  FSdbFileName <- strFSdbFileLocation[length(strFSdbFileLocation)]
  FSdbFileLocation <- paste0(strFSdbFileLocation, collapse = "/")
  ##
  FSdb <- FSA_loadRdata(FSdbFileLocation)
  ##
  massWindowIntegration <- FSdb[["logFSdb"]][["massWindowIntegration"]]
  allowedWeightedSpectralEntropy <- FSdb[["logFSdb"]][["allowedWeightedSpectralEntropy"]]
  noiseRemovalPercentage <- FSdb[["logFSdb"]][["noiseRemovalPercentage"]]
  ##
  LNumPeaks <- length(FSdb[["Num Peaks"]])
  ##
  rawIntensityValues <- FSdb[["MSPLibraryParameters"]][["basepeakintensity"]]
  if (is.null(rawIntensityValues)) {
    rawIntensityValues <- rep(1000, LNumPeaks) # 1000 is a standard NIST normalization absolute value
  } else {
    rawIntensityValues <- as.numeric(rawIntensityValues)
  }
  ##
  mspProprties <- colnames(FSdb[["MSPLibraryParameters"]])
  x_name <- which(mspProprties == "name")
  x_precursormz <- which(mspProprties == "precursormz")
  x_ex <- setdiff(seq(1, length(mspProprties), 1), c(x_name, x_precursormz))
  mspProprties_colon <- paste0(mspProprties, ": ")
  ##
  ##############################################################################
  ##
  if (allowedWeightedSpectralEntropy & UnweightMSP) {
    allowedWeightedSpectralEntropy <- FALSE
    ##
    call_S_FragmentList <- function(i) {
      FragList <- FSdb[["FragmentList"]][[i]]
      S <- FSdb[["Spectral Entropy"]][i]
      SFL <- list(S, FragList)
      if (S > 3) {
        w <- 0.25 + S * 0.25
        ##
        FragList[, 2] <- FragList[, 2]^(1/w)
        SFL <- spectral_entropy_calculator(FragList, allowedWeightedSpectralEntropy = FALSE, noiseRemovalRatio = 0)
      }
      ##
      return(SFL)
    }
  }
  ##
  call_MSP <- function(i) {
    mspBlockName <- paste0("Name: ", FSdb[["MSPLibraryParameters"]][i, x_name], "\n")
    ##
    if (!is.infinite(FSdb[["PrecursorMZ"]][i])) {
      mspBlockPrecursorMZ <- paste0("PrecursorMZ: ", FSdb[["PrecursorMZ"]][i], "\n")
    } else {
      mspBlockPrecursorMZ <- NULL
    }
    ##
    mspBlockMSPproprties <- do.call(c, lapply(x_ex, function(j) {
      paste0(mspProprties_colon[j], FSdb[["MSPLibraryParameters"]][i, j])
    }))
    mspBlockMSPproprties <- paste0(mspBlockMSPproprties, collapse = "\n")
    ##
    mspBlockSE1 <- paste0("spectral_entropy: ", FSdb[["Spectral Entropy"]][i], "\n")
    ##
    mspBlockSE2 <- paste0("Weighted_Spectral_Entropy_Transformation: ", allowedWeightedSpectralEntropy, "\n")
    ##
    mspBlocklogFSdb <- paste0("massWindowIntegration: ", massWindowIntegration, " Da\n", "noiseRemovalPercentage: ", noiseRemovalPercentage, " percent\n")
    ##
    mspBlockNumPeak <- paste0("Num Peaks: ", FSdb[["Num Peaks"]][i], "\n")
    ##
    FList <- FSdb[["FragmentList"]][[i]]
    mspBlockFList <- paste0(FList[, 1], " ", round(FList[, 2]/max(FList[, 2])*rawIntensityValues[i], digits = 0), collapse = "\n")
    ##
    paste0(mspBlockName, mspBlockPrecursorMZ, mspBlockMSPproprties, "\n", mspBlockSE1, mspBlockSE2, mspBlocklogFSdb, mspBlockNumPeak, mspBlockFList, "\n\n")
  }
  ##
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    if (UnweightMSP) {
      S_FragmentList <- lapply(1:LNumPeaks, function(i) {
        call_S_FragmentList(i)
      })
      ##
      SpectralEntropy <- do.call(c, lapply(1:LNumPeaks, function(i) {
        S_FragmentList[[i]][[1]]
      }))
      ##
      FragmentList <- lapply(1:LNumPeaks, function(i) {
        S_FragmentList[[i]][[3]]
      })
    }
    ##
    MSP <- do.call(c, lapply(1:LNumPeaks, function(i) {
      call_MSP(i)
    }))
  } else {
    osType <- Sys.info()[['sysname']]
    if (osType == "Linux") {
      if (UnweightMSP) {
        S_FragmentList <- mclapply(1:LNumPeaks, function(i) {
          call_S_FragmentList(i)
        }, mc.cores = number_processing_threads)
        ##
        SpectralEntropy <- do.call(c, mclapply(1:LNumPeaks, function(i) {
          S_FragmentList[[i]][[1]]
        }, mc.cores = number_processing_threads))
        ##
        FragmentList <- lapply(1:LNumPeaks, function(i) {
          S_FragmentList[[i]][[3]]
        }, mc.cores = number_processing_threads)
      }
      ##
      MSP <- do.call(c, mclapply(1:LNumPeaks, function(i) {
        call_MSP(i)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      if (UnweightMSP) {
        S_FragmentList <- foreach(i = 1:LNumPeaks, .verbose = FALSE) %dopar% {
          call_S_FragmentList(i)
        }
        ##
        SpectralEntropy <- foreach(i = 1:LNumPeaks, .combine = 'c', .verbose = FALSE) %dopar% {
          S_FragmentList[[i]][[1]]
        }
        ##
        FragmentList <- foreach(i = 1:LNumPeaks, .verbose = FALSE) %dopar% {
          S_FragmentList[[i]][[3]]
        }
      }
      ##
      MSP <- foreach(i = 1:LNumPeaks, .combine = 'c', .verbose = FALSE) %dopar% {
        call_MSP(i)
      }
      ##
      stopCluster(clust)
    }
  }
  ##
  ##############################################################################
  ##
  mspFileName <- gsub("[.]Rdata$", ".msp", FSdbFileLocation, ignore.case = TRUE)
  ##
  write.table(MSP, file = mspFileName, quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
  ##
  FSA_logRecorder(paste0("`", FSdbFileName, "` was convereted into the .msp format and stored as `", mspFileName, "`!"))
  ##
  return()
}