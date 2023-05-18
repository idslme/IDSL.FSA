FSA_uniqueMSPblockTaggerUntargeted <- function(path, MSPfile_vector, minCSAdetectionFrequency = 20, minEntropySimilarity = 0.75,
                                               massError = 0.01, massErrorPrecursor = 0.01, RTtolerance = 0.1, noiseRemovalRatio = 0.01,
                                               allowedNominalMass = FALSE, allowedWeightedSpectralEntropy = TRUE, plotSpectra = FALSE,
                                               number_processing_threads = 1) {
  ##
  ##############################################################################
  ##
  if (is.na(massErrorPrecursor)) {
    precursorMZcheck <- FALSE
  } else {
    precursorMZcheck <- TRUE
  }
  ##
  ##############################################################################
  ##
  FSdb <- msp2FSdb(path, MSPfile_vector, massIntegrationWindow = massError, allowedNominalMass, allowedWeightedSpectralEntropy, noiseRemovalRatio, number_processing_threads)
  ##
  basePeakMZ <- as.numeric(FSdb[["MSPLibraryParameters"]][["basepeakmz"]])
  if (!is.null(basePeakMZ)) {
    ##
    ############################################################################
    ##
    basePeakIntensity <- as.numeric(FSdb[["MSPLibraryParameters"]][["basepeakintensity"]])
    orderBasePeakIntensity <- order(basePeakIntensity, decreasing = TRUE)
    ##
    nSpectra <- length(basePeakIntensity)
    spectraIDs <- seq(1, nSpectra, 1)
    res_list <- vector(mode = "list", length = nSpectra)
    res_list_len <- rep(0, nSpectra)
    ##
    progressBARboundaries <- txtProgressBar(min = 0, max = nSpectra, initial = 0, style = 3)
    ##
    i <- 0
    for (k in orderBasePeakIntensity) {
      ##
      if (spectraIDs[k] > 0) {
        if (allowedNominalMass) {
          if (precursorMZcheck) {
            ID <- which(abs(FSdb[["Retention Time"]] - FSdb[["Retention Time"]][k]) <= RTtolerance &
                          abs(FSdb[["Spectral Entropy"]] - FSdb[["Spectral Entropy"]][k]) <= 2 &
                          (basePeakMZ == basePeakMZ[k]) &
                          (FSdb[["PrecursorMZ"]] == FSdb[["PrecursorMZ"]][k]))
          } else {
            ID <- which(abs(FSdb[["Retention Time"]] - FSdb[["Retention Time"]][k]) <= RTtolerance &
                          abs(FSdb[["Spectral Entropy"]] - FSdb[["Spectral Entropy"]][k]) <= 2 &
                          (basePeakMZ == basePeakMZ[k]))
          }
        } else {
          if (precursorMZcheck) {
            ID <- which(abs(FSdb[["Retention Time"]] - FSdb[["Retention Time"]][k]) <= RTtolerance &
                          abs(FSdb[["Spectral Entropy"]] - FSdb[["Spectral Entropy"]][k]) <= 2 &
                          abs(basePeakMZ - basePeakMZ[k]) <= massError &
                          abs(FSdb[["PrecursorMZ"]] - FSdb[["PrecursorMZ"]][k]) <= massErrorPrecursor)
          } else {
            ID <- which(abs(FSdb[["Retention Time"]] - FSdb[["Retention Time"]][k]) <= RTtolerance &
                          abs(FSdb[["Spectral Entropy"]] - FSdb[["Spectral Entropy"]][k]) <= 2 &
                          abs(basePeakMZ - basePeakMZ[k]) <= massError)
          }
        }
        ##
        LID <- length(ID)
        if (LID == 1) {
          uniqueMSPblock <- ID
        } else {
          ##
          S_PEAK_B <- FSdb[["Spectral Entropy"]][k]
          PEAK_B <- FSdb[["FragmentList"]][[k]]
          IDspectralSimilarity <- do.call(c, lapply(2:LID, function(j) {
            S_PEAK_A <- FSdb[["Spectral Entropy"]][ID[j]]
            PEAK_A <- FSdb[["FragmentList"]][[ID[j]]]
            ##
            SESS <- spectral_entropy_similarity_score(PEAK_A, S_PEAK_A, PEAK_B, S_PEAK_B, massError, allowedNominalMass)
            ##
            if (SESS >= minEntropySimilarity) {
              ID[j]
            }
          }))
          uniqueMSPblock <- c(k, IDspectralSimilarity)
        }
        ##
        res_list[[k]] <- uniqueMSPblock # The order is same as the FSDB
        spectraIDs[spectraIDs %in% uniqueMSPblock] <- 0
        ##
        res_list_len[k] <- length(uniqueMSPblock)
      }
      ##
      i <- i + 1
      setTxtProgressBar(progressBARboundaries, i)
    }
    ##
    close(progressBARboundaries)
    ##
    ############################################################################
    ##
    selectedFSdbIDs <- which(res_list_len >= minCSAdetectionFrequency)
    if (length(selectedFSdbIDs) > 0) {
      ##
      ##########################################################################
      ##########################################################################
      ##
      if (plotSpectra) {
        ##
        outputUniqueSpectra <- paste0(path, "/UNIQUETAGS/uniqueTagsSpectra/")
        FSA_dir.create(outputUniqueSpectra, allowedUnlink = TRUE)
        FSA_logRecorder(paste0("Unique tag spectra figures are stored in the `", outputUniqueSpectra,"` folder!"))
        ##
        dev.offCheck <- TRUE
        while (dev.offCheck) {
          dev.offCheck <- tryCatch(dev.off(), error = function(e) {FALSE})
        }
        ##
        call_plotUniqueTag <- function(i) {
          ##
          variantFolder <- paste0(outputUniqueSpectra, "/variant_", i, "_RT_", FSdb[["Retention Time"]][i])
          ##
          FSA_dir.create(variantFolder, allowedUnlink = TRUE)
          ##
          for (j in res_list[[i]]) {
            filenameUniqueTag <- paste0(variantFolder, "/", j, "_", FSdb[["MSPLibraryParameters"]][["MSPfilename"]][j], "_.png")
            fileCreateRCheck <- file.create(file = filenameUniqueTag, showWarnings = FALSE)
            if (fileCreateRCheck) {
              png(filenameUniqueTag, width = 16, height = 8, units = "in", res = 100)
              ##
              plotFSdb2SpectraCore(FSdb, index = j)
              ##
              dev.off()
            } else {
              FSA_logRecorder(paste0("WARNING!!! Figure can not be created for `", filenameUniqueTag, "` due to character length limit in the `", variantFolder, "`!"))
            }
          }
          return()
        }
        ##
        ########################################################################
        ##
        if (number_processing_threads == 1) {
          ##
          progressBARboundaries <- txtProgressBar(min = 0, max = nSpectra, initial = 0, style = 3)
          ##
          for (i in selectedFSdbIDs) {
            setTxtProgressBar(progressBARboundaries, i)
            ##
            null_variable <- call_plotUniqueTag(i)
          }
          ##
          close(progressBARboundaries)
          ##
        } else {
          ## Processing OS
          osType <- Sys.info()[['sysname']]
          ##
          ######################################################################
          ##
          if (osType == "Windows") {
            clust <- makeCluster(number_processing_threads)
            clusterExport(clust, setdiff(ls(), c("clust", "selectedFSdbIDs")), envir = environment())
            ##
            null_variable <- parLapply(clust, selectedFSdbIDs, function(i) {
              call_plotUniqueTag(i)
            })
            ##
            stopCluster(clust)
            ##
            ####################################################################
            ##
          } else {
            ##
            null_variable <- mclapply(selectedFSdbIDs, function(i) {
              call_plotUniqueTag(i)
            }, mc.cores = number_processing_threads)
            ##
            closeAllConnections()
            ##
            ####################################################################
            ##
          }
        }
      }
      ##
      ##########################################################################
      ##########################################################################
      ##
      FSdb[["MSPLibraryParameters"]][["CSAvariantDetectionFreq"]] <- res_list_len
      ##
      uniqueMSPvariants <- FSdb_subsetter(FSdb, inclusionIDs = selectedFSdbIDs)
      ##
      FSA_dir.create(paste0(path, "/UNIQUETAGS/"), allowedUnlink = FALSE)
      ##
      save(uniqueMSPvariants, file = paste0(path, "/UNIQUETAGS/uniqueMSPtagsUntargeted.Rdata"))
      ##
      FSdb2msp(path = paste0(path,"/UNIQUETAGS"), FSdbFileName = "uniqueMSPtagsUntargeted.Rdata", UnweightMSP = FALSE, number_processing_threads)
    } else {
      FSA_logRecorder(paste0("No common msp blocks was found using absolute frequency of detection >= `", minCSAdetectionFrequency,"` in the untargeted unique tag analysis!"))
    }
  } else {
    FSA_logRecorder("The `basePeakIntensity` meta-variable was not found in the .msp files!")
  }
  return()
}
