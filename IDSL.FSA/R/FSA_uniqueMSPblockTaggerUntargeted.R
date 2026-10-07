FSA_uniqueMSPblockTaggerUntargeted <- function(path, MSPfile_vector = "", peak_alignment_folder = NA, minCSAdetectionFrequency = 20, massError = 0.01, massErrorPrecursor = 0.01,
                                               RTtolerance = 0.05, minEntropySimilarity = 0.75, noiseRemovalRatio = 0.01, minCosineSimilarity = 0.75, allowedNominalMass = FALSE,
                                               allowedWeightedSpectralEntropy = TRUE, plotSpectra = FALSE, number_processing_threads = 1) {
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
    nFSdb <- nrow(FSdb[["MSPLibraryParameters"]])
    ##
    ############################################################################
    ## Updating Retention Times
    if (!is.na(peak_alignment_folder)) {
      msp_mode = tolower(FSdb[["MSPLibraryParameters"]][["msp_mode"]][1])
      if (is.character(msp_mode) && !is.na(msp_mode) && nzchar(msp_mode)) {
        ##
        peakType = NA
        if (msp_mode == "csa") {
          peakType = "idsl.ipa_collective_peakids"
        } else if ((msp_mode == "dda") || (msp_mode == "dia")) {
          peakType = "idsl.ipa_peakid"
        }
        if (is.character(peakType) && !is.na(peakType) && nzchar(peakType)) {
          if (any(colnames(FSdb[["MSPLibraryParameters"]]) == peakType)) {
            ##
            listCorrectedRTpeaklists <- FSA_loadRdata(paste0(peak_alignment_folder, "/listCorrectedRTpeaklists.Rdata"))
            names(listCorrectedRTpeaklists) <- gsub("^peaklist_|.Rdata$", "", names(listCorrectedRTpeaklists), ignore.case = TRUE)
            ##
            if (all(grepl("uniqueMSPtagsUntargeted.msp", FSdb[["MSPLibraryParameters"]][["MSPfilename"]]))) {
              mspfilename <- "mspfilename"
            } else {
              mspfilename <- "MSPfilename"
            }
            mzMLfilename <- gsub("^DDA_MSP_|^DDA_REF_MSP_|^DIA_MSP_|^DIA_REF_MSP_|^CSA_MSP_|^CSA_REF_MSP_|.msp$", "", FSdb[["MSPLibraryParameters"]][[mspfilename]], ignore.case = TRUE)
            tmzMLfilename <- base::tapply(seq(1, nFSdb, 1), mzMLfilename, FUN = 'c', simplify = FALSE)
            ##
            ####################################################################
            ## Processing OS
            osType <- Sys.info()[['sysname']]
            ##
            if ((number_processing_threads == 1) || (osType == "Windows")) {
              IPA_collective_peakids <- lapply(strsplit(FSdb[["MSPLibraryParameters"]][[peakType]], ","), function(k) {as.numeric(k[k != "0"])})
            } else {
              IPA_collective_peakids <- mclapply(strsplit(FSdb[["MSPLibraryParameters"]][[peakType]], ","), function(k) {as.numeric(k[k != "0"])}, mc.cores = number_processing_threads)
            }
            ##
            ####################################################################
            ##
            for (mzML in names(tmzMLfilename)) {
              CorrectedRTs <- listCorrectedRTpeaklists[[mzML]]
              for (k in tmzMLfilename[[mzML]]) {
                if (length(IPA_collective_peakids[[k]]) > 0) {
                  FSdb[["Retention Time"]][k] <- median(CorrectedRTs[IPA_collective_peakids[[k]]])
                }
              }
            }
          }
        }
      }
    }
    ##
    ############################################################################
    ##
    basePeakIntensity <- as.numeric(FSdb[["MSPLibraryParameters"]][["basepeakintensity"]])
    orderBasePeakIntensity <- order(basePeakIntensity, decreasing = TRUE)
    ##
    nSpectra <- length(basePeakIntensity)
    spectraIDs <- seq(1, nSpectra, 1)
    clusterIndices <- vector(mode = "list", length = nSpectra)
    lengthClusterIndices <- rep(0, nSpectra)
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
        if (length(ID) == 1) {
          uniqueMSPblock <- ID
        } else {
          ##
          ID = setdiff(ID, k)
          ##
          S_PEAK_A <- FSdb[["Spectral Entropy"]][k]
          PEAK_A <- FSdb[["FragmentList"]][[k]]
          L_A <- FSdb[["Num Peaks"]][k]
          sumPEAK_A2 <- sum(PEAK_A[, 2]^2)
          ##
          IDspectralSimilarity <- do.call(c, lapply(ID, function(j) {
            S_PEAK_B <- FSdb[["Spectral Entropy"]][j]
            PEAK_B <- FSdb[["FragmentList"]][[j]]
            ##
            SESS <- spectral_entropy_similarity_score(PEAK_A, S_PEAK_A, PEAK_B, S_PEAK_B, massError, allowedNominalMass)
            ##
            if (SESS >= minEntropySimilarity) {
              ##
              matchedPeak_B <- matrix(0, nrow = L_A, ncol = 2)
              for (f in 1:L_A) {
                x_f <- which(abs(PEAK_A[f, 1] - PEAK_B[, 1]) <= massError)
                L_x_f <- length(x_f)
                if (L_x_f > 0) {
                  if (L_x_f > 1) {
                    x_min <- which.min(abs(PEAK_A[f, 1] - PEAK_B[x_f, 1]))
                    x_f <- x_f[x_min[1]]
                  }
                  matchedPeak_B[f, ] <- PEAK_B[x_f, ]
                }
              }
              CS <- sum(matchedPeak_B[, 2]*PEAK_A[, 2])/sqrt(sum(matchedPeak_B[, 2]^2)*sumPEAK_A2)
              ##
              if (CS >= minCosineSimilarity) {
                j
              }
            }
          }))
          ##
          uniqueMSPblock <- c(k, IDspectralSimilarity)
        }
        ##
        clusterIndices[[k]] <- uniqueMSPblock # The order is same as the FSDB
        spectraIDs[spectraIDs %in% uniqueMSPblock] <- 0
        ##
        lengthClusterIndices[k] <- length(uniqueMSPblock)
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
    if (minCSAdetectionFrequency == 0) {
      minCSAdetectionFrequency <- 1e-16
    }
    ##
    selectedFSdbIDs <- which(lengthClusterIndices >= minCSAdetectionFrequency)
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
          for (j in clusterIndices[[i]]) {
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
      FSdb[["MSPLibraryParameters"]][["CSAvariantDetectionFreq"]] <- lengthClusterIndices
      FSdb[["MSPLibraryParameters"]][["fsa_unique_id"]] <- seq(1, nFSdb, 1)
      ##
      uniqueMSPvariants <- FSdb_subsetter(FSdb, inclusionIDs = selectedFSdbIDs)
      rownames(uniqueMSPvariants[["MSPLibraryParameters"]]) <- NULL
      ##
      FSA_dir.create(paste0(path, "/UNIQUETAGS/"), allowedUnlink = FALSE)
      ##
      peakType = NA
      mzMLfilename <- NA
      ##
      msp_mode = tolower(FSdb[["MSPLibraryParameters"]][["msp_mode"]][1])
      if (is.character(msp_mode) && !is.na(msp_mode) && nzchar(msp_mode)) {
        ##
        if (msp_mode == "csa") {
          peakType = "idsl.ipa_collective_peakids"
        } else if ((msp_mode == "dda") || (msp_mode == "dia")) {
          peakType = "idsl.ipa_peakid"
        }
        if (is.character(peakType) && !is.na(peakType) && nzchar(peakType)) {
          if (any(colnames(FSdb[["MSPLibraryParameters"]]) == peakType)) {
            ##
            mzMLfilename <- gsub("^DDA_MSP_|^DDA_REF_MSP_|^DIA_MSP_|^DIA_REF_MSP_|^CSA_MSP_|^CSA_REF_MSP_|.msp$", "", FSdb[["MSPLibraryParameters"]][["MSPfilename"]], ignore.case = TRUE)
          }
        } else {
          FSA_logRecorder(paste0("`", peakType,"` was not found in the msp files!"))
        }
      }
      ##
      ##########################################################################
      ##
      uniqueMSPvariants[["uniqueClusterInfo"]] <- list(minCSAdetectionFrequency = minCSAdetectionFrequency,
                                                       ClusterIndices = clusterIndices,
                                                       lengthClusterIndices = lengthClusterIndices,
                                                       mzMLfilename = mzMLfilename,
                                                       CollectivePeakIDs = FSdb[["MSPLibraryParameters"]][[peakType]])
      ##
      save(uniqueMSPvariants, file = paste0(path, "/UNIQUETAGS/uniqueMSPtagsUntargeted.Rdata"))
      ##
      FSdb2msp(path = paste0(path,"/UNIQUETAGS"), FSdbFileName = "uniqueMSPtagsUntargeted.Rdata", UnweightMSP = FALSE, number_processing_threads)
      ##
      ##
      ##########################################################################
      ##########################################################################
      ##
    } else {
      FSA_logRecorder(paste0("No common msp blocks was found using absolute frequency of detection >= `", minCSAdetectionFrequency,"` in the untargeted unique tag analysis!"))
    }
  } else {
    FSA_logRecorder("The `basePeakIntensity` meta-variable was not found in the .msp files!")
  }
  return(NULL)
}
