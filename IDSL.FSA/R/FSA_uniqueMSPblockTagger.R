FSA_uniqueMSPblockTagger <- function(path, MSPfile = "", aggregateBy = "Name", massError = 0.01, RTtolerance = NA, minEntropySimilarity = 0.75,
                                     noiseRemovalRatio = 0.01, allowedNominalMass = FALSE, allowedWeightedSpectralEntropy = TRUE, plotSpectra = FALSE,
                                     number_processing_threads = 1) {
  ##
  listSimilarMSPvariants <- NULL
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
  strmspFileName <- strsplit(mspFileName, "[.]")[[1]]
  formatFileName <- tolower(strmspFileName[length(strmspFileName)])
  if (formatFileName == "msp") {
    FSdb <- msp2FSdb(path = mspFilePath, MSPfile_vector = mspFileName, massIntegrationWindow = massError,
                     allowedNominalMass, allowedWeightedSpectralEntropy, noiseRemovalRatio, number_processing_threads)
  } else if (formatFileName == "rdata") {
    FSdb <- FSA_loadRdata(paste0(mspFilePath, "/", mspFileName))
  } else {
    stop("Allowed inputs : FSDB in `.Rdata` or MSP files in `.msp` formats!")
  }
  ##
  precursorINT <- as.numeric(FSdb[["MSPLibraryParameters"]][["precursor_intensity"]])
  ##
  if (!is.null(precursorINT)) {
    ##
    colMSPLibraryParameters <- colnames(FSdb[["MSPLibraryParameters"]])
    xColAgg <- which(colMSPLibraryParameters == tolower(aggregateBy))
    if (length(xColAgg) > 0) {
      ##
      listID <- FSA_R.aggregate(FSdb[["MSPLibraryParameters"]][[tolower(aggregateBy)]])
      uniqueAggCompound <- names(listID)
      ##
      if (!is.na(RTtolerance)) {
        RTcheck <- TRUE
      } else {
        RTcheck <- FALSE
      }
      ##
      ##########################################################################
      ##
      if (plotSpectra) {
        outputUniqueSpectra <- paste0(path, "/uniqueTagsSpectra/")
        FSA_dir.create(outputUniqueSpectra, allowedUnlink = TRUE)
        FSA_logRecorder(paste0("Unique tag spectra figures are stored in the `", outputUniqueSpectra,"` folder!"))
        ##
        dev.offCheck <- TRUE
        while (dev.offCheck) {
          dev.offCheck <- tryCatch(dev.off(), error = function(e) {FALSE})
        }
      }
      ##
      ##########################################################################
      ##
      call_listSimilarMSPvariants <- function(i) {
        ##
        ID <- listID[[i]]
        LID <- length(ID)
        ##
        if (LID == 1) {
          uniqueMSPblock <- list(ID)
        } else {
          ##
          IDspectralSimilarity <- do.call(rbind, lapply(1:(LID - 1), function(j) {
            ##
            RTj <- FSdb[["Retention Time"]][ID[j]]
            if (RTcheck) {
              RTjInfCheck <- if (!is.infinite(RTj)) {TRUE} else {FALSE}
            }
            ##
            S_PEAK_A <- FSdb[["Spectral Entropy"]][ID[j]]
            PEAK_A <- FSdb[["FragmentList"]][[ID[j]]]
            ##
            do.call(rbind, lapply((j + 1):LID, function(k) {
              S_PEAK_B <- FSdb[["Spectral Entropy"]][ID[k]]
              ##
              if (abs(S_PEAK_B - S_PEAK_A) <= 2) {
                ##
                RTk <- FSdb[["Retention Time"]][ID[k]]
                ##
                passLoopCheck <- TRUE
                if (RTcheck) {
                  if (RTjInfCheck) {
                    if (!is.infinite(RTk)) {
                      if (abs(RTj - RTk) > RTtolerance) {
                        passLoopCheck <- FALSE
                      }
                    } else {
                      passLoopCheck <- FALSE
                    }
                  } else {
                    passLoopCheck <- FALSE
                  }
                }
                ##
                if (passLoopCheck) {
                  ##
                  PEAK_B <- FSdb[["FragmentList"]][[ID[k]]]
                  ##
                  SESS <- spectral_entropy_similarity_score(PEAK_A, S_PEAK_A, PEAK_B, S_PEAK_B, massError, allowedNominalMass)
                  ##
                  if (SESS >= minEntropySimilarity) {
                    c(ID[j], ID[k])
                  }
                }
              }
            }))
          }))
          ##
          IDspectralSimilarity <- rbind(IDspectralSimilarity, cbind(ID, ID))
          ##
          orderINT <- order(precursorINT[ID], decreasing = TRUE)
          ##
          uniqueMSPblock <- vector(mode = "list", LID)
          counter <- 0
          for (j in orderINT) {
            xj <- which((IDspectralSimilarity[, 1] == ID[j]) | (IDspectralSimilarity[, 2] == ID[j]))
            if (length(xj) > 0) {
              counter <- counter + 1
              uniqueMSPblock[[counter]] <- unique(c(IDspectralSimilarity[xj, 1], IDspectralSimilarity[xj, 2]))
              ##
              IDspectralSimilarity[IDspectralSimilarity[, 1] %in% uniqueMSPblock[[counter]], ] <- 0
              IDspectralSimilarity[IDspectralSimilarity[, 2] %in% uniqueMSPblock[[counter]], ] <- 0
            }
          }
          uniqueMSPblock <- uniqueMSPblock[1:counter]
        }
        ##
        ########################################################################
        ##
        if (plotSpectra) {
          uniqueTagAggSpectra <- gsub('/|[\\]|\t|\n|:|[*]|[?]|<|>|[|]|[.][.][.]|[.][.]|^[.]', '_', i, fixed = FALSE) # To replace invalid characters
          uniqueTagAggSpectra <- paste0(outputUniqueSpectra, "/", uniqueTagAggSpectra)
          ##
          FSA_dir.create(uniqueTagAggSpectra, allowedUnlink = TRUE)
          ##
          jCounter <- 0
          for (q in uniqueMSPblock) {
            jCounter <- jCounter + 1
            ##
            variantFolder <- paste0(uniqueTagAggSpectra, "/variant_", jCounter)
            ##
            FSA_dir.create(variantFolder, allowedUnlink = TRUE)
            ##
            for (j in q) {
              filenameUniqueTag <- paste0(variantFolder, "/", j, "_", FSdb[["MSPLibraryParameters"]][["filename"]][j], "_.png")
              fileCreateRCheck <- file.create(file = filenameUniqueTag, showWarnings = FALSE)
              if (fileCreateRCheck) {
                png(filenameUniqueTag, width = 16, height = 8, units = "in", res = 100)
                ##
                plotFSdb2SpectraCore(FSdb, index = j)
                ##
                dev.off()
              }
            }
          }
        }
        ##
        ########################################################################
        ##
        return(uniqueMSPblock)
      }
      ##
      ##########################################################################
      ##
      call_indexUniqueMSPvariants <- function(i) {
        uniqueMSPblock <- listSimilarMSPvariants[[i]]
        ##
        do.call(c, lapply(uniqueMSPblock, function(j) {
          j[1]
        }))
      }
      ##
      ##########################################################################
      ##
      if (number_processing_threads == 1) {
        ##
        listSimilarMSPvariants <- lapply(uniqueAggCompound, function(i) {
          call_listSimilarMSPvariants(i)
        })
        ##
        names(listSimilarMSPvariants) <- uniqueAggCompound
        ##
        indexUniqueMSPvariants <- do.call(c, lapply(uniqueAggCompound, function(i) {
          call_indexUniqueMSPvariants(i)
        }))
        ##
      } else {
        ##
        ########################################################################
        ##
        osType <- Sys.info()[['sysname']]
        ##
        if (osType == "Windows") {
          ####
          clust <- makeCluster(number_processing_threads)
          clusterExport(clust, setdiff(ls(), c("clust", "uniqueAggCompound")), envir = environment())
          ##
          listSimilarMSPvariants <- parLapply(clust, uniqueAggCompound, function(i) {
            call_listSimilarMSPvariants(i)
          })
          ##
          stopCluster(clust)
          ####
          names(listSimilarMSPvariants) <- uniqueAggCompound
          ####
          clust <- makeCluster(number_processing_threads)
          clusterExport(clust, setdiff(ls(), c("clust", "uniqueAggCompound")), envir = environment())
          ##
          indexUniqueMSPvariants <- do.call(c, parLapply(clust, uniqueAggCompound, function(i) {
            call_indexUniqueMSPvariants(i)
          }))
          ##
          stopCluster(clust)
          ##
          ######################################################################
          ##
        } else {
          ##
          listSimilarMSPvariants <- mclapply(uniqueAggCompound, function(i) {
            call_listSimilarMSPvariants(i)
          }, mc.cores = number_processing_threads)
          ##
          names(listSimilarMSPvariants) <- uniqueAggCompound
          ##
          indexUniqueMSPvariants <- do.call(c, mclapply(uniqueAggCompound, function(i) {
            call_indexUniqueMSPvariants(i)
          }, mc.cores = number_processing_threads))
          ##
          closeAllConnections()
          ##
          ######################################################################
          ##
        }
      }
      ##
      ##########################################################################
      ##
      if (length(indexUniqueMSPvariants) > 0) {
        ##
        LnumPeaks <- length(FSdb[["Num Peaks"]])
        CSAvariantIDs <- rep("", LnumPeaks)
        CSAvariantFreq <- rep(0, LnumPeaks)
        ##
        for (i in listSimilarMSPvariants) {
          for (j in i) {
            indicesCSA <- paste0(j, collapse = ",")
            nIndicesCSA <- length(j)
            for (k in j) {
              CSAvariantIDs[k] <- indicesCSA
              CSAvariantFreq[k] <- nIndicesCSA
            }
          }
        }
        FSdb[["MSPLibraryParameters"]] <- cbind(CSAvariantIDs, CSAvariantFreq, FSdb[["MSPLibraryParameters"]])
        ##
        uniqueMSPvariants <- FSdb_subsetter(FSdb, inclusionIDs = indexUniqueMSPvariants)
        FSdbFileName <- paste0("uniqueMSPtags_", gsub("[.]msp$|[.]Rdata$", ".Rdata", mspFileName, ignore.case = TRUE))
        save(uniqueMSPvariants, file = paste0(mspFilePath, "/", FSdbFileName))
        ##
        FSdb2msp(path = mspFilePath, FSdbFileName = FSdbFileName, UnweightMSP = FALSE, number_processing_threads)
        ##
      } else {
        FSA_logRecorder(paste0("No common msp blocks was found in the unique tag analysis!"))
      }
    } else {
      FSA_logRecorder(paste0("The `aggregateBy = ", aggregateBy, "` variable is not available in the msp file!!!"))
    }
  } else {
    FSA_logRecorder("NOTICE: exclusion lists can not be generated because precursor intensity values were not detected in the .msp file!!!")
  }
  return(listSimilarMSPvariants)
}