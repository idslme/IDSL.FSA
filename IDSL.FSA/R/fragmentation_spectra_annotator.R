fragmentation_spectra_annotator <- function(path, MSPfile = "", libFSdb, libFSdbIDlist, targetedPrecursorType = NA, ratio2basePeak4nSpectraMarkers = 0,
                                            allowedNominalMass = FALSE, allowedWeightedSpectralEntropy = TRUE, noiseRemovalRatio = 0.01, roundingDigitPrefiltering = 1,
                                            minMatchedNumPeaks = 1, massError = 0, maxNEME = 0, minIonRangeDifference = 0, minCosineSimilarity, minEntropySimilarity,
                                            minRatioMatchedNspectraMarkers, spectralEntropyDeviationPrefiltering, massErrorPrecursor = NA, RTtolerance = NA,
                                            exportSpectraParameters = NULL, number_processing_threads = 1) {
  ##
  mspMatch <- NULL
  ##
  ##############################################################################
  ##############################################################################
  ##
  sampleFSdb <- msp2FSdb(path, MSPfile, massError, allowedNominalMass, allowedWeightedSpectralEntropy, noiseRemovalRatio, number_processing_threads)
  ##
  sampleFSdbIDlist <- FSA_spectra_marker_generator(sampleFSdb, ratio2basePeak4nSpectraMarkers, aggregationLevel = NA)
  ##
  ##############################################################################
  ##############################################################################
  ##
  ############################ Precursor m/z check #############################
  if (is.na(massErrorPrecursor)) {
    precursorMZcheck <- FALSE
  } else {
    precursorMZcheck <- TRUE
  }
  ############################ retention time check ############################
  if (is.na(RTtolerance)) {
    precursorRTcheck <- FALSE
  } else {
    precursorRTcheck <- TRUE
  }
  ##############################################################################
  if (minCosineSimilarity > 0) {
    cosineSimilarityCheck <- TRUE
  } else {
    cosineSimilarityCheck <- FALSE
  }
  ##
  if (minEntropySimilarity > 0) {
    entropySimilarityCheck <- TRUE
  } else {
    entropySimilarityCheck <- FALSE
  }
  ##############################################################################
  if (!is.na(targetedPrecursorType[1])) {
    refPrecursorTypeCheck <- TRUE
  } else {
    refPrecursorTypeCheck <- FALSE
  }
  ##############################################################################
  if (is.null(exportSpectraParameters)) {
    exportSpectraCheck <- FALSE
  } else {
    exportSpectraCheck <- TRUE
    maxAllowedNumberHits <- as.numeric(exportSpectraParameters[1])
    SpectraDevice <- tolower(exportSpectraParameters[2])
    file_name_msp <- exportSpectraParameters[3]
    outputMSMSspectra <- exportSpectraParameters[4]
    ##
    ############################################################################
    ## Graphical device adjustment
    colors <- c("blue", "red")
    ##
    if (SpectraDevice == "png") {
      pngDeviceCheck <- TRUE
    } else {
      pngDeviceCheck <- FALSE
    }
    ##
    if (SpectraDevice == "pdf") {
      pdfDeviceCheck <- TRUE
    } else {
      pdfDeviceCheck <- FALSE
    }
    ##
    if (SpectraDevice == "svg") {
      svgDeviceCheck <- TRUE
    } else {
      svgDeviceCheck <- FALSE
    }
    ##
    FSA_dir.create(outputMSMSspectra, allowedUnlink = TRUE)
    ##
  }
  ##############################################################################
  sampleAncillaryTable <- sampleFSdb[["MSPLibraryParameters"]]
  ##
  if (exportSpectraCheck) {
    sampleCompoundName <- sampleAncillaryTable$name
    if (is.null(sampleCompoundName)) {
      sampleCompoundName <- paste0("NA_", seq(1, nrow(sampleAncillaryTable), 1))
    }
  }
  ##
  libAncillaryTable <- libFSdb[["MSPLibraryParameters"]]
  ##
  if (exportSpectraCheck) {
    refMSPfilename <- libAncillaryTable$MSPfilename
    ##
    refCompoundName <- libAncillaryTable$name
    if (is.null(refCompoundName)) {
      refCompoundName <- paste0("NA_", seq(1, nrow(libAncillaryTable), 1))
    }
  }
  ##############################################################################
  ##################### Aggregated Library Spectra Markers #####################
  libSpectraMarkersIndex <- libFSdbIDlist[[1]]
  libNspectraMarkers <- libFSdbIDlist[[2]]
  ##############################################################################
  ################## Matching against the aggregated library ###################
  sampleIndexListSpectraMarkers <- sampleFSdbIDlist[[1]]
  sampleNspectraMarkers <- sampleFSdbIDlist[[2]]
  ####
  sampleSpectraMarkers <- as.character(round(sampleIndexListSpectraMarkers[, 1], digits = roundingDigitPrefiltering))
  if (length(sampleSpectraMarkers) > 0) {
    libIDspectraMarkersList <- libSpectraMarkersIndex[sampleSpectraMarkers]
    ####
    x_Index <- c(0, which(abs(diff(sampleIndexListSpectraMarkers[, 2])) > 0), dim(sampleIndexListSpectraMarkers)[1])
    L_x_Index <- length(x_Index) - 1
    ############################################################################
    ########################## Pre-filtering function ##########################
    ############################################################################
    call_listFSApeakIDlibID <- function(i) {
      ##
      FSApeakID <- sampleIndexListSpectraMarkers[(x_Index[i] + 1), 2]
      ##
      FSdbLibID <- do.call(c, lapply((x_Index[i] + 1):x_Index[i + 1], function(j) {
        libIDspectraMarkersList[[sampleSpectraMarkers[j]]]
      }))
      if (!is.null(FSdbLibID)) {
        ######################## Match precursor types #########################
        if (refPrecursorTypeCheck) {
          indexRefPrecursorType <- do.call(c, lapply(targetedPrecursorType, function(j) {
            which(libFSdb[["Precursor Type"]][FSdbLibID] == j)
          }))
          ##
          if (length(indexRefPrecursorType) > 0) {
            FSdbLibID <- FSdbLibID[indexRefPrecursorType]
          } else {
            FSdbLibID <- NULL
          } 
        }
        if (!is.null(FSdbLibID)) {
          ################ Checking number of spectra markers ##################
          FSApeakIDnSpectraMarkers <- sampleNspectraMarkers[FSApeakID]
          ##
          tFSdbLibID <- table(FSdbLibID)
          xT <- which((tFSdbLibID <= ceiling(FSApeakIDnSpectraMarkers/minRatioMatchedNspectraMarkers)) & (tFSdbLibID >= floor(FSApeakIDnSpectraMarkers*minRatioMatchedNspectraMarkers)))
          if (length(xT) > 0) {
            FSdbLibID <- as.numeric(names(tFSdbLibID[xT]))
            #################### Checking spectral entropy #####################
            FSApeakID_SE <- sampleFSdb[["Spectral Entropy"]][FSApeakID]
            FSdbLibID_SE <- libFSdb[["Spectral Entropy"]][FSdbLibID]
            ##
            xSE <- which(abs(FSApeakID_SE - FSdbLibID_SE) <= spectralEntropyDeviationPrefiltering)
            if (length(xSE) > 0) {
              FSdbLibID <- FSdbLibID[xSE]
              #################### Checking precursor m/z ######################
              if (precursorMZcheck) {
                FSApeakIDprecursorMZ <- sampleFSdb[["PrecursorMZ"]][FSApeakID]
                FSdbLibIDprecursorMZ <- libFSdb[["PrecursorMZ"]][FSdbLibID]
                ##
                if (allowedNominalMass) {
                  xPrecursor <- which(FSdbLibIDprecursorMZ == FSApeakIDprecursorMZ)
                } else {
                  xPrecursor <- which(abs(FSdbLibIDprecursorMZ - FSApeakIDprecursorMZ) <= massErrorPrecursor)
                }
                ##
                if (length(xPrecursor) > 0) {
                  FSdbLibID <- FSdbLibID[xPrecursor]
                } else {
                  FSdbLibID <- NULL
                }
              }
              #################### Checking retention time #####################
              if (precursorRTcheck) {
                if (!is.null(FSdbLibID)) {
                  FSApeakIDrt <- sampleFSdb[["Retention Time"]][FSApeakID]
                  FSdbLibIDrt <- libFSdb[["Retention Time"]][FSdbLibID]
                  ##
                  xRT <- which(abs(FSdbLibIDrt - FSApeakIDrt) <= RTtolerance)
                  if (length(xRT) > 0) {
                    FSdbLibID <- FSdbLibID[xRT]
                  } else {
                    FSdbLibID <- NULL
                  }
                }
              }
              ##################################################################
            } else {
              FSdbLibID <- NULL
            }
          } else {
            FSdbLibID <- NULL
          }
        } else {
          FSdbLibID <- NULL
        }
      } else {
        FSdbLibID <- NULL
      }
      ##
      list(FSApeakID, FSdbLibID)
    }
    ############################################################################
    call_fragment_matcher <- function(i) {
      matchedFSAlist <- NULL
      ####
      FSdbLibID <- listFSApeakIDlibID[[i]][[2]]
      ####
      if (!is.null(FSdbLibID)) {
        FSApeakID <- listFSApeakIDlibID[[i]][[1]]
        ####
        sampleFSdb5i <- sampleFSdb[["FragmentList"]][[FSApeakID]]
        sampleNumPeaks <- sampleFSdb[["Num Peaks"]][FSApeakID]
        sampleSpectralEntropy <- sampleFSdb[["Spectral Entropy"]][FSApeakID]
        ####
        matchedFSAlist <- do.call(rbind, lapply(FSdbLibID, function(IDj) {
          ##
          libNumPeaks <- libFSdb[["Num Peaks"]][IDj]
          if (libNumPeaks >= minMatchedNumPeaks) {
            ##
            passNspectraMarkers <- TRUE
            ##
            sampleFragmentList <- sampleFSdb5i
            ##
            libFragmentList <- libFSdb[["FragmentList"]][[IDj]]
            ##
            nSpectraMarkersIDj <- libNspectraMarkers[IDj]
            ##
            matchedSampleFragmentList <- matrix(rep(0, libNumPeaks*2), ncol = 2)
            indexMatchedLib <- rep(0, libNumPeaks)
            matchedNumPeaks <- 0
            f <- 1
            passNspectraMarkersWhile <- TRUE
            while (f <= libNumPeaks) {
              ##
              if (allowedNominalMass) {
                x_f <- which(sampleFragmentList[, 1] == libFragmentList[f, 1])
              } else {
                x_f <- which(abs(sampleFragmentList[, 1] - libFragmentList[f, 1]) <= massError)
              }
              ##
              L_x_f <- length(x_f)
              if (L_x_f > 0) {
                matchedNumPeaks <- matchedNumPeaks + 1
                ##
                if (L_x_f > 1) {
                  x_f_min <- which.min(abs(sampleFragmentList[x_f, 1] - libFragmentList[f, 1]))
                  x_f <- x_f[x_f_min[1]]
                }
                ##
                indexMatchedLib[matchedNumPeaks] <- f
                matchedSampleFragmentList[f, ] <- sampleFragmentList[x_f, ]
                sampleFragmentList[x_f, 1] <- 0
              }
              ##
              if (passNspectraMarkersWhile) {
                if (f == nSpectraMarkersIDj) {
                  ratioMatchedSpectraMarkers <- matchedNumPeaks/nSpectraMarkersIDj
                  if (ratioMatchedSpectraMarkers < minRatioMatchedNspectraMarkers) {
                    f <- libNumPeaks
                    passNspectraMarkers <- FALSE
                  } else {
                    numberMatchedSpectraMarkers <- matchedNumPeaks
                    passNspectraMarkersWhile <-  FALSE
                  }
                }
              }
              ##
              f <- f + 1
            }
            ##
            if (passNspectraMarkers) {
              if (matchedNumPeaks >= minMatchedNumPeaks) {
                ###################### Entropy Similarity ######################
                libSpectralEntropy <- libFSdb[["Spectral Entropy"]][IDj]
                ##
                entropySimilarity <- spectral_entropy_similarity_score(sampleFSdb5i, sampleSpectralEntropy, libFragmentList, libSpectralEntropy, massError, allowedNominalMass)
                if (entropySimilarity >= minEntropySimilarity) {
                  ##################### Cosine Similarity ######################
                  cosineSimilarity <- sum(matchedSampleFragmentList[, 2]*libFragmentList[, 2])/sqrt(sum(matchedSampleFragmentList[, 2]^2)*sum(libFragmentList[, 2]^2))
                  if (cosineSimilarity >= minCosineSimilarity) {
                    #################### Matched Ion Range #####################
                    indexMatchedLib <- indexMatchedLib[seq(1, matchedNumPeaks, 1)]
                    matchedSampleFragmentList <- matchedSampleFragmentList[indexMatchedLib, ]
                    if (minMatchedNumPeaks == 1) {
                      matchedSampleFragmentList <- matrix(matchedSampleFragmentList, nrow = 1)
                    }
                    ionRangeDifference <- max(matchedSampleFragmentList[, 1]) - min(matchedSampleFragmentList[, 1])
                    if (ionRangeDifference >= minIonRangeDifference) {
                      ############ Normalized Euclidean Mass Error #############
                      if (allowedNominalMass) {
                        NEME <- NA
                      } else {
                        matchedLibFragmentList <- libFragmentList[indexMatchedLib, ]
                        if (minMatchedNumPeaks == 1) {
                          matchedLibFragmentList <- matrix(matchedLibFragmentList, nrow = 1)
                        }
                        NEME <- sqrt(sum((matchedSampleFragmentList[, 1] - matchedLibFragmentList[, 1])^2)/matchedNumPeaks)*1000 # in mDa
                      }
                      if (NEME <= maxNEME) {
                        ################ Calculate sorting score ###############
                        rankingScore <- cosineSimilarity*(entropySimilarity^2)
                        ##
                        c(FSApeakID, matchedNumPeaks, numberMatchedSpectraMarkers, ratioMatchedSpectraMarkers, sampleSpectralEntropy, NEME, cosineSimilarity, entropySimilarity, libSpectralEntropy, rankingScore, IDj)
                      }
                    }
                  }
                }
              }
            }
          }
        }))
        ##
        ########################################################################
        ##
        if (!is.null(matchedFSAlist)) {
          ##
          numberMatchedFSAspectra <- nrow(matchedFSAlist)
          if (numberMatchedFSAspectra > 1) {
            ##
            if (entropySimilarityCheck == cosineSimilarityCheck) {
              matchedFSAlist <- matchedFSAlist[order(matchedFSAlist[, 10], decreasing = TRUE), ] # To sort by ranking scores
            } else if (entropySimilarityCheck) {
              matchedFSAlist <- matchedFSAlist[order(matchedFSAlist[, 8], decreasing = TRUE), ]
            } else if (cosineSimilarityCheck) {
              matchedFSAlist <- matchedFSAlist[order(matchedFSAlist[, 7], decreasing = TRUE), ]
            }
            ##
            matchedFSAlist <- cbind(matchedFSAlist, seq(1, numberMatchedFSAspectra, 1))
          } else {
            matchedFSAlist <- matrix(c(matchedFSAlist, 1), nrow = 1)
          }
          ##
          ######################################################################
          ## To sort spectra plots
          if (exportSpectraCheck) {
            ##
            NameSampleCompoundFSApeakID <- gsub('/|[\\]|\t|\n|:|[*]|[?]|<|>|[|]|[.][.][.]|[.][.]|^[.]', '_', sampleCompoundName[FSApeakID], fixed = FALSE) # To replace invalid characters
            folderSampleCompoundFSApeakID <- paste0(outputMSMSspectra, "/", FSApeakID, "_", NameSampleCompoundFSApeakID)
            ##
            dirCreationRCheck <- FSA_dir.create(folderSampleCompoundFSApeakID, allowedUnlink = TRUE)
            if (dirCreationRCheck) {
              ##
              ##################################################################
              ##
              strSampleSpectralEntropy <- paste0("S = ", sampleSpectralEntropy)
              ##
              maxNumberHits <- min(c(maxAllowedNumberHits, numberMatchedFSAspectra))
              ##
              for (counterRank in 1:maxNumberHits) {
                ##
                matchedNumPeaks <- matchedFSAlist[counterRank, 2]
                cosineSimilarity <- matchedFSAlist[counterRank, 7]
                entropySimilarity <- matchedFSAlist[counterRank, 8]
                libSpectralEntropy <- matchedFSAlist[counterRank, 9]
                IDj <- matchedFSAlist[counterRank, 11]
                ##
                filenameMSMSmatch <- paste0(folderSampleCompoundFSApeakID, "/", counterRank, "_ID_", IDj, ".", SpectraDevice)
                fileCreateRCheck <- file.create(file = filenameMSMSmatch, showWarnings = FALSE)
                if (fileCreateRCheck) {
                  ##
                  if (allowedNominalMass) {
                    samplePlotText <- paste0("Entropy Similarity = ", round(entropySimilarity, 3), " | Cosine Similarity = ", round(cosineSimilarity, 3), " | Num Matched Peaks = ", matchedNumPeaks)
                  } else {
                    NEME <- matchedFSAlist[counterRank, 6]
                    samplePlotText <- paste0("Entropy Similarity = ", round(entropySimilarity, 3), " | Cosine Similarity = ", round(cosineSimilarity, 3), " | NEME = ", round(NEME, 2), " mDa | Num Matched Peaks = ", matchedNumPeaks)
                  }
                  ##
                  if (precursorRTcheck) {
                    rtDevSamplePlotText <- round((sampleFSdb[["Retention Time"]][FSApeakID] - libFSdb[["Retention Time"]][IDj]), digits = 3)
                    signRtDevSamplePlotText <- if (rtDevSamplePlotText > 0) {"+"} else {""}
                    samplePlotText <- paste0(samplePlotText, " | RT dev = ", signRtDevSamplePlotText, rtDevSamplePlotText, " min")
                  }
                  ##
                  ##############################################################
                  ##
                  libPlotText <- ""
                  ##
                  if (!is.infinite(libFSdb[["PrecursorMZ"]][IDj])) {
                    libPlotText <- paste0("PrecursorMZ = ", libFSdb[["PrecursorMZ"]][IDj])
                  }
                  ##
                  if (libFSdb[["Precursor Type"]][IDj] != "") {
                    refPrecursorTypeIDj <- libFSdb[["MSPLibraryParameters"]][["precursor_type"]][IDj]
                    if (refPrecursorTypeIDj != "") {
                      if (libPlotText == "") {
                        libPlotText <- paste0("PrecursorType = ", refPrecursorTypeIDj)
                      } else {
                        libPlotText <- paste0(libPlotText, " | PrecursorType = ", refPrecursorTypeIDj)
                      }
                    }
                  }
                  ##
                  if (precursorRTcheck) {
                    if (!is.infinite(libFSdb[["Retention Time"]][IDj])) {
                      if (libPlotText == "") {
                        libPlotText <- paste0("RT = ", libFSdb[["Retention Time"]][IDj], " min")
                      } else {
                        libPlotText <- paste0(libPlotText, " | RT = ", libFSdb[["Retention Time"]][IDj], " min")
                      }
                    }
                  }
                  ##
                  ##############################################################
                  ##
                  sampleFragmentList <- sampleFSdb5i
                  sampleFragmentList[, 2] <- sampleFragmentList[, 2]/max(sampleFragmentList[, 2])*100
                  libFragmentList <- libFSdb[["FragmentList"]][[IDj]]
                  libFragmentList[, 2] <- -libFragmentList[, 2]/max(libFragmentList[, 2])*100
                  ##
                  combinedFragmentList <- data.frame(rbind(sampleFragmentList, libFragmentList))
                  typeCombinedFragmentList <- c(rep("Sample", sampleNumPeaks), rep("Library", libFSdb[["Num Peaks"]][IDj]))
                  combinedFragmentList <- cbind(combinedFragmentList, typeCombinedFragmentList)
                  colnames(combinedFragmentList) <- c("mz", "int", "type")
                  combinedFragmentList$type <- factor(combinedFragmentList$type, levels = c("Sample", "Library"))
                  ##
                  xMinLimPlot <- min(combinedFragmentList$mz)
                  xMaxLimPlot <- max(combinedFragmentList$mz)
                  x_position <- (xMinLimPlot + xMaxLimPlot)/2
                  ##
                  yPositionLabel <- c((sampleFragmentList[, 2] + 5), (libFragmentList[, 2] - 5))
                  ##
                  plotlabelsMZ <- FSA_annotation_text_repel(FSAspectra = combinedFragmentList[, 1:2], nGridX = 12, nGridY = 6)
                  ## Spectra Markers Signs
                  nSpectraMarkersIDj <- libNspectraMarkers[IDj]
                  xSpectraMarkerFSApeakID <- c(sampleFragmentList[seq(1, sampleNspectraMarkers[FSApeakID], 1), 1], libFragmentList[seq(1, nSpectraMarkersIDj, 1), 1])
                  ySpectraMarkerFSApeakID <- c((sampleFragmentList[seq(1, sampleNspectraMarkers[FSApeakID], 1), 2] + 10.5), (libFragmentList[seq(1, nSpectraMarkersIDj, 1), 2] - 10.5))
                  ##
                  ##############################################################
                  ##
                  if (pngDeviceCheck) {
                    png(filenameMSMSmatch, width = 16, height = 8, units = "in", res = 100)
                  }
                  ##
                  if (pdfDeviceCheck) {
                    pdf(filenameMSMSmatch, width = 16, height = 8, colormodel = "cmyk")
                  }
                  ##
                  if (svgDeviceCheck) {
                    svg(filenameMSMSmatch, width = 16, height = 8)
                  }
                  ##
                  ##############################################################
                  ##
                  plot(combinedFragmentList$mz, combinedFragmentList$int, type = "h", lend = 2, xlim = c((xMinLimPlot - 10), (xMaxLimPlot + 10)), ylim = c(-125, 125), 
                       lwd = 2, lend = 2, col = colors[combinedFragmentList$type], cex = 4, xaxt = "n", xlab = "", ylab = "", yaxt = "n")
                  ##
                  points(xSpectraMarkerFSApeakID, ySpectraMarkerFSApeakID, type = "p", col = "purple", cex = 1.5, pch = 8)
                  ##
                  abline(h = 0, col = "black")
                  ##
                  text(x = combinedFragmentList$mz, y = yPositionLabel, label = plotlabelsMZ)
                  text(x = x_position , y = -125, cex = 1.125, label = libPlotText)
                  text(x = x_position , y = 125, cex = 1.125, label = samplePlotText)
                  ## y-axis
                  mtext("Adjusted intensity (%)", side = 2, adj = 0.5, line = 0.35, cex = 1.35)
                  ## sample info
                  mtext(strSampleSpectralEntropy, side = 3, adj = 1, line = 0.225, cex = 1.0, col = colors[1])
                  mtext(file_name_msp, side = 3, adj = 0, line = 0.25, cex = 1.15, col = colors[1])
                  mtext(sampleCompoundName[FSApeakID], side = 3, adj = 1, line = 1.5, cex = 1.10, col = colors[1])
                  ## sample info
                  mtext(paste0("S = ", libSpectralEntropy), side = 1, adj = 1, line = 0.225, cex = 1.0, col = colors[2])
                  mtext(refMSPfilename[IDj], side = 1, adj = 0, line = 0.25, cex = 1.15, col = colors[2])
                  mtext(refCompoundName[IDj], side = 1, adj = 1, line = 1.5, cex = 1.10, col = colors[2])
                  ## legend
                  if (precursorMZcheck) {
                    points(libFSdb[["PrecursorMZ"]][IDj], -6, type = "p", col = "green", cex = 2, pch = 17)
                    ##
                    legend(x = "topleft", legend = c("Sample", "Library", "Marker", "Precursor"), 
                           lty = c(1, 1, NA, NA),
                           col = c(colors, "purple", "green"),
                           lwd = c(2, 2, NA, NA),
                           pch = c(NA, NA, 8, 17),
                           pt.cex = c(NA, NA, 1.5, 2),
                           cex = 1.125, seg.len = 1, x.intersp = 0.5, y.intersp = 0.9)
                  } else {
                    legend(x = "topleft", legend = c("Sample", "Library", "Marker"),
                           lty = c(1, 1, NA),
                           col = c(colors, "purple"),
                           lwd = c(2, 2, NA),
                           pch = c(NA, NA, 8),
                           pt.cex = c(NA, NA, 1.5),
                           cex = 1.125, seg.len = 1, x.intersp = 0.5, y.intersp = 0.9)
                  }
                  ##
                  dev.off()
                } else {
                  FSA_logRecorder(paste0("WARNING!!! Figure can not be created for `", NameSampleCompoundFSApeakID, "` due to character length limit for `", MSPfile, "`!"))
                }
              }
            } else {
              FSA_logRecorder(paste0("WARNING!!! Directory can not be created for `", folderSampleCompoundFSApeakID, "` due to character length limit for `", MSPfile, "`!"))
            }
          }
        }
      }
      ##
      return(matchedFSAlist)
    }
    ############################################################################
    if (number_processing_threads == 1) {
      ####
      listFSApeakIDlibID <- lapply(1:L_x_Index, function(i) {
        call_listFSApeakIDlibID(i)
      })
      ####
      libSpectraMarkersIndex <- NULL
      sampleSpectraMarkers <- NULL
      sampleIndexListSpectraMarkers <- NULL
      libIDspectraMarkersList <- NULL
      ####
      FSAannotationList <- do.call(rbind, lapply(1:L_x_Index, function(i) {
        call_fragment_matcher(i)
      }))
      ##
    } else {
      osType <- Sys.info()[['sysname']]
      ##
      if (osType == "Linux") {
        ####
        listFSApeakIDlibID <- mclapply(1:L_x_Index, function(i) {
          call_listFSApeakIDlibID(i)
        }, mc.cores = number_processing_threads)
        ####
        libSpectraMarkersIndex <- NULL
        sampleSpectraMarkers <- NULL
        sampleIndexListSpectraMarkers <- NULL
        libIDspectraMarkersList <- NULL
        ####
        FSAannotationList <- do.call(rbind, mclapply(1:L_x_Index, function(i) {
          call_fragment_matcher(i)
        }, mc.cores = number_processing_threads))
        ####
        closeAllConnections()
        ##
      } else if (osType == "Windows") {
        ##
        clust <- makeCluster(number_processing_threads)
        registerDoParallel(clust)
        ####
        listFSApeakIDlibID <- foreach(i = 1:L_x_Index, .verbose = FALSE) %dopar% {
          call_listFSApeakIDlibID(i)
        }
        ####
        libSpectraMarkersIndex <- NULL
        sampleSpectraMarkers <- NULL
        sampleIndexListSpectraMarkers <- NULL
        libIDspectraMarkersList <- NULL
        ####
        FSAannotationList <- foreach(i = 1:L_x_Index, .combine = 'rbind', .verbose = FALSE) %dopar% {
          call_fragment_matcher(i)
        }
        ####
        stopCluster(clust)
      }
    }
    ############################################################################
    if (!is.null(FSAannotationList)) {
      FSAannotationList <- matrix(FSAannotationList, ncol = 12)
      FSAannotationList <- matrix(FSAannotationList[, -10], ncol = 11)
      FSAannotationList[, 4] <- round(FSAannotationList[, 4]*100, 2)
      FSAannotationList[, 5] <- round(FSAannotationList[, 5], 5)
      FSAannotationList[, 6] <- round(FSAannotationList[, 6], 2)
      FSAannotationList[, 7] <- round(FSAannotationList[, 7], 3)
      FSAannotationList[, 8] <- round(FSAannotationList[, 8], 3)
      FSAannotationList[, 9] <- round(FSAannotationList[, 9], 5)
      ##########################################################################
      analyteAncillary <- data.frame(sampleAncillaryTable[FSAannotationList[, 1], ])
      colnames(analyteAncillary) <- paste0("analyte_", colnames(analyteAncillary))
      ##########################################################################
      matched_sample <- data.frame(cbind(matrix(FSAannotationList[, 2:9], ncol = 8), matrix(FSAannotationList[, 11], ncol = 1), matrix(FSAannotationList[, 10], ncol = 1)))
      colnames(matched_sample) <- c("matchedTotalNumPeaks", "matchedNumberSpectraMarkers", "matchedPercentageSpectraMarkers", "sampleSpectralEntropy", "matchedNEME_mDa", "matchedCosineSimilarity", "matchedEntropySimilarity", "libSpectralEntropy", "matchedRank", "IDSL.FSA_FSDBreferenceID")
      ##########################################################################
      matched_ref <- data.frame(libAncillaryTable[FSAannotationList[, 10], ])
      colnames(matched_ref) <- paste0("FSDB_", colnames(libAncillaryTable))
      ##
      ##########################################################################
      ##
      if (allowedNominalMass) {
        matched_sample$'matchedNEME_mDa' <- NULL
      }
      ##
      ##########################################################################
      ##
      mspMatch <- data.frame(cbind(analyteAncillary, matched_sample, matched_ref))
      ##
      ##########################################################################
      ##########################################################################
      ##########################################################################
      ## Special and specific option for annotating .msp files created by the IDSL.CSA package ##
      ##
      colx1Check <- grep("analyte_idsl.ipa_collective_peakids", colnames(analyteAncillary), ignore.case = TRUE)
      if (length(colx1Check) > 0) {
        ## TO skip correcting MSP of unique tag spectra
        colx2Check <- grep("analyte_csavariantdetectionfreq", colnames(analyteAncillary), ignore.case = TRUE)
        if (length(colx2Check) == 0) {
          ##
          mspMatch$analyte_idsl.ipa_peakid <- NULL
          ##
          listIPApeakIDname <- FSA_R.aggregate(mspMatch$analyte_name)
          uniqueAnalyteName <- names(listIPApeakIDname)
          ##
          call_collective_peakIDs <- function(i) {
            x_i <- listIPApeakIDname[[i]]
            subsetMSPmatch <- mspMatch[x_i, ]
            L_x_i <- length(x_i)
            ##
            collectiveIPApeakID <- eval(parse(text = paste0("c(", subsetMSPmatch$analyte_idsl.ipa_collective_peakids[1], ")")))
            CSApeakID <- which(sampleFSdb[["MSPLibraryParameters"]][["csapeakgrouping_id"]] == subsetMSPmatch$analyte_csapeakgrouping_id[1])
            if (length(CSApeakID) > 1) {
              xCSApeakID <- which(sampleFSdb[["MSPLibraryParameters"]][["idsl.ipa_collective_peakids"]][CSApeakID] == subsetMSPmatch$analyte_idsl.ipa_collective_peakids[1])
              CSApeakID <- CSApeakID[xCSApeakID[1]]
            }
            ##
            iCollectiveIPApeakID <- collectiveIPApeakID[1:sampleFSdb[["Num Peaks"]][CSApeakID]]
            ##
            xiIPApeakID <- which(iCollectiveIPApeakID != 0)
            IPApeakID <- iCollectiveIPApeakID[xiIPApeakID]
            ##
            IPA12CMZ <- sampleFSdb[["FragmentList"]][[CSApeakID]][xiIPApeakID, 1]
            IPART <- sampleFSdb[["Retention Time"]][CSApeakID]
            ##
            do.call(rbind, lapply(1:length(xiIPApeakID), function(j) {
              tempSubsetMSPmatch <- subsetMSPmatch
              tempSubsetMSPmatch$analyte_idsl.ipa_collective_peakids <- rep(IPApeakID[j], L_x_i)
              ##
              analyte_name <- paste0("IDSL.IPA_PeakID_", IPApeakID[j], "_mz_", IPA12CMZ[j], "_RT_", IPART)
              cbind(rep(analyte_name, L_x_i), rep(IPA12CMZ[j], L_x_i), rep(IPART, L_x_i), tempSubsetMSPmatch)
            }))
          }
          ##
          if (number_processing_threads == 1) {
            ##
            mspMatch <- do.call(rbind, lapply(uniqueAnalyteName, function(i) {
              call_collective_peakIDs(i)
            }))
            ##
          } else {
            ##
            if (osType == "Linux") {
              ##
              mspMatch <- do.call(rbind, mclapply(uniqueAnalyteName, function(i) {
                call_collective_peakIDs(i)
              }, mc.cores = number_processing_threads))
              ##
              closeAllConnections()
              ##
            } else if (osType == "Windows") {
              ##
              clust <- makeCluster(number_processing_threads)
              registerDoParallel(clust)
              ##
              mspMatch <- foreach(i = uniqueAnalyteName, .combine = 'rbind', .verbose = FALSE) %dopar% {
                call_collective_peakIDs(i)
              }
              ##
              stopCluster(clust)
            }
          }
          ##
          colnames(mspMatch)[colnames(mspMatch) == 'analyte_name'] <- 'analyte_csa_name'
          colnames(mspMatch)[colnames(mspMatch) == 'analyte_idsl.ipa_collective_peakids'] <- 'analyte_idsl.ipa_peakid'
          colnames(mspMatch)[1:3] <- c("analyte_name", "analyte_mz", "analyte_rt")
        }
      }
      ##########################################################################
      ##########################################################################
      ##########################################################################
      ##
      rownames(mspMatch) <- NULL
    }
  }
  ##
  return(mspMatch)
}