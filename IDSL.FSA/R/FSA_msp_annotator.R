FSA_msp_annotator <- function(PARAM_SPEC, libFSdb, address_input_msp, output_path, allowedVerbose = TRUE) {
  ##
  if (allowedVerbose) {
    ##
    ############################################################################
    ## To create log record for IDSL.FSA
    initiation_time <- Sys.time()
    timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
    FSA_opendir(output_path)
    FSA_logRecorder(paste0("MSP: ", address_input_msp))
    FSA_logRecorder(paste0("OUTPUT: ", output_path))
    FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
    FSA_logRecorder("Initiated `.msp` annotation workflow!")
    FSA_logRecorder(paste0(as.character(initiation_time), " ", timeZone))
    FSA_logRecorder("", allowedPrinting = FALSE)
    FSA_logRecorder("", allowedPrinting = FALSE)
    FSA_logRecorder(paste0(PARAM_SPEC[, 1], "\t", PARAM_SPEC[, 2]),  allowedPrinting = FALSE)
    FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
    ##
    ############################################################################
    ##
  }
  file_name_sample_msp <- dir(path = address_input_msp)
  file_name_sample_msp <- file_name_sample_msp[grepl(".msp$", file_name_sample_msp, ignore.case = TRUE)]
  if (length(file_name_sample_msp) == 0) {
    stop(FSA_logRecorder("ERROR!!! EMPTY MSP FOLDER!!!"))
  }
  ##
  output_path_spectra_tables <- paste0(output_path, "/annotated_spectra_tables")
  FSA_dir.create(output_path_spectra_tables, allowedUnlink = FALSE)
  FSA_opendir(output_path_spectra_tables)
  NPT <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0001'), 2])
  parallelizationMode <- tolower(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0002'), 2])
  ##
  x0003 <- which(PARAM_SPEC[, 1] == 'SPEC0003')
  if (length(x0003) > 0) {
    allowedNominalMass <- eval(parse(text = (PARAM_SPEC[x0003, 2])))
  } else {
    allowedNominalMass <- FALSE
  }
  ##
  targetedPrecursorType <- tryCatch(eval(parse(text = paste0("c(", PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0004'), 2], ")"))), error = function(e){NA})
  if (!is.na(targetedPrecursorType[1])) {
    targetedPrecursorType <- unique(UFSA_precursorType_corrector(precursorType = targetedPrecursorType, ionMode = NULL))
    targetedPrecursorType <- setdiff(targetedPrecursorType, "")
    if (allowedVerbose) {FSA_logRecorder(paste0("The following ", length(targetedPrecursorType), " precursor type(s) are selected for screening: ", paste0(targetedPrecursorType, collapse = ", "), " !"))}
  } else {
    if (allowedVerbose) {FSA_logRecorder("No precursor types was selected for screening!")}
  }
  ##
  massErrorPrecursor <- tryCatch(as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0005'), 2]), warning = function(w){NA})
  RTtolerance <- tryCatch(as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0006'), 2]), warning = function(w){NA})
  spectralEntropyDeviationPrefiltering <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0007'), 2])
  noiseRemovalRatio <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0008'), 2])/100
  ratio2basePeak4nSpectraMarkers <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0009'), 2])/100
  minRatioMatchedNspectraMarkers <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0010'), 2])/100
  ##
  if (allowedNominalMass) {
    if (is.na(massErrorPrecursor)) {
      if (allowedVerbose) {FSA_logRecorder("NOTICE: Precursor m/z match (SPEC0005) was not selected for spectra annotations!",
                                         allowedPrinting = FALSE)}
    } else {
      if (allowedVerbose) {FSA_logRecorder("NOTICE: Precursor m/z match (SPEC0005) was selected for spectra annotations! Empty annotation tables will be generated when .msp files do not contain 'PrecursorMZ' information!",
                                         allowedPrinting = FALSE)}
    }
    ##
    roundingDigitPrefiltering <- 0
    massError <- 0
    maxNEME <- 0
    ##
  } else {
    ##
    if (is.na(massErrorPrecursor)) {
      if (allowedVerbose) {FSA_logRecorder("NOTICE: Precursor m/z match (SPEC0005) was not selected for spectra annotations!",
                                         allowedPrinting = FALSE)}
    } else {
      if (allowedVerbose) {FSA_logRecorder(paste0("NOTICE: Mass accuracy = '" , massErrorPrecursor, " Da' for precursor m/z (SPEC0005) was selected for spectra annotations! Empty annotation tables will be generated when .msp files do not contain 'PrecursorMZ' information!"),
                                         allowedPrinting = FALSE)}
    }
    ##
    roundingDigitPrefiltering <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0011'), 2])
    massError <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0012'), 2])
    maxNEME <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0018'), 2])
  }
  ##
  if (is.na(RTtolerance)) {
    if (allowedVerbose) {FSA_logRecorder("NOTICE: Retention time tolerance (SPEC0006) was not selected and retention time values will not be used for spectra annotations!",
                                       allowedPrinting = FALSE)}
  } else {
    if (allowedVerbose) {FSA_logRecorder(paste0("NOTICE: Retention time tolerance = '" , RTtolerance, " min' (SPEC0006) was selected for spectra annotations! Empty annotation tables will be generated when .msp files do not contain 'Retention Time' information!"),
                                       allowedPrinting = FALSE)}
  }
  ##
  allowedWeightedSpectralEntropy <- eval(parse(text = (PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0013'), 2])))
  minIonRangeDifference = as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0014'), 2])
  minMatchedNumPeaks <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0015'), 2])
  minEntropySimilarity <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0016'), 2])
  minCosineSimilarity <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0017'), 2])
  maxAllowedNumberHits <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0019'), 2])
  ##
  if (maxAllowedNumberHits > 0) {
    dev.offCheck <- TRUE
    while (dev.offCheck) {
      dev.offCheck <- tryCatch(dev.off(), error = function(e) {FALSE})
    }
    ##
    exportSpectraCheck <- TRUE
    outputMSMSspectra <- paste0(output_path, "/plot")
    FSA_dir.create(outputMSMSspectra, allowedUnlink = FALSE)
    FSA_opendir(outputMSMSspectra)
    if (allowedVerbose) {FSA_logRecorder("FSA spectra comparison plots with the FSDB library are stored in the `plot` folder!")}
    ##
    SpectraDevice <- tolower(gsub(" ", "", PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0020'), 2]))
    if ((SpectraDevice != "pdf") & (SpectraDevice != "svg")) {
      SpectraDevice <- "png"
    }
    exportSpectraParameters <- c(maxAllowedNumberHits, SpectraDevice, NA, NA)
  } else {
    exportSpectraCheck <- FALSE
    exportSpectraParameters <- NULL
  }
  ############################## Details FSDB ##################################
  detailsFSdb <- libFSdb[["logFSdb"]]
  FSdbMassIntegrationWindow <- detailsFSdb$massWindowIntegration
  FSdbAllowedMergeNominalMass <- detailsFSdb$allowedNominalMass
  FSdbWeightedSpectralEntropyMode <- detailsFSdb$allowedWeightedSpectralEntropy
  FSdbNoiseRemovalRatio <- detailsFSdb$noiseRemovalPercentage/100
  ##
  if (allowedVerbose) {
    if ((FSdbMassIntegrationWindow != massError) | (FSdbAllowedMergeNominalMass != allowedNominalMass) | (FSdbWeightedSpectralEntropyMode != allowedWeightedSpectralEntropy) | (FSdbNoiseRemovalRatio != noiseRemovalRatio)) {
      FSA_logRecorder("WARNING!!! Parameters used to generate FSDB shown below are different from the search parameters!!!")
      FSA_logRecorder(detailsFSdb)
    }
  }
  ##############################################################################
  ########################## Spectra marker generator ##########################
  ##
  if (allowedVerbose) {FSA_logRecorder("Initiated aggregating spectra markers from the FSDB for `.msp` prefiltering!")}
  ##
  libFSdbIDlist <- FSA_spectra_marker_generator(libFSdb, ratio2basePeak4nSpectraMarkers, aggregationLevel = roundingDigitPrefiltering)
  ##
  if (allowedVerbose) {FSA_logRecorder("Completed aggregating spectra markers!")}
  ##############################################################################
  ##
  if (allowedVerbose) {FSA_logRecorder("Initiated spectra annotation on individual `.msp` files!")}
  if (allowedVerbose) {FSA_logRecorder("Individual annotated spectra tables are stored in `.Rdata` and `.csv` formats in the `annotated_spectra_tables` folder!")}
  ##
  ##############################################################################
  ##############################################################################
  ##
  call_msp_annotator <- function(address_input_msp, iFileNameMSP, exportSpectraCheck, exportSpectraParameters, outputMSMSspectra,
                                 output_path_spectra_tables, libFSdb, libFSdbIDlist, targetedPrecursorType, ratio2basePeak4nSpectraMarkers,
                                 allowedNominalMass, allowedWeightedSpectralEntropy, noiseRemovalRatio, roundingDigitPrefiltering,
                                 minMatchedNumPeaks, massError, maxNEME, minIonRangeDifference, minCosineSimilarity, minEntropySimilarity,
                                 minRatioMatchedNspectraMarkers, spectralEntropyDeviationPrefiltering, massErrorPrecursor, RTtolerance, NPT) {
    ##
    if (exportSpectraCheck) {
      exportSpectraParameters[3] <- iFileNameMSP
      exportSpectraParameters[4] <- paste0(outputMSMSspectra, "/", iFileNameMSP)
    }
    ##
    SpectraAnnotationTable <- fragmentation_spectra_annotator(path = address_input_msp, MSPfile = iFileNameMSP, libFSdb, libFSdbIDlist, targetedPrecursorType,
                                                              ratio2basePeak4nSpectraMarkers, allowedNominalMass, allowedWeightedSpectralEntropy,
                                                              noiseRemovalRatio, roundingDigitPrefiltering, minMatchedNumPeaks, massError, maxNEME,
                                                              minIonRangeDifference, minCosineSimilarity, minEntropySimilarity, minRatioMatchedNspectraMarkers,
                                                              spectralEntropyDeviationPrefiltering, massErrorPrecursor, RTtolerance, exportSpectraParameters,
                                                              number_processing_threads = NPT)
    ##
    if (!is.null(SpectraAnnotationTable)) {
      save(SpectraAnnotationTable, file = paste0(output_path_spectra_tables, "/SpectraAnnotationTable_", iFileNameMSP, ".Rdata"))
      write.csv(SpectraAnnotationTable, file = paste0(output_path_spectra_tables, "/SpectraAnnotationTable_", iFileNameMSP, ".csv"), row.names = TRUE)
    } else {
      FSA_logRecorder(paste0("No annotation was made for `", iFileNameMSP, "`!"))
    }
    ##
    return()
  }
  ##############################################################################
  if (NPT == 1 | parallelizationMode == "peakmode") {
    L_MSP <- length(file_name_sample_msp)
    ##
    if (allowedVerbose) {
      iCounter <- 0
      progressBARboundaries <- txtProgressBar(min = 0, max = L_MSP, initial = 0, style = 3)
    }
    ##
    for (i in file_name_sample_msp) {
      ##
      null_variable <- tryCatch(call_msp_annotator(address_input_msp, iFileNameMSP = i, exportSpectraCheck, exportSpectraParameters, outputMSMSspectra,
                                                   output_path_spectra_tables, libFSdb, libFSdbIDlist, targetedPrecursorType, ratio2basePeak4nSpectraMarkers,
                                                   allowedNominalMass, allowedWeightedSpectralEntropy, noiseRemovalRatio, roundingDigitPrefiltering,
                                                   minMatchedNumPeaks, massError, maxNEME, minIonRangeDifference, minCosineSimilarity, minEntropySimilarity,
                                                   minRatioMatchedNspectraMarkers, spectralEntropyDeviationPrefiltering, massErrorPrecursor, RTtolerance, NPT),
                                error = function(e) {FSA_logRecorder(paste0("Problem with `", i,"`!"))})
      ##
      if (allowedVerbose) {
        iCounter <- iCounter + 1
        setTxtProgressBar(progressBARboundaries, iCounter)
      }
    }
    if (allowedVerbose) {close(progressBARboundaries)}
    ##
  } else if (parallelizationMode == "samplemode") {
    NPT0 <- NPT
    NPT <- 1
    ##
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Linux") {
      ##
      null_variable <- mclapply(file_name_sample_msp, function(i) {
        ##
        tryCatch(call_msp_annotator(address_input_msp, iFileNameMSP = i, exportSpectraCheck, exportSpectraParameters, outputMSMSspectra,
                                    output_path_spectra_tables, libFSdb, libFSdbIDlist, targetedPrecursorType, ratio2basePeak4nSpectraMarkers,
                                    allowedNominalMass, allowedWeightedSpectralEntropy, noiseRemovalRatio, roundingDigitPrefiltering,
                                    minMatchedNumPeaks, massError, maxNEME, minIonRangeDifference, minCosineSimilarity, minEntropySimilarity,
                                    minRatioMatchedNspectraMarkers, spectralEntropyDeviationPrefiltering, massErrorPrecursor, RTtolerance, NPT),
                 error = function(e) {FSA_logRecorder(paste0("Problem with `", i,"`!"))})
      }, mc.cores = NPT0)
      ##
      closeAllConnections()
      ##
    } else if (osType == "Windows") {
      ##
      clust <- makeCluster(NPT0)
      registerDoParallel(clust)
      ##
      null_variable <- foreach(i = file_name_sample_msp, .verbose = FALSE) %dopar% {
        ##
        tryCatch(call_msp_annotator(address_input_msp, iFileNameMSP = i, exportSpectraCheck, exportSpectraParameters, outputMSMSspectra,
                                    output_path_spectra_tables, libFSdb, libFSdbIDlist, targetedPrecursorType, ratio2basePeak4nSpectraMarkers,
                                    allowedNominalMass, allowedWeightedSpectralEntropy, noiseRemovalRatio, roundingDigitPrefiltering,
                                    minMatchedNumPeaks, massError, maxNEME, minIonRangeDifference, minCosineSimilarity, minEntropySimilarity,
                                    minRatioMatchedNspectraMarkers, spectralEntropyDeviationPrefiltering, massErrorPrecursor, RTtolerance, NPT),
                 error = function(e) {FSA_logRecorder(paste0("Problem with `", i,"`!"))})
      }
      ##
      stopCluster(clust)
    }
  }
  if (allowedVerbose) {
    ##
    ############################################################################
    ##
    completion_time <- Sys.time()
    FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
    required_time <- completion_time - initiation_time
    FSA_logRecorder(paste0("The required processing time was `", required_time, " ", attributes(required_time)$units, "`"))
    FSA_logRecorder(paste0(as.character(completion_time), " ", timeZone), allowedPrinting = FALSE)
    FSA_logRecorder("", allowedPrinting = FALSE)
    FSA_logRecorder("", allowedPrinting = FALSE)
    FSA_logRecorder("Successfully completed spectra annotation on individual `.msp` files!")
    FSA_logRecorder(paste0(rep("", 100), collapse = "="), allowedPrinting = FALSE)
    ##
    ############################################################################
    ##
  }
  ##
  gc()
  closeAllConnections()
  ##
  return()
}