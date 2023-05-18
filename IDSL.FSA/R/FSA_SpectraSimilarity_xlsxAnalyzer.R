FSA_SpectraSimilarity_xlsxAnalyzer <- function(spreadsheet) {
  ##
  checkpoint_parameter <- FALSE
  #
  if (typeof(spreadsheet) == "list") {
    if (ncol(spreadsheet) >= 4) {
      PARAM_SPEC <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
      ##
    } else if (ncol(spreadsheet) == 2) {
      PARAM_SPEC <- spreadsheet
      checkpoint_parameter <- TRUE
      ##
    } else {
      FSA_message("The `SpectraSimilarity` spreadsheet tab was not produced properly!")
    }
  } else if (typeof(spreadsheet) == "character") {
    if (length(spreadsheet) == 1) {
      if (file.exists(spreadsheet)) {
        ##
        readxlPackageCheck <- tryCatch(requireNamespace('readxl', quietly = TRUE), error = function(e) {FALSE})
        if (!readxlPackageCheck) {
          warning("IDSL.FSA requires the 'readxl' package of R to read Excel spreadsheets!")
          stop(" <<< install.packages('readxl') >>> ")
        }
        ##
        PARAM_SPEC <- readxl::read_xlsx(spreadsheet, sheet = "SpectraSimilarity")
        PARAM_SPEC <- cbind(PARAM_SPEC[, 2], PARAM_SPEC[, 4])
        checkpoint_parameter <- TRUE
      } else {
        FSA_message("The `SpectraSimilarity` spreadsheet tab not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      FSA_message("The `SpectraSimilarity` spreadsheet tab was not produced properly!")
    }
  } else {
    FSA_message("The `SpectraSimilarity` spreadsheet tab was not produced properly!")
  }
  ##
  if (checkpoint_parameter) {
    ############################################################################
    number_processing_threads <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0001'), 2])
    if (length(number_processing_threads) == 0) {
      FSA_message("ERROR!!! Problem with SPEC0001! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (number_processing_threads >= 1) {
        if ((number_processing_threads %% 1) != 0) {
          FSA_message("ERROR!!! Problem with SPEC0001! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        FSA_message("ERROR!!! Problem with SPEC0001! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (number_processing_threads > 1) {
      x0002 <- which(PARAM_SPEC[, 1] == 'SPEC0002')
      SPEC0002 <- PARAM_SPEC[x0002, 2]
      if (is.na(SPEC0002)) {
        FSA_message("ERROR!!! Problem with SPEC0002!")
        checkpoint_parameter <- FALSE
      } else {
        SPEC0002 <- gsub(" ", "", tolower(SPEC0002))
        if (SPEC0002 == "samplemode" | SPEC0002 == "peakmode") {
          PARAM_SPEC[x0002, 2] <- SPEC0002
        } else {
          FSA_message("ERROR!!! Problem with SPEC0002!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x0003 <- which(PARAM_SPEC[, 1] == 'SPEC0003')
    if (length(x0003) > 0) {
      allowedNominalMass <- tolower(gsub(" ", "", PARAM_SPEC[x0003, 2]))
      if (allowedNominalMass == "1" | allowedNominalMass == "t" | allowedNominalMass == "true") {
        allowedNominalMass <- TRUE
        FSA_message("NOTICE: Nominal masses will be utilized for spectra matching!")
      } else {
        allowedNominalMass <- FALSE
      }
    } else {
      allowedNominalMass <- FALSE
    }
    PARAM_SPEC[x0003, 2] <- allowedNominalMass
    ##
    targetedPrecursorType <- tryCatch(eval(parse(text = paste0("c(", PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0004'), 2], ")"))), error = function(e){NULL})
    if (is.null(targetedPrecursorType)) {
      FSA_message("ERROR!!! Problem with SPEC0004!")
      checkpoint_parameter <- FALSE
    }
    ##
    x0005 <- which(PARAM_SPEC[, 1] == 'SPEC0005')
    massErrorPrecursor <- tryCatch(as.numeric(PARAM_SPEC[x0005, 2]), warning = function(w){NA})
    PARAM_SPEC[x0005, 2] <- massErrorPrecursor
    ##
    x0006 <- which(PARAM_SPEC[, 1] == 'SPEC0006')
    RTtolerance <- tryCatch(as.numeric(PARAM_SPEC[x0006, 2]), warning = function(w){NA})
    PARAM_SPEC[x0006, 2] <- RTtolerance
    ##
    spectralEntropyDeviationPrefiltering <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0007'), 2])
    if (length(spectralEntropyDeviationPrefiltering) == 0) {
      FSA_message("ERROR!!! Problem with SPEC0007! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (spectralEntropyDeviationPrefiltering < 0) {
        FSA_message("ERROR!!! Problem with SPEC0007! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    noiseRemovalRatio <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0008'), 2])
    if (length(noiseRemovalRatio) == 0) {
      FSA_message("ERROR!!! Problem with SPEC0008! This parameter should be a positive number between 0 - 100!")
      checkpoint_parameter <- FALSE
    } else {
      if (!((noiseRemovalRatio >= 0) & (noiseRemovalRatio <= 100))) {
        FSA_message("ERROR!!! Problem with SPEC0008! This parameter should be a positive number between 0 - 100!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    percentage2basePeak4nSpectraMarkers <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0009'), 2])
    if (length(percentage2basePeak4nSpectraMarkers) == 0) {
      FSA_message("ERROR!!! Problem with SPEC0009! This parameter should be a positive number between 0 - 100!")
      checkpoint_parameter <- FALSE
    } else {
      if (!((percentage2basePeak4nSpectraMarkers >= 0) & (percentage2basePeak4nSpectraMarkers <= 100))) {
        FSA_message("ERROR!!! Problem with SPEC0009! This parameter should be a positive number between 0 - 100!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minPercentageMatchedNspectraMarkers <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0010'), 2])
    if (length(minPercentageMatchedNspectraMarkers) == 0) {
      FSA_message("ERROR!!! Problem with SPEC0010! This parameter should be a positive number greater than 0 and less than 100!")
      checkpoint_parameter <- FALSE
    } else {
      if (!((minPercentageMatchedNspectraMarkers > 0) & (minPercentageMatchedNspectraMarkers <= 100))) {
        FSA_message("ERROR!!! Problem with SPEC0010! This parameter should be a positive number greater than 0 and less than 100!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (!allowedNominalMass) {
      roundingDigitPrefiltering <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0011'), 2])
      if (length(roundingDigitPrefiltering) == 0) {
        FSA_message("ERROR!!! Problem with SPEC0011! This parameter should be either 0, 1 or 2!")
        checkpoint_parameter <- FALSE
      } else {
        if (!((roundingDigitPrefiltering == 0) | (roundingDigitPrefiltering == 1) | (roundingDigitPrefiltering == 2))) {
          FSA_message("ERROR!!! Problem with SPEC0011! This parameter should be either 0, 1 or 2!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      massError <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0012'), 2])
      if (length(massError) == 0) {
        FSA_message("ERROR!!! Problem with SPEC0012! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      } else {
        if (massError <= 0) {
          FSA_message("ERROR!!! Problem with SPEC0012! This parameter should be greater than `0 Da`!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x0013 <- which(PARAM_SPEC[, 1] == 'SPEC0013')
    allowedWeightedSpectralEntropy <- tolower(gsub(" ", "", PARAM_SPEC[x0013, 2]))
    if (allowedWeightedSpectralEntropy == "1" | allowedWeightedSpectralEntropy == "t" | allowedWeightedSpectralEntropy == "true") {
      allowedWeightedSpectralEntropy <- TRUE
    } else {
      allowedWeightedSpectralEntropy <- FALSE
    }
    PARAM_SPEC[x0013, 2] <- allowedWeightedSpectralEntropy
    ##
    minIonRangeDifference <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0014'), 2])
    if (length(minIonRangeDifference) == 0) {
      FSA_message("ERROR!!! Problem with SPEC0014! This parameter should be >= 0!")
      checkpoint_parameter <- FALSE
    } else {
      if (minIonRangeDifference < 0) {
        FSA_message("ERROR!!! Problem with SPEC0014! This parameter should be >= 0!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minMatchedNumPeaks <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0015'), 2])
    if (length(minMatchedNumPeaks) == 0) {
      FSA_message("ERROR!!! Problem with SPEC0015! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (minMatchedNumPeaks >= 1) {
        if ((minMatchedNumPeaks %% 1) != 0) {
          FSA_message("ERROR!!! Problem with SPEC0015! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        FSA_message("ERROR!!! Problem with SPEC0015! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minEntropySimilarity <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0016'), 2])
    if (length(minEntropySimilarity) == 0) {
      FSA_message("ERROR!!! Problem with SPEC0016! This parameter should be a positive number between 0 - 1!")
      checkpoint_parameter <- FALSE
    } else {
      if (!((minEntropySimilarity >= 0) & (minEntropySimilarity <= 1))) {
        FSA_message("ERROR!!! Problem with SPEC0016! This parameter should be a positive number between 0 - 1!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    minCosineSimilarity <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0017'), 2])
    if (length(minCosineSimilarity) == 0) {
      FSA_message("ERROR!!! Problem with SPEC0017! This parameter should be a positive number between 0 - 1!")
      checkpoint_parameter <- FALSE
    } else {
      if (!((minCosineSimilarity >= 0) & (minCosineSimilarity <= 1))) {
        FSA_message("ERROR!!! Problem with SPEC0017! This parameter should be a positive number between 0 - 1!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (!allowedNominalMass) {
      maxNEME <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0018'), 2])
      if (length(maxNEME) == 0) {
        FSA_message("ERROR!!! Problem with SPEC0018! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      } else {
        if (maxNEME < 0) {
          FSA_message("ERROR!!! Problem with SPEC0018! This parameter should be a positive number!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    maxAllowedNumberHits <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0019'), 2])
    if (length(maxAllowedNumberHits) == 0) {
      FSA_message("ERROR!!! Problem with SPEC0019! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (maxAllowedNumberHits < 0) {
        FSA_message("ERROR!!! Problem with SPEC0019! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (maxAllowedNumberHits > 0) {
      x0020 <- which(PARAM_SPEC[, 1] == 'SPEC0020')
      SPEC0020 <- PARAM_SPEC[x0020, 2]
      if (is.na(SPEC0020)) {
        FSA_message("ERROR!!! Problem with SPEC0020!")
        checkpoint_parameter <- FALSE
      } else {
        SPEC0020 <- gsub(" ", "", tolower(SPEC0020))
        PARAM_SPEC[x0020, 2] <- SPEC0020
      }
    }
  }
  ##############################################################################
  if (!checkpoint_parameter) {
    PARAM_SPEC <- NULL
  }
  ##
  return(PARAM_SPEC)
}