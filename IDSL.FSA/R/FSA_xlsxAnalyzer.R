FSA_xlsxAnalyzer <- function(spreadsheet) {
  ##
  FSA_message("Initiated testing the IDSL.FSA workflow spreadsheet consistency!", failedMessage = FALSE)
  ##
  checkpoint_parameter <- FALSE
  ##
  if (length(spreadsheet) >= 4) {
    if (typeof(spreadsheet) == "list") {
      PARAM_FSA <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
    } else {
      FSA_message("The IDSL.FSA workflow spreadsheet tab was not produced properly!")
      checkpoint_parameter <- FALSE
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)) {
        ##
        readxlPackageCheck <- tryCatch(requireNamespace('readxl', quietly = TRUE), error = function(e) {FALSE})
        if (!readxlPackageCheck) {
          warning("IDSL.FSA requires the 'readxl' package of R to read Excel spreadsheets!")
          stop(" <<< install.packages('readxl') >>> ")
        }
        ##
        PARAM_FSA <- readxl::read_xlsx(spreadsheet, sheet = "Start")
        PARAM_FSA <- cbind(PARAM_FSA[, 2], PARAM_FSA[, 4])
        checkpoint_parameter <- TRUE
      } else {
        FSA_message("The IDSL.FSA workflow spreadsheet not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      FSA_message("The IDSL.FSA workflow spreadsheet was not produced properly!")
    }
  } else {
    FSA_message("The IDSL.FSA workflow spreadsheet was not produced properly!")
  }
  ##
  if (checkpoint_parameter) {
    ############################################################################
    x0001 <- which(PARAM_FSA[, 1] == 'FSA0001')
    if (length(x0001) == 0) {
      FSA_message("ERROR!!! Problem with FSA0001!")
      checkpoint_parameter <- FALSE
    } else {
      FSA0001 <- tolower(PARAM_FSA[x0001, 2])
      if (FSA0001 == "yes" | FSA0001 == "no") {
        PARAM_FSA[x0001, 2] <- FSA0001
      } else {
        FSA_message("ERROR!!! Problem with FSA0001!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (FSA0001 == "yes") {
      FSA_message("Initiated testing the `FSDB` spreadsheet tab consistency!", failedMessage = FALSE)
      ##
      PARAM_FSdb <- FSA_FSdb_xlsxAnalyzer(spreadsheet)
      if (is.null(PARAM_FSdb)) {
        FSA_message("ERROR!!! Problem with the `FSDB` spreadsheet tab!")
        checkpoint_parameter <- FALSE
      }
    } else if (FSA0001 == "no") {
      PARAM_FSdb <- NULL
      ##
      x0002 <- which(PARAM_FSA[, 1] == 'FSA0002')
      ##
      if (length(x0002) == 0) {
        FSA_message("ERROR!!! Problem with FSA0002!")
        checkpoint_parameter <- FALSE
      } else {
        ##
        FSdb_file <- PARAM_FSA[x0002, 2]
        FSdb_file <- gsub("\\", "/", FSdb_file, fixed = TRUE)
        PARAM_FSA[x0002, 2] <- FSdb_file
        if (!file.exists(FSdb_file)) {
          FSA_message("ERROR!!! Problem with FSA0002! Please ensure the full path is provided for the FSDB in .Rdata format OR select 'YES' for FSA0004!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    FSA_message("Initiated testing the `SpectraSimilarity` spreadsheet tab consistency!", failedMessage = FALSE)
    ##
    PARAM_SPEC <- FSA_SpectraSimilarity_xlsxAnalyzer(spreadsheet)
    if (is.null(PARAM_SPEC)) {
      FSA_message("ERROR!!! Problem with the `SpectraSimilarity` spreadsheet tab!")
      checkpoint_parameter <- FALSE
    }
    ##
    ############################################################################
    ##
    if (!is.null(PARAM_FSdb) & !is.null(PARAM_SPEC)) {
      ##
      FSdb0006 <- eval(parse(text = PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0006'), 2]))
      SPEC0003 <- eval(parse(text = PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0003'), 2]))
      if (FSdb0006 != SPEC0003) {
        FSA_message("Inconsistency between `FSdb0006` and `SPEC0003` parameters in the `FSDB` and `SpectraSimilarity` tabs!")
        checkpoint_parameter <- FALSE
      }
      ##
      FSdb0007 <- as.numeric(PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0007'), 2])
      SPEC0008 <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0008'), 2])
      if (FSdb0007 != SPEC0008) {
        FSA_message("Inconsistency between `FSdb0007` and `SPEC0008` parameters in the `FSDB` and `SpectraSimilarity` tabs!")
        checkpoint_parameter <- FALSE
      }
      ##
      if (FSdb0006 & SPEC0003) {
        FSdb0008 <- as.numeric(PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0008'), 2])
        SPEC0012 <- as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0012'), 2])
        if (FSdb0008 != SPEC0012) {
          FSA_message("Inconsistency between `FSdb0008` and `SPEC0012` parameters in the `FSDB` and `SpectraSimilarity` tabs!")
          checkpoint_parameter <- FALSE
        }
      }
      ##
      FSdb0009 <- eval(parse(text = PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0009'), 2]))
      SPEC0013 <- eval(parse(text = PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0013'), 2]))
      if (FSdb0009 != SPEC0013) {
        FSA_message("Inconsistency between `FSdb0009` and `SPEC0013` parameters in the `FSDB` and `SpectraSimilarity` tabs!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    ############################################################################
    ##
    x0003 <- which(PARAM_FSA[, 1] == 'FSA0003')
    if (length(x0003) == 0) {
      FSA_message("ERROR!!! Problem with FSA0003!")
      checkpoint_parameter <- FALSE
    } else {
      address_input_msp <- PARAM_FSA[x0003, 2]
      address_input_msp <- gsub("\\", "/", address_input_msp, fixed = TRUE)
      PARAM_FSA[x0003, 2] <- address_input_msp
      if (!dir.exists(address_input_msp)) {
        FSA_message("ERROR!!! Problem with FSA0003! Please make sure the full path is provided!")
        checkpoint_parameter <- FALSE
      } else {
        file_name_sample_msp <- dir(path = address_input_msp)
        file_name_sample_msp <- file_name_sample_msp[c(grepl(".msp$", file_name_sample_msp, ignore.case = TRUE), grepl(".mgf$", file_name_sample_msp, ignore.case = TRUE))]
        if (length(file_name_sample_msp) == 0) {
          FSA_message("ERROR!!! Problem with FSA0003! No MSP or MGF file was detected in the designated folder!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x0004 <- which(PARAM_FSA[, 1] == 'FSA0004')
    if (length(x0004) == 0) {
      FSA_message("ERROR!!! Problem with FSA0004!")
      checkpoint_parameter <- FALSE
    } else {
      output_sample <- PARAM_FSA[x0004, 2]
      output_sample <- gsub("\\", "/", output_sample, fixed = TRUE)
      PARAM_FSA[x0004, 2] <- output_sample
      if (!dir.exists(output_sample)) {
        tryCatch(dir.create(output_sample, recursive = TRUE), warning = function(w){warning("Problem with FSA0004! R cannot create the folder!")})
        if (!dir.exists(output_sample)) {
          checkpoint_parameter <- FALSE
        }
      }
    }
  }
  ##
  ##############################################################################
  ##
  if (!checkpoint_parameter) {
    PARAM_total <- vector(mode = "list", 3)
    ##
  } else {
    ##
    PARAM_total <- list(PARAM_FSA, PARAM_FSdb, PARAM_SPEC)
    ##
    if (FSA0001 == "yes") {
      FSA_message("The `FSDB` spreadsheet tab is consistent with the IDSL.FSA workflow!", failedMessage = FALSE)
    }
    ##
    allowedNominalMass <- eval(parse(text = (PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0003'), 2])))
    massErrorPrecursor <- tryCatch(as.numeric(PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0005'), 2]), warning = function(w){NA})
    ##
    if (allowedNominalMass) {
      ##
      if (!is.na(massErrorPrecursor)) {
        FSA_message("NOTICE: Precursor m/z match (SPEC0005) was selected for spectra annotations! Empty annotation tables will be generated when .msp files do not contain 'PrecursorMZ' information!", failedMessage = FALSE)
      } else {
        FSA_message("NOTICE: Precursor m/z match (SPEC0005) was not selected for spectra annotations!", failedMessage = FALSE)
      }
      ##
    } else {
      ##
      if (!is.na(massErrorPrecursor)) {
        FSA_message(paste0("NOTICE: Mass accuracy = '" , massErrorPrecursor, " Da' for precursor m/z (SPEC0005) was selected for spectra annotations! Empty annotation tables will be generated when .msp files do not contain 'PrecursorMZ' information!"), failedMessage = FALSE)
      } else {
        FSA_message("NOTICE: Mass accuracy for precursor m/z (SPEC0005) was not selected and precursor values will not be used for spectra annotations!", failedMessage = FALSE)
      }
      ##
    }
    ##
    RTtolerance <- PARAM_SPEC[which(PARAM_SPEC[, 1] == 'SPEC0006'), 2]
    if (!is.na(RTtolerance)) {
      FSA_message(paste0("NOTICE: Retention time tolerance = '" , RTtolerance, " min' (SPEC0006) was selected for spectra annotations! Empty annotation tables will be generated when .msp files do not contain 'Retention Time' information!"), failedMessage = FALSE)
    } else {
      FSA_message("NOTICE: Retention time tolerance (SPEC0006) was not selected and retention time values will not be used for spectra annotations!", failedMessage = FALSE)
    }
    ##
    FSA_message("The `SpectraSimilarity` spreadsheet tab is consistent with the IDSL.FSA workflow!", failedMessage = FALSE)
    ##
    FSA_message("The FSA spreadsheet is consistent with the IDSL.FSA workflow!", failedMessage = FALSE)
  }
  ##
  names(PARAM_total) <- c("PARAM_FSA", "PARAM_FSdb", "PARAM_SPEC")
  ##############################################################################
  return(PARAM_total)
}