FSA_FSdb_xlsxAnalyzer <- function(spreadsheet) {
  ##
  checkpoint_parameter <- FALSE
  #
  if (typeof(spreadsheet) == "list") {
    if (ncol(spreadsheet) >= 4) {
      PARAM_FSdb <- cbind(spreadsheet[, 2], spreadsheet[, 4])
      checkpoint_parameter <- TRUE
      ##
    } else if (ncol(spreadsheet) == 2) {
      PARAM_FSdb <- spreadsheet
      checkpoint_parameter <- TRUE
      ##
    } else {
      FSA_message("The `FSDB` spreadsheet tab was not produced properly!")
    }
  } else if (length(spreadsheet) == 1) {
    if (typeof(spreadsheet) == "character") {
      if (file.exists(spreadsheet)) {
        PARAM_FSdb <- readxl::read_xlsx(spreadsheet, sheet = "FSDB")
        PARAM_FSdb <- cbind(PARAM_FSdb[, 2], PARAM_FSdb[, 4])
        checkpoint_parameter <- TRUE
      } else {
        FSA_message("The `FSDB` spreadsheet tab not found! It should be an Excel file with .xlsx extention!")
      }
    } else {
      FSA_message("The `FSDB` spreadsheet tab was not produced properly!")
    }
  } else {
    FSA_message("The `FSDB` spreadsheet tab was not produced properly!")
  }
  if (checkpoint_parameter) {
    ################### MS/MS library import and export ########################
    x0001 <- which(PARAM_FSdb[, 1] == 'FSdb0001')
    if (length(x0001) == 0) {
      FSA_message("ERROR!!! Problem with FSdb0001!")
      checkpoint_parameter <- FALSE
    } else {
      input_path_library <- PARAM_FSdb[x0001, 2]
      input_path_library <- gsub("\\", "/", input_path_library, fixed = TRUE)
      PARAM_FSdb[x0001, 2] <- input_path_library
      if (!dir.exists(input_path_library)) {
        FSA_message("ERROR!!! Problem with FSdb0001! Please make sure the full path is provided!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0002 <- which(PARAM_FSdb[, 1] == 'FSdb0002')
    if (length(x0002) == 0) {
      FSA_message("ERROR!!! Problem with FSdb0002!")
      checkpoint_parameter <- FALSE
    } else {
      library_string <- PARAM_FSdb[x0002, 2]
    }
    ##
    if (tolower(library_string) == "all") {
      PARAM_FSdb[x0002, 2] <- "all"
      file_name_library_msp <- dir(path = input_path_library)
      file_name_library_msp <- file_name_library_msp[grepl(".msp$", file_name_library_msp, ignore.case = TRUE)]
      if (length(file_name_library_msp) > 0) {
        if (is.na(file_name_library_msp[1])) {
          FSA_message("ERROR!!! Problem with FSdb0002! No .msp file was detected in the designated folder!")
          checkpoint_parameter <- FALSE
        }
      } else {
        FSA_message("ERROR!!! Problem with FSdb0002! No .msp file was detected in the designated folder!")
        checkpoint_parameter <- FALSE
      }
    } else if (grepl(".msp", library_string, ignore.case = TRUE)) {
      file_name_library_msp <- strsplit(library_string, ";")[[1]]
      ID <- do.call(c, lapply(file_name_library_msp, function(i) {
        if (!file.exists(paste0(input_path_library, "/", i))) {
          i
        }
      }))
      if (length(ID) > 0) {
        FSA_message("ERROR!!! Problem with FSdb0002! not detected the following file(s) (case sensitive even for file extensions):")
        for (i in ID) {
          i
        }
        checkpoint_parameter <- FALSE
      }
    } else {
      FSA_message("ERROR!!! Problem with FSdb0002!!!")
      checkpoint_parameter <- FALSE
    }
    ##
    address_FSDB <- PARAM_FSdb[which(PARAM_FSdb[, 1] == "FSdb0003"), 2]
    exportFSdbCheck <- if (tolower(address_FSDB) == "na") {FALSE} else {TRUE}
    if (exportFSdbCheck) {
      tryCatch(dir.create(address_FSDB, recursive = TRUE), error = function(e){stop("Problem with FSdb0003! R cannot create the folder!")}, warning = function(w){NULL})
      if (!dir.exists(address_FSDB)) {
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0004 <- which(PARAM_FSdb[, 1] == 'FSdb0004')
    if (length(x0004) == 0) {
      FSA_message("ERROR!!! Problem with FSdb0004!")
      checkpoint_parameter <- FALSE
    }
    ##
    number_processing_threads <- as.numeric(PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0005'), 2])
    if (length(number_processing_threads) == 0) {
      FSA_message("ERROR!!! Problem with FSdb0005! This parameter should be a positive integer!")
      checkpoint_parameter <- FALSE
    } else {
      if (number_processing_threads >= 1) {
        if ((number_processing_threads %% 1) != 0) {
          FSA_message("ERROR!!! Problem with FSdb0005! This parameter should be a positive integer!")
          checkpoint_parameter <- FALSE
        }
      } else {
        FSA_message("ERROR!!! Problem with FSdb0005! This parameter should be at least 1 !")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    x0006 <- which(PARAM_FSdb[, 1] == 'FSdb0006')
    if (length(x0006) > 0) {
      allowedNominalMass <- tolower(gsub(" ", "", PARAM_FSdb[x0006, 2]))
      if (allowedNominalMass == "1" | allowedNominalMass == "t" | allowedNominalMass == "true") {
        allowedNominalMass <- TRUE
      } else {
        allowedNominalMass <- FALSE
      }
    } else {
      allowedNominalMass <- FALSE
    }
    PARAM_FSdb[x0006, 2] <- allowedNominalMass
    ##
    noiseRemovalRatio <- as.numeric(PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0007'), 2])
    if (length(noiseRemovalRatio) == 0) {
      FSA_message("ERROR!!! Problem with FSdb0007! This parameter should be a positive number!")
      checkpoint_parameter <- FALSE
    } else {
      if (noiseRemovalRatio < 0 | noiseRemovalRatio > 100) {
        FSA_message("ERROR!!! Problem with FSdb0007! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      }
    }
    ##
    if (!allowedNominalMass) {
      massIntegrationWindow <- as.numeric(PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0008'), 2])
      if (length(massIntegrationWindow) == 0) {
        FSA_message("ERROR!!! Problem with FSdb0008! This parameter should be a positive number!")
        checkpoint_parameter <- FALSE
      } else {
        if (massIntegrationWindow <= 0) {
          FSA_message("ERROR!!! Problem with FSdb0008! This parameter should be a positive number!")
          checkpoint_parameter <- FALSE
        }
      }
    }
    ##
    x0009 <- which(PARAM_FSdb[, 1] == 'FSdb0009')
    allowedWeightedSpectralEntropy <- tolower(gsub(" ", "", PARAM_FSdb[x0009, 2]))
    if (allowedWeightedSpectralEntropy == "1" | allowedWeightedSpectralEntropy == "t" | allowedWeightedSpectralEntropy == "true") {
      allowedWeightedSpectralEntropy <- TRUE
    } else {
      allowedWeightedSpectralEntropy <- FALSE
    }
    PARAM_FSdb[x0009, 2] <- allowedWeightedSpectralEntropy
  }
  ##############################################################################
  if (checkpoint_parameter == FALSE) {
    PARAM_FSdb <- NULL
  }
  ##
  return(PARAM_FSdb)
}