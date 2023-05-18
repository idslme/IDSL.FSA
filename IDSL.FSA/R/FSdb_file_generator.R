FSdb_file_generator <- function(PARAM_FSdb, output_path = NULL) {
  ##
  ##############################################################################
  ## To create log record for IDSL.FSA
  initiation_time_FSdb <- Sys.time()
  timeZone <- tryCatch(Sys.timezone(), warning = function(w) {"UTC"}, error = function(e) {"UTC"})
  input_path_library <- PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0001'), 2]
  address_FSDB <- PARAM_FSdb[which(PARAM_FSdb[, 1] == "FSdb0003"), 2]
  exportFSdbCheck <- if (tolower(address_FSDB) == "na") {FALSE} else {TRUE}
  ##
  if (!exportFSdbCheck) {
    address_FSDB <- output_path
  }
  .logFSA <- NULL
  .logFSA <<- paste0(address_FSDB, "/logFSA_FSDB.txt")
  FSA_logRecorder(paste0(rep("", 100), collapse = "="))
  FSA_logRecorder("Type <<< citation('IDSL.FSA') >>> for citing this R package in publications.")
  FSA_logRecorder(paste0("Reference MSP libraries: ", input_path_library))
  if (exportFSdbCheck) {
    FSA_logRecorder(paste0("OUTPUT: ", address_FSDB))
  } else {
    FSA_logRecorder("OUTPUT: Not saving the FSDB!")
  }
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  FSA_logRecorder("Initiated generating the FSDB workflow!")
  FSA_logRecorder(paste0(as.character(initiation_time_FSdb), " ", timeZone))
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder(paste0(PARAM_FSdb[, 1], "\t", PARAM_FSdb[, 2]),  allowedPrinting = FALSE)
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  ##
  ##############################################################################
  ##
  FSA_logRecorder("Initiated deconvoluting the `.msp` reference libraries!")
  ##
  library_string <- PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0002'), 2]
  ##
  if (tolower(library_string) == "all") {
    file_name_library_msp <- dir(path = input_path_library)
    file_name_library_msp <- file_name_library_msp[grepl(".msp$", file_name_library_msp, ignore.case = TRUE)]
  } else if (grepl(".msp$", library_string, ignore.case = TRUE)) {
    file_name_library_msp <- strsplit(library_string, ";")[[1]]
  }
  ##
  NPT <- as.numeric(PARAM_FSdb[which(PARAM_FSdb[, 1] == "FSdb0005"), 2])
  ##
  x0006 <- which(PARAM_FSdb[, 1] == 'FSdb0006')
  if (length(x0006) > 0) {
    allowedNominalMass <- eval(parse(text = PARAM_FSdb[x0006, 2]))
  } else {
    allowedNominalMass <- FALSE
  }
  ##
  noiseRemovalRatio <- as.numeric(PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0007'), 2])/100
  ##
  if (allowedNominalMass) {
    massIntegrationWindow_lib <- 0
  } else {
    massIntegrationWindow_lib <- as.numeric(PARAM_FSdb[which(PARAM_FSdb[, 1] == "FSdb0008"), 2])
  }
  ##
  allowedWeightedSpectralEntropy <- eval(parse(text = PARAM_FSdb[which(PARAM_FSdb[, 1] == 'FSdb0009'), 2]))
  ##
  libFSdb <- msp2FSdb(path = input_path_library, file_name_library_msp, massIntegrationWindow_lib, allowedNominalMass, allowedWeightedSpectralEntropy, noiseRemovalRatio, number_processing_threads = NPT)
  ##
  FSA_logRecorder(paste0("Completed deconvoluting ", length(libFSdb[["Num Peaks"]])," MSP reference blocks!"))
  ##
  if (exportFSdbCheck) {
    name_FSDB <- PARAM_FSdb[which(PARAM_FSdb[, 1] == "FSdb0004"), 2]
    FSdb_file <- paste0(address_FSDB, "/", name_FSDB, ".Rdata")
    FSdb_file <- gsub("\\", "/", FSdb_file, fixed = TRUE)
    FSA_logRecorder("Saving the FSDB library!")
    save(libFSdb, file = FSdb_file)
  }
  ##
  ##############################################################################
  ##
  completion_time_FSdb <- Sys.time()
  FSA_logRecorder(paste0(rep("", 100), collapse = "-"))
  required_time <- completion_time_FSdb - initiation_time_FSdb
  FSA_logRecorder(paste0("The required processing time was `", required_time, " ", attributes(required_time)$units, "`"))
  FSA_logRecorder(paste0(as.character(completion_time_FSdb), " ", timeZone), allowedPrinting = FALSE)
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder("", allowedPrinting = FALSE)
  FSA_logRecorder("Completed generating the FSDB workflow!")
  FSA_logRecorder(paste0(rep("", 100), collapse = "="), allowedPrinting = FALSE)
  ##
  ##############################################################################
  ##
  return(libFSdb)
}