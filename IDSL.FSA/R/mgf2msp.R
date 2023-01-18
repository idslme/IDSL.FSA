mgf2msp <- function(path, MGFfile = "") {
  ##
  MGFFileLocation <- paste0(path, "/", MGFfile)
  MGFFileLocation <- gsub("\\", "/", MGFFileLocation, fixed = TRUE)
  strMGFFileLocation <- strsplit(MGFFileLocation, "/")[[1]]
  MGFfile <- strMGFFileLocation[length(strMGFFileLocation)]
  MGFFileLocation <- paste0(strMGFFileLocation, collapse = "/")
  ##
  mgf <- tryCatch(readLines(MGFFileLocation, warn = FALSE),
                  error = function(e){FSA_logRecorder(paste0("Problem with deconvoluting `.mgf` file --> `", MGFfile, "`!"))})
  ##
  ##############################################################################
  ##
  xPEPMASS <- grep("PEPMASS=", mgf, ignore.case = TRUE)
  if (length(xPEPMASS) > 0) {
    precursorCheck <- TRUE
    PrecPep <- c("PrecursorMZ: ", "PrecursorIntensity: ", paste0("PREC", seq(1, 10, 1), ": "))
    ##
    strPEPMASS <- strsplit(gsub("PEPMASS=", "", mgf[xPEPMASS]), " ")
    ##
    PEPMASS <- do.call(rbind, lapply(strPEPMASS, function(i) {i}))
    Lpep <- dim(PEPMASS)[2]
    ##
    PrecPep <- PrecPep[seq(1, Lpep, 1)]
    ##
    xRemove <- grep("PEPMASS=", mgf, ignore.case = TRUE)
    mgf <- mgf[-xRemove]
  } else {
    precursorCheck <- FALSE
  }
  ##
  ##############################################################################
  ##
  scan <- NULL
  xName <- grep("NAME=", mgf, ignore.case = TRUE)
  ##
  if (length(xName) == 0) {
    ##
    mgf <- gsub('"', '', mgf)
    ##
    xTITLE <- grep("TITLE=", mgf, ignore.case = TRUE)
    ##
    if (length(xTITLE) > 0) {
      xSCAN <- grep("SCAN=", mgf, ignore.case = TRUE)
      ##
      if (length(xSCAN) == 0) {
        scan <- do.call(c, lapply(xTITLE, function(i) {
          xScan <- FSA_locate_regex(mgf[i], "scan=", ignore.case = TRUE)[2]
          substr(mgf[i], (xScan + 1), nchar(mgf[i]))
        }))
      }
      ##
      mgf <- gsub("TITLE=", "NAME=", mgf, ignore.case = TRUE)
    }
  }
  ##
  if (is.null(scan)) {
    scanCheck <- FALSE
  } else {
    scanCheck <- TRUE
  }
  ##
  ##############################################################################
  ##
  xBeginIons <- grep("BEGIN IONS", mgf, ignore.case = TRUE)
  xNameGrep <- grep("NAME=", mgf, ignore.case = TRUE)
  mgf <- gsub("NAME=", "Name=", mgf, ignore.case = TRUE)
  ##
  xName <- do.call(c, lapply(xNameGrep, function(j) {
    if (substr(mgf[j], 1, 5) == "Name=") {
      j
    }
  }))
  ##
  xEqualSign <- grep("=", mgf)
  xNum <- grep("^[[:digit:]]", mgf)
  xNum <- setdiff(xNum, xEqualSign)
  xDiff <- which(diff(xNum) > 1)
  xNum1 <- c(xNum[1], xNum[xDiff + 1])
  xNum2 <- c(xNum[xDiff], xNum[length(xNum)])
  ##
  ##############################################################################
  ##
  mgf <- gsub("=", ": ", mgf)
  ##
  MSP <- do.call(c, lapply(1:length(xBeginIons), function(i) {
    pasteMSP <- mgf[xName[i]]
    ##
    if (precursorCheck) {
      pastePrecursor <- do.call(c, lapply(1:Lpep, function(j) {
        paste0(PrecPep[j], PEPMASS[i, j])
      }))
      pastePrecursor <- paste0(pastePrecursor, collapse = "\n")
      pasteMSP <- paste0(pasteMSP, "\n", pastePrecursor)
    }    
    ##
    xMGFblock <- setdiff(seq((xBeginIons[i] + 1), (xNum1[i] - 1), 1), xName[i])
    pasteMSP <- paste0(pasteMSP, "\n", paste0(mgf[xMGFblock], collapse = "\n"))
    ##
    if (scanCheck) {
      pasteMSP <- paste0(pasteMSP, "\nscan: ", scan[i])
    }
    ##
    pasteMSP <- paste0(pasteMSP, "\nNum Peaks: ", (xNum2[i] - xNum1[i] + 1), "\n")
    ##
    paste0(pasteMSP, paste0(mgf[xNum1[i]:xNum2[i]], collapse = "\n"), "\n\n")
  }))
  ##
  mspFileName <- gsub("[.]MGF$", ".msp", MGFFileLocation, ignore.case = TRUE)
  ##
  tryCatch({write.table(MSP, file = mspFileName, quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)
    FSA_logRecorder(paste0("'", MGFfile, "' was converted into `.msp` format and stored in the same folder!"))},
    ##
    error = function(e) {FSA_logRecorder(paste0("The saving address is problematic for `", mspFileName, "`!!!"))})
  ##
  return()
}