mspPosNegSplitter <- function(path, MSPfile = "", number_processing_threads = 1) {
  ##
  mspFileLocation <- paste0(path, "/", MSPfile)
  mspFileLocation <- gsub("\\", "/", mspFileLocation, fixed = TRUE)
  strmspFileLocation <- strsplit(mspFileLocation, "/")[[1]]
  mspFileName <- strmspFileLocation[length(strmspFileLocation)]
  mspFileLocation <- paste0(strmspFileLocation, collapse = "/")
  ##
  msp <- tryCatch(readLines(mspFileLocation, warn = FALSE),
                  error = function(e){stop(paste0("Problem with deconvoluting .msp file --> ", mspFileName))})
  msp <- c("", msp, "")
  ##############################################################################
  x_msp <- which((msp == "") | (msp == " ") | (msp == "  ") | (msp == "\t"))
  xdiff <- which(diff(x_msp) > 1)
  x1msp <- x_msp[xdiff] + 1
  x2msp <- x_msp[xdiff + 1] - 1
  ##############################################################################
  xIonModeNeg <- c(grep("IONMODE: N", msp, ignore.case = TRUE), grep("ION_MODE: N", msp, ignore.case = TRUE),
                   grep("IONMODE:  N", msp, ignore.case = TRUE), grep("ION_MODE:  N", msp, ignore.case = TRUE))
  ##
  xIonModeNeg <- unique(xIonModeNeg)
  ##
  if (length(xIonModeNeg) > 0) {
    ionModeNegCheck <- TRUE
  } else {
    ionModeNegCheck <- FALSE
    FSA_message("'Negative' mode blocks were not detected!", failedMessage = TRUE)
  }
  ##############################################################################
  xIonModePos <- c(grep("IONMODE: P", msp, ignore.case = TRUE), grep("ION_MODE: P", msp, ignore.case = TRUE),
                   grep("IONMODE:  P", msp, ignore.case = TRUE), grep("ION_MODE:  P", msp, ignore.case = TRUE))
  ##
  xIonModePos <- unique(xIonModePos)
  ##
  if (length(xIonModePos) > 0) {
    ionModePosCheck <- TRUE
  } else {
    ionModePosCheck <- FALSE
    FSA_message("'Positive' mode blocks were not detected!", failedMessage = TRUE)
  }
  ##############################################################################
  call_MSPblocks <- function(i) {
    x1 <- which(i > x1msp)
    x1 <- x1msp[length(x1)]
    x2 <- x2msp[which(i < x2msp)[1]]
    ##
    c(msp[x1:x2], "", "")
  }
  ##############################################################################
  if (number_processing_threads == 1) {
    ##
    if (ionModeNegCheck) {
      negativeMSPblocks <- do.call(c, lapply(xIonModeNeg, function(i) {
        call_MSPblocks(i)
      }))
    }
    ##
    if (ionModePosCheck) {
      positiveMSPblocks <- do.call(c, lapply(xIonModePos, function(i) {
        call_MSPblocks(i)
      }))
    }
    ##
  } else {
    ##
    osType <- Sys.info()[['sysname']]
    if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      if (ionModeNegCheck) {
        negativeMSPblocks <- foreach(i = xIonModeNeg, .combine = 'c', .verbose = FALSE) %dopar% {
          call_MSPblocks(i)
        }
      }
      ##
      if (ionModePosCheck) {
        positiveMSPblocks <- foreach(i = xIonModePos, .combine = 'c', .verbose = FALSE) %dopar% {
          call_MSPblocks(i)
        }
      }
      ##
      stopCluster(clust)
      ##
    } else if (osType == "Linux") {
      ##
      if (ionModeNegCheck) {
        negativeMSPblocks <- do.call(c, mclapply(xIonModeNeg, function(i) {
          call_MSPblocks(i)
        }, mc.cores = number_processing_threads))
      }
      ##
      if (ionModePosCheck) {
        positiveMSPblocks <- do.call(c, mclapply(xIonModePos, function(i) {
          call_MSPblocks(i)
        }, mc.cores = number_processing_threads))
      }
      ##
      closeAllConnections()
    }
  }
  ##############################################################################
  ##
  filename <- gsub("[.]msp$", "", mspFileLocation, ignore.case = TRUE)
  ##
  ##############################################################################
  if (ionModeNegCheck) {
    ##
    tryCatch(write.table(negativeMSPblocks, file = paste0(filename, "_Neg.msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE),
             error = function(e) {warning("The saving address is problematic!!!")})
  }
  ##
  if (ionModePosCheck) {
    ##
    tryCatch(write.table(positiveMSPblocks, file = paste0(filename, "_Pos.msp"), quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE),
             error = function(e) {warning("The saving address is problematic!!!")})
  }
  ##
  if ((ionModeNegCheck) | (ionModePosCheck)) {
    FSA_message(paste0("'", mspFileName, "' was splitted in positive and negative msp blocks and stored in the same folder!"), failedMessage = FALSE)
  }
  ##
  return()
}