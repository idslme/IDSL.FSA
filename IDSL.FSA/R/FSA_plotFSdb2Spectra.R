FSA_plotFSdb2Spectra <- function(path, allowedUnlink = TRUE, annexName = "", FSdb, selectedFSdbIDs = NULL, number_processing_threads = 1, allowedVerbose = TRUE) {
  ##
  dev.offCheck <- TRUE
  while (dev.offCheck) {
    dev.offCheck <- tryCatch(dev.off(), error = function(e) {FALSE})
  }
  ##
  FSA_dir.create(path, allowedUnlink)
  ##
  call_plotFSdb2SpectraCore <- function(i) {
    filenameSpectraPNG <- paste0(path, "/", annexName, "_", i, "_.png")
    fileCreateRCheck <- file.create(file = filenameSpectraPNG, showWarnings = FALSE)
    if (fileCreateRCheck) {
      png(filenameSpectraPNG, width = 16, height = 8, units = "in", res = 100)
      ##
      plotFSdb2SpectraCore(FSdb, index = i)
      ##
      dev.off()
    }
    return()
  }
  ##
  nSpectra <- length(FSdb[["Num Peaks"]])
  if (is.null(selectedFSdbIDs)) {
    selectedFSdbIDs <- seq(1, nSpectra, 1)
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  if (number_processing_threads == 1) {
    ##
    if (allowedVerbose) {progressBARboundaries <- txtProgressBar(min = 0, max = nSpectra, initial = 0, style = 3)}
    ##
    for (i in selectedFSdbIDs) {
      if (allowedVerbose) {setTxtProgressBar(progressBARboundaries, i)}
      ##
      null_variable <- call_plotFSdb2SpectraCore(i)
    }
    ##
    if (allowedVerbose) {close(progressBARboundaries)}
    ##
  } else {
    ## Processing OS
    osType <- Sys.info()[['sysname']]
    ##
    ############################################################################
    ##
    if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust", "selectedFSdbIDs")), envir = environment())
      ##
      null_variable <- parLapply(clust, selectedFSdbIDs, function(i) {
        call_plotFSdb2SpectraCore(i)
      })
      ##
      stopCluster(clust)
      ##
      ##########################################################################
      ##
    } else {
      ##
      null_variable <- mclapply(selectedFSdbIDs, function(i) {
        call_plotFSdb2SpectraCore(i)
      }, mc.cores = number_processing_threads)
      ##
      closeAllConnections()
      ##
      ##########################################################################
      ##
    }
  }
  ##
  ##############################################################################
  ##############################################################################
  ##
  return()
}