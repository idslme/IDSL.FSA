FSdb2PeakXcolSubsetter <- function(FSdb_address, peak_alignment_folder, metavariable = "idsl.ipa_collective_peakids", number_processing_threads = 1) {
  ##
  FSdb <- FSA_loadRdata(FSdb_address)
  if (!is.null(FSdb)) {
    file_name_hrms <- FSdb[["MSPLibraryParameters"]][["filename"]]
    if (is.null(file_name_hrms)) {
      file_name_hrms <- gsub("^CSA_MSP_|^DDA_MSP_|^DIA_MSP_|.msp$", "", FSdb[["MSPLibraryParameters"]][["MSPfilename"]], ignore.case = TRUE)
    }
    if (is.null(file_name_hrms)) {
      stop(FSA_logRecorder("No metavariable representing file names was detected to subset aligned tables!!!"))
    }
    listFileNameHRMS <- FSA_R.aggregate(file_name_hrms)
    uniqueFileNameHRMS <- names(listFileNameHRMS)
    ##
    metavariable <- tolower(metavariable)
    ##
    peakXcol <- FSA_loadRdata(paste0(peak_alignment_folder, "/peakXcol.Rdata"))
    ##
    ############################################################################
    ############################################################################
    ##
    call_FSdb2PeakXcolSubsetter <- function(i) {
      ##
      peakIDj <- do.call(c, lapply(listFileNameHRMS[[i]], function(j) {
        strIDj <- FSdb[["MSPLibraryParameters"]][[metavariable]][j]
        IDj <- eval(parse(text = paste0("c(", strIDj, ")")))
        IDj <- IDj[IDj != 0]
        IDj[1]
      }))
      ##
      peakXcolID <- peakXcol[, i]
      ##
      detectedIDs <- which(peakXcolID %in% peakIDj)
      ##
      return(detectedIDs)
    }
    ##
    ############################################################################
    ############################################################################
    ##
    if (number_processing_threads == 1) {
      ##
      subsetAlignedPeakIDs <- do.call(c, lapply(uniqueFileNameHRMS, function(i) {
        call_FSdb2PeakXcolSubsetter(i)
      }))
      ##
    } else {
      ## Processing OS
      osType <- Sys.info()[['sysname']]
      ##
      ##########################################################################
      ##
      if (osType == "Windows") {
        clust <- makeCluster(number_processing_threads)
        clusterExport(clust, setdiff(ls(), c("clust", "uniqueFileNameHRMS")), envir = environment())
        ##
        subsetAlignedPeakIDs <- do.call(c, parLapply(clust, uniqueFileNameHRMS, function(i) {
          call_FSdb2PeakXcolSubsetter(i)
        }))
        ##
        stopCluster(clust)
        ##
        ########################################################################
        ##
      } else {
        ##
        subsetAlignedPeakIDs <- do.call(c, mclapply(uniqueFileNameHRMS, function(i) {
          call_FSdb2PeakXcolSubsetter(i)
        }, mc.cores = number_processing_threads))
        ##
        closeAllConnections()
        ##
        ########################################################################
        ##
      }
    }
    ##
    ############################################################################
    ############################################################################
    ##
    subsetAlignedPeakIDs <- sort(unique(subsetAlignedPeakIDs), decreasing = FALSE)
    peakXcol <- peakXcol[subsetAlignedPeakIDs, ]
    rownames(peakXcol) <- subsetAlignedPeakIDs
    ##
    peak_height <- FSA_loadRdata(paste0(peak_alignment_folder, "/peak_height.Rdata"))
    peak_height <- peak_height[subsetAlignedPeakIDs, ]
    rownames(peak_height) <- subsetAlignedPeakIDs
    ##
    peak_area <- FSA_loadRdata(paste0(peak_alignment_folder, "/peak_area.Rdata"))
    peak_area <- peak_area[subsetAlignedPeakIDs, ]
    rownames(peak_area) <- subsetAlignedPeakIDs
    ##
    peak_R13C <- FSA_loadRdata(paste0(peak_alignment_folder, "/peak_R13C.Rdata"))
    peak_R13C <- peak_R13C[subsetAlignedPeakIDs, ]
    rownames(peak_R13C) <- subsetAlignedPeakIDs
    ##
    listXcolHeightAreaR13C <- list(peakXcol, peak_height, peak_area, peak_R13C)
    names(listXcolHeightAreaR13C) <- c("peakXcol", "peak_height", "peak_area", "peak_R13C")
  } else {
    listXcolHeightAreaR13C <- NULL
    FSA_logRecorder("NULL FSDB!")
  }
  ##
  return(listXcolHeightAreaR13C)
}