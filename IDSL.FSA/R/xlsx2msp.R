xlsx2msp <- function(path, xlsxFileName = "", number_processing_threads = 1) {
  ##
  readxlPackageCheck <- tryCatch(requireNamespace('readxl', quietly = TRUE), error = function(e) {FALSE})
  if (!readxlPackageCheck) {
    warning("IDSL.FSA requires the 'readxl' package of R to read Excel spreadsheets!")
    stop(" <<< install.packages('readxl') >>> ")
  }
  ##
  xlsxFileLocation <- paste0(path, "/", xlsxFileName)
  xlsxFileLocation <- gsub("\\", "/", xlsxFileLocation, fixed = TRUE)
  strxlsxFileLocation <- strsplit(xlsxFileLocation, "/")[[1]]
  xlsxFileName <- strxlsxFileLocation[length(strxlsxFileLocation)]
  xlsxFileLocation <- paste0(strxlsxFileLocation, collapse = "/")
  ##
  msp_xlsx <- readxl::read_xlsx(xlsxFileLocation)
  msp_xlsx <- data.frame(msp_xlsx)
  xlsxColNames <- colnames(msp_xlsx)
  ##############################################################################
  x_IDcol <- which(xlsxColNames == "ID")
  x_mzFcol <- which(xlsxColNames == "mz_fragment")
  x_intFcol <- which(xlsxColNames == "int_fragment")
  x_namecol <- which(xlsxColNames == "Name")
  ##
  if ((length(x_IDcol) != 1) | (length(x_mzFcol) != 1) | (length(x_intFcol) != 1) | (length(x_namecol) != 1)) {
    stop("The xlsx file should have only one column for the following headers (case-senstive):
    'ID'
    'mz_fragment'
    'int_fragment'
    'Name'")
  }
  ##
  xlsxCol <- setdiff(xlsxColNames, c("ID", "mz_fragment", "int_fragment","Name"))
  if (length(xlsxCol) > 0) {
    xlsxColCheck <- TRUE
    x_xlsxCol <- do.call(c, lapply(xlsxCol, function(i) {
      which(xlsxColNames == i)
    }))
  } else {
    xlsxColCheck <- FALSE
  }
  ##
  msp_xlsx$ID <- as.numeric(msp_xlsx$ID)
  msp_xlsx <- msp_xlsx[order(msp_xlsx$ID), ]
  x_diffID <- c(0, which(abs(diff(msp_xlsx$ID)) > 0), length(msp_xlsx$ID))
  ##
  call_xlsx2msp <- function(i, msp_xlsx, x_diffID, x_namecol, xlsxColCheck,
                            xlsxColNames, x_xlsxCol, x_mzFcol, x_intFcol) {
    x_ID <- seq((x_diffID[i] + 1), x_diffID[i + 1], 1)
    ID_i <-  x_diffID[i] + 1
    ##
    MSPname <- paste0("Name: ", msp_xlsx[ID_i, x_namecol], "\n")
    ##
    if (xlsxColCheck) {
      MSPid <- do.call(paste0, lapply(x_xlsxCol, function(j) {
        paste0(xlsxColNames[j], ": ", msp_xlsx[ID_i, j], "\n")
      }))
    } else {
      MSPid <- ""
    }
    ##
    MSPid <- paste0(MSPid, "Num Peaks: ", length(x_ID), "\n")
    ##
    MSPid_mz_int <- paste0(msp_xlsx[x_ID, x_mzFcol], " ", msp_xlsx[x_ID, x_intFcol], "\n", collapse = "")
    ##
    paste0(MSPname, MSPid, MSPid_mz_int, "\n")
  }
  ##############################################################################
  if (number_processing_threads == 1) {
    MSPvector <- do.call(c, lapply(1:(length(x_diffID) - 1), function(i) {
      ##
      call_xlsx2msp(i, msp_xlsx, x_diffID, x_namecol, xlsxColCheck,
                    xlsxColNames, x_xlsxCol, x_mzFcol, x_intFcol)
    }))
  } else {
    ##
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Windows") {
      ##
      clust <- makeCluster(number_processing_threads)
      clusterExport(clust, setdiff(ls(), c("clust")), envir = environment())
      ##
      MSPvector <- do.call(c, parLapply(clust, 1:(length(x_diffID) - 1), function(i) {
        ##
        call_xlsx2msp(i, msp_xlsx, x_diffID, x_namecol, xlsxColCheck,
                      xlsxColNames, x_xlsxCol, x_mzFcol, x_intFcol)
      }))
      ##
      stopCluster(clust)
      ##
    } else {
      ##
      MSPvector <- do.call(c, mclapply(1:(length(x_diffID) - 1), function(i) {
        ##
        call_xlsx2msp(i, msp_xlsx, x_diffID, x_namecol, xlsxColCheck,
                      xlsxColNames, x_xlsxCol, x_mzFcol, x_intFcol)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
      ##
    }
  }
  ##
  mspFileName <- gsub("[.]xlsx$", ".msp", xlsxFileLocation, ignore.case = TRUE)
  tryCatch({
    write.table(MSPvector, file = mspFileName, quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
    FSA_message(paste0("`", xlsxFileName, "` was converted into .msp format and stored in the same folder!"), failedMessage = FALSE)
    ##
  }, error = function(e) {
    stop(paste0("The saving address is problematic for `", mspFileName, "`!"))
    }
  )
  ##
  return()
}