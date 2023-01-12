FSdb2precursorType <- function(InChIKeyVector, libFSdb, tableIndicator = "Frequency", number_processing_threads = 1) {
  ##
  precursor_type_frequency <- NULL
  #
  AncillaryTable <- libFSdb[["MSPLibraryParameters"]]
  colAT <- colnames(AncillaryTable)
  x_inchikey <- grep("inchikey", colAT, ignore.case = TRUE)
  if (length(x_inchikey) == 0) {
    stop("There is no `inchikey` information in the FSDB!")
  }
  ##
  libPrecursorType <- libFSdb[["Precursor Type"]]
  uLibPrecursorType <- unique(libPrecursorType)
  uLibPrecursorType <- setdiff(uLibPrecursorType, "")
  if (length(uLibPrecursorType) == 0) {
    stop("There is no `precursor_type` information in the FSDB!")
  }
  ##
  libInChIKey <- AncillaryTable[, x_inchikey]
  strLibInChIKey <- strsplit(libInChIKey, "-")
  libInChIKey14 <- do.call(c, lapply(strLibInChIKey, function(i) {i[1]}))
  ##
  listLibInChIKeyID <- FSA_R.aggregate(libInChIKey14)
  ##
  L_InChIKey <- length(InChIKeyVector)
  strInChIKeyVector <- strsplit(InChIKeyVector, "-")
  InChIKeyVector14 <- do.call(c, lapply(strInChIKeyVector, function(i) {i[1]}))
  ##
  if (tolower(tableIndicator) == "frequency") {
    precursor_type_frequency_call <- function(i) {
      ##
      tPrecursorType <- table(libPrecursorType[listLibInChIKeyID[[i]]])
      if (length(tPrecursorType) > 0) {
        nPrecursorType <- names(tPrecursorType)
        xNULL <- which(nPrecursorType == "")
        if (length(xNULL) > 0) {
          nPrecursorType <- setdiff(nPrecursorType, "")
          if (length(nPrecursorType) > 0) {
            tPrecursorType <- do.call(c, lapply(nPrecursorType, function(j) {
              tPrecursorType[[j]]
            }))
            names(tPrecursorType) <- nPrecursorType
          } else {
            nPrecursorType <- NULL
          }
        }
      } else {
        nPrecursorType <- NULL
      }
      ##
      do.call(cbind, lapply(uLibPrecursorType, function(j) {
        jCheck <- j %in% nPrecursorType
        if (jCheck) {
          tPrecursorType[[j]]
        } else {
          0
        }
      }))
    }
  } else if (tolower(tableIndicator) == "precursormz") {
    precursormz <- libFSdb[["PrecursorMZ"]]
    ##
    precursor_type_frequency_call <- function(i) {
      ##
      libID <- listLibInChIKeyID[[i]]
      if (!is.null(libID)) {
        precursorMZtype <- data.frame(PrecursorMZ = precursormz[libID], PrecursorType = libPrecursorType[libID])
        precursorMZtype.agg <- aggregate(precursorMZtype$PrecursorMZ, by = list(precursorMZtype$PrecursorType), FUN = c) ## aggregation function
        ##
        xNULL <- which(precursorMZtype.agg$Group.1 == "")
        if (length(xNULL) > 0) {
          precursorMZtype.agg <- precursorMZtype.agg[-xNULL, ]
        }
        ##
        if (length(precursorMZtype.agg) > 0) {
          listPrecursorMZtype <- precursorMZtype.agg$x
          namesListPrecursorMZtype <- precursorMZtype.agg$Group.1
          names(listPrecursorMZtype) <- namesListPrecursorMZtype
        } else {
          namesListPrecursorMZtype <- NULL
        }
        ##
      } else {
        namesListPrecursorMZtype <- NULL
      }
      ##
      do.call(cbind, lapply(uLibPrecursorType, function(j) {
        jCheck <- j %in% namesListPrecursorMZtype
        if (jCheck) {
          precMZ <- listPrecursorMZtype[[j]]
          precMZ <- precMZ[!is.infinite(precMZ)]
          if (length(precMZ) > 0) {
            round(median(precMZ), digits = 5)
          } else {
            0
          }
        } else {
          0
        }
      }))
    }
  } else {
    stop("tableIndicator should be 'Frequency' or 'PrecursorMZ'!")
  }
  ##############################################################################
  if (number_processing_threads == 1) {
    ##
    precursor_type_frequency <- do.call(rbind, lapply(InChIKeyVector14, function(i) {
      precursor_type_frequency_call(i)
    }))
    ##
  } else {
    osType <- Sys.info()[['sysname']]
    if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      precursor_type_frequency <- foreach(i = InChIKeyVector14, .combine = 'rbind', .verbose = FALSE) %dopar% {
        precursor_type_frequency_call(i)
      }
      ##
      stopCluster(clust)
      ##
    } else if (osType == "Linux") {
      ##
      precursor_type_frequency <- do.call(rbind, mclapply(InChIKeyVector14, function(i) {
        precursor_type_frequency_call(i)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
    }
  }
  ##
  if (length(precursor_type_frequency) > 0) {
    precursor_type_frequency <- data.frame(cbind(InChIKeyVector, precursor_type_frequency))
    colnames(precursor_type_frequency) <- c("InChIKey", uLibPrecursorType)
    rownames(precursor_type_frequency) <- NULL
    ##
    x_non0 <- do.call(c, lapply(1:(length(uLibPrecursorType) + 1), function(i) {
      x_pt <- which(precursor_type_frequency[, i] == 0)
      if (length(x_pt) != L_InChIKey) {
        i
      }
    }))
    ##
    precursor_type_frequency <- precursor_type_frequency[, x_non0]
  }
  return(precursor_type_frequency)
}