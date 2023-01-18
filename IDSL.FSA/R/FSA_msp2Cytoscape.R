FSA_msp2Cytoscape <- function(path, MSPfile = "", mspVariableVector = NULL, mspNodeID = NULL, massError = 0.01, RTtolerance = NA, minEntropySimilarity = 0.75,
                              allowedNominalMass = FALSE, allowedWeightedSpectralEntropy = TRUE, noiseRemovalRatio = 0.01, number_processing_threads = 1) {
  ##
  mspFileLocation <- paste0(path, "/", MSPfile)
  mspFileLocation <- gsub("\\", "/", mspFileLocation, fixed = TRUE)
  strmspFileLocation <- strsplit(mspFileLocation, "/")[[1]]
  LstrmspFileLocation <- length(strmspFileLocation)
  mspFileName <- strmspFileLocation[LstrmspFileLocation]
  mspFilePath <- paste0(strmspFileLocation[1:(LstrmspFileLocation - 1)], collapse = "/")
  ##
  ##############################################################################
  ##
  strmspFileName <- strsplit(mspFileName, "[.]")[[1]]
  formatFileName <- tolower(strmspFileName[length(strmspFileName)])
  if (formatFileName == "msp") {
    FSdb <- msp2FSdb(path = mspFilePath, MSPfile_vector = mspFileName, massIntegrationWindow = massError,
                     allowedNominalMass, allowedWeightedSpectralEntropy, noiseRemovalRatio, number_processing_threads)
  } else if (formatFileName == "rdata") {
    FSdb <- FSA_loadRdata(paste0(mspFilePath, "/", mspFileName))
  } else {
    stop("Allowed inputs : FSDB in `.Rdata` or MSP files in `.msp` formats!")
  }
  ##
  if (!is.na(RTtolerance)) {
    RTcheck <- TRUE
  } else {
    RTcheck <- FALSE
  }
  ##
  lengthFSdb <- length(FSdb[["Spectral Entropy"]])
  if (lengthFSdb > 2) {
    ##
    xIDfsdb <- which(colnames(FSdb[["MSPLibraryParameters"]]) == tolower(mspNodeID))
    if ((is.null(mspNodeID)) | (length(xIDfsdb) == 0)) {
      mspNodeID <- "IDcytoscape"
      FSdb[["MSPLibraryParameters"]][["idcytoscape"]] <- seq(1, lengthFSdb, 1)
    }
    ##
    variableVector <- tolower(c(mspNodeID, mspVariableVector))
    ##
    call_cytoscape_df <- function(i) {
      ##
      RTi <- FSdb[["Retention Time"]][i]
      if (RTcheck) {
        RTiInfCheck <- if (!is.infinite(RTi)) {TRUE} else {FALSE}
      }
      ##
      S_PEAK_A <- FSdb[["Spectral Entropy"]][i]
      PEAK_A <- FSdb[["FragmentList"]][[i]]
      ##
      do.call(rbind, lapply((i + 1):lengthFSdb, function(j) {
        ##
        S_PEAK_B <- FSdb[["Spectral Entropy"]][j]
        ##
        if (abs(S_PEAK_B - S_PEAK_A) <= 2) {
          ##
          RTj <- FSdb[["Retention Time"]][j]
          ##
          passLoop <- TRUE
          if (RTcheck) {
            if (RTiInfCheck) {
              if (!is.infinite(RTj)) {
                if (abs(RTi - RTj) > RTtolerance) {
                  passLoop <- FALSE
                }
              } else {
                passLoop <- FALSE
              }
            } else {
              passLoop <- FALSE
            }
          }
          ##
          if (passLoop) {
            ##
            PEAK_B <- FSdb[["FragmentList"]][[j]]
            ##
            SESS <- spectral_entropy_similarity_score(PEAK_A, S_PEAK_A, PEAK_B, S_PEAK_B, massError, allowedNominalMass)
            ##
            if (SESS >= minEntropySimilarity) {
              ##
              clapplyVar <- do.call(c, lapply(variableVector, function(k) {
                c(FSdb[["MSPLibraryParameters"]][[k]][[i]], FSdb[["MSPLibraryParameters"]][[k]][[j]])
              }))
              ##
              c(SESS, clapplyVar)
            }
          }
        }
      }))
    }
    ##
    if (number_processing_threads == 1) {
      ##
      ##########################################################################
      ##
      cytoscape_df <- do.call(rbind, lapply(1:(lengthFSdb - 1), function(i) {
        call_cytoscape_df(i)
      }))
      ##
      ##########################################################################
      ##
    } else {
      ##
      osType <- Sys.info()[['sysname']]
      ##
      if (osType == "Linux") {
        ##
        ########################################################################
        ##
        cytoscape_df <- do.call(rbind, mclapply(1:(lengthFSdb - 1), function(i) {
          call_cytoscape_df(i)
        }, mc.cores = number_processing_threads))
        ##
        ########################################################################
        ##
        closeAllConnections()
        ##
      } else if (osType == "Windows") {
        ##
        clust <- makeCluster(number_processing_threads)
        registerDoParallel(clust)
        ##
        ########################################################################
        ##
        cytoscape_df <- foreach(i = 1:(lengthFSdb - 1), .combine = 'rbind', .verbose = FALSE) %dopar% {
          call_cytoscape_df(i)
        }
        ##
        ########################################################################
        ##
        stopCluster(clust)
      }
    }
    ##
    ############################################################################
    ############################################################################
    ##
    if (!is.null(cytoscape_df)) {
      ##
      cytoscape_df <- data.frame(cytoscape_df, stringsAsFactors = FALSE)
      cytoscape_df[, 1] <- as.numeric(cytoscape_df[, 1])
      cytoscape_df <- cytoscape_df[order(cytoscape_df[, 1], decreasing = TRUE), ]
      colnames(cytoscape_df) <- NULL
      rownames(cytoscape_df) <- NULL
      ##
      edge_dataFrame <- data.frame(edgeid = paste0(cytoscape_df[, 2], " (specsim) ", cytoscape_df[, 3]), spectral_entropy_similarity_score = cytoscape_df[, 1])
      colnames(edge_dataFrame) <- c("edgeid", "Spectral Entropy")
      rownames(edge_dataFrame) <- NULL
      ##
      correlation_network <- data.frame(cytoscape_df[, 2], rep("specsim", dim(edge_dataFrame)[1]), cytoscape_df[, 3], stringsAsFactors = FALSE)
      colnames(correlation_network) <- c("node1", "link", "node2")
      ##
      mspnodeid <- as.character(FSdb[["MSPLibraryParameters"]][[tolower(mspNodeID)]])
      xUnlinked <- which(!(mspnodeid %in% c(correlation_network[, 1], correlation_network[, 3])))
      t1df <- data.frame(mspnodeid[xUnlinked], rep("specsim", length(xUnlinked)), rep("", length(xUnlinked)), stringsAsFactors = FALSE)
      colnames(t1df) <-  c("node1", "link", "node2")
      ##
      correlation_network <- rbind(correlation_network, t1df)
      colnames(correlation_network) <- NULL
      rownames(correlation_network) <- NULL
      ##
      node_attributes_dataFrame <- do.call(cbind, lapply(variableVector , function(i) {FSdb[["MSPLibraryParameters"]][[tolower(i)]]}))
      ##
      rownames(node_attributes_dataFrame) <- NULL
      node_attributes_dataFrame <- unique(as.matrix(node_attributes_dataFrame)) # To remove redundant rows
      ##
      seqFSdb <- seq(1, lengthFSdb, 1)
      names(seqFSdb) <- mspnodeid
      #
      FragmentList <- do.call(c, lapply(node_attributes_dataFrame[, 1], function(i) {
        x <- seqFSdb[[i]]
        FL <- FSdb[["FragmentList"]][[x]]
        FLcolon <- paste0(round(FL[, 1], digits = 3), ":", round(FL[, 2]/max(FL[, 2]), digits = 2))
        paste0(FLcolon, collapse = ",")
      }))
      ##
      node_attributes_dataFrame <- cbind(node_attributes_dataFrame, FragmentList)
      node_attributes_dataFrame <- data.frame(node_attributes_dataFrame)
      rownames(node_attributes_dataFrame) <- NULL
      ##########################################################################
      colNameNodeAttributes <- c("nodeid", mspVariableVector, "FragmentList")
      ## NOTICE: Cytoscape doesn't allow column with a`Name` header
      xName <- which(colNameNodeAttributes == 'Name')
      if (length(xName) > 0) {
        colNameNodeAttributes[xName] <- 'MSP_Name'
      }
      ##
      colnames(node_attributes_dataFrame) <- colNameNodeAttributes
      ##########################################################################
      if (length(xIDfsdb) == 0) {
        FSdb[["MSPLibraryParameters"]][["idcytoscape"]] <- NULL
      }
      ##
      ##########################################################################
      listCytoscape <- list(node_attributes_dataFrame, edge_dataFrame, correlation_network, FSdb)
      names(listCytoscape) <- c("node_attributes_dataFrame", "edge_dataFrame", "correlation_network", "FSDB")
      ##
      if (RTcheck) { ## we find the unique MSP blocks using an exclusion approach
        ##
        intensityHeightVec <- as.numeric(FSdb[["MSPLibraryParameters"]][["precursor_intensity"]])
        if (length(intensityHeightVec) > 0) {
          orderINT <- order(intensityHeightVec, decreasing = FALSE)
          ##
          retentionTimeVec <- FSdb[["Retention Time"]]
          ##
          netdf <- correlation_network[, c(1, 3)]
          exclusionList <- rep(0, lengthFSdb)
          for (i in orderINT) {
            x_netdf <- which((netdf[, 1] == mspnodeid[i]) | (netdf[, 2] == mspnodeid[i]))
            if (length(x_netdf) > 1) {
              x_rt <- which(abs(retentionTimeVec[x_netdf] - retentionTimeVec[i]) <= RTtolerance)
              if (length(x_rt) > 0) {
                iExclusion <- c(i, x_netdf[x_rt])
                xMaxInt <- which.max(intensityHeightVec[iExclusion])
                iExclusion <- iExclusion[-xMaxInt]
                exclusionList[iExclusion] <- iExclusion
              }
            }
          }
          exclusionList <- exclusionList[exclusionList != 0]
          exclusionMSPnoideid <- mspnodeid[exclusionList]
          ##
          netdf <- netdf[!(netdf[, 1] %in% exclusionMSPnoideid), ]
          netdf <- netdf[!(netdf[, 2] %in% exclusionMSPnoideid), ]
          ##
          colnames(netdf) <- c("node1", "node2")
          t1df <- data.frame(node_attributes_dataFrame$nodeid, rep("", length(node_attributes_dataFrame$nodeid)), stringsAsFactors = FALSE)
          colnames(t1df) <-  c("node1", "node2")
          ##
          netdf_final <- rbind(netdf, t1df)
          netdf_final_sif <- paste0(netdf_final[,1],"\tspecsim\t",netdf_final[,2])
          listCytoscape$filteredNetworkSIF <- netdf_final_sif
          listCytoscape$exclusionMSPnoideid <- exclusionMSPnoideid # this list can be used to subset the FSDB to unique MSP blocks
        } else {
          FSA_logRecorder("NOTICE: exclusion list can not be generated because precursor intensity values were not detected in the .msp file!!!")
        }
      } 
    } else {
      FSA_logRecorder("No pairwise correlation was detected!")
      listCytoscape <- NULL
    }
  } else {
    FSA_logRecorder("No pairwise correlation was detected!")
    listCytoscape <- NULL
  }
  ##
  return(listCytoscape)
}