## This scope of function is to standardize precursor types from various sources
## to be consistent with the `ionization_pathway_deconvoluter` module of the
## IDSL.UFA and IDSL.SUFA packages. (Sadjad)
##
UFSA_precursorType_corrector <- function(precursorType, ionMode = NULL) {
  ##
  precursorType <- gsub("*", "", precursorType, fixed = TRUE)
  precursorType <- gsub("[.]", "", precursorType)
  precursorType <- gsub(" ", "", precursorType)
  precursorType <- gsub("Cat", "M+H", precursorType, ignore.case = TRUE)
  ##
  ##############################################################################
  ##############################################################################
  ## To adjust adduct types from https://fiehnlab.ucdavis.edu/staff/kind/metabolomics/ms-adduct-calculator
  precursorType <- gsub("ACN", "C2H3N", precursorType, ignore.case = TRUE)
  precursorType <- gsub("TFA", "C2HF3O2", precursorType, ignore.case = TRUE)
  precursorType <- gsub("FA", "CH2O2", precursorType, ignore.case = TRUE)
  precursorType <- gsub("formate", "CHO2", precursorType, ignore.case = TRUE)
  precursorType <- gsub("MeOH", "CH3OH", precursorType, ignore.case = TRUE)
  precursorType <- gsub("Hac", "C2H3O", precursorType, ignore.case = TRUE)
  precursorType <- gsub("IsoProp", "C3H8O", precursorType, ignore.case = TRUE)
  precursorType <- gsub("DMSO", "C2H6OS", precursorType, ignore.case = TRUE)
  ##
  precursorType <- gsub("[]]--$", "]2-", precursorType, ignore.case = TRUE)
  precursorType <- gsub("[]]1-$", "]-", precursorType, ignore.case = TRUE)
  ##
  precursorType <- gsub("[]][+][+]$", "]2+", precursorType, ignore.case = TRUE)
  precursorType <- gsub("[]]1[+]$", "]+", precursorType, ignore.case = TRUE)
  ##
  precursorType <- gsub("-[[:digit:]]i|-[[:digit:]]e|[+][[:digit:]]i", "", precursorType)
  precursorType <- gsub("-i|-e|[+]i", "", precursorType)
  ##
  ##############################################################################
  ##############################################################################
  ##
  x_blank <- which((precursorType == "") | (precursorType == "\t"))
  ##
  x_bracket2 <- which(!grepl("[[:punct:]]$", precursorType))
  if (length(x_bracket2) > 0) {
    precursorType[x_bracket2] <- paste0(precursorType[x_bracket2], "]")
  }
  ##
  x_bracket1 <- which(!grepl("^[[]", precursorType))
  if (length(x_bracket1) > 0) {
    precursorType[x_bracket1] <- paste0("[", precursorType[x_bracket1])
  }
  ##
  if (!is.null(ionMode)) {
    xPolePos <- which(!grepl("[+]$", precursorType) & grepl("^P", ionMode, ignore.case = TRUE))
    if (length(xPolePos) > 0) {
      precursorType[xPolePos] <- gsub("-$", "", precursorType[xPolePos], fixed = FALSE)
      precursorType[xPolePos] <- paste0(precursorType[xPolePos], "+")
    }
    ##
    xPoleNeg <- which(!grepl("-$", precursorType) & grepl("^N", ionMode, ignore.case = TRUE))
    if (length(xPoleNeg) > 0) {
      precursorType[xPoleNeg] <- gsub("[+]$", "", precursorType[xPoleNeg], fixed = FALSE)
      precursorType[xPoleNeg] <- paste0(precursorType[xPoleNeg], "-")
    }
  }
  ##
  precursorType[x_blank] <- ""
  ##
  ##############################################################################
  ##############################################################################
  ##
  LprecursorType <- length(precursorType)
  listIDprecursorType <- FSA_R.aggregate(precursorType)
  uniquePrecursorType <- names(listIDprecursorType)
  LuniquePrecursorType <- length(uniquePrecursorType)
  ##
  Elements <- UFSA_element_sorter()
  LElements <- 84 ## length(Elements) # Number of elements
  ##
  listDeconvolutedPathways <- UFSA_ionization_pathway_deconvoluter(uniquePrecursorType, Elements, LElements)
  ##
  coeffMvec <- do.call(c, lapply(listDeconvolutedPathways, function(i) {
    if (!is.null(i)) {
      i[[1]]
    } else {
      0
    }
  }))
  ##
  correctedPrecursorType <- rep("", LprecursorType)
  x0 <- which(coeffMvec == 0)
  if (length(x0) < LuniquePrecursorType) {
    ##
    deductMolVecMat <- do.call(rbind, lapply(listDeconvolutedPathways, function(i) {
      rep0Elements <- rep(0, LElements)
      if (!is.null(i)) {
        iVec <- i[[2]]
        ##
        xNeg <- which(iVec < 0)
        if (length(xNeg) > 0) {
          rep0Elements[xNeg] <- -iVec[xNeg]
        }
      }
      rep0Elements
    }))
    deductFormula <- UFSA_hill_molecular_formula_printer(deductMolVecMat, Elements, LElements)
    ##
    adductMolVecMat <- do.call(rbind, lapply(listDeconvolutedPathways, function(i) {
      rep0Elements <- rep(0, LElements)
      if (!is.null(i)) {
        iVec <- i[[2]]
        ##
        xPos <- which(iVec > 0)
        if (length(xPos) > 0) {
          rep0Elements[xPos] <- iVec[xPos]
        }
      }
      rep0Elements
    }))
    adductFormula <- UFSA_hill_molecular_formula_printer(adductMolVecMat, Elements, LElements)
    ##
    chargeAdduct <- do.call(c, lapply(listDeconvolutedPathways, function(i) {
      if (!is.null(i)) {
        i[[3]]
      } else {
        ""
      }
    }))
    ##
    correctedUniquePrecursorType <- do.call(c, lapply(1:LuniquePrecursorType, function(i) {
      if (coeffMvec[i] > 0) {
        if (coeffMvec[i] == 1) {
          coeff <- ""
        } else {
          coeff <- coeffMvec[i]
        }
        ##
        if (is.na(adductFormula[i])) {
          addForm <- ""
        } else {
          addForm <- paste0("+", adductFormula[i])
        }
        ##
        if (is.na(deductFormula[i])) {
          dedForm <- ""
        } else {
          dedForm <- paste0("-", deductFormula[i])
        }
        paste0("[", coeff, "M", dedForm, addForm, "]", chargeAdduct[i])
      } else {
        NA
      }
    }))
    ##
    for (i in 1:LuniquePrecursorType) {
      if (!is.na(correctedUniquePrecursorType[i])) {
        indexPrecursorType <- listIDprecursorType[[uniquePrecursorType[i]]]
        correctedPrecursorType[indexPrecursorType] <- correctedUniquePrecursorType[i]
      }
    }
  }
  ##
  return(correctedPrecursorType)
}