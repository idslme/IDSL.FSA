UFSA_hill_molecular_formula_printer <- function(MolVecMat, Elements, LElements = length(Elements)) {
  HillElements <- sort(Elements, decreasing = FALSE)
  x_c <- which(HillElements == "C")
  x_h <- which(HillElements == "H")
  if ((length(x_c) == 1) & (length(x_h) == 1)) {
    x_c_h <- c(x_c, x_h)
    HillElements <- HillElements[-x_c_h]
    HillElements <- c("C", "H", HillElements)
  } else if ((length(x_c) == 1) & (length(x_h) == 0)) {
    HillElements <- HillElements[-x_c]
    HillElements <- c("C", HillElements)
  } else if ((length(x_c) == 0) & (length(x_h) == 1)) {
    HillElements <- HillElements[-x_h]
    HillElements <- c(HillElements, "H")
  }
  ##
  ##############################################################################
  ##
  orderHillElements <- do.call('c', lapply(HillElements, function(i) {
    which(Elements == i)
  }))
  ##
  ##############################################################################
  ##
  molecular_formula_printer <- function(orderHillElements, Elements, MolVec) {
    ##
    strMolecularFormula <- ""
    ##
    for (i in orderHillElements) {
      if (MolVec[i] > 0) {
        if (MolVec[i] == 1) {
          strMolecularFormula <- paste0(strMolecularFormula, Elements[i])
        } else {
          strMolecularFormula <- paste0(strMolecularFormula, Elements[i], MolVec[i])
        }
      }
    }
    ##
    if (strMolecularFormula == "") {
      strMolecularFormula <- NA
    }
    return(strMolecularFormula)
  }
  ##
  ##############################################################################
  ##
  MolVecMat <- matrix(MolVecMat, ncol = LElements)
  ##
  MolFormList <- do.call('c', lapply(1:dim(MolVecMat)[1], function(i) {
    molecular_formula_printer(orderHillElements, Elements, MolVecMat[i, ])
  }))
  ##
  return(MolFormList)
}