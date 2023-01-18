FSA_locate_regex <- function(string, pattern, ignore.case = FALSE, perl = FALSE,
                             fixed = FALSE, useBytes = FALSE) {
  ##
  FSAgreg <- gregexpr(pattern, string, ignore.case, perl, fixed, useBytes)[[1]]
  if (FSAgreg[1] > 0) {
    FSAgreg_lengthchar <- attributes(FSAgreg)$match.length
    #
    loc_mat <- do.call(rbind, lapply(1:length(FSAgreg_lengthchar), function(i) {
      c(FSAgreg[i], (FSAgreg[i] + FSAgreg_lengthchar[i] - 1))
    }))
  } else {
    loc_mat <- NULL
  }
  #
  return(loc_mat)
}