FSA_R.aggregate <- function(FSAvec) {
  ##
  listIDFSAvec <- base::tapply(seq(1, length(FSAvec), 1), FSAvec, FUN = 'c', simplify = FALSE)
  ##
  return(listIDFSAvec)
}