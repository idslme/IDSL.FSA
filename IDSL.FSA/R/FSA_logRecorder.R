FSA_logRecorder <- function(messageQuote, allowedPrinting = TRUE) {
  if (allowedPrinting) {
    cat(paste0("\n", messageQuote))
  }
  ##
  if (exists('.logFSA')) {
    if (typeof(messageQuote) == "list") {
      namesMessageQuote <- names(messageQuote)
      if (length(messageQuote) > 0) {
        for (i in 1:length(messageQuote)) {
          write(paste0(i, ": In ", paste0(deparse(messageQuote[[i]]), collapse = "\n")), file = .logFSA, append = TRUE, sep = "\n")
          write(paste0("  ", namesMessageQuote[i]), file = .logFSA, append = TRUE, sep = "\n")
        }
      }
      ##
    } else {
      write(messageQuote, file = .logFSA, append = TRUE, sep = "\n")
    }
  }
}
