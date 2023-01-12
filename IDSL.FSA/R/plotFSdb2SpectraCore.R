plotFSdb2SpectraCore <- function(FSdb, index) {
  ##
  FragmentList <- FSdb[["FragmentList"]][[index]]
  ##
  minFragmentMass <- min(FragmentList[, 1])
  maxFragmentMass <- max(FragmentList[, 1])
  maxFragmentInt <- max(FragmentList[, 2])
  nudge_y = maxFragmentInt/80
  ##
  plotlabelsMZ <- FSA_annotation_text_repel(FSAspectra = FragmentList, nGridX = 12, nGridY = 8)
  ##
  plot(FragmentList[, 1], FragmentList[, 2], type = "h", xlim = c((minFragmentMass - 10), (maxFragmentMass + 10)),
       ylim = c(0, maxFragmentInt*1.1), xaxs = "i", yaxs = "i", lwd = 2, cex = 4, yaxt = "n", xlab = "", ylab = "")
  text(x = FragmentList[, 1], y = (FragmentList[, 2] + nudge_y), label = plotlabelsMZ, col = "red")
  mtext("m/z", side = 1, adj = 0.5, line = 2, cex = 1.35)
  if (!is.null(FSdb[["MSPLibraryParameters"]][["filename"]])) {
    mtext(FSdb[["MSPLibraryParameters"]][["filename"]][index], side = 3, adj = 0, line = 0.25, cex = 1.15)
  } else {
    mtext(FSdb[["MSPLibraryParameters"]][["MSPfilename"]][index], side = 3, adj = 0, line = 0.25, cex = 1.15)
  }
  mtext(FSdb[["MSPLibraryParameters"]][["name"]][index], side = 3, adj = 1, line = 1.25, cex = 1.00)
  mtext(paste0("S = ", round(FSdb[["Spectral Entropy"]][index], 5)), side = 3, adj = 1, line = 0.20, cex = 0.90)
}