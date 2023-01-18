FSA_dir.create <- function(folder, allowedUnlink = FALSE) {
  dirCreationRCheck <- FALSE
  dirCreationCheck <- dir.exists(folder)
  if (allowedUnlink) {
    if (dirCreationCheck) {
      tryCatch({
        unlink(folder, recursive = TRUE)
        dirCreationCheck <- FALSE
      },
      warning = function(w) {
        dirCreationCheck <- TRUE
        dirCreationRCheck <- FALSE
        FSA_logRecorder(paste0("Can't delete `", folder, "`!"))
      },
      error = function(e) {
        dirCreationCheck <- TRUE
        dirCreationRCheck <- FALSE
        FSA_logRecorder(paste0("Can't delete `", folder, "`!"))
      })
    }
  }
  if (!dirCreationCheck) {
    tryCatch({
      dir.create(folder, recursive = TRUE)
      dirCreationRCheck <- TRUE
    },
    warning = function(w) {
      dirCreationRCheck <- FALSE
      FSA_logRecorder(paste0("Can't create `", folder, "`!"))
    },
    error = function(e) {
      dirCreationRCheck <- FALSE
      FSA_logRecorder(paste0("Can't create `", folder, "`!"))
    })
  }
  return(dirCreationRCheck)
}