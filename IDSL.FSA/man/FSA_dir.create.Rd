\name{FSA_dir.create}
\alias{FSA_dir.create}
\title{
FSA_dir.create
}
\description{
A module to create directories after removing the existing directory with the same name to prevent data interferences.
}
\usage{
FSA_dir.create(folder, allowedUnlink = FALSE)
}
\arguments{
  \item{folder}{
folder
}
  \item{allowedUnlink}{
allowedUnlink
}
}
\value{
when the original folder was deleted and recreated successfully, `TRUE` is returned by this function.
}