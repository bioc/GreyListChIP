\name{greyListBS}
\alias{greyListBS}
\title{Construct a grey list with default arguments}
\description{
This function is a convenience function, wrapping several steps of grey list
construction into one step.  If you are content to accept the package's
defaults, and are using a \code{\linkS4class{BSgenome}} object to supply
the karyotype, this function might be of use to you.
}
\usage{
greyListBS(genome, bam)
}
\arguments{
  \item{genome}{a BSgenome object, for the relevant genome.}
  \item{bam}{a BAM file to use for making the grey list.}
}
\value{
An object of class \code{GreyList}.
}
\author{
Gord Brown
}
\seealso{
  \code{\link{GreyList}}
}
\examples{
# If you want to accept the defaults for everything, you can create the
# GreyList in one step using a BSgenome object:
library(BSgenome.Hsapiens.UCSC.hg19)
path <- system.file("extra", package="GreyListChIP")
\dontrun{fn <- file.path(path,"sample_chr21.bam")}
\dontrun{gl <- greyListBS(BSgenome.Hsapiens.UCSC.hg19,fn)}
}
\keyword{grey list}
