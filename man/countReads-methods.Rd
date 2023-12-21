\name{countReads-methods}
\alias{countReads}
\alias{countReads,GreyList-method}
\title{Count reads from a \code{\linkS4class{BamFile}}}
\description{
  Given the tiling of the genome created when the \code{\linkS4class{GreyList}}
  object was created (or replaced via \code{getKaryotype}), count reads
  overlapping the bins, in preparation for estimating the threshold for
  grey-listing bins.
}
\usage{
countReads(obj,bamFile,yieldSize=NA_integer_)
}
\arguments{
  \item{obj}{A \code{GreyList} object on which to count reads.}
  \item{bamFile}{A \code{\linkS4class{BamFile}} from which to count reads.}
  \item{yieldSize}{Number of records to yield each time the BAM file is read.}
}
\details{
  This method counts reads contained within the bins that make up the genome
  tiling.  Bins are overlapping (by default 1Kb bins at 512b intervals) so
  reads are counted once for each bin that wholly contains them.

  Setting the \code{yieldSize} can help control the memory usage. If unset, a
  value of 1,000,000 is used by the \code{summarizeOverlaps} function
  from the \code{GenomicAlignments} package that this method calls. Setting the
  \code{yieldSize} to a lower value will reduce the memory requirement at the
  expense of longer run times.
}
\value{The modified \code{GreyList} object, with added counts.}
\author{Gord Brown}
\seealso{
  \linkS4class{GreyList}, \linkS4class{BamFile}
}
\examples{
# Load a pre-built GreyList object.
data(greyList)

path <- system.file("extra", package="GreyListChIP")
\dontrun{fn <- file.path(path,"sample_chr21.bam")}
\dontrun{gl <- countReads(greyList,fn)}
}
\keyword{methods}
