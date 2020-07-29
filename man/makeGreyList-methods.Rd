\name{makeGreyList-methods}
\alias{makeGreyList}
\alias{makeGreyList,GreyList-method}
\title{Generate a grey list from a \code{\linkS4class{GreyList}} object}
\description{
  Create the actual grey list, based on the threshold calculated by
  \code{calcThreshold} and the read counts from \code{countReads}.
}
\usage{
makeGreyList(obj,maxGap=16384)
}
\arguments{
  \item{obj}{The \code{\linkS4class{GreyList}} object to create the list
             for.}
  \item{maxGap}{If the distance between neighbouring grey regions is less than
                or equal to \code{maxGap}, the regions will be merged into
                one big region.}
}
\details{
  Create the grey list as a \code{\linkS4class{GRanges}} object.  Merge grey
  regions if they are separated by up to \code{maxGap} bp.
}
\value{The modified \code{\linkS4class{GreyList}}, with the regions added.}
\author{Gord Brown}
\seealso{
  \code{\linkS4class{GreyList}}, \code{\linkS4class{GRanges}}
}
\examples{
# load a pre-built GreyList object:
data(greyList)

# calculate the actual regions:
gl <- makeGreyList(greyList)
}
\keyword{methods}
