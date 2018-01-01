\name{setRegions-methods}
\alias{setRegions}
\alias{setRegions,GreyList-method}
\title{Replace the set of regions of a \code{\linkS4class{GreyList}} object}
\description{
  Though a \code{\linkS4class{BSgenome}} object (or a karyotype file, or a
  \code{\linkS4class{GRanges}})  is
  supplied when the 
  \code{\linkS4class{GreyList}} object is created, it is conceivable that
  the user might want to replace it.  This method allows that.
}
\usage{
setRegions(obj, regions, tileSize=1024)
}
\arguments{
  \item{obj}{A \code{GreyList} object.}
  \item{regions}{A \code{\linkS4class{GRanges}} object, from which to take
  the karyotype.}
  \item{tileSize}{The size in nucleotides of each tile.  Overlapping tiles
  will be generated, spaced at 1/2 the width of the tiles.}
}
\value{
  Returns the \code{\linkS4class{GreyList}} object with a new
  tiling.
}
\author{Gord Brown}
\seealso{
  \code{\linkS4class{GreyList}}, \code{\linkS4class{GRanges}}
}
\examples{
# Load a pre-built GreyList objct.
data(greyList)

# Replace the karyotype, updating the genome tiling.
regions=GRanges(seqnames=Rle(c('chr21','chr21')),ranges=IRanges(c(1,20000),end=c(10000,30000)))
gl <- setRegions(gl,regions)
}
\keyword{methods}
