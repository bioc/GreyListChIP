\name{getKaryotype-methods}
\alias{getKaryotype}
\alias{getKaryotype,GreyList-method}
\title{Replace the karyotype of a \code{\linkS4class{GreyList}} object}
\description{
  Though a \code{\linkS4class{BSgenome}} object (or a karyotype file)  is
  supplied when the 
  \code{\linkS4class{GreyList}} object is created, it is conceivable that
  the user might want to replace it.  This method allows that.
}
\usage{
getKaryotype(obj, genome, tileSize=1024)
}
\arguments{
  \item{obj}{A \code{GreyList} object.}
  \item{genome}{A \code{\linkS4class{BSgenome}} object, from which to take
  the karyotype.}
  \item{tileSize}{The size in nucleotides of each tile.  Overlapping tiles
  will be generated, spaced at 1/2 the width of the tiles.}
}
\value{
  Returns the \code{\linkS4class{GreyList}} object with a new genome and
  tiling.
}
\author{Gord Brown}
\seealso{
  \code{\linkS4class{GreyList}}, \code{\linkS4class{BSgenome}}
}
\examples{
# Load a pre-built GreyList objct.
data(greyList)
library(BSgenome.Hsapiens.UCSC.hg19)

# Replace the karyotype, updating the genome tiling.
\dontrun{gl <- getKaryotype(greyList,BSgenome.Hsapiens.UCSC.hg19)}
}
\keyword{methods}
