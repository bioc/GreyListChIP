\name{loadKaryotype-methods}
\alias{loadKaryotype}
\alias{loadKaryotype,GreyList-method}
\title{Load a karyotype from a file}
\description{
  Load a karyotype from a file, (re)generate the tiling, for a
  \code{\linkS4class{GreyList}} object.
}
\usage{
loadKaryotype(obj, karyoFile, tileSize=1024)
}
\arguments{
  \item{obj}{A \code{\linkS4class{GreyList}} object to set a new
             karyotype for.}
  \item{karyoFile}{A text file describing a genome's karyotype.  The format
                   is one line per chromosome (or contig or whatever), with
                   the name of the chromosome, some white space, and an integer
                   giving the length of the chromosome in nucleotides.}
  \item{tileSize}{The width of tiles on which to count.  Tiles will be placed
                  every tileSize/2 nucleotides, to catch regions of high
                  signal that might otherwise be split across (non-overlapping)
                  tiles and hence missed.}
}
\value{
  Returns the \code{\linkS4class{GreyList}} object with a new genome and
  tiling, loaded from the provided file.
}
\author{Gord Brown}
\seealso{
  \code{\linkS4class{GreyList}}
}
\examples{
# load a pre-built GreyList object:
data(greyList)

# Get a karyotype file:
path <- system.file("extra", package="GreyListChIP")
fn <- file.path(path,"karyotype_chr21.txt")

# Replace the karyotype in the GreyList:
gl <- loadKaryotype(greyList,fn)
}
\keyword{methods}
