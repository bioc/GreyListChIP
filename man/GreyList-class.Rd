\name{GreyList-class}
\Rdversion{1.1}
\docType{class}
\alias{GreyList-class}
\alias{GreyList}
\alias{export,GreyList,character,missing-method}
\alias{show,GreyList-method}
\alias{initialize,GreyList,BSgenome-method}

\title{\code{"GreyList"} Objects}
\description{
Regions of high signal in the input samples of a ChIP experiment can lead to
artefacts in peak calling.  This class generates "grey lists" of such regions,
for use in filtering reads before peak calling (or filtering peaks after peak
calling, though it is generally safer to filter first).
}
\section{Objects from the Class}{
Objects can be created by calls of the
form \code{new("GreyList", genome, ...)},
where \code{genome} is a \code{"BSgenome"} object describing a genome, such as
\code{BSgenome.Hsapiens.UCSC.hg19}.  Alternatively, a karyotype
file can be provided explictly: \code{new("GreyList", karyoFile=fn, ...)}.
Either \code{genome} or \code{karyoFile} must be provided; if both are present,
the \code{\linkS4class{BSgenome}} object takes precedence.  Alternatively,
an explicit list of regions may be provided as a \code{\linkS4class{GRanges}}
object.
}
\section{Slots}{
  \describe{
    \item{\code{genome}:}{The \code{\linkS4class{BSgenome}} object
          corresponding with the genome the reads are aligned to}
    \item{\code{karyotype}:}{The \code{\linkS4class{Seqinfo}} object from
          the \code{\linkS4class{BSgenome}} object, or made from the
          \code{karyo_file}}
    \item{\code{karyo_file}:}{The name of a file containing chromosome sizes
          for the reference genome of interest, one per line, as
          "chromName chromLength" pairs.}
    \item{\code{genomicRegions}:}{a \code{\linkS4class{GRanges}} object
          that defines which regions of the genome should be used to build
          the grey list.  This is to allow the list to be built on just part
          of the genome.}
    \item{\code{tiles}:}{A \code{\linkS4class{GRanges}} object with an
          overlapping tiling of the genome (by default 1Kb tiles every 512b). }
    \item{\code{counts}:}{A numeric vector holding the counts corresponding to
          the tiling and the BAM file provided.}
    \item{\code{files}:}{A vector of BAM filenames that were used to generate
          the counts (currently only accepts one).}
    \item{\code{size_param}:}{The computed estimates of the "size" parameter
          of the negative binomial distribution, estimated by
          \code{MASS::fitdistr} from repeated sampling from the counts.}
    \item{\code{size_stderr}:}{The standard errors of the "size" parameters, as
          estimated by \code{MASS::fitdistr}.}
    \item{\code{size_mean}:}{The mean of the "size" estimates.}
    \item{\code{mu_param}:}{Computed estimates of the "mu" parameter of the
          negative binomial distribution, estimated by \code{MASS::fitdistr}
          from repeated sampling from the counts.}
    \item{\code{mu_stderr}:}{The standard errors of the "mu" parameter.}
    \item{\code{mu_mean}:}{The mean of the "mu" estimates.}
    \item{\code{reps}:}{How many samples from the counts were taken.}
    \item{\code{sample_size}:}{How many values were sampled
          from the counts, for each estimate of "size" and "mu".}
    \item{\code{pvalue}:}{The requested p-value threshold.}
    \item{\code{threshold}:}{The calculated threshold, based on the p-value.}
    \item{\code{max_gap}:}{The largest gap to consider when merging nearby
          regions (i.e. if there are "grey" regions up to this many nucleotides
          apart, merge them into one long region).}
    \item{\code{regions}:}{A \code{\linkS4class{GRanges}} object defining the
          final grey list regions.}
    \item{\code{coverage}:}{The percentage of the genome covered by the grey
          list regions.}
  }
}
\section{Methods}{
  \describe{
    \item{calcThreshold}{\code{signature(obj = "GreyList")}: Calculate the
          cutoff for reads in bins, based on fitting the counts to a negative
          binomial distribution.}
    \item{countReads}{\code{signature(obj = "GreyList")}: Count reads in bins
          across the genome.}
    \item{export}{\code{signature(object = "GreyList",
                  con = "character", format = "missing")}: Write the grey list
          to a file.}
    \item{initialize}{\code{signature(.Object = "GreyList")}: Create an initial
          object (invoked automatically by \code{new("GreyList",...)}).}
    \item{loadKaryotype}{\code{signature(obj = "GreyList")}: Load a genome
          description from a file.  The file format is one line per chromosome,
          with the name of the chromosome followed by white space followed by
          an integer indicating the length of the chromosome.}
    \item{getKaryotype}{\code{signature(obj = "GreyList")}: Get the karyotype
          of a genome from a \code{\linkS4class{BSgenome}} object.}
    \item{setRegions}{\code{signature(obj = "GreyList")}: Set the region(s)
          of a genome to use in making the \code{\linkS4class{GreyList}} object.}
    \item{makeGreyList}{\code{signature(obj = "GreyList")}: Compute the actual
          grey list, after calculating the threshold.}
    \item{show}{\code{signature(object = "GreyList")}: Display the grey list.}
  }
}
\author{
Gord Brown (\email{gdbzork@gmail.com})
}
\seealso{
  \code{\linkS4class{BSgenome}},
  \code{\linkS4class{Seqinfo}},
  \code{\linkS4class{GRanges}}
}
\examples{
showClass("GreyList")

# Load a karyotype file:
path <- system.file("extra", package="GreyListChIP")
fn <- file.path(path,"karyotype_chr21.txt")

# Create a GreyList object:
gl <- new("GreyList",karyoFile=fn)
}
\keyword{classes}
