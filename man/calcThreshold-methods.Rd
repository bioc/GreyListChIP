\name{calcThreshold-methods}
\alias{calcThreshold}
\alias{calcThreshold,GreyList-method}
\title{Calculate Read Count Threshold}
\description{
  Based on the counts from \code{countReads}, sample counts from the set
  several times, estimate the parameters of the negative binomial distribution
  for each sample, then calculate the mean of the parameters (\emph{size} and
  \emph{mu}).  Use these values to calculate the read count threshold, given
  the specified p-value threshold.
}
\usage{
calcThreshold(obj, reps=100, sampleSize=30000, p=0.99,cores=1)
}
\arguments{
  \item{obj}{A \code{GreyList} object for which to calculate
             the threshold.}
  \item{reps}{The number of times to sample bins and estimate the parameters
              of the negative binomial distribution.}
  \item{sampleSize}{The number of bins to sample on each repetition.}
  \item{p}{The p-value threshold for marking bins as ``grey''.}
  \item{cores}{The number of CPU cores (parallel threads) to use when
               sampling repeatedly from the set of counts}
}
\details{

  This method samples from the set of counts generated during the
  \code{countReads} step.  Each sample is fitted to the negative binomial
  distribution, and the parameters estimated.  The means of the \code{mu}
  and \code{size} parameters is calculated, then used to choose a read count
  threshold, given the p-value cutoff provided.  If \code{cores} is given,
  the process will use that many cores to parallelize the parameter estimation.
}
\value{The modified \code{GreyList} object, with the threshold added.}
\author{Gord Brown}
\examples{
# Load a pre-built R object with counts.
data(greyList)

# Calculate the threshold:
gl <- calcThreshold(gl,reps=10,sampleSize=1000,p=0.99,cores=1)
}
\references{
Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S.
Fourth edition. Springer.
}
\keyword{methods}
