\name{ce10.blacklist}
\docType{data}
\alias{ce10.blacklist}
\title{A \code{\linkS4class{GRanges}} object representing ENCODE signal artefact
or "black list" regions.}
\description{These objects are made from the so-called black list or signal
artefact regions
defined by Anshul Kundaje Alan Boyle as part of the ENCODE and modENCODE projects.
They are regions which show a signal in essentially any ChIP or similar experiment.
For this reason it is useful to remove reads aligning to these regions before
carrying out any sort of functional genomics analysis.

These objects can be merged with \code{\linkS4class{GreyList}} objects to perform
grey list and black list filtering in one step.  They are included in the package
as a convenience, with the permission of the creators.
}

\usage{ce10.blacklist}
\value{\code{\linkS4class{GRanges}} objects defining the black list regions.}
\format{An S4 \code{\linkS4class{GRanges}} object.}
\source{ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the
human genome. Nature. 2012 Sep 6;489(7414):57-74. doi: 10.1038/nature11247.
https://sites.google.com/site/anshulkundaje/projects/blacklists}
\keyword{datasets}
