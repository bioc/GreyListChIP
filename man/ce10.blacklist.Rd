\name{ce10.blacklist}
\docType{data}
\alias{ce10.blacklist}
\alias{ce11.blacklist}
\alias{dm3.blacklist}
\alias{dm6.blacklist}
\alias{grch37.blacklist}
\alias{hg19.blacklist}
\alias{grch38.blacklist}
\alias{hg38.blacklist}
\alias{mm9.blacklist}
\alias{mm10.blacklist}
\title{A \code{\linkS4class{GRanges}} object representing ENCODE signal artefact
or "black list" regions.}
\description{These objects are made from the so-called black list or signal
artefact regions defined by Anshul Kundaje and Alan Boyle as part of the
ENCODE and modENCODE projects.  They are regions which show a signal in
essentially any ChIP or similar experiment.  For this reason it is useful to
remove reads aligning to these regions before carrying out any sort of
functional genomics analysis.

These objects can be merged with \code{\linkS4class{GreyList}} objects to
perform grey list and black list filtering in one step.  They are included
in the package as a convenience, with the permission of the creators.

The package includes black lists for worm (ce10, ce11), fly (dm3, dm6), 
human (hg19, grch37, hg38, grch38), and mouse (mm9, mm10).  (Note that
hg19 and grch37 differ only in chromosome naming, and likewise for hg38 and
grch38.)
}

\usage{data(ce10.blacklist,package=GreyListChIP)
data(ce11.blacklist,package=GreyListChIP)
data(dm3.blacklist,package=GreyListChIP)
data(dm6.blacklist,package=GreyListChIP)
data(grch37.blacklist,package=GreyListChIP)
data(hg19.blacklist,package=GreyListChIP)
data(grch38.blacklist,package=GreyListChIP)
data(hg38.blacklist,package=GreyListChIP)
data(mm9.blacklist,package=GreyListChIP)
data(mm10.blacklist,package=GreyListChIP)
}
\value{A set of intervals defining the black list regions
for the named genome.}
\format{An S4 \code{\linkS4class{GRanges}} object.}
\source{https://sites.google.com/site/anshulkundaje/projects/blacklists}
\references{Amemiya HM, Kundaje A, Boyle AP. The ENCODE blacklist:
identification of problematic regions of the genome.
Sci Rep. 2019 Dec; 9(1) 9354 DOI: 10.1038/s41598-019-45839-z. 

ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the
human genome. Nature. 2012 Sep 6;489(7414):57-74. doi: 10.1038/nature11247.}
\keyword{datasets}
