setClass("GreyList",
         representation(genome="BSgenome",
                        karyotype="Seqinfo",
                        karyo_file="character",
                        genomeRegions="GRanges",
                        tiles="GRanges",
                        counts="numeric",
                        files="character",
                        size_param="numeric",
                        size_stderr="numeric",
                        size_mean="numeric",
                        mu_param="numeric",
                        mu_stderr="numeric",
                        mu_mean="numeric",
                        reps="numeric",
                        sample_size="numeric",
                        pvalue="numeric",
                        threshold="numeric",
                        max_gap="numeric",
                        regions="GRanges",
                        coverage="numeric"),
         prototype(genome=new("BSgenome"),
                   karyotype=new("Seqinfo"),
                   karyo_file=NA_character_,
                   genomeRegions=new("GRanges"),
                   tiles=new("GRanges"),
                   counts=NA_integer_,
                   files=vector(mode="character"),
                   size_param=NA_real_,
                   size_stderr=NA_real_,
                   size_mean=NA_real_,
                   mu_param=NA_real_,
                   mu_stderr=NA_real_,
                   mu_mean=NA_real_,
                   reps=NA_integer_,
                   sample_size=NA_integer_,
                   pvalue=NA_real_,
                   threshold=NA_integer_,
                   max_gap=NA_integer_,
                   regions=new("GRanges"),
                   coverage=NA_real_),
         validity=function(object) { return(TRUE) }
)

initialize.GreyList <- function(.Object,
                                genome=NA,
                                karyotype=NA,
                                karyoFile=NA,
                                regions=NA, ...) {
  if (!missing(genome) && isS4(genome)) {
    .Object <- getKaryotype(.Object, genome)
  } else if (!missing(karyotype)) {
    .Object <- setKaryotype(.Object, karyotype)
  } else if (!is.na(karyoFile)) {
    .Object <- loadKaryotype(.Object, karyoFile)
  } else if (!missing(regions)) {
    .Object <- setRegions(.Object,regions)
  } 
  callNextMethod(.Object, ...)
}

loadKaryotype.GreyList <- function(obj,karyoFile,tileSize) {
  tbl <- read.table(karyoFile,header=FALSE,stringsAsFactors=FALSE)
  colnames(tbl) <- c("Chrom","Length")
  kInfo <- Seqinfo(tbl$Chrom,
                   seqlengths=tbl$Length,
                   isCircular=rep(FALSE,nrow(tbl)),
                   genome=NA)
  obj@karyo_file <- karyoFile
  obj@karyotype <- kInfo
  x <- tileGenome(obj@karyotype,tilewidth=tileSize/2)
  tiles_half <- unlist(x)
  obj@tiles <- suppressWarnings(resize(tiles_half,tileSize))
  obj@tiles <- trim(obj@tiles)
  return(obj)
}

setKaryotype.GreyList <- function(obj,karyo,tileSize) {
  obj@karyotype <- karyo
  x <- tileGenome(obj@karyotype,tilewidth=tileSize/2)
  tiles_half <- unlist(x)
  obj@tiles <- suppressWarnings(resize(tiles_half,tileSize))
  obj@tiles <- trim(obj@tiles)
  return(obj)
}

getKaryotype.GreyList <- function(obj,genome,tileSize) {
  if (missing(genome)) {
    stop("Please supply a BSGenome object ('genome=...').")
  } else {
    obj@genome <- genome
    obj@karyotype <- seqinfo(genome)
    x <- tileGenome(obj@karyotype,tilewidth=tileSize/2)
    tiles_half <- unlist(x)
    obj@tiles <- suppressWarnings(resize(tiles_half,tileSize))
    obj@tiles <- trim(obj@tiles)
  }
  return(obj)
}

#tileRegions <- function(regions,tileSize) {
#  tiles <- slidingWindows(regions,width=tileSize,step=tileSize/2)
#  obj@tiles = do.call("c",tiles)
#}

setRegions.GreyList <- function(obj,regions,tileSize) {
  obj@genomeRegions <- regions
  x <- slidingWindows(regions,width=tileSize,step=tileSize/2)
  if (length(x) > 1) {
#    obj@tiles <- do.call("c",x)
    obj@tiles <- unlist(x)
  } else {
    obj@tiles <- x
  }
  obj@tiles <- trim(obj@tiles)
  return(obj)
}

countReads.GreyList <- function(obj,bamFile,yieldSize=NA_integer_) {
  fd <- BamFile(bamFile,yieldSize=yieldSize)
  counts <- summarizeOverlaps(obj@tiles,fd,inter.feature=FALSE,
                              ignore.strand=TRUE,
                              mode="IntersectionStrict")
  obj@counts <- assays(counts)[[1]][,1]
  obj@files <- c(obj@files,bamFile)
  return(obj)
}

fitDist <- function(n,x,size) {
  s <- sample(x,size=size)
  p <- tryCatch(
    suppressWarnings(fitdistr(s,"negative binomial")),
    error=function(x) NA
  )
  return(p)
}

estimateParams <- function(obj,reps,size,cores) {
  stuff <- mclapply(seq_len(reps),fitDist,obj@counts,size,mc.cores=cores)
  stuff <- stuff[!is.na(stuff)]
  reps <- length(stuff)
  obj@size_param <- vapply(seq_len(reps),function(x) stuff[[x]]$estimate['size'],numeric(1))
  obj@size_stderr <- vapply(seq_len(reps),function(x) stuff[[x]]$sd['size'],numeric(1))
  obj@mu_param <- vapply(seq_len(reps),function(x) stuff[[x]]$estimate['mu'],numeric(1))
  obj@mu_stderr <- vapply(seq_len(reps),function(x) stuff[[x]]$sd['mu'],numeric(1))
  obj@size_mean <- mean(obj@size_param)
  obj@mu_mean <- mean(obj@mu_param)
  return(obj)
}

calcThreshold.GreyList <- function(obj,reps,sampleSize,p,cores) {
  obj@reps <- reps
  obj@sample_size <- sampleSize
  obj@pvalue <- p
  obj <- estimateParams(obj,reps,sampleSize,cores)
  obj@threshold <- qnbinom(p,size=obj@size_mean,mu=obj@mu_mean)
  return(obj)
}

makeGreyList.GreyList <- function(obj,maxGap) {
  obj@max_gap <- maxGap
  regions <- obj@tiles[obj@counts > obj@threshold]
  obj@regions <- reduce(regions,min.gapwidth=maxGap)
  bp <- sum(width(ranges(obj@regions)))
  total <- sum(as.numeric(seqlengths(obj@karyotype)))
  obj@coverage <- (bp / total) * 100
  cat(sprintf("coverage: %.0f bp (%.2f%%)\n",bp,obj@coverage))
  return(obj)
}

export.GreyList <- function(object, con, ...) {
  ir <- ranges(object@regions)
  names <- seqnames(object@regions)
  df <- data.frame(as.character(names), start(ir), end(ir))
  write.table(df,file=con,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}

show.GreyList <- function(object) {
  words <- character(10)
  if (is.na(object@karyo_file)) {
    if (length(organism(object@genome)) == 0) {
      sinf <- object@karyotype
      g <- genome(sinf)
      if (length(g) == 0) {
        words[1] <- "GreyList on unknown genome\n"
      } else {
        words[1] <- sprintf("GreyList on %s\n",g[1])
      }
    } else {
      words[1] <- sprintf("GreyList on %s (%s %s)\n",
                          organism(object@genome),
                          provider(object@genome),
                          providerVersion(object@genome))
    }
  } else {
    words[1] <- sprintf("GreyList on karyotype file %s\n",
                        basename(object@karyo_file))
  }
  current <- 2
  if (length(object@tiles) > 0) {
    words[current] <- sprintf("  tiles: %s\n",length(object@tiles))
    current <- current + 1
  }
  if (length(object@files) > 0) {
    words[current] <- sprintf("  files: %s\n",
                              paste0(object@files,collapse=", "))
    current <- current + 1
  }
  if (!is.na(object@size_mean)) {
    words[current] <- sprintf("  size (mean): %s\n",object@size_mean)
    current <- current + 1
  }
  if (!is.na(object@mu_mean)) {
    words[current] <- sprintf("  mu (mean): %s\n",object@mu_mean)
    current <- current + 1
  }
  if (!is.na(object@reps)) {
    words[current] <-
                  sprintf("  params: reps=%d, sample size=%d, p-value=%.2f\n",
                          object@reps,object@sample_size,object@pvalue)
    current <- current + 1
  }
  if (!is.na(object@threshold)) {
    words[current] <- sprintf("  threshold: %s\n",object@threshold)
    current <- current + 1
  }
  if (length(object@regions) > 0) {
    words[current] <- sprintf("  regions: %d\n",length(object@regions))
    current <- current + 1
  }
  if (!is.na(object@coverage)) {
    words[current] <- sprintf("  coverage: %.2f%%\n",object@coverage)
  }
  cat(paste0(words,collapse=""))
}
  
greyListBS <- function(genome,bam) {
  gl <- new("GreyList",genome)
  gl <- countReads(gl,bam)
  gl <- calcThreshold(gl)
  gl <- makeGreyList(gl)
  return(gl)
}

setMethod("initialize","GreyList",initialize.GreyList)

setGeneric("getKaryotype",def=function(obj,genome,tileSize=1024) {
  standardGeneric("getKaryotype")})
setMethod(getKaryotype,"GreyList",getKaryotype.GreyList)

setGeneric("setKaryotype",def=function(obj,karyo,tileSize=1024) {
  standardGeneric("setKaryotype")})
setMethod(setKaryotype,"GreyList",setKaryotype.GreyList)

setGeneric("loadKaryotype",def=function(obj,karyoFile,tileSize=1024) {
  standardGeneric("loadKaryotype")})
setMethod(loadKaryotype,"GreyList",loadKaryotype.GreyList)

setGeneric("setRegions",def=function(obj,regions,tileSize=1024) {
  standardGeneric("setRegions")})
setMethod(setRegions,"GreyList",setRegions.GreyList)

setGeneric("countReads",def=function(obj,bamFile,yieldSize=NA_integer_) {
  standardGeneric("countReads")})
setMethod(countReads,"GreyList",countReads.GreyList)
setGeneric("calcThreshold",
           def=function(obj,reps=100,sampleSize=30000,p=0.99,cores=1) {
                        standardGeneric("calcThreshold")}
)
setMethod(calcThreshold,"GreyList",calcThreshold.GreyList)
setGeneric("makeGreyList",def=function(obj,maxGap=16384) {
  standardGeneric("makeGreyList")})
setMethod(makeGreyList,"GreyList",makeGreyList.GreyList)
setMethod(export,c("GreyList","character","missing"),export.GreyList)
setMethod(show,"GreyList",show.GreyList)
