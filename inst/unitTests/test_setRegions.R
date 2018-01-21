test_setRegions <- function() {
  path <- system.file("extra", package="GreyListChIP")
  fn <- file.path(path,"karyotype_chr21.txt")
  regions=GRanges(seqnames=Rle(c('chr21','chr21')),ranges=IRanges(c(1,20000),end=c(10000,30000)))
  x <- new("GreyList",regions=regions)
  checkEquals(length(x@tiles),38)
}
