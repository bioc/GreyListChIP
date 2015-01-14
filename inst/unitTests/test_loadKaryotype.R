test_loadKaryo <- function() {
  path <- system.file("extra", package="GreyListChIP")
  fn <- file.path(path,"karyotype_chr21.txt")
  x <- new("GreyList",karyoFile=fn)
  checkEquals(length(x@karyotype),1)
  checkEquals(seqnames(x@karyotype)[1],"chr21")
}
