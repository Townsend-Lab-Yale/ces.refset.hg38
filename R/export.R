# Adapted from BSgenome.Hsapiens.UCSC.hg19's zzz.R
.onLoad <- function(libname, pkgname)
{
  data_path <- system.file("refset", package=pkgname, lib.loc=libname, mustWork=TRUE)
  refset = cancereffectsizeR:::preload_ref_data(data_path)
  ns <- asNamespace(pkgname)
  objname <- pkgname
  assign(objname, refset, envir=ns)
  namespaceExport(ns, objname)
}

