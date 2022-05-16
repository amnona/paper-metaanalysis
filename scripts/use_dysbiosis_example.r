library(reticulate)
ca <- import("calour")
expc <- ca$read_qiime2(biompath, metadatapath, min_reads=2000, normalize=10000)
source_python("UniDI.py")
nsf <- py$import_nsf()
dbi <- py$dbi_ranks(expc, nsf)