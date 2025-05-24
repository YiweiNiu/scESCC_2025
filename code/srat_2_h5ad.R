#!/usr/bin/env Rscript

# usage: Rscript srat_2_h5ad.R xxx.rds
# output: xxx.h5ad

require(Seurat)

# args
args = commandArgs(trailingOnly=TRUE)

# test args
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
} else {
  fin = args[1]
  fout = stringr::str_replace(fin, '.rds', '.h5ad')
}

srat = readRDS(fin)
sceasy::convertFormat(srat, from="seurat", to="anndata",
                      main_layer = 'counts',
                      outFile=fout)

