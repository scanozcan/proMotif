## Builds inst/extdata/chipatlas_pulldown_{hg38,mm10}.tsv.gz: the ChIP-Atlas
## experiments that pulled down an epitope-tagged TF (FLAG, HA, V5, GFP, HALO,
## biotin ...), used by chip_pulldown. Reads ChIP-Atlas experimentList.tab
## (about 400 MB, streamed). Run from the package folder.
devtools::load_all(".")
for (g in c("hg38", "mm10"))
  build_chip_pulldown_table(genome = g, out = sprintf("inst/extdata/chipatlas_pulldown_%s.tsv.gz", g))
