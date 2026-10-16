## Split the FANTOM5 promoter TPM table into one small .rds per chromosome and
## strand (fantom5_promoter_tpm_hg38/chr17_+.rds ...), so proMotif only loads
## the piece a gene needs. TPM stored as integer x 100 (0.01 TPM resolution).
## Usage: Rscript make_rds.R fantom5_promoter_tpm_hg38.tsv.gz out_dir
args <- commandArgs(trailingOnly = TRUE)
x <- read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE)
dir.create(args[2], showWarnings = FALSE, recursive = TRUE)
scols <- names(x)[-(1:6)]
samples <- data.frame(label = scols,
                      sample_type = sub(" \\[[^]]*\\]$", "", scols),
                      category = sub("^.* \\[([^]]*)\\]$", "\\1", scols),
                      stringsAsFactors = FALSE)
saveRDS(samples, file.path(args[2], "samples.rds"))
key <- paste0(x$chr, "_", x$strand)
for (k in unique(key)) {
  i <- which(key == k)
  m <- as.matrix(x[i, scols]); storage.mode(m) <- "double"
  m <- matrix(as.integer(round(m * 100)), nrow = length(i))
  saveRDS(list(peaks = x[i, 1:6], tpm100 = m), file.path(args[2], paste0(k, ".rds")),
          compress = "xz")
}
cat("files:", length(unique(key)) + 1, "\n")
