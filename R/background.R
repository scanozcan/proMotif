## ===========================================================================
## 9. PROMOTER BACKGROUND FOR DISCOVERY MODE
## ---------------------------------------------------------------------------
## Discovery mode ranks TFs by how unusual this gene's best site is compared with
## the best site of the same motif in typical human promoters - not by the raw
## per-site p-value, which favours long / GC-rich motifs that score well in
## almost any promoter.
## Background (built once by data-raw/make_promoter_background.py): 2,000 random
## human protein-coding genes (Ensembl canonical TSS, -2500/+500 bp, hg38); for
## every JASPAR 2026 CORE vertebrate matrix (converted exactly as in
## .pfm_to_logodds) the best-site score in each promoter. Stored per matrix as
## upper-tail quantiles q_p (P(best site >= q_p) = p, dense in the tail) plus an
## exponential tail (peaks over threshold) used only beyond the strongest
## background promoter. bg_pvalue = P(a typical promoter's best site >= this
## gene's best site); its floor is ~1/2000 within the data, smaller only by the
## tail extrapolation. For other window sizes all quantiles are shifted by
## beta * log(L / L0) (extreme-value max-stability, beta from the moments).
## Human / JASPAR 2026 only; otherwise discovery falls back to FIMO-p ranking.
## ===========================================================================
.promoter_background <- function(release = "2026", species = "homo_sapiens",
                                 upstream = 2500, downstream = 500,
                                 custom_db = FALSE) {
  if (isTRUE(custom_db)) {
    message("  promoter background: not available for a user-supplied jaspar_db.")
    return(NULL)
  }
  if (!identical(species, "homo_sapiens")) {
    message("  promoter background: human only.")
    return(NULL)
  }
  fn   <- sprintf("promoter_background_JASPAR%s_hg38.tsv", release %||% "NA")
  cand <- c(system.file("extdata", fn, package = "proMotif"),
            file.path(getOption("promotif.bg_dir", "."), fn))
  f <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(f)) {
    message("  promoter background file ", fn, " not found.")
    return(NULL)
  }
  bg <- utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
  qcols <- grep("^q_", names(bg), value = TRUE)
  L0 <- 2500 + 500 + 1; L <- upstream + downstream + 1
  shift <- bg$beta * log(L / L0)
  Q <- as.matrix(bg[qcols]) + shift                 # rows = matrices
  attr(bg, "Q") <- Q
  attr(bg, "probs") <- as.numeric(sub("^q_", "", qcols))
  bg$mean <- bg$mean + shift
  message(sprintf("  promoter background: JASPAR %s, %d motifs x %d human promoters%s",
                  release, nrow(bg), bg$n_bg[1],
                  if (L != L0) sprintf(" (built for a 3,001-bp window; shifted for %s bp)",
                                       .fmt_bp(L)) else ""))
  bg
}

## P(best site of a typical promoter >= score) for scores of matrices bg[i, ].
.bg_pvalue <- function(score, bg, i) {
  Q <- attr(bg, "Q"); probs <- attr(bg, "probs")
  vapply(seq_along(score), function(k) {
    r <- i[k]
    if (is.na(r) || is.na(score[k])) return(NA_real_)
    q <- Q[r, ]; s <- score[k]; last <- length(q)
    if (s > q[last])                                  # beyond the data: exponential tail
      return(probs[last] * exp(-(s - q[last]) / bg$tail_sigma[r]))
    if (length(unique(q)) < 2) return(if (s <= q[1]) 1 else probs[last])
    stats::approx(q, probs, xout = s, ties = max, rule = 2)$y
  }, numeric(1))
}
