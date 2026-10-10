#' proMotif: promoter transcription-factor binding-site scanner
#'
#' For a gene symbol and a TF, proMotif enumerates the gene's transcripts from
#' Ensembl, fetches each promoter window, scans both strands for JASPAR and
#' HOCOMOCO motif matches, and reports predicted sites in TSS-relative and
#' absolute genomic coordinates, with per-transcript and genomic plots. A Shiny
#' front-end is available via \code{\link{launch_app}}.
#'
#' Key functions: \code{\link{scan_binding_sites}} (data),
#' \code{\link{plot_binding_sites}} (plots), \code{\link{run_binding_site_analysis}}
#' (scan + plot + write files), \code{\link{launch_app}} (interactive UI).
#'
#' @keywords internal
#' @import ggplot2
#' @importFrom stats embed setNames
#' @importFrom utils write.csv head
"_PACKAGE"

## ggplot2 aesthetics reference data-frame columns by bare name; declare them
## so R CMD check does not flag "no visible binding for global variable".
utils::globalVariables(c("position", "score_frac", "database", "strand",
                         "genomic_position", "is_canonical", "tss",
                         "transcript_label", "x", "xend", "y", "yend",
                         "neglogq", "neglogp", "tf_disp", "family_lab"))

## ---- internal helpers ------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

.fmt_bp <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)

.palette_pool <- c(
  "#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D",
  "#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF",
  "#1F78B4","#33A02C","#6A3D9A","#B15928","#008080","#CC6677","#332288",
  "#117733","#882255","#44AA99")

## Vectorised sliding-window log-odds PWM scorer.
.score_pwm_positions <- function(seq_chars, pwm) {
  L <- ncol(pwm); n <- length(seq_chars); ns <- n - L + 1
  if (ns <= 0) return(numeric(0))
  base_idx <- c(A = 1L, C = 2L, G = 3L, T = 4L)
  idx <- base_idx[seq_chars]
  win <- embed(idx, L)[, L:1, drop = FALSE]
  pwm_aug <- rbind(pwm, rep(0, L)); win[is.na(win)] <- 5L
  score_mat <- vapply(seq_len(L), function(j) pwm_aug[win[, j], j], numeric(ns))
  if (!is.matrix(score_mat)) score_mat <- matrix(score_mat, ncol = L)
  rowSums(score_mat)
}

## Scan one promoter window for one PWM; returns a data.frame of TSS-relative hits.
.scan_promoter <- function(promoter_info, pwm_obj, threshold_frac) {
  if (is.null(pwm_obj)) return(NULL)
  pwm <- pwm_obj$pwm; L <- pwm_obj$len
  seq <- promoter_info$seq; tss_i <- promoter_info$tss_offset + 1
  n <- nchar(seq)
  max_score <- sum(apply(pwm, 2, max)); min_score <- sum(apply(pwm, 2, min))
  threshold <- min_score + threshold_frac * (max_score - min_score)
  fwd <- strsplit(seq, "")[[1]]
  rev_map <- c(A = "T", T = "A", G = "C", C = "G", N = "N")
  rev <- rev(rev_map[fwd])
  fwd_scores <- .score_pwm_positions(fwd, pwm)
  rev_scores <- .score_pwm_positions(rev, pwm)
  fi <- which(fwd_scores >= threshold); ri <- which(rev_scores >= threshold)
  rev_start <- n - (ri + L - 2)
  hits <- rbind(
    if (length(fi)) data.frame(start = fi,        strand = "+", score = fwd_scores[fi]),
    if (length(ri)) data.frame(start = rev_start, strand = "-", score = rev_scores[ri])
  )
  if (is.null(hits) || nrow(hits) == 0)
    return(data.frame(position = integer(0), strand = character(0),
                      score = numeric(0), score_frac = numeric(0), site = character(0)))
  hits <- hits[order(hits$start, -hits$score), ]
  hits <- hits[!duplicated(hits$start), ]
  hits$position   <- hits$start - tss_i
  hits$score_frac <- (hits$score - min_score) / (max_score - min_score)
  hits$site       <- substring(seq, hits$start, hits$start + L - 1)
  hits <- hits[order(hits$position), ]; rownames(hits) <- NULL
  hits[, c("position", "strand", "score", "score_frac", "site")]
}

## Ensembl symbol lookup (expanded with transcripts). Reuses GLproxScape's
## ensembl_get when available; else a direct REST call.
.ensembl_lookup_expand <- function(gene, species) {
  path <- paste0("/lookup/symbol/", species, "/", gene)
  ns <- asNamespace("GLproxScape")
  if (exists("ensembl_get", envir = ns, inherits = FALSE)) {
    js <- tryCatch(get("ensembl_get", envir = ns)(
      path, params = list(expand = 1, `content-type` = "application/json")),
      error = function(e) NULL)
    if (!is.null(js)) return(js)
  }
  res <- httr::GET(paste0("https://rest.ensembl.org", path),
                   query = list(expand = 1), httr::accept_json(), httr::timeout(25))
  if (httr::status_code(res) != 200)
    stop("Ensembl lookup failed for ", gene, " (", species, ")")
  httr::content(res, as = "parsed", type = "application/json")
}

## Ensembl current-assembly -> UCSC genome name (for the sequence fallback).
.ucsc_assembly <- c(homo_sapiens = "hg38", mus_musculus = "mm39",
                    rattus_norvegicus = "rn7", danio_rerio = "danRer11",
                    drosophila_melanogaster = "dm6",
                    caenorhabditis_elegans = "ce11", gallus_gallus = "galGal6",
                    saccharomyces_cerevisiae = "sacCer3")

## Fallback promoter-sequence fetch from the UCSC REST API.
.fetch_seq_ucsc <- function(chr, tss, strand, upstream, downstream, species) {
  asm <- .ucsc_assembly[[species]]
  if (is.null(asm)) return(NULL)
  if (strand == 1) { ss <- max(1, tss - upstream);   se <- tss + downstream }
  else             { ss <- max(1, tss - downstream); se <- tss + upstream }
  chrom <- if (grepl("^chr", chr)) chr else paste0("chr", chr)
  url <- sprintf(
    "https://api.genome.ucsc.edu/getData/sequence?genome=%s;chrom=%s;start=%d;end=%d",
    asm, chrom, ss - 1L, se)
  res <- tryCatch(httr::GET(url, httr::timeout(30)), error = function(e) NULL)
  if (is.null(res) || httr::status_code(res) != 200) return(NULL)
  js  <- httr::content(res, as = "parsed", type = "application/json")
  dna <- toupper(js$dna %||% "")
  if (nchar(dna) < 10) return(NULL)
  if (strand == -1) {
    comp <- chartr("ACGT", "TGCA", dna)
    dna  <- paste(rev(strsplit(comp, "")[[1]]), collapse = "")
  }
  list(seq = dna, tss_offset = upstream, seq_start = ss, seq_end = se, strand = strand)
}

## Retry a flaky network call; returns the value or NULL after exhausting tries.
.with_retry <- function(fn, tries = 4, waits = c(3, 8, 15, 30), what = "request") {
  for (k in seq_len(tries)) {
    res <- tryCatch(fn(),
                    error = function(e) structure(list(msg = conditionMessage(e)),
                                                  class = "retry_fail"))
    if (!inherits(res, "retry_fail")) return(res)
    if (k < tries) {
      message("  ", what, " failed (", res$msg, "); retry ", k, "/", tries - 1,
              " in ", waits[k], "s ...")
      Sys.sleep(waits[k])
    }
  }
  message("  ", what, " failed after ", tries, " attempts - skipping.")
  NULL
}

## Mononucleotide background composition of a sequence (named A,C,G,T).
.seq_bg <- function(seq) {
  ch <- strsplit(toupper(seq), "")[[1]]
  ch <- ch[ch %in% c("A", "C", "G", "T")]
  tb <- table(factor(ch, levels = c("A", "C", "G", "T")))
  p  <- as.numeric(tb); p <- p / sum(p)
  setNames(p, c("A", "C", "G", "T"))
}

## Per-site PWM p-values via TFMPvalue (score -> P[random >= score] under `bg`).
## Returns NA for every site if TFMPvalue is not installed. Scores are
## deduplicated so the (relatively slow) exact calculation runs once per value.
.site_pvalues <- function(scores, pwm, bg) {
  if (!requireNamespace("TFMPvalue", quietly = TRUE))
    return(rep(NA_real_, length(scores)))
  m <- pwm; rownames(m) <- c("A", "C", "G", "T")
  uq <- unique(scores)
  pv <- vapply(uq, function(s)
    tryCatch(TFMPvalue::TFMsc2pv(m, s, bg, type = "PWM"),
             error = function(e) NA_real_), numeric(1))
  pv[match(scores, uq)]
}

## FIMO-style Benjamini-Hochberg q-values (FDR) for the reported sites.
## `p` are the p-values of sites kept above the score threshold; because the
## per-site p-value is monotone in the PWM score, these are the smallest
## p-values among the `m` positions tested (both strands), occupying ranks
## 1..length(p). Thus q_(i) = p_(i) * m / i, enforced monotone non-decreasing.
## This matches FIMO's p -> q conversion. Returns NA if p-values are
## unavailable (TFMPvalue not installed).
.bh_qvalue <- function(p, m) {
  if (all(is.na(p))) return(rep(NA_real_, length(p)))
  o   <- order(p)
  ps  <- p[o]
  q   <- rev(cummin(rev(ps * m / seq_along(ps))))
  res <- numeric(length(p)); res[o] <- pmin(q, 1)
  res
}

## ---- discovery-mode helpers ------------------------------------------------

## PFM counts (4 x L, rows A,C,G,T) -> log2-odds PWM against a flat background.
.pfm_to_logodds <- function(counts, pseudo = 0.8,
                            bg = c(A = .25, C = .25, G = .25, T = .25)) {
  counts <- counts[c("A", "C", "G", "T"), , drop = FALSE]
  L  <- ncol(counts); cs <- colSums(counts)
  ppm <- vapply(seq_len(L),
                function(j) (counts[, j] + pseudo * bg) / (cs[j] + pseudo),
                numeric(4))
  rownames(ppm) <- c("A", "C", "G", "T")
  log2(ppm / bg)
}

## Memory-bounded FIMO p-value. TFMPvalue::TFMsc2pv on a continuous real PWM can
## refine its score lattice until it needs tens of GB. FIMO avoids this by
## discretising the PWM to a fixed-resolution integer lattice; we do the same
## (bins steps across the matrix min..max score), capping memory per matrix.
## Returns P[random site scores >= observed] under background bg, or NA.
.fimo_pvalue <- function(pwm, score, bg, bins = 1000) {
  if (!requireNamespace("TFMPvalue", quietly = TRUE)) return(NA_real_)
  rng <- sum(apply(pwm, 2, max)) - sum(apply(pwm, 2, min))
  if (!is.finite(rng) || rng <= 0) return(NA_real_)
  gran <- rng / bins
  m <- round(pwm / gran); rownames(m) <- c("A", "C", "G", "T")
  ## score the site on the SAME integer matrix and clamp to its achievable range,
  ## else a perfect match can exceed the rounded max and TFMsc2pv returns 0.
  s <- round(score / gran)
  s <- min(max(s, sum(apply(m, 2, min))), sum(apply(m, 2, max)))
  tryCatch(TFMPvalue::TFMsc2pv(m, s, bg, type = "PWM"),
           error = function(e) NA_real_)
}

## Best single site (both strands) for one motif in a promoter sequence.
.best_site <- function(seq, tss_i, pwm_obj) {
  pwm <- pwm_obj$pwm; L <- pwm_obj$len; n <- nchar(seq)
  if (n < L) return(NULL)
  max_s <- sum(apply(pwm, 2, max)); min_s <- sum(apply(pwm, 2, min))
  fwd <- strsplit(seq, "")[[1]]
  rmap <- c(A = "T", T = "A", G = "C", C = "G", N = "N"); rv <- rev(rmap[fwd])
  fs <- .score_pwm_positions(fwd, pwm); rs <- .score_pwm_positions(rv, pwm)
  bf <- if (length(fs)) which.max(fs) else NA_integer_
  br <- if (length(rs)) which.max(rs) else NA_integer_
  sf <- if (!is.na(bf)) fs[bf] else -Inf
  sr <- if (!is.na(br)) rs[br] else -Inf
  if (sf >= sr) { s <- sf; start <- bf;               strand <- "+" }
  else          { s <- sr; start <- n - (br + L - 2); strand <- "-" }
  list(score = s, start = start, strand = strand,
       score_frac = (s - min_s) / (max_s - min_s),
       position = start - tss_i,
       site = substring(seq, start, start + L - 1))
}

## Fetch the canonical promoter for a gene (Ensembl, UCSC fallback).
.canonical_promoter <- function(gene, species, upstream, downstream) {
  js <- .ensembl_lookup_expand(gene, species)
  gene_chr <- js$seq_region_name; gene_strand <- js$strand
  txs <- js$Transcript %||% list()
  if (length(txs) == 0) stop("No transcripts returned for ", gene)
  txdf <- do.call(rbind, lapply(txs, function(tx) data.frame(
    transcript_id = tx$id %||% NA_character_,
    is_canonical  = isTRUE((tx$is_canonical %||% 0) == 1),
    start = as.numeric(tx$start), end = as.numeric(tx$end),
    stringsAsFactors = FALSE)))
  txdf$tss <- if (gene_strand == 1) txdf$start else txdf$end
  rp <- if (any(txdf$is_canonical)) txdf[txdf$is_canonical, ][1, ] else txdf[1, ]
  gene_info <- list(name = gene, chr = gene_chr, strand = gene_strand,
                    tss = rp$tss, start = rp$start, end = rp$end,
                    species = species, transcript_id = rp$transcript_id)
  pinfo <- .with_retry(
    function() GLproxScape::fetch_promoter_seq(gene_info, upstream = upstream,
                                               downstream = downstream),
    tries = 2, waits = c(3, 8), what = paste("Ensembl seq", rp$transcript_id))
  if (is.null(pinfo))
    pinfo <- .fetch_seq_ucsc(gene_chr, rp$tss, gene_strand, upstream, downstream, species)
  if (is.null(pinfo)) stop("Could not fetch canonical promoter for ", gene, ".")
  list(pinfo = pinfo, chr = gene_chr, strand = gene_strand,
       tss = rp$tss, transcript_id = rp$transcript_id)
}

## Fetch the JASPAR CORE motif collection as a list of pwm objects
## (pwm = log-odds matrix, len, id, name, family).
.jaspar_core_pwms <- function(collection = "CORE", tax_group = "vertebrates",
                              jaspar_db = NULL) {
  if (!requireNamespace("TFBSTools", quietly = TRUE))
    stop("Discovery mode needs TFBSTools: BiocManager::install('TFBSTools').")
  con <- NULL
  if (is.null(jaspar_db)) {
    cand <- c("JASPAR2024", "JASPAR2022", "JASPAR2020")
    have <- cand[vapply(cand, requireNamespace, logical(1), quietly = TRUE)]
    if (!length(have))
      stop("Install a JASPAR data package, e.g. BiocManager::install('JASPAR2024').")
    jpkg <- have[1]; message("  using motif database: ", jpkg)
    if (jpkg == "JASPAR2024") {
      if (!requireNamespace("RSQLite", quietly = TRUE))
        stop("JASPAR2024 needs RSQLite: install.packages('RSQLite').")
      jobj  <- get("JASPAR2024", asNamespace("JASPAR2024"))()
      dbfun <- get("db", asNamespace("JASPAR2024"))
      con   <- RSQLite::dbConnect(RSQLite::SQLite(), dbfun(jobj))
      jaspar_db <- con
    } else {
      jaspar_db <- get(jpkg, asNamespace(jpkg))
    }
  }
  pfms <- TFBSTools::getMatrixSet(
    jaspar_db, list(collection = collection, tax_group = tax_group,
                    matrixtype = "PFM"))
  if (!is.null(con)) try(RSQLite::dbDisconnect(con), silent = TRUE)
  if (!length(pfms)) stop("No JASPAR matrices returned for the chosen options.")
  message("  fetched ", length(pfms), " ", tax_group, " ", collection, " motifs.")
  lapply(as.list(pfms), function(pf) {
    cm <- as.matrix(TFBSTools::Matrix(pf)); rownames(cm) <- c("A", "C", "G", "T")
    fam <- tryCatch(TFBSTools::tags(pf)$family, error = function(e) NULL)
    fam <- fam[nzchar(fam)]
    fam <- if (length(fam)) paste(unique(fam), collapse = "/") else NA_character_
    list(pwm = .pfm_to_logodds(cm), len = ncol(cm),
         id = TFBSTools::ID(pf), name = TFBSTools::name(pf), family = fam)
  })
}
