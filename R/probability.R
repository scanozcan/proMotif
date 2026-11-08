## proMotif score of predicted sites, ranked site tables, top-site plot
## and motif logo. See ?score_binding_sites for the model.

#' proMotif score of predicted sites
#'
#' Combines the evidence for every predicted site of a scan into one score
#' between 0 and 1:
#' \deqn{B = 1 / (1 + e^{-L}), \quad L = \beta_0 + \beta_m m + \beta_a a + \beta_S S + \beta_R R}
#' with defaults \eqn{\beta_0 = -3, \beta_m = 4, \beta_a = 0.5, \beta_S = 2.5,
#' \beta_R = 2}.
#' \itemize{
#'   \item `motif_strength` (m): \eqn{\log_{10} p / \log_{10} p^*}, where
#'     \eqn{p^*} is the p-value of a perfect match to the same matrix (0-1).
#'   \item `db_agreement` (a): 1 if the other database (JASPAR / HOCOMOCO) calls
#'     a site within 10 bp on the same strand.
#'   \item `chip_support` (S): summit-weighted share of the ChIP-Atlas
#'     experiments binding the locus that have a peak centred on the site; the
#'     credit of a peak is shared among its motifs by relative motif strength
#'     (see [add_chip_peaks()]).
#'   \item `chip_replication` (R): grows with the number of experiments with a
#'     peak over the site (Hill curve, 0.5 at 20 experiments), scaled by the
#'     site's relative motif strength.
#' }
#' A perfect motif without ChIP-seq data reaches B = 0.73 (0.82 if both
#' databases agree). `motif_only_score` is the same score without the ChIP terms.
#' Sites are `high` (B >= 0.8), `medium` (>= 0.5) or `low` confidence, and
#' `significant` when their robust z-score of L against all binding loci
#' of the scanned windows is >= 2 and B >= 0.5. The weights are set by design,
#' not fitted: the proMotif score is a transparent evidence score for ranking
#' sites, not a probability that a site is bound.
#'
#' Binding loci: overlapping hits of the TF (the two strands of a palindrome,
#' the JASPAR and HOCOMOCO matrices, shifted frames of one match) are one
#' binding event. Hits are grouped greedily, best first; a hit joins a locus when
#' it overlaps the locus' best hit by at least `options(promotif.locus_overlap)`
#' (default 0.8) of the longer of the two spans, so a short core inside a longer,
#' different match stays separate. A locus takes the proMotif score of its
#' best hit.
#'
#' Called automatically by [scan_binding_sites()] and [add_chip_peaks()].
#'
#' @param scan A scan from [scan_binding_sites()].
#' @param weights Optional named list overriding `intercept`, `motif`,
#'   `databases`, `chip`, `replication`; default
#'   `getOption("promotif.score_weights")`. Significance thresholds:
#'   `options(promotif.sig_z = 2, promotif.sig_min_score = 0.5)`.
#' @return The scan with `motif_strength`, `db_agreement`, `motif_only_score`,
#'   `promotif_score`, `confidence`, `score_z`, `significant` and `locus` added to
#'   `scan$result`; the ranked binding loci in `scan$sites` (best hit of each
#'   locus, with `locus_start`, `locus_end`, `n_hits`, `hit_motifs`, `databases`)
#'   and every distinct hit, with its locus (`site_id`), in `scan$hits`.
#' @export
#' @examples
#' \dontrun{
#' s <- scan_binding_sites("TFF1", "ESR1")
#' s <- add_chip_peaks(s, cell_pattern = "MCF-7", max_mb = 1000)  # re-scores
#' head(s$sites)
#' }
score_binding_sites <- function(scan, weights = getOption("promotif.score_weights")) {
  scan <- .compat_scores(scan)   # drop columns of the old name (binding_prob, motif_prob)
  res <- scan$result
  if (is.null(res) || !nrow(res)) { scan$sites <- NULL; return(scan) }
  w <- utils::modifyList(list(intercept = -3, motif = 4, databases = 0.5, chip = 2.5,
                              replication = 2),
                         as.list(weights %||% list()))
  ## motif strength: -log10 p relative to a perfect match of the same matrix (0..1)
  if (all(c("pvalue", "pvalue_best") %in% names(res)) && any(is.finite(res$pvalue))) {
    res$motif_strength <- pmin(1, pmax(0, log10(res$pvalue) / log10(res$pvalue_best)))
  } else {
    message("  proMotif score uses PWM score fractions (install.packages('TFMPvalue') for p-values).")
    res$motif_strength <- pmin(1, pmax(0, (res$score_frac - (scan$meta$threshold_frac %||% 0.8)) /
                                         (1 - (scan$meta$threshold_frac %||% 0.8))))
  }
  ## database agreement: the other database calls a site within 10 bp on the same strand
  ndb <- length(unique(res$database))
  res$db_agreement <- if (ndb < 2 || !"site_center" %in% names(res)) 0 else
    vapply(seq_len(nrow(res)), function(i) {
      any(res$database != res$database[i] & abs(res$site_center - res$site_center[i]) <= 10 &
          res$motif_strand_genomic == res$motif_strand_genomic[i]) }, logical(1)) * 1
  cs <- if ("chip_support" %in% names(res)) res$chip_support else rep(0, nrow(res))
  cs[is.na(cs)] <- 0
  cr <- if ("chip_replication" %in% names(res)) res$chip_replication else rep(0, nrow(res))
  cr[is.na(cr)] <- 0
  res$motif_only_score   <- stats::plogis(w$intercept + w$motif * res$motif_strength +
                                    w$databases * res$db_agreement)
  res$promotif_score <- stats::plogis(w$intercept + w$motif * res$motif_strength +
                                    w$databases * res$db_agreement + w$chip * cs +
                                    w$replication * cr)
  res$confidence <- as.character(cut(res$promotif_score, c(-Inf, 0.5, 0.8, Inf),
                                     labels = c("low", "medium", "high"), right = FALSE))
  ## binding loci: overlapping hits of the TF (both strands, both databases, shifted
  ## frames of one match) are one binding event. Distinct sites (one per span and
  ## strand) are grouped greedily, best first: a hit joins a locus when it overlaps
  ## the locus' best hit by >= locus_overlap (default 0.8) of the longer of the two.
  ## A short core inside a longer, different match (e.g. 5 bp in 9 bp) stays apart.
  lp   <- stats::qlogis(pmin(pmax(res$promotif_score, 1e-6), 1 - 1e-6))
  if ("genomic_start" %in% names(res)) {
    key <- paste(res$genomic_start, res$genomic_end, res$motif_strand_genomic)
    o   <- order(-lp); first <- o[!duplicated(key[o])]
    loc <- .site_loci(res$genomic_start[first], res$genomic_end[first], lp[first],
                      getOption("promotif.locus_overlap", 0.8))
    res$locus <- loc[match(key, key[first])]
  } else res$locus <- seq_len(nrow(res))
  ## significance against this gene's own background of motif matches: robust
  ## z of L (logit of the score) over the binding loci of the scanned windows (median, MAD);
  ## z >= sig_z (default 2) marks a locus that stands out from the gene
  best <- tapply(lp, res$locus, max)                 # one value per binding locus
  sc   <- stats::mad(best); if (!is.finite(sc) || sc <= 0) sc <- stats::sd(best)
  if (!is.finite(sc) || sc <= 0) sc <- 1
  sig_z <- getOption("promotif.sig_z", 2)
  res$score_z     <- (lp - stats::median(best)) / sc
  ## ... and a score of at least 0.5: with a tight background (small MAD) the
  ## z-cutoff alone falls too low
  min_p <- getOption("promotif.sig_min_score", getOption("promotif.sig_min_prob", 0.5))
  res$score_cutoff <- max(min_p, stats::plogis(stats::median(best) + sig_z * sc))
  res$significant <- res$score_z >= sig_z & res$promotif_score >= min_p
  scan$result <- res
  scan$meta$score_weights <- w
  scan$meta$score_cutoff  <- res$score_cutoff[1]; scan$meta$sig_z <- sig_z
  st <- .site_table(scan)
  scan$hits  <- attr(st, "hits"); if (!is.null(st)) attr(st, "hits") <- NULL
  scan$sites <- st
  scan
}

## Greedy grouping of hits into binding loci (best first): returns a locus number
## per hit (1 = the locus of the best hit). A hit joins the current seed's locus
## when their overlap is >= frac of the longer of the two spans.
.site_loci <- function(start, end, score, frac = 0.8) {
  n <- length(start); lab <- integer(n); len <- end - start + 1; k <- 0L
  for (i in order(-score)) {
    if (lab[i]) next
    k <- k + 1L
    ov <- pmin(end, end[i]) - pmax(start, start[i]) + 1
    lab[lab == 0L & ov >= frac * pmax(len, len[i])] <- k
  }
  lab
}

## One row per binding locus, ranked by proMotif score. A locus is the best
## hit of a group of overlapping hits (see score_binding_sites); its coordinates,
## motif and scores are those of the best hit, locus_start / locus_end span all of
## its hits (useful for primer design), databases / matrix_ids / hit_motifs list
## them, motif_only_score is the best motif-only score among them. The hit-level
## table (every distinct span and strand, with its locus) is attr(, "hits").
.site_table <- function(scan) {
  res <- scan$result; m <- scan$meta
  if (is.null(res) || !nrow(res) || !"genomic_start" %in% names(res)) return(NULL)
  if (is.null(res$locus)) res$locus <- seq_len(nrow(res))
  key <- paste(res$genomic_start, res$genomic_end, res$motif_strand_genomic)
  o   <- order(key, -res$promotif_score, res$pvalue, -res$score_frac)
  r   <- res[o, , drop = FALSE]; key <- key[o]
  hit <- r[!duplicated(key), , drop = FALSE]          # one row per distinct span + strand
  agg <- function(col) tapply(r[[col]], key, function(v) paste(sort(unique(v)), collapse = ";"))
  hk  <- key[!duplicated(key)]
  hit$databases  <- agg("database")[hk]
  hit$matrix_ids <- agg("matrix_id")[hk]
  hit$in_promoter_windows <- if ("in_promoter_windows" %in% names(r))
    tapply(r$in_promoter_windows, key, function(v)
      paste(sort(unique(unlist(strsplit(v, ";")))), collapse = ";"))[hk] else hit$promoter_id
  ## other database's motif within 10 bp on the same genomic strand
  hit$n_databases_nearby <- vapply(seq_len(nrow(hit)), function(i) {
    near <- abs(res$site_center - hit$site_center[i]) <= 10 &
            res$motif_strand_genomic == hit$motif_strand_genomic[i]
    length(unique(res$database[near])) }, integer(1))
  hit$gene <- m$gene; hit$tf <- m$tf
  hit$chrom <- paste0("chr", sub("^chr", "", m$chr))
  hit$position_rel_tss <- as.integer(hit$position)
  for (cc in c("genomic_start", "genomic_end", "genomic_position", "chip_n_exp"))
    if (cc %in% names(hit)) hit[[cc]] <- as.integer(round(hit[[cc]]))
  hit <- hit[order(-hit$promotif_score, hit$pvalue, -hit$score_frac), , drop = FALSE]

  ## loci: best hit of each locus, with the evidence of all its hits
  lo <- hit[!duplicated(hit$locus), , drop = FALSE]
  by_l <- split(seq_len(nrow(hit)), hit$locus)[as.character(lo$locus)]
  cat_u <- function(v) paste(sort(unique(unlist(strsplit(as.character(v), ";")))), collapse = ";")
  lo$locus_start <- vapply(by_l, function(i) min(hit$genomic_start[i]), integer(1))
  lo$locus_end   <- vapply(by_l, function(i) max(hit$genomic_end[i]), integer(1))
  lo$n_hits      <- lengths(by_l)
  lo$hit_strands <- vapply(by_l, function(i) cat_u(hit$motif_strand_genomic[i]), "")
  lo$hit_motifs  <- vapply(by_l, function(i) paste(sprintf("%s(%s)", hit$site_motif[i],
                                                           hit$motif_strand_genomic[i]), collapse = ";"), "")
  lo$databases   <- vapply(by_l, function(i) cat_u(hit$databases[i]), "")
  lo$matrix_ids  <- vapply(by_l, function(i) cat_u(hit$matrix_ids[i]), "")
  lo$in_promoter_windows <- vapply(by_l, function(i) cat_u(hit$in_promoter_windows[i]), "")
  lo$n_databases_nearby  <- vapply(by_l, function(i) max(hit$n_databases_nearby[i]), integer(1))
  if ("motif_only_score" %in% names(hit))
    lo$motif_only_score <- vapply(by_l, function(i) max(hit$motif_only_score[i]), numeric(1))
  lo$rank    <- seq_len(nrow(lo))
  lo$site_id <- paste0("S", lo$rank)
  hit$site_id <- lo$site_id[match(hit$locus, lo$locus)]        # locus each hit belongs to
  hit$hit_id  <- paste0("H", seq_len(nrow(hit)))

  keep <- c("site_id", "rank", "gene", "tf", "chrom", "genomic_start", "genomic_end",
            "locus_start", "locus_end", "n_hits", "hit_strands", "hit_motifs",
            "motif_strand_genomic", "site_center", "genomic_position", "promoter_id",
            "in_promoter_windows", "transcript_id", "position_rel_tss",
            "databases", "matrix_ids", "n_databases_nearby",
            "site_motif", "site_genomic", "context_genomic",
            "score_frac", "pvalue", "qvalue", "motif_strength", "db_agreement", "motif_only_score",
            "chip_n_exp", "chip_support", "chip_support_all", "chip_replication", "chip_summit_dist", "chip_cell_types",
            "promotif_score", "confidence", "score_z", "significant")
  out <- lo[, intersect(keep, names(lo)), drop = FALSE]; rownames(out) <- NULL
  hk2 <- c("hit_id", "site_id", setdiff(keep, c("site_id", "rank", "locus_start", "locus_end",
                                                "n_hits", "hit_strands", "hit_motifs")))
  hits <- hit[, intersect(hk2, names(hit)), drop = FALSE]; rownames(hits) <- NULL
  attr(out, "hits") <- hits
  out
}

## Results saved before the rename (binding_prob, motif_prob) still load:
## the old columns are renamed to promotif_score and motif_only_score.
.compat_scores <- function(scan) {
  fix <- function(d) {
    if (!is.data.frame(d)) return(d)
    old_new <- c(binding_prob = "promotif_score", motif_prob = "motif_only_score")
    for (o in names(old_new)) {
      n <- old_new[[o]]
      if (!o %in% names(d)) next
      if (n %in% names(d)) d[[o]] <- NULL else names(d)[names(d) == o] <- n
    }
    d
  }
  for (k in c("result", "sites", "hits")) if (!is.null(scan[[k]])) scan[[k]] <- fix(scan[[k]])
  scan
}

#' Top binding sites of a scan
#'
#' The likely binding sites: proMotif score >= `min_score`, and at least
#' the `min_n` best sites. One row per distinct site with exact genomic
#' coordinates, sequences, p- and q-values, evidence layers and proMotif score.
#'
#' @param scan A scan from [scan_binding_sites()] (optionally with [add_chip_peaks()]).
#' @param min_score Minimum proMotif score (default 0.5).
#' @param min_n Minimum number of sites returned (default 10).
#' @return A data frame (see the README for the columns).
#' @export
top_binding_sites <- function(scan, min_score = 0.5, min_n = 10) {
  scan <- .compat_scores(scan)
  s <- scan$sites %||% .site_table(scan)
  if (is.null(s) || !nrow(s)) return(s)
  k <- max(sum(s$promotif_score >= min_score, na.rm = TRUE), min(min_n, nrow(s)))
  s[seq_len(k), , drop = FALSE]
}

#' Plot the top binding sites
#'
#' Genomic map of a scan with the proMotif score on the y axis, the
#' ChIP-Atlas peak coverage underneath (if added) and the `top_n` best sites
#' labelled (site id, motif-oriented sequence, score). Rings mark sites
#' that are significant for the locus: black = inside a ChIP-Atlas peak,
#' grey = motif evidence only.
#'
#' @inheritParams plot_promoter_map
#' @param top_n Number of sites to label (default 10; 0 = no labels).
#' @param chip_track Draw the ChIP coverage panel (default `TRUE`).
#' @return A ggplot (or patchwork with the ChIP panel).
#' @export
plot_top_sites <- function(scan, top_n = 10, base_size = 12, point_size = 2.6,
                           color_seed = NULL, chip_ring_p = 1, chip_track = TRUE) {
  plot_promoter_map(scan, base_size = base_size, point_size = point_size,
                    color_seed = color_seed, chip_ring_p = chip_ring_p,
                    y_metric = "score", label_top = top_n, chip_track = chip_track)
}

## Block-letter glyphs in a unit box for the sequence logo.
.logo_glyph <- function(letter) {
  arc <- function(a0, a1, r0 = 0.5, r1 = 0.31) {
    a <- seq(a0, a1, length.out = 60) * pi / 180
    data.frame(x = c(0.5 + r0 * cos(a), rev(0.5 + r1 * cos(a))),
               y = c(0.5 + r0 * sin(a), rev(0.5 + r1 * sin(a))), piece = 1L, sub = 1L)
  }
  switch(letter,
    A = rbind(data.frame(x = c(0, 0.4, 0.6, 1, 0.78, 0.68, 0.32, 0.22),
                         y = c(0, 1, 1, 0, 0, 0.3, 0.3, 0), piece = 1L, sub = 1L),
              data.frame(x = c(0.37, 0.63, 0.5), y = c(0.45, 0.45, 0.8), piece = 1L, sub = 2L)),
    C = arc(40, 320),
    G = rbind(arc(40, 330),
              data.frame(x = c(0.55, 1, 1, 0.55), y = c(0.40, 0.40, 0.54, 0.54), piece = 2L, sub = 1L),
              data.frame(x = c(0.82, 1, 1, 0.82), y = c(0.22, 0.22, 0.54, 0.54), piece = 3L, sub = 1L)),
    T = rbind(data.frame(x = c(0, 1, 1, 0), y = c(0.84, 0.84, 1, 1), piece = 1L, sub = 1L),
              data.frame(x = c(0.39, 0.61, 0.61, 0.39), y = c(0, 0, 0.84, 0.84), piece = 2L, sub = 1L)))
}

## probability matrix (4 x L) of a scan's PWM (log-odds vs flat background)
.pwm_to_ppm <- function(pwm) {
  p <- 0.25 * 2^pwm; p <- sweep(p, 2, colSums(p), "/"); rownames(p) <- c("A", "C", "G", "T"); p
}

#' Motif matrices of the queried TF
#'
#' The JASPAR and HOCOMOCO matrices used by a scan, as probabilities per
#' position with information content and consensus base.
#'
#' @param scan A scan from [scan_binding_sites()].
#' @return A long data frame: `database`, `matrix_id`, `position`, `A`, `C`,
#'   `G`, `T`, `information_content`, `consensus`.
#' @export
motif_matrices <- function(scan) {
  if (!length(scan$pwms)) return(NULL)
  do.call(rbind, lapply(names(scan$pwms), function(db) {
    p  <- .pwm_to_ppm(scan$pwms[[db]]$pwm)
    ic <- 2 + colSums(ifelse(p > 0, p * log2(p), 0))
    data.frame(database = db, matrix_id = scan$pwms[[db]]$id %||% NA_character_,
               position = seq_len(ncol(p)), A = p["A", ], C = p["C", ], G = p["G", ], T = p["T", ],
               information_content = ic, consensus = c("A", "C", "G", "T")[apply(p, 2, which.max)],
               stringsAsFactors = FALSE)
  }))
}

#' Sequence logo of the queried TF
#'
#' Logos of the JASPAR and HOCOMOCO matrices used by a scan, in forward and
#' reverse-complement orientation (sites on the - strand read like the reverse
#' complement).
#'
#' @param scan A scan from [scan_binding_sites()].
#' @param base_size Base font size (default 12).
#' @param revcomp Also draw the reverse complement (default `TRUE`).
#' @return A ggplot.
#' @export
plot_motif_logo <- function(scan, base_size = 12, revcomp = TRUE) {
  if (!length(scan$pwms)) return(NULL)
  cols <- c(A = "#109648", C = "#255C99", G = "#F7B32B", T = "#D62839")
  rows <- list()
  for (db in names(scan$pwms)) {
    p <- .pwm_to_ppm(scan$pwms[[db]]$pwm)
    orient <- list(forward = p)
    if (revcomp) { rc <- p[c("T", "G", "C", "A"), ncol(p):1, drop = FALSE]
                   rownames(rc) <- c("A", "C", "G", "T"); orient$`reverse complement` <- rc }
    for (o in names(orient)) {
      q  <- orient[[o]]
      ic <- 2 + colSums(ifelse(q > 0, q * log2(q), 0))
      panel <- sprintf("%s %s (%s)", db, scan$pwms[[db]]$id %||% "", o)
      for (j in seq_len(ncol(q))) {
        h <- q[, j] * ic[j]; ord <- order(h); y0 <- 0
        for (b in names(h)[ord]) {
          if (h[b] > 0.005) {
            gl <- .logo_glyph(b)
            rows[[length(rows) + 1]] <- data.frame(
              x = j - 0.45 + 0.9 * gl$x, y = y0 + h[b] * gl$y,
              grp = paste(panel, j, b, gl$piece), sub = gl$sub, base = b, panel = panel)
          }
          y0 <- y0 + h[b]
        }
      }
    }
  }
  d <- do.call(rbind, rows)
  d$panel <- factor(d$panel, levels = unique(d$panel))
  ggplot(d, aes(x, y, group = grp, subgroup = sub, fill = base)) +
    geom_polygon(color = NA) +
    facet_wrap(~ panel, ncol = 1, scales = "free_x") +
    scale_fill_manual(values = cols, name = NULL) +
    scale_x_continuous(breaks = function(l) seq(1, floor(l[2]))) +
    coord_cartesian(ylim = c(0, 2)) +
    labs(x = "Position", y = "Bits",
         title = sprintf("%s motif (the matrices used for this scan)", scan$meta$tf)) +
    theme_bw(base_size = .tsz(base_size)) +
    theme(panel.grid = element_blank(), legend.position = "none",
          strip.background = element_rect(fill = "grey92", color = NA))
}
