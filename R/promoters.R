## ===========================================================================
## 5. ALTERNATIVE PROMOTERS
## ---------------------------------------------------------------------------
## Annotated transcripts often start a few bp to a few hundred bp apart; their
## promoter windows then overlap almost completely and carry the same motifs.
## So a PROMOTER (P1, P2, ... numbered 5'->3' along the gene) = a group of TSSs
## whose scan windows overlap by >= `max_overlap` (default 20%), i.e. TSSs no
## more than (1 - max_overlap) * window length apart (2,401 bp for -2500/+500).
## `merge_dist` overrides that distance. Different promoters therefore have
## (largely) independent sequence and can be compared for differential TF content.
## ===========================================================================

## Add promoter_id (P1, P2, ... 5'->3') to a transcript table with a `tss` column.
.cluster_tss <- function(txdf, gene_strand, merge_dist = 500) {
  txdf <- txdf[order(txdf$tss), , drop = FALSE]
  cl <- cumsum(c(TRUE, diff(txdf$tss) > merge_dist))
  if (gene_strand == -1) cl <- max(cl) + 1L - cl
  txdf$promoter_id <- paste0("P", cl)
  txdf
}

## Genomic scan window(s) for TSS(s) on a strand (matrix: start, end).
.promoter_window <- function(tss, strand, upstream, downstream) {
  if (strand == 1) cbind(start = tss - upstream,   end = tss + downstream)
  else             cbind(start = tss - downstream, end = tss + upstream)
}

## One row per promoter: representative TSS (canonical, else the TSS used by
## most transcripts), spread of TSSs, window, region, overlap with the others.
.promoter_table <- function(txdf, gene_strand, upstream, downstream,
                            max_overlap = 0.2) {
  prom <- do.call(rbind, lapply(split(txdf, txdf$promoter_id), function(d) {
    rep_tss <- if (any(d$is_canonical)) d$tss[d$is_canonical][1] else {
      tt <- table(d$tss); as.numeric(names(tt)[which.max(tt)]) }
    data.frame(promoter_id = d$promoter_id[1], tss = rep_tss,
               tss_min = min(d$tss), tss_max = max(d$tss),
               span_bp = max(d$tss) - min(d$tss),
               n_transcripts = nrow(d), n_tss = length(unique(d$tss)),
               has_canonical = any(d$is_canonical),
               transcripts = paste(d$transcript_id, collapse = ";"),
               stringsAsFactors = FALSE)
  }))
  prom <- prom[order(as.integer(sub("^P", "", prom$promoter_id))), , drop = FALSE]
  rownames(prom) <- NULL
  n <- nrow(prom); len <- upstream + downstream + 1
  ## overlap between representative windows (fraction of window length)
  w  <- .promoter_window(prom$tss, gene_strand, upstream, downstream)
  ov <- diag(1, n); dimnames(ov) <- list(prom$promoter_id, prom$promoter_id)
  for (i in seq_len(n)) for (j in seq_len(n)) if (i != j)
    ov[i, j] <- max(0, min(w[i, 2], w[j, 2]) - max(w[i, 1], w[j, 1]) + 1) / len
  ## union of the windows of every TSS in a promoter (used to assign sites)
  wmin <- .promoter_window(prom$tss_min, gene_strand, upstream, downstream)
  wmax <- .promoter_window(prom$tss_max, gene_strand, upstream, downstream)
  prom$win_start <- pmin(wmin[, "start"], wmax[, "start"])
  prom$win_end   <- pmax(wmin[, "end"],   wmax[, "end"])
  ## regions: consecutive promoters whose equal-length windows overlap by
  ## >= max_overlap, i.e. TSSs no more than (1 - max_overlap) * window apart
  prom$region <- prom$promoter_id   # one level only: each promoter is its own window
  ref <- if (any(prom$has_canonical)) which(prom$has_canonical)[1] else 1L
  prom$dist_to_canonical_bp <- abs(prom$tss - prom$tss[ref])
  prom$max_overlap_other <- vapply(seq_len(n), function(i)
    if (n > 1) max(ov[i, -i]) else 0, numeric(1))
  list(promoters = prom, overlap = ov)
}

.report_promoters <- function(prom, chr, merge_dist, max_overlap, scan_unit) {
  nr <- nrow(prom)
  message(sprintf("  promoter structure: %d promoter(s) (TSSs within %s bp, i.e. windows overlapping >= %.0f%%, form one promoter); scanning per %s",
                  nr, .fmt_bp(merge_dist), 100 * max_overlap, scan_unit))
  for (i in seq_len(nr))
    message(sprintf("    %s %s| %s:%s | %d transcript(s), %d TSS, span %s bp | %s bp from canonical",
                    prom$promoter_id[i], if (prom$has_canonical[i]) "[canonical] " else "",
                    chr, .fmt_bp(prom$tss[i]), prom$n_transcripts[i], prom$n_tss[i],
                    .fmt_bp(prom$span_bp[i]), .fmt_bp(prom$dist_to_canonical_bp[i])))
  if (nr == 1)
    message("    -> one promoter: all TSSs share almost the same sequence, so a differential promoter comparison is not informative for this gene.")
  invisible(nr)
}

## Which promoters' windows contain each genomic position ("P1;P2", or "").
.covering_promoters <- function(gpos, prom) {
  vapply(gpos, function(g)
    paste(prom$promoter_id[g >= prom$win_start & g <= prom$win_end], collapse = ";"),
    character(1))
}

## Per promoter x database summary of a confirm-mode scan: distinct genomic sites
## in each promoter's window, the strongest one, and how many are specific to
## that promoter (inside no other promoter's window).
#' Per-promoter summary of a confirm-mode scan
#'
#' For each promoter (P1, P2, ... numbered 5' to 3') and database: the number of
#' distinct genomic sites in the promoter's window, how many are specific to it
#' (inside no other promoter's window), the strongest p/q-value, and, if ChIP
#' annotation was added, how many sites lie in a ChIP-Atlas peak.
#'
#' @param scan A list returned by \code{\link{scan_binding_sites}} (optionally
#'   annotated by \code{\link{add_chip_peaks}}).
#' @return A data.frame, or \code{NULL} if there are no sites.
#' @export
summarise_promoters <- function(scan) {
  res <- scan$result; prom <- scan$promoters
  if (is.null(prom) || is.null(res) || !nrow(res)) return(NULL)
  g <- res[!duplicated(res[c("genomic_position", "strand", "database", "matrix_id")]), ,
           drop = FALSE]
  out <- list()
  for (i in seq_len(nrow(prom))) for (db in unique(g$database)) {
    d <- g[g$database == db & g$genomic_position >= prom$win_start[i] &
             g$genomic_position <= prom$win_end[i], , drop = FALSE]
    spec <- if (nrow(d)) d$in_promoter_windows == prom$promoter_id[i] else logical(0)
    out[[length(out) + 1]] <- data.frame(
      promoter_id = prom$promoter_id[i],
      has_canonical = prom$has_canonical[i], tss = prom$tss[i], database = db,
      n_sites = nrow(d), n_specific_sites = sum(spec),
      best_pvalue = if (nrow(d) && any(is.finite(d$pvalue))) min(d$pvalue, na.rm = TRUE) else NA_real_,
      best_qvalue = if (nrow(d) && any(is.finite(d$qvalue))) min(d$qvalue, na.rm = TRUE) else NA_real_,
      best_score_frac = if (nrow(d)) max(d$score_frac) else NA_real_,
      n_chip_supported = if ("chip_supported" %in% names(d)) sum(d$chip_supported, na.rm = TRUE)
                         else NA_integer_,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

## NEW plot: genomic view with each promoter's window shaded (P1, P2, ...;
## * = canonical; one panel per promoter). Points = distinct sites
## (y = -log10 p if available, else score fraction). If ChIP-Atlas peaks were
## added (add_chip_peaks), sites inside a peak get a black ring and a peak
## coverage track (number of experiments with a peak) is drawn underneath.
## chip_ring_p / chip_ring_min_exp: a site gets the black ring only if it lies
## in a ChIP-Atlas peak AND its per-site p <= chip_ring_p (falls back to all
## in-peak sites when p-values are unavailable) AND >= chip_ring_min_exp
## experiments support it. Set chip_ring_p = 1 to ring every in-peak site.
## The CSV keeps chip_supported for every site regardless.
#' Promoter map: sites across all promoters of a gene
#'
#' Genomic view with each promoter's window shaded and labelled with the number
#' of transcripts grouped there (\code{*} = canonical), every TSS drawn
#' (canonical, alternative with its own window, alternative merged / not
#' scanned), and the transcription direction. Promoters far apart get their own
#' panel. If \code{\link{add_chip_peaks}} was run, a ChIP-Atlas peak-coverage
#' track is added and qualifying in-peak sites get a black ring.
#'
#' @inheritParams summarise_promoters
#' @param base_size,point_size Plot font and marker size.
#' @param color_seed Optional integer to fix the database colours.
#' @param chip_ring_p Ring only in-peak sites with per-site p at or below this
#'   (default 1e-3; 1 rings every in-peak site).
#' @param chip_ring_min_exp Ring only sites supported by at least this many
#'   ChIP experiments (default 1).
#' @param facet_regions One panel per promoter when there are several
#'   (default \code{TRUE}); \code{FALSE} draws one continuous axis.
#' @return A ggplot (or patchwork) object, or \code{NULL} if there are no sites.
#' @export
plot_promoter_map <- function(scan, base_size = 12, point_size = 2.6,
                              color_seed = NULL, chip_ring_p = 1e-3,
                              chip_ring_min_exp = 1, facet_regions = TRUE) {
  res <- scan$result; prom <- scan$promoters; m <- scan$meta
  if (is.null(prom) || is.null(res) || !nrow(res)) return(NULL)
  g <- res[!duplicated(res[c("genomic_position", "strand", "database", "matrix_id")]), ,
           drop = FALSE]
  use_p <- "pvalue" %in% names(g) && any(is.finite(g$pvalue))
  g$y  <- if (use_p) -log10(pmax(g$pvalue, 1e-300)) else g$score_frac
  ylab <- if (use_p) "-log10 per-site p" else "PWM score (fraction of max)"
  if (!is.null(color_seed)) set.seed(as.integer(color_seed))
  db_cols <- setNames(sample(.palette_pool, 2), c("JASPAR", "HOCOMOCO"))
  ybot   <- if (use_p) 0 else m$threshold_frac
  ylab_y <- ybot + (max(g$y, na.rm = TRUE) - ybot) * 1.25     # promoter labels
  ytop   <- ybot + (ylab_y - ybot) * 0.97                      # top of shading

  ## Independent regions, ordered left -> right on the chromosome. When there is
  ## more than one, each region gets its own panel with its own x axis
  ## (facet_regions = TRUE), so promoters tens of kb apart stay readable.
  reg_rng <- do.call(rbind, lapply(split(prom, prom$region), function(d) data.frame(
    region = d$region[1], xmin = min(d$win_start), xmax = max(d$win_end),
    members = sprintf("%d transcript%s", sum(d$n_transcripts),
                      if (sum(d$n_transcripts) == 1) "" else "s"),
    stringsAsFactors = FALSE)))
  reg_rng <- reg_rng[order(reg_rng$xmin), , drop = FALSE]
  reg_rng$panel <- sprintf("%s (%s)", reg_rng$region, reg_rng$members)
  lev   <- reg_rng$panel
  facet <- isTRUE(facet_regions) && nrow(reg_rng) > 1
  panel_of_region <- setNames(reg_rng$panel, reg_rng$region)
  .panel_at <- function(x) {                 # region panel of genomic positions
    i <- vapply(x, function(v) {
      inside <- which(v >= reg_rng$xmin & v <= reg_rng$xmax)
      if (length(inside)) inside[1]
      else which.min(pmin(abs(v - reg_rng$xmin), abs(v - reg_rng$xmax)))
    }, integer(1))
    factor(reg_rng$panel[i], levels = lev)
  }
  xr <- range(c(reg_rng$xmin, reg_rng$xmax))

  ## promoter windows + labels; labels of nearby promoters are staggered
  po <- prom[order(prom$region, prom$tss), , drop = FALSE]
  k  <- stats::ave(seq_len(nrow(po)), po$region, FUN = seq_along)
  rects <- data.frame(xmin = po$win_start, xmax = po$win_end, ymin = ybot, ymax = ytop,
                      lab_y = ylab_y - ((k - 1) %% 3) * (ylab_y - ybot) * 0.08,
                      region = po$region, tss = po$tss,
                      lab = sprintf("%d tx%s", as.integer(po$n_transcripts),
                                    ifelse(po$has_canonical, "*", "")))
  rects$panel <- factor(panel_of_region[rects$region], levels = lev)

  ## all TSSs: canonical, alternative with its own scanned window, and
  ## alternative TSSs merged into a promoter / not scanned on their own
  ta <- scan$tss_all %||% data.frame(tss = prom$tss, is_canonical = prom$has_canonical,
                                     scanned = TRUE)
  tss_kinds <- c("canonical TSS", "alt TSS (own window)", "alt TSS (merged / not scanned)")
  tl <- data.frame(tss = ta$tss, kind = factor(
    ifelse(ta$is_canonical, tss_kinds[1], ifelse(ta$scanned, tss_kinds[2], tss_kinds[3])),
    levels = tss_kinds))
  tl$panel <- .panel_at(tl$tss)
  g$panel  <- .panel_at(g$genomic_position)

  ## transcription-direction arrow from the canonical TSS
  ci <- which(prom$has_canonical)[1]; if (is.na(ci)) ci <- 1L
  can_panel <- panel_of_region[[prom$region[ci]]]
  span <- if (facet) with(reg_rng[reg_rng$panel == can_panel, ], xmax - xmin) else diff(xr)
  y_arr <- ybot + (ylab_y - ybot) * 0.86
  arr <- data.frame(x = prom$tss[ci],
                    xend = prom$tss[ci] + (if (m$strand == 1) 1 else -1) *
                           (if (facet) 0.08 else 0.05) * span,
                    y = y_arr, yend = y_arr, panel = factor(can_panel, levels = lev))

  base_theme <- theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = rel(0.95)),
          plot.subtitle = element_text(size = rel(0.78)))
  has_chip <- "chip_supported" %in% names(g)
  p <- ggplot() +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                                fill = region), alpha = 0.12, color = NA) +
    geom_vline(data = tl, aes(xintercept = tss, linetype = kind), color = "grey40") +
    geom_point(data = g, aes(genomic_position, y, color = database, shape = strand),
               size = point_size) +
    geom_segment(data = arr, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"),
                 color = "grey30", linewidth = 0.6, inherit.aes = FALSE) +
    geom_label(data = rects, aes(x = tss, y = lab_y, label = lab), size = base_size / 3.8,
               label.size = 0, fill = "white", alpha = 0.85)
  if (use_p) p <- p + geom_hline(yintercept = -log10(0.05), linetype = "dotted",
                                 color = "grey50")
  ring <- if (has_chip) g$chip_supported %in% TRUE &
                        g$chip_n_exp >= chip_ring_min_exp &
                        (if (use_p) g$pvalue <= chip_ring_p else TRUE)
          else logical(0)
  ring[is.na(ring)] <- FALSE
  if (has_chip && any(ring))
    p <- p + geom_point(data = g[ring, , drop = FALSE],
                        aes(genomic_position, y), shape = 21, size = point_size * 2.1,
                        stroke = 0.8, color = "black", fill = NA)
  sub <- sprintf("shaded = promoter windows (-%d/+%d bp); label = transcripts grouped at that promoter (* = canonical) | %d promoter(s)%s%s",
                 m$upstream, m$downstream, nrow(prom),
                 if (facet) "\none panel per promoter (the x axis is not continuous between panels)" else "",
                 if (!has_chip) "" else
                   sprintf("\nblack ring = in a ChIP-Atlas peak%s%s (%d of %d in-peak sites)",
                           if (use_p) sprintf(" & p <= %g", chip_ring_p) else "",
                           if (chip_ring_min_exp > 1) sprintf(" & >= %d experiments",
                                                              as.integer(chip_ring_min_exp)) else "",
                           sum(ring), sum(g$chip_supported %in% TRUE)))
  xscale <- scale_x_continuous(labels = .fmt_bp, n.breaks = if (facet) 3 else 5)
  layout <- if (facet) facet_grid(. ~ panel, scales = "free_x", space = "free_x")
            else coord_cartesian(xlim = xr)
  p <- p + scale_shape_manual(values = c("+" = 16, "-" = 17)) +
    scale_color_manual(values = db_cols, name = "Database") +
    scale_fill_brewer(palette = "Set2", name = "Promoter") +
    scale_linetype_manual(values = setNames(c("dashed", "dotted", "dotdash"), tss_kinds),
                          name = "TSS", drop = TRUE) +
    xscale + layout +
    labs(x = sprintf("Absolute position on chromosome %s (bp)", m$chr), y = ylab,
         shape = "Motif strand",
         title = sprintf("Predicted %s sites across the %s promoters", m$tf, m$gene),
         subtitle = sub) +
    base_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                       strip.background = element_rect(fill = "grey92", color = NA),
                       legend.text = element_text(size = rel(0.75)),
                       legend.title = element_text(size = rel(0.85)),
                       legend.key.height = grid::unit(0.9, "lines"),
                       legend.spacing.y = grid::unit(0.2, "lines"))
  if (!is.null(scan$chip) && nrow(scan$chip$peaks)) {
    cov <- if (facet) do.call(rbind, lapply(seq_len(nrow(reg_rng)), function(i)
             cbind(.peak_coverage(scan$chip$peaks, reg_rng$xmin[i], reg_rng$xmax[i]),
                   panel = reg_rng$panel[i])))
           else cbind(.peak_coverage(scan$chip$peaks, xr[1], xr[2]), panel = lev[1])
    cov$panel <- factor(cov$panel, levels = lev)
    cp <- scan$chip$cell_pattern
    pc <- ggplot(cov, aes(x, n)) +
      geom_area(fill = "grey35", alpha = 0.7) +
      xscale + layout +
      labs(x = sprintf("Absolute position on chromosome %s (bp)", m$chr),
           y = "ChIP experiments\nwith a peak",
           subtitle = sprintf("ChIP-Atlas %s peaks (MACS2 q < 1e-%s)%s: %d experiment(s) in these regions",
                              m$tf, scan$chip$threshold,
                              if (!is.null(cp)) paste0(", cells ~ '", cp, "'") else "",
                              scan$chip$n_experiments)) +
      base_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1),
                         strip.text = element_blank(), strip.background = element_blank())
    p <- patchwork::wrap_plots(
      p + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), pc,
      ncol = 1, heights = c(3, 1))
  }
  p
}


## ===========================================================================
## 7. DIFFERENTIAL PROMOTER REGULATORS (discovery across independent promoters)
## ---------------------------------------------------------------------------
## For genes with >= 2 independent promoter regions (section 5), scan each
## region's representative promoter against JASPAR CORE and ask which TFs have a
## strong site in one promoter but not the other(s).
## Statistic: each TF's best-site FIMO p in a window becomes a WINDOW p =
## 1 - (1 - p)^m (m = positions tested, both strands): the chance that a random
## window of this size has a site at least this strong. All windows have the
## same length, so window p-values are comparable between promoters.
## present = window p < present_p; absent = window p > absent_p. One matrix per
## TF (the best one overall) is used in every promoter. Like discovery mode this
## is PRIORITISATION, not proof of binding: about present_p of all motifs pass
## by chance in each promoter. BH q-values over all motifs are reported too.
## ===========================================================================

.gene_transcripts <- function(gene, species, transcripts = "protein_coding") {
  js  <- .ensembl_lookup_expand(gene, species)
  txs <- js$Transcript %||% list()
  if (length(txs) == 0) stop("No transcripts returned for ", gene)
  txdf <- do.call(rbind, lapply(txs, function(tx) data.frame(
    transcript_id = tx$id %||% NA_character_,
    biotype       = tx$biotype %||% NA_character_,
    is_canonical  = isTRUE((tx$is_canonical %||% 0) == 1),
    start = as.numeric(tx$start), end = as.numeric(tx$end),
    stringsAsFactors = FALSE)))
  txdf$tss <- if (js$strand == 1) txdf$start else txdf$end
  if (transcripts == "protein_coding" &&
      any(txdf$biotype == "protein_coding", na.rm = TRUE))
    txdf <- txdf[txdf$biotype == "protein_coding", , drop = FALSE]
  list(txdf = txdf, chr = js$seq_region_name, strand = js$strand)
}

#' TFs that distinguish a gene's promoters
#'
#' For genes with two or more promoters (TSS groups whose windows overlap by
#' less than \code{max_overlap}), scans each promoter against JASPAR CORE and
#' compares, per TF, the window-level p-value of its best site
#' (1 - (1 - p)^m, m = positions tested). Prioritisation, not proof of binding.
#'
#' @inheritParams scan_binding_sites
#' @inheritParams discover_regulators
#' @param present_p Window p below which a TF counts as present (default 0.01).
#' @param absent_p Window p above which a TF counts as absent (default 0.2).
#' @return A list with \code{comparison} (one row per TF), \code{long},
#'   \code{promoters} and \code{meta}; \code{comparison} is \code{NULL} for a
#'   single-promoter gene.
#' @export
compare_promoter_regulators <- function(gene, species = "homo_sapiens",
                                        upstream = 2500, downstream = 500,
                                        transcripts = c("protein_coding", "all"),
                                        merge_dist = NULL, max_overlap = 0.2,
                                        present_p = 0.01, absent_p = 0.2,
                                        collection = "CORE", tax_group = "vertebrates",
                                        use_motifmatchr = TRUE, jaspar_db = NULL,
                                jaspar_release = "2026") {
  transcripts <- match.arg(transcripts)
  message("== differential promoter regulators: ", gene, " ==")
  if (!requireNamespace("TFMPvalue", quietly = TRUE))
    stop("Needs TFMPvalue: install.packages('TFMPvalue').")
  gt   <- .gene_transcripts(gene, species, transcripts)
  merge_dist <- merge_dist %||% round((1 - max_overlap) * (upstream + downstream + 1))
  txdf <- .cluster_tss(gt$txdf, gt$strand, merge_dist)
  prom <- .promoter_table(txdf, gt$strand, upstream, downstream, max_overlap)$promoters
  .report_promoters(prom, gt$chr, merge_dist, max_overlap, "region")
  if (length(unique(prom$region)) < 2) {
    message("  only one promoter - nothing to compare. ",
            "Try transcripts = \"all\", or a lower max_overlap.")
    return(invisible(list(promoters = prom, comparison = NULL, long = NULL,
                          meta = list(gene = gene, n_regions = 1L))))
  }
  ## representative promoter per region: canonical, else most transcripts
  rp <- do.call(rbind, lapply(split(prom, prom$region), function(d)
    d[order(!d$has_canonical, -d$n_transcripts), ][1, ]))
  rp <- rp[order(as.integer(sub("^R", "", rp$region))), ]

  seqs <- list()
  for (i in seq_len(nrow(rp))) {
    d  <- txdf[txdf$promoter_id == rp$promoter_id[i], ]
    tx <- if (any(d$is_canonical)) d[d$is_canonical, ][1, ] else d[d$tss == rp$tss[i], ][1, ]
    gene_info <- list(name = gene, chr = gt$chr, strand = gt$strand, tss = rp$tss[i],
                      start = tx$start, end = tx$end, species = species,
                      transcript_id = tx$transcript_id)
    pinfo <- .with_retry(
      function() GLproxScape::fetch_promoter_seq(gene_info, upstream = upstream, downstream = downstream),
      tries = 2, waits = c(3, 8), what = paste("Ensembl seq", tx$transcript_id))
    if (is.null(pinfo))
      pinfo <- .fetch_seq_ucsc(gt$chr, rp$tss[i], gt$strand, upstream, downstream, species)
    if (is.null(pinfo)) { message("  could not fetch ", rp$promoter_id[i], " - dropped."); next }
    seqs[[rp$region[i]]] <- list(pinfo = pinfo, promoter_id = rp$promoter_id[i],
                                 tss = rp$tss[i], transcript_id = tx$transcript_id,
                                 has_canonical = rp$has_canonical[i],
                                 n_tx = rp$n_transcripts[i])
  }
  if (length(seqs) < 2) stop("Fewer than two promoter sequences could be fetched.")
  pwms <- .jaspar_core_pwms(collection, tax_group, jaspar_db, release = jaspar_release)
  n_motifs <- length(pwms)

  ## optional speed-up: keep motifs with a decent match (p < 5e-4) in ANY promoter
  ## (the rest have window p ~ 1 everywhere, so they cannot be "present")
  keep_idx <- seq_along(pwms)
  if (use_motifmatchr && all(vapply(c("motifmatchr", "Biostrings", "TFBSTools"),
                                    requireNamespace, logical(1), quietly = TRUE))) {
    allseq <- vapply(seqs, function(s) s$pinfo$seq, "")
    sl <- tryCatch({
      pwmlist <- do.call(TFBSTools::PWMatrixList, lapply(pwms, function(p) {
        mm <- p$pwm; rownames(mm) <- c("A", "C", "G", "T")
        TFBSTools::PWMatrix(ID = p$id %||% "NA", name = p$name %||% "NA", profileMatrix = mm)
      }))
      mm <- motifmatchr::matchMotifs(pwmlist, Biostrings::DNAStringSet(allseq),
                                     bg = .seq_bg(paste(allseq, collapse = "")),
                                     out = "matches", p.cutoff = 5e-4)
      which(colSums(as.matrix(motifmatchr::motifMatches(mm))) > 0)
    }, error = function(e) { message("  motifmatchr shortlist skipped: ", conditionMessage(e)); NULL })
    if (length(sl)) {
      keep_idx <- sl
      message("  motifmatchr shortlisted ", length(sl), " motifs (match in any promoter).")
    }
  }

  long <- list()
  for (r in names(seqs)) {
    s <- seqs[[r]]; sq <- s$pinfo$seq; tss_i <- s$pinfo$tss_offset + 1; bg <- .seq_bg(sq)
    message("  scanning ", r, " (", s$promoter_id, ") against ", length(keep_idx), " motifs ...")
    rows <- lapply(keep_idx, function(k) {
      p  <- pwms[[k]]; bs <- .best_site(sq, tss_i, p)
      if (is.null(bs)) return(NULL)
      pv   <- .fimo_pvalue(p$pwm, bs$score, bg)
      mpos <- 2 * max(1, nchar(sq) - p$len + 1)
      data.frame(tf = p$name %||% p$id, matrix_id = p$id %||% NA_character_,
                 family = p$family, region = r, promoter_id = s$promoter_id,
                 position = bs$position, strand = bs$strand, score_frac = bs$score_frac,
                 site = bs$site,
                 genomic_position = if (gt$strand == 1) s$tss + bs$position
                                    else s$tss - bs$position,
                 p_site = pv,
                 p_window = if (is.finite(pv)) -expm1(mpos * log1p(-pv)) else NA_real_,
                 stringsAsFactors = FALSE)
    })
    d <- do.call(rbind, rows)
    d$q_window <- stats::p.adjust(d$p_window, "BH",
                                  n = max(n_motifs, sum(!is.na(d$p_window))))
    long[[r]] <- d
  }
  long <- do.call(rbind, long); rownames(long) <- NULL

  ## one matrix per TF: the one with the smallest window p in any promoter
  best_mat <- suppressWarnings(tapply(seq_len(nrow(long)), long$tf, function(ix) {
    mn <- tapply(long$p_window[ix], long$matrix_id[ix], min, na.rm = TRUE)
    names(mn)[which.min(mn)]
  }))
  long <- long[long$matrix_id == as.character(best_mat[long$tf]), , drop = FALSE]

  regs <- names(seqs); tfs <- unique(long$tf)
  first <- long[!duplicated(long$tf), ]
  wide  <- data.frame(tf = tfs, matrix_id = first$matrix_id[match(tfs, first$tf)],
                      family = first$family[match(tfs, first$tf)], stringsAsFactors = FALSE)
  for (r in regs) {
    d <- long[long$region == r, ]; ix <- match(tfs, d$tf)
    wide[[paste0("p_window_", r)]] <- d$p_window[ix]
    wide[[paste0("q_window_", r)]] <- d$q_window[ix]
    wide[[paste0("p_site_", r)]]   <- d$p_site[ix]
    wide[[paste0("position_", r)]] <- d$position[ix]
  }
  P <- as.matrix(wide[paste0("p_window_", regs)]); P[is.na(P)] <- 1
  present <- P < present_p; absent <- P > absent_p
  wide$n_present   <- rowSums(present)
  wide$best_region <- regs[apply(P, 1, which.min)]
  srt <- t(apply(P, 1, sort))
  wide$specificity <- -log10(pmax(srt[, 1], 1e-300)) + log10(pmax(srt[, 2], 1e-300))
  wide$class <- ifelse(wide$n_present == 0, "none",
                ifelse(wide$n_present == 1 & rowSums(absent) == length(regs) - 1,
                       paste0(wide$best_region, "-specific"),
                ifelse(wide$n_present >= 2, "shared", "weak/ambiguous")))
  rank_cls <- ifelse(grepl("-specific$", wide$class), 1L,
              ifelse(wide$class == "shared", 2L, ifelse(wide$class == "none", 4L, 3L)))
  wide <- wide[order(rank_cls, -wide$specificity), ]; rownames(wide) <- NULL

  region_labels <- vapply(regs, function(r) {
    s <- seqs[[r]]
    sprintf("%s%s (%d tx)\n%s:%s", s$promoter_id, if (s$has_canonical) "*" else "",
            as.integer(s$n_tx), gt$chr, .fmt_bp(s$tss))
  }, "")
  cnt <- table(factor(wide$class))
  message("  calls: ", paste(sprintf("%s %d", names(cnt), as.integer(cnt)), collapse = " | "))
  list(comparison = wide, long = long, promoters = prom,
       meta = list(gene = gene, species = species, chr = gt$chr, strand = gt$strand,
                   upstream = upstream, downstream = downstream, regions = regs,
                   region_labels = region_labels, present_p = present_p,
                   absent_p = absent_p, n_motifs = n_motifs,
                   n_tested = length(keep_idx), n_regions = length(regs)))
}

## Dot plot: TFs (rows) x independent promoters (columns); size = -log10 window
## p, opacity = present / weak / absent, colour = call (specific / shared).
#' Dot plot of promoter-specific and shared TFs
#'
#' @param cmp A list returned by \code{\link{compare_promoter_regulators}}.
#' @param top_n Number of TFs to show (default 30).
#' @param base_size Base font size.
#' @return A ggplot, or \code{NULL}.
#' @export
plot_promoter_regulators <- function(cmp, top_n = 30, base_size = 12) {
  w <- cmp$comparison; m <- cmp$meta
  if (is.null(w)) return(NULL)
  sel <- head(w[w$class != "none", , drop = FALSE], top_n)
  if (!nrow(sel)) return(NULL)
  regs <- m$regions
  d <- do.call(rbind, lapply(regs, function(r) data.frame(
    tf = sel$tf, class = sel$class, region = m$region_labels[[r]],
    p = sel[[paste0("p_window_", r)]], stringsAsFactors = FALSE)))
  d$p[is.na(d$p)] <- 1
  d$neglogp <- -log10(pmax(d$p, 1e-300))
  d$status  <- factor(ifelse(d$p < m$present_p, "present",
                      ifelse(d$p > m$absent_p, "absent", "weak")),
                      levels = c("present", "weak", "absent"))
  d$tf     <- factor(d$tf, levels = rev(unique(sel$tf)))
  d$region <- factor(d$region, levels = unname(m$region_labels[regs]))
  ggplot(d, aes(region, tf)) +
    geom_point(aes(size = neglogp, color = class, alpha = status)) +
    scale_alpha_manual(values = c(present = 1, weak = 0.45, absent = 0.15),
                       name = "Site in window", drop = FALSE) +
    scale_size_continuous(name = "-log10 window p", range = c(0.6, 7)) +
    labs(x = NULL, y = NULL, color = "Call",
         title = sprintf("TFs that distinguish the %s promoters", m$gene),
         subtitle = sprintf("window -%d/+%d bp; window p = chance a random window this size has a site this strong\npresent < %g, absent > %g | best matrix per TF | %d of %d motifs scanned | * = canonical | prioritisation, not binding",
                            m$upstream, m$downstream, m$present_p, m$absent_p,
                            m$n_tested, m$n_motifs)) +
    theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = rel(0.95)),
          plot.subtitle = element_text(size = rel(0.72)))
}

## Convenience wrapper: compare + plot + write files to results/<GENE>_promoters/.
#' Compare a gene's promoters and write the results
#'
#' Runs \code{\link{compare_promoter_regulators}} and
#' \code{\link{plot_promoter_regulators}}, writing CSVs and a PDF to
#' \code{out_root/<GENE>_promoters/}.
#'
#' @inheritParams compare_promoter_regulators
#' @inheritParams plot_promoter_regulators
#' @param out_root Output folder (default \code{"results"}).
#' @param write_files Write the files (default \code{TRUE}).
#' @param ... Passed to \code{\link{compare_promoter_regulators}}.
#' @return Invisibly, the comparison list plus \code{plot}.
#' @examples
#' \dontrun{
#' run_promoter_comparison("TP73", transcripts = "all")
#' }
#' @export
run_promoter_comparison <- function(gene, top_n = 30, out_root = "results",
                                    write_files = TRUE, base_size = 12, ...) {
  cmp <- compare_promoter_regulators(gene, ...)
  if (is.null(cmp$comparison)) return(invisible(cmp))
  plot <- plot_promoter_regulators(cmp, top_n = top_n, base_size = base_size)
  if (write_files) {
    out_dir <- file.path(out_root, paste0(gene, "_promoters"))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(cmp$promoters, file.path(out_dir, sprintf("%s_promoter_structure.csv", gene)),
              row.names = FALSE)
    write.csv(cmp$comparison, file.path(out_dir, sprintf("%s_promoter_regulators.csv", gene)),
              row.names = FALSE)
    write.csv(cmp$long, file.path(out_dir, sprintf("%s_promoter_regulators_long.csv", gene)),
              row.names = FALSE)
    if (!is.null(plot))
      ggsave(file.path(out_dir, sprintf("%s_promoter_regulators.pdf", gene)), plot,
             width = 4 + 2.2 * cmp$meta$n_regions,
             height = max(4, 0.28 * min(top_n, sum(cmp$comparison$class != "none")) + 2.2),
             limitsize = FALSE)
    message("Wrote: ", out_dir)
  }
  message("Done (promoter comparison).")
  invisible(c(cmp, list(plot = plot)))
}
