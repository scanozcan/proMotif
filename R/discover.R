## ===========================================================================
## Discovery mode: "what are the top candidate regulators of this gene?"
## ---------------------------------------------------------------------------
## Scan the canonical promoter against the whole JASPAR CORE collection and rank
## TFs by their single strongest match, scored with the SAME FIMO per-site
## p-value as confirm mode. This is motif-match PRIORITISATION, not enrichment
## and not evidence of binding - see the README "Discovery mode" section.
## ===========================================================================

#' Rank candidate regulators of a gene's promoter (discovery mode)
#'
#' Scans every promoter of \code{gene} (TSS groups, as in
#' \code{\link{scan_binding_sites}}) against the whole JASPAR CORE
#' motif collection and ranks transcription factors by their single strongest
#' predicted site, scored with the same FIMO per-site p-value used by
#' \code{\link{scan_binding_sites}}. TFs whose best site falls at the same locus
#' (within \code{site_window} bp) are grouped so one element is not reported
#' many times; the full per-TF table is kept, with the shared members listed.
#'
#' This is a motif-match \emph{prioritisation}, not a motif-enrichment test
#' (there is no background set of sequences) and not evidence of binding
#' (chromatin, concentration, cofactors and cell type decide that). Treat the
#' output as "which TFs are worth following up".
#'
#' Requires the Bioconductor packages \pkg{TFBSTools} and a JASPAR data package
#' (\pkg{JASPAR2024}, with \pkg{RSQLite}), plus \pkg{TFMPvalue} for the p-value.
#' \pkg{motifmatchr} + \pkg{Biostrings} are optional and only pre-shortlist
#' motifs for speed.
#'
#' @param gene HGNC gene symbol.
#' @param species Ensembl species token (default \code{"homo_sapiens"}).
#' @param upstream,downstream Promoter window around the canonical TSS (bp).
#' @param collection,tax_group JASPAR collection and taxonomic group
#'   (defaults \code{"CORE"} / \code{"vertebrates"}).
#' @param top_n How many representative TFs to return in \code{top} / plot.
#' @param collapse_site Group TFs whose best site is the same locus (default TRUE).
#' @param site_window Window (bp) within which best sites count as the same
#'   locus (default 20).
#' @param use_motifmatchr Use \pkg{motifmatchr} to pre-shortlist motifs if it is
#'   installed (speed only; the ranking p-value is always the exact value).
#' @param jaspar_release JASPAR release year downloaded from jaspar.elixir.no
#'   (default \code{"2026"}; cached). Ignored when \code{jaspar_db} is given; set
#'   \code{NULL} to use a Bioconductor JASPAR package instead.
#' @param rank_by \code{"background"} (default): rank TFs by how unusual this
#'   gene's best site is compared with the best site of the same motif in all
#'   20,094 human protein-coding promoters (empirical; column \code{bg_pvalue}). Removes the
#'   advantage of long or GC-rich motifs that match well in almost any promoter.
#'   Human and JASPAR 2026 only; otherwise falls back to \code{"fimo"} (the
#'   best site's FIMO per-site p-value).
#' @param promoters \code{"all"} (default; one ranking per promoter) or
#'   \code{"canonical"} (canonical promoter only).
#' @inheritParams scan_binding_sites
#' @param jaspar_db Optional pre-opened JASPAR database/connection to use
#'   instead of auto-detecting an installed JASPAR package.
#'
#' @return A list with \code{result} (the full per-TF table: best-site position,
#'   strand, score, \code{score_frac}, \code{pvalue}, \code{qvalue},
#'   \code{genomic_position}, \code{family}, \code{site_group},
#'   \code{is_representative}, \code{n_at_site}, \code{shares_site_with}),
#'   \code{top} (the \code{top_n} representatives per promoter; column
#'   \code{promoter_id}), \code{promoters} and \code{meta}.
#' @examples
#' \dontrun{
#' d <- discover_regulators("SNCA")
#' head(d$top[, c("tf", "pvalue", "n_at_site", "shares_site_with")])
#' }
#' @export
discover_regulators <- function(gene, species = "homo_sapiens",
                                upstream = 2500, downstream = 500,
                                collection = "CORE", tax_group = "vertebrates",
                                top_n = 20, collapse_site = TRUE, site_window = 20,
                                use_motifmatchr = TRUE, jaspar_db = NULL,
                                jaspar_release = "2026",
                                rank_by = c("background", "fimo"),
                                promoters = c("all", "canonical"),
                                transcripts = c("protein_coding", "all"),
                                merge_dist = NULL, max_overlap = 0.2) {
  rank_by <- match.arg(rank_by); promoters <- match.arg(promoters)
  transcripts <- match.arg(transcripts)
  message("== discovery: candidate regulators of ", gene, " ==")
  if (!requireNamespace("TFMPvalue", quietly = TRUE))
    stop("Discovery mode needs TFMPvalue: install.packages('TFMPvalue').")

  ## the gene's promoters (TSS groups, as in confirm mode); one window each,
  ## representative TSS = canonical, else the TSS used by most transcripts
  gt <- .gene_transcripts(gene, species, transcripts)
  merge_dist <- merge_dist %||% round((1 - max_overlap) * (upstream + downstream + 1))
  txdf <- .cluster_tss(gt$txdf, gt$strand, merge_dist)
  prom <- .promoter_table(txdf, gt$strand, upstream, downstream, max_overlap)$promoters
  .report_promoters(prom, gt$chr, merge_dist, max_overlap, "promoter")
  ci <- which(prom$has_canonical)[1]; if (is.na(ci)) ci <- 1L
  prom_scan <- if (promoters == "canonical") prom[ci, , drop = FALSE] else prom

  pwms <- .jaspar_core_pwms(collection, tax_group, jaspar_db, release = jaspar_release)
  pbg  <- if (rank_by == "background")
    .promoter_background(jaspar_release, species, upstream, downstream,
                         custom_db = !is.null(jaspar_db)) else NULL
  if (rank_by == "background" && is.null(pbg)) {
    message("  -> ranking by FIMO per-site p instead.")
    rank_by <- "fimo"
  }

  res_all <- list(); n_tested <- integer(0)
  for (k in seq_len(nrow(prom_scan))) {
    pr <- prom_scan[k, ]
    d  <- txdf[txdf$promoter_id == pr$promoter_id, ]
    tx <- if (any(d$is_canonical)) d[d$is_canonical, ][1, ] else d[d$tss == pr$tss, ][1, ]
    gene_info <- list(name = gene, chr = gt$chr, strand = gt$strand, tss = pr$tss,
                      start = tx$start, end = tx$end, species = species,
                      transcript_id = tx$transcript_id)
    pinfo <- .with_retry(
      function() GLproxScape::fetch_promoter_seq(gene_info, upstream = upstream, downstream = downstream),
      tries = 2, waits = c(3, 8), what = paste("Ensembl seq", tx$transcript_id))
    if (is.null(pinfo)) {
      message("  falling back to UCSC for ", pr$promoter_id, " ...")
      pinfo <- .fetch_seq_ucsc(gt$chr, pr$tss, gt$strand, upstream, downstream, species)
      if (!is.null(pinfo)) message("  Got ", nchar(pinfo$seq), " bp (UCSC)")
    }
    if (is.null(pinfo)) { message("  could not fetch ", pr$promoter_id, " - skipped."); next }
    if (nrow(prom_scan) > 1)
      message("  scanning ", pr$promoter_id, if (pr$has_canonical) " (canonical)" else "",
              " against ", length(pwms), " motifs ...")
    one <- .discover_window(pinfo$seq, pinfo$tss_offset + 1, pr$tss, gt$strand, pwms,
                            pbg, rank_by, use_motifmatchr, collapse_site, site_window)
    if (is.null(one)) next
    one$res$promoter_id        <- pr$promoter_id
    one$res$promoter_tss       <- pr$tss
    one$res$canonical_promoter <- pr$has_canonical
    res_all[[pr$promoter_id]]  <- one$res
    n_tested[pr$promoter_id]   <- one$n_tested
  }
  if (!length(res_all)) stop("No motif matches found for ", gene, ".")
  res <- do.call(rbind, res_all); rownames(res) <- NULL
  top <- do.call(rbind, lapply(res_all, function(r)
    utils::head(r[r$is_representative, , drop = FALSE], top_n)))
  rownames(top) <- NULL

  list(result = res, top = top, promoters = prom,
       meta = list(gene = gene, species = species, chr = gt$chr, strand = gt$strand,
                   tss = prom$tss[ci],
                   transcript_id = txdf$transcript_id[txdf$is_canonical][1],
                   upstream = upstream, downstream = downstream,
                   n_tested = max(n_tested), n_tested_by_promoter = n_tested,
                   top_n = top_n, collection = collection, tax_group = tax_group,
                   jaspar_release = jaspar_release, rank_by = rank_by,
                   transcripts = transcripts, n_promoters = length(res_all),
                   bg_n = if (!is.null(pbg)) pbg$n_bg[1] else NA_integer_,
                   promoters_scanned = names(res_all)))
}

.discover_window <- function(seq, tss_i, tss, strand, pwms, pbg, rank_by,
                             use_motifmatchr = TRUE, collapse_site = TRUE,
                             site_window = 20) {
  bg <- .seq_bg(seq)
  ## optional speed-up: keep motifs with any match at p < 5e-4
  keep_idx <- seq_along(pwms)
  if (use_motifmatchr && requireNamespace("motifmatchr", quietly = TRUE) &&
      requireNamespace("Biostrings", quietly = TRUE) &&
      requireNamespace("TFBSTools", quietly = TRUE)) {
    sl <- tryCatch({
      pwmlist <- do.call(TFBSTools::PWMatrixList, lapply(pwms, function(p) {
        m <- p$pwm; rownames(m) <- c("A", "C", "G", "T")
        TFBSTools::PWMatrix(ID = p$id %||% "NA", name = p$name %||% "NA",
                            profileMatrix = m)
      }))
      mm <- motifmatchr::matchMotifs(pwmlist, Biostrings::DNAStringSet(seq), bg = bg,
                                     out = "matches", p.cutoff = 5e-4)
      which(as.logical(motifmatchr::motifMatches(mm)[1, ]))
    }, error = function(e) {
      message("  motifmatchr shortlist skipped: ", conditionMessage(e)); NULL })
    if (length(sl)) keep_idx <- sl
  }
  n_tested <- length(keep_idx)
  rows <- lapply(keep_idx, function(k) {
    p  <- pwms[[k]]
    bs <- .best_site(seq, tss_i, p)
    if (is.null(bs)) return(NULL)
    data.frame(tf = p$name %||% p$id, matrix_id = p$id %||% NA_character_,
               family = p$family, strand = bs$strand,
               position = bs$position, score = bs$score, score_frac = bs$score_frac,
               pvalue = .fimo_pvalue(p$pwm, bs$score, bg), site = bs$site,
               genomic_position = if (strand == 1) tss + bs$position else tss - bs$position,
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)
  if (is.null(res) || !nrow(res)) return(NULL)

  ## promoter background: how unusual is this best site for the motif?
  if (!is.null(pbg)) {
    i <- match(res$matrix_id, pbg$matrix_id)
    res$bg_pvalue <- .bg_pvalue(res$score, pbg, i)
    res$bg_zscore <- (res$score - pbg$mean[i]) / pbg$sd[i]
  }
  if (rank_by == "background") {
    ## background p-values within 2% of the floor (1 / number of background
    ## promoters) are ties - their differences are interpolation noise - so
    ## they are ordered by the FIMO per-site p-value instead
    fl <- min(attr(pbg, "probs"))
    bp <- ifelse(is.na(res$bg_pvalue), 1, res$bg_pvalue)
    res$bg_at_floor <- bp <= 1.02 * fl
    res$rank_p <- ifelse(abs(bp - fl) <= 0.02 * fl, fl, bp)
  } else res$rank_p <- res$pvalue

  ## one row per TF name (best matrix), best-first
  res <- res[order(res$rank_p, res$pvalue, -res$score_frac), ]
  res <- res[!duplicated(res$tf), ]
  res$n_at_site <- 1L; res$shares_site_with <- ""
  res$is_representative <- TRUE; res$site_group <- seq_len(nrow(res))

  ## group TFs whose best site is the same locus (within site_window bp);
  ## nothing is dropped, the representative lists the folded members
  if (collapse_site && nrow(res) > 1) {
    grp <- integer(nrow(res)); reps_pos <- numeric(0)
    for (i in seq_len(nrow(res))) {
      hit <- which(abs(reps_pos - res$genomic_position[i]) <= site_window)
      if (length(hit)) grp[i] <- hit[1]
      else { reps_pos <- c(reps_pos, res$genomic_position[i]); grp[i] <- length(reps_pos) }
    }
    mem <- split(res$tf, grp)
    res$site_group        <- grp
    res$n_at_site         <- lengths(mem)[as.character(grp)]
    res$is_representative <- !duplicated(grp)
    res$shares_site_with  <- vapply(seq_len(nrow(res)), function(i) {
      others <- setdiff(mem[[as.character(grp[i])]], res$tf[i])
      if (length(others)) paste(others, collapse = ", ") else ""
    }, character(1))
  }
  res$qvalue <- .bh_qvalue(res$pvalue, n_tested)
  if (!is.null(res$bg_pvalue))
    res$bg_qvalue <- stats::p.adjust(res$bg_pvalue, "BH",
                                     n = max(length(pwms), sum(!is.na(res$bg_pvalue))))
  res <- res[order(res$rank_p, res$pvalue, -res$score_frac), ]
  res$rank_p <- NULL; rownames(res) <- NULL
  list(res = res, n_tested = n_tested)
}

#' Targeted follow-up of discovery hits
#'
#' For the top TFs of a \code{\link{discover_regulators}} result (or the TFs
#' given), runs the confirm-mode scan (\code{\link{scan_binding_sites}}) on the
#' same gene and window - all promoters, JASPAR and HOCOMOCO, per-site p/q - and,
#' with \code{chip = TRUE}, adds ChIP-Atlas peaks (\code{\link{add_chip_peaks}}).
#' Plot each element with \code{\link{plot_promoter_map}} or
#' \code{\link{plot_binding_sites}}.
#'
#' @param disc A list returned by \code{\link{discover_regulators}}.
#' @param tfs Optional character vector of TFs; default the top \code{n}.
#' @param n Number of top TFs when \code{tfs} is NULL.
#' @param by Which top TFs: \code{"auto"} (default; \code{"chip"} if the
#'   discovery result has ChIP-Atlas annotation, else \code{"rank"}),
#'   \code{"rank"} (motif ranking) or \code{"chip"} (TFs with ChIP-Atlas peaks
#'   at the gene first, by share of experiments bound, then the rest by rank).
#' @param chip Add ChIP-Atlas peaks of each TF (default FALSE).
#' @param threshold_frac,hocomoco_version,scan_unit Passed to
#'   \code{\link{scan_binding_sites}}.
#' @param chip_cell_pattern,chip_max_mb Passed to \code{\link{add_chip_peaks}}.
#' @param ... Further arguments to \code{\link{scan_binding_sites}}.
#' @return A named list of scan objects (one per TF that could be scanned).
#' @examples
#' \dontrun{
#' d  <- discover_regulators("ABCB1")
#' ex <- expand_regulators(d, n = 3, chip = TRUE)
#' plot_promoter_map(ex[[1]])
#' }
#' @export
expand_regulators <- function(disc, tfs = NULL, n = 5, by = c("auto", "rank", "chip"),
                              chip = FALSE,
                              threshold_frac = 0.80, hocomoco_version = "v14",
                              scan_unit = "promoter", chip_cell_pattern = NULL,
                              chip_max_mb = 300, ...) {
  m <- disc$meta; by <- match.arg(by)
  if (is.null(tfs)) {
    t   <- disc$top
    key <- if (identical(m$rank_by, "background") && !is.null(t$bg_pvalue)) t$bg_pvalue else t$pvalue
    t   <- t[order(key, t$pvalue), , drop = FALSE]
    has_chip <- "chip_n_bound" %in% names(t)
    if (by == "auto") by <- if (has_chip) "chip" else "rank"
    if (by == "chip" && !has_chip) {
      message("  no ChIP-Atlas annotation (use chip = TRUE / add_chip_targets) - expanding by rank.")
      by <- "rank"
    }
    if (by == "chip") {
      ## TFs with ChIP-Atlas peaks at this gene first (share of experiments
      ## bound, then number bound), then the rest in motif-rank order
      bound <- !is.na(t$chip_n_bound) & t$chip_n_bound > 0
      tb <- t[bound, , drop = FALSE]
      tb <- tb[order(-tb$chip_frac_exp, -tb$chip_n_bound), , drop = FALSE]
      t  <- rbind(tb, t[!bound, , drop = FALSE])
    }
    tfs <- utils::head(unique(t$tf), n)
  }
  message("== expanding ", length(tfs), " top TF(s) of ", m$gene, ": ",
          paste(tfs, collapse = ", "), " ==")
  out <- list()
  for (tf in tfs) {
    s <- tryCatch(
      scan_binding_sites(m$gene, tf, species = m$species, upstream = m$upstream,
                         downstream = m$downstream, threshold_frac = threshold_frac,
                         hocomoco_version = hocomoco_version,
                         transcripts = m$transcripts %||% "protein_coding",
                         scan_unit = scan_unit, jaspar_release = m$jaspar_release, ...),
      error = function(e) { message("  ", tf, ": ", conditionMessage(e)); NULL })
    if (is.null(s)) next
    if (!nrow(s$result))
      message("  ", tf, ": no sites at threshold ", threshold_frac,
              " in any promoter window - lower the threshold to see it.")
    if (isTRUE(chip) && nrow(s$result))
      s <- add_chip_peaks(s, cell_pattern = chip_cell_pattern, max_mb = chip_max_mb)
    out[[tf]] <- s
  }
  out
}

#' Plots for discovery mode
#'
#' Builds three ggplots from a \code{\link{discover_regulators}} result: a
#' ranked lollipop of the top TFs, a rank-vs-position map, and a promoter needle
#' map (position vs significance, showing where candidate sites cluster). Points
#' are coloured by TF family and shaped by strand; \code{(+N)} on a label means N
#' further TFs share that site. \pkg{ggrepel}, if installed, is used for tidy
#' needle-map labels.
#'
#' @param disc A list returned by \code{\link{discover_regulators}}.
#' @param base_size Base font size.
#' @param point_size Marker size.
#' @return A list of ggplot objects: \code{ranking}, \code{map}, \code{needle}.
#' @examples
#' \dontrun{
#' p <- plot_regulators(discover_regulators("SNCA"))
#' p$needle
#' }
#' @export
plot_regulators <- function(disc, base_size = 12, point_size = 3) {
  top <- disc$top; m <- disc$meta
  if (is.null(top) || !nrow(top)) return(list(ranking = NULL, map = NULL, needle = NULL))
  use_bg <- identical(m$rank_by, "background") && !is.null(top$bg_pvalue)
  top$neglogp <- -log10(pmax(if (use_bg) top$bg_pvalue else top$pvalue, 1e-300))
  plab <- if (use_bg) "-log10 promoter-background p (best site)" else
                      "-log10 FIMO per-site p (best site)"
  nshare <- if (!is.null(top$n_at_site)) top$n_at_site else rep(1L, nrow(top))
  disp   <- ifelse(nshare > 1, sprintf("%s (+%d)", top$tf, nshare - 1L), as.character(top$tf))
  top$label <- disp
  famlab <- ifelse(is.na(top$family) | top$family == "", "unassigned", top$family)
  top$family_lab <- factor(famlab)

  ## one panel per scanned promoter (labelled with its transcripts / TSS)
  if (is.null(top$promoter_id)) top$promoter_id <- "P1"
  pr <- disc$promoters
  pan_lab <- function(pid) {
    i <- if (is.null(pr)) NA else match(pid, pr$promoter_id)
    ifelse(is.na(i), pid,
           sprintf("%s%s | %d tx | %s:%s", pid, ifelse(pr$has_canonical[i], " (canonical)", ""),
                   pr$n_transcripts[i], m$chr, .fmt_bp(pr$tss[i])))
  }
  pids <- unique(top$promoter_id)
  pids <- pids[order(as.integer(sub("^P", "", pids)))]
  top$panel <- factor(pan_lab(top$promoter_id), levels = pan_lab(pids))
  multi <- length(pids) > 1
  ## y keys unique per panel; most significant TF at the top of each panel
  top$ykey <- paste(disp, top$promoter_id, sep = "___")
  o <- order(top$panel, top$neglogp)
  top$tf_disp <- factor(top$ykey, levels = unique(top$ykey[o]))
  ylabs <- function(x) sub("___.*$", "", x)
  facet <- if (multi) facet_wrap(~ panel, nrow = 1, scales = "free_y") else NULL

  base_theme <- theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = rel(0.95)),
          plot.subtitle = element_text(size = rel(0.78)),
          strip.background = element_rect(fill = "grey92", color = NA))
  sub <- sprintf("%s -%d/+%d bp | %d motifs tested | %s | match-priority, not binding | (+N) = more TFs share that site",
                 if (multi) sprintf("each of %d promoters", length(pids)) else "promoter",
                 m$upstream, m$downstream, m$n_tested,
                 if (use_bg) sprintf("ranked vs best sites in %s human promoters",
                                     format(m$bg_n %||% 20094, big.mark = ",")) else "FIMO per-site p")

  ranking <- ggplot(top, aes(neglogp, tf_disp)) +
    geom_segment(aes(x = 0, xend = neglogp, yend = tf_disp), color = "grey70") +
    geom_point(aes(color = family_lab), size = point_size) +
    scale_y_discrete(labels = ylabs) + facet +
    labs(x = plab, y = NULL, color = "TF family",
         title = sprintf("Top candidate regulators of %s", m$gene), subtitle = sub) +
    base_theme

  map <- ggplot(top, aes(position, tf_disp)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    geom_point(aes(color = family_lab, shape = strand), size = point_size) +
    scale_shape_manual(values = c("+" = 16, "-" = 17)) +
    scale_y_discrete(labels = ylabs) + facet +
    labs(x = "Best-site position relative to TSS (bp)", y = NULL,
         color = "TF family", shape = "Motif strand",
         title = sprintf("Where the top %s regulators' best sites sit", m$gene),
         subtitle = sub) + base_theme

  ## needle map: position vs significance per promoter; top hits labelled
  needle <- ggplot(top, aes(position, neglogp)) +
    geom_segment(aes(xend = position, yend = 0), color = "grey75", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(aes(color = family_lab, shape = strand), size = point_size) +
    scale_shape_manual(values = c("+" = 16, "-" = 17)) +
    (if (multi) facet_wrap(~ panel, nrow = 1) else NULL) +
    labs(x = "Position relative to TSS (bp)", y = sub(" \\(best site\\)", "", plab),
         color = "TF family", shape = "Motif strand",
         title = sprintf("Where %s's candidate regulators sit", m$gene),
         subtitle = sub) + base_theme
  lab_df <- do.call(rbind, lapply(split(top, top$panel, drop = TRUE), function(d)
    utils::head(d[order(-d$neglogp), , drop = FALSE], 10)))
  needle <- if (requireNamespace("ggrepel", quietly = TRUE))
    needle + ggrepel::geom_text_repel(
      data = lab_df, aes(label = label), size = base_size / 4,
      min.segment.length = 0, max.overlaps = Inf, seed = 1)
  else
    needle + geom_text(data = lab_df, aes(label = label), vjust = -0.6,
                       size = base_size / 4, check_overlap = TRUE)

  list(ranking = ranking, map = map, needle = needle)
}

#' Discover candidate regulators, plot, and write files
#'
#' Convenience wrapper: runs \code{\link{discover_regulators}}, builds the plots
#' with \code{\link{plot_regulators}}, and (by default) writes a CSV plus the
#' three PDFs to \code{results/<GENE>_regulators/}.
#'
#' @inheritParams discover_regulators
#' @param out_root Output folder (default \code{"results"}).
#' @param write_files Write the CSV and PDFs (default TRUE).
#' @param base_size,point_size Plot font and marker size.
#' @param chip Add ChIP-Atlas Target Genes support for the top TFs
#'   (\code{\link{add_chip_targets}}; default \code{FALSE}).
#' @param chip_distance 1, 5 or 10 kb around the TSS.
#' @param chip_cell_pattern Optional regular expression restricting ChIP
#'   experiments to cell types.
#' @param chip_max_mb Skip ChIP-Atlas files larger than this (MB).
#' @param expand_top Run the targeted scan (\code{\link{expand_regulators}})
#'   for this many top TFs and write their promoter maps (with ChIP-Atlas track
#'   and rings when \code{chip = TRUE}) to \code{top_TFs/} (default 5; 0 = off).
#' @param expand_threshold PWM threshold for those targeted scans (default 0.80).
#' @param expand_tfs Optional TFs to expand instead of the top ones.
#' @param expand_by How to pick the top TFs: \code{"auto"}, \code{"rank"} or
#'   \code{"chip"} (see \code{\link{expand_regulators}}).
#' @param ... Passed to \code{\link{discover_regulators}}.
#' @return Invisibly, a list with \code{result}, \code{top}, \code{plots}, \code{meta}.
#' @examples
#' \dontrun{
#' run_regulator_discovery("SNCA", top_n = 25)
#' }
#' @export
run_regulator_discovery <- function(gene, species = "homo_sapiens",
                                    upstream = 2500, downstream = 500,
                                    top_n = 20, out_root = "results",
                                    write_files = TRUE,
                                    base_size = 12, point_size = 3,
                                    chip = FALSE, chip_distance = 1,
                                    chip_cell_pattern = NULL, chip_max_mb = 300,
                                    expand_top = 5, expand_threshold = 0.80,
                                    expand_tfs = NULL, expand_by = "auto", ...) {
  disc  <- discover_regulators(gene, species = species, upstream = upstream,
                               downstream = downstream, top_n = top_n, ...)
  if (isTRUE(chip))
    disc <- add_chip_targets(disc, distance = chip_distance,
                             cell_pattern = chip_cell_pattern, max_mb = chip_max_mb)
  plots <- plot_regulators(disc, base_size = base_size, point_size = point_size)
  if (isTRUE(chip))
    plots$chip <- plot_regulators_chip(disc, base_size = base_size, point_size = point_size)
  expanded <- if (expand_top > 0 || length(expand_tfs))
    expand_regulators(disc, tfs = expand_tfs, n = expand_top, by = expand_by,
                      chip = chip, threshold_frac = expand_threshold,
                      chip_cell_pattern = chip_cell_pattern, chip_max_mb = chip_max_mb)
    else list()

  if (write_files) {
    out_dir <- file.path(out_root, paste0(gene, "_regulators"))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    csv_path <- file.path(out_dir, sprintf("%s_candidate_regulators.csv", gene))
    write.csv(disc$result, csv_path, row.names = FALSE)
    message("Wrote: ", csv_path, "  (", nrow(disc$result), " TF x promoter rows)")
    npan <- length(unique(disc$top$promoter_id))
    h <- max(4, 0.32 * max(table(disc$top$promoter_id)) + 1.8)
    w <- if (npan > 1) min(6.5 * npan, 26) else 9
    if (!is.null(plots$ranking))
      ggsave(file.path(out_dir, sprintf("%s_regulators_ranking.pdf", gene)),
             plots$ranking, width = w, height = h, limitsize = FALSE)
    if (!is.null(plots$map))
      ggsave(file.path(out_dir, sprintf("%s_regulators_map.pdf", gene)),
             plots$map, width = w, height = h, limitsize = FALSE)
    if (!is.null(plots$needle))
      ggsave(file.path(out_dir, sprintf("%s_regulators_needle.pdf", gene)),
             plots$needle, width = if (npan > 1) min(7 * npan, 26) else 10, height = 5,
             limitsize = FALSE)
    if (!is.null(plots$chip))
      ggsave(file.path(out_dir, sprintf("%s_regulators_chip.pdf", gene)),
             plots$chip, width = 10, height = 6.5)
    ## targeted scans of the top TFs: one folder each
    for (tf in names(expanded)) {
      s <- expanded[[tf]]
      if (is.null(s) || !nrow(s$result)) {
        message("  top_TFs: ", tf, " skipped (no sites at threshold ", expand_threshold, ").")
        next
      }
      d <- file.path(out_dir, "top_TFs", gsub("[^A-Za-z0-9_.-]+", "-", tf))
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      stem <- file.path(d, sprintf("%s_%s", gene, gsub("[^A-Za-z0-9_.-]+", "-", tf)))
      write.csv(s$result, paste0(stem, "_binding_sites.csv"), row.names = FALSE)
      with_chip <- !is.null(s$chip) && nrow(s$chip$peaks) > 0
      pm <- plot_promoter_map(s, base_size = base_size)
      if (!is.null(pm))
        ggsave(paste0(stem, if (with_chip) "_promoter_map_chip.pdf" else "_promoter_map.pdf"),
               pm, width = 11, height = if (with_chip) 7.5 else 6)
      pb <- plot_binding_sites(s, base_size = base_size)
      if (!is.null(pb$genomic_sig))
        ggsave(paste0(stem, "_genomic_significance.pdf"), pb$genomic_sig, width = 10, height = 5)
    }
    message("Wrote plots to: ", out_dir,
            if (length(expanded)) sprintf(" (+ top_TFs/ for %d TF(s))", length(expanded)) else "")
  }
  message("Done (discovery).")
  invisible(list(result = disc$result, top = disc$top, plots = plots, meta = disc$meta,
                 promoters = disc$promoters, expanded = expanded))
}
