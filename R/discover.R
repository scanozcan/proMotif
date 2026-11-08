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
#' @param tss_source `"ensembl"` (default) or `"epd"` (EPDnew mapped TSSs).
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
                                merge_dist = NULL, max_overlap = 0.2,
                                tss_source = c("ensembl", "epd")) {
  rank_by <- match.arg(rank_by); promoters <- match.arg(promoters)
  transcripts <- match.arg(transcripts); tss_source <- match.arg(tss_source)
  message("== discovery: candidate regulators of ", gene, " ==")
  if (!requireNamespace("TFMPvalue", quietly = TRUE))
    stop("Discovery mode needs TFMPvalue: install.packages('TFMPvalue').")

  ## the gene's promoters (TSS groups, as in confirm mode); one window each,
  ## representative TSS = canonical, else the TSS used by most transcripts
  gt <- .gene_transcripts(gene, species, transcripts, tss_source)
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
    pinfo <- .promoter_seq(gene_info, upstream, downstream, species,
                           ensembl = gt$source != "epd", label = pr$promoter_id)
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
                   transcripts = transcripts, tss_source = gt$source,
                   n_promoters = length(res_all),
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
                         tss_source = m$tss_source %||% "ensembl",
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
  if (!is.null(pr)) pr$n_ensembl_tx <- .ensembl_tx_count(pr, m$gene, m$species %||% "homo_sapiens")
  pan_lab <- function(pid) {
    i <- if (is.null(pr)) NA else match(pid, pr$promoter_id)
    ifelse(is.na(i), pid,
           sprintf("%s%s | %s | %s:%s", pid, ifelse(pr$has_canonical[i], " (canonical)", ""),
                   .tss_tx_label(pr[i, , drop = FALSE], pr$n_ensembl_tx[i], short = TRUE), m$chr, .fmt_bp(pr$tss[i])))
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

  base_theme <- theme_bw(base_size = .tsz(base_size)) +
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
      data = lab_df, aes(label = label), size = .tsz(base_size) / 4,
      min.segment.length = 0, max.overlaps = Inf, seed = 1)
  else
    needle + geom_text(data = lab_df, aes(label = label), vjust = -0.6,
                       size = .tsz(base_size) / 4, check_overlap = TRUE)

  list(ranking = ranking, map = map, needle = needle)
}

## ---- Discovery: the most likely regulators of each promoter -------------------
## A TF that regulates many genes can regulate this one too, so regulators are
## not required to be unusual for this promoter. Three steps:
##   1. motif shortlist, per promoter: the TF's best site there is at least as
##      good as in most promoters (background p <= shortlist_bg_p, 0.2) and a
##      real match (FIMO per-site p <= shortlist_p, 1e-4);
##   2. ChIP screen: ChIP-Atlas Target Genes - how many of the TF's experiments
##      have a peak within +/- screen_kb of the gene's TSSs (cached per TF);
##   3. the best n_chip ChIP-bound TFs per promoter (most experiments bound,
##      then background p) and the n_motif best motif-only TFs (no ChIP data or
##      no peak; by background p) are scored with the targeted scan and
##      promoter-resolved ChIP-Atlas peaks; regulators are ranked per promoter
##      by their best proMotif score there.
## Returns one row per TF per promoter (its best site there) with the screen
## and the selection; summarise_candidates() adds the proMotif scores.
.chip_screen <- function(tfs, gene, species, distance = 5, cell_pattern = NULL, max_mb = 300) {
  genome <- .chipatlas_genome(species, coordinates = FALSE)
  if (is.null(genome) || !length(tfs)) return(NULL)
  fix_case <- if (species == "homo_sapiens") toupper else
    function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
  message(sprintf("  ChIP screen: ChIP-Atlas Target Genes (+/- %d kb of the %s TSSs) for %d TF(s); files are cached ...",
                  as.integer(distance), gene, length(tfs)))
  out <- do.call(rbind, lapply(tfs, function(tf) {
    parts <- fix_case(unique(trimws(strsplit(tf, "::", fixed = TRUE)[[1]])))
    rows  <- do.call(rbind, lapply(parts, .chip_target_row, gene = gene, genome = genome,
                                   distance = distance, cell_pattern = cell_pattern,
                                   max_mb = max_mb))
    ok   <- which(rows$chip_available %in% TRUE)
    best <- if (length(ok)) ok[which.max(rows$chip_n_bound[ok])] else 1L
    cbind(tf = tf, rows[best, , drop = FALSE], stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' Select the most likely regulators of each promoter (discovery mode)
#'
#' A TF that regulates many genes can regulate this one too, so regulators are
#' not required to be unusual for the promoter. Three steps: (1) motif shortlist
#' per promoter - the TF's best site there has promoter-background
#' \eqn{p \le} \code{shortlist_bg_p} and FIMO per-site \eqn{p \le}
#' \code{shortlist_p}; (2) ChIP screen - ChIP-Atlas Target Genes, the number of
#' the TF's experiments with a peak within \code{screen_kb} kb of the gene's TSSs
#' (one cached file per TF; TFs screened with the same ChIP-Atlas data, e.g. JUN
#' and FOS::JUN, count once); (3) per promoter the \code{n_chip} TFs bound in
#' most experiments (at least \code{min_bound}) and the \code{n_motif} best
#' motif-only TFs are marked for scoring with the proMotif score.
#'
#' @param disc Result of \code{\link{discover_regulators}}.
#' @param shortlist_bg_p,shortlist_p Motif shortlist thresholds (0.2, 1e-4).
#' @param n_chip,n_motif TFs per promoter selected by ChIP screen and by motif (5, 2).
#' @param min_bound Minimum ChIP-Atlas experiments bound near the gene (2).
#' @param screen Run the ChIP screen (default TRUE).
#' @param screen_kb 1, 5 or 10 kb around the TSSs (Target Genes distance; 5).
#' @param chip_cell_pattern Optional regular expression restricting ChIP
#'   experiments to cell types.
#' @param screen_max_mb Skip Target Genes files larger than this (MB).
#' @param extra_tfs TFs to score in addition ("requested").
#' @return One row per TF per promoter (its best site there) with the motif
#'   evidence, the screen (\code{screen_n_bound}, \code{screen_n_exp}), the
#'   selection (\code{selected}, \code{route}, \code{status}) and the
#'   specificity against the other promoters (\code{pref}, reported only).
#' @seealso \code{\link{summarise_candidates}}, \code{\link{run_regulator_discovery}}
#' @export
select_regulators <- function(disc, shortlist_bg_p = 0.2, shortlist_p = 1e-4,
                              n_chip = 5, n_motif = 2, min_bound = 2,
                              screen = TRUE, screen_kb = 5, chip_cell_pattern = NULL,
                              screen_max_mb = 300, extra_tfs = NULL) {
  r <- disc$result
  if (is.null(r) || !nrow(r)) return(NULL)
  if (is.null(r$promoter_id)) r$promoter_id <- "P1"
  use_bg <- !is.null(r$bg_pvalue) && any(is.finite(r$bg_pvalue))
  r$p_rank <- if (use_bg) ifelse(is.na(r$bg_pvalue), 1, r$bg_pvalue) else r$pvalue
  r$lp <- -log10(pmax(r$p_rank, 1e-300))
  ## specificity is reported, not used: -log10 p here minus the best elsewhere
  best <- tapply(r$lp, list(r$tf, r$promoter_id), max)
  pids <- colnames(best); multi <- length(pids) > 1
  r$lp_other <- NA_real_; r$other_promoter <- NA_character_
  if (multi) for (i in seq_len(nrow(r))) {
    oth <- setdiff(pids, r$promoter_id[i])
    o <- best[r$tf[i], oth]; names(o) <- oth; o[is.na(o)] <- 0
    r$lp_other[i] <- max(o); r$other_promoter[i] <- names(o)[which.max(o)]
  }
  r$pref <- if (multi) r$lp - r$lp_other else NA_real_
  ## one row per TF per promoter (its best site there); every TF keeps its own
  ## row, also when another TF leads the same site
  r <- r[order(r$promoter_id, r$p_rank, r$pvalue), , drop = FALSE]
  r <- r[!duplicated(r[c("promoter_id", "tf")]), , drop = FALSE]
  r$shortlisted <- r$p_rank <= shortlist_bg_p & r$pvalue <= shortlist_p
  ## ChIP screen of the shortlisted TFs (gene-level, ChIP-Atlas Target Genes)
  scr_cols <- c("chip_tf_queried", "chip_available", "chip_n_exp", "chip_n_bound", "chip_frac_exp")
  for (cc in scr_cols) r[[paste0("screen_", sub("^chip_", "", cc))]] <- NA
  if (isTRUE(screen)) {
    sc <- .chip_screen(unique(c(r$tf[r$shortlisted], extra_tfs)), disc$meta$gene,
                       disc$meta$species %||% "homo_sapiens", distance = screen_kb,
                       cell_pattern = chip_cell_pattern, max_mb = screen_max_mb)
    if (!is.null(sc)) {
      i <- match(r$tf, sc$tf)
      r$screen_tf_queried <- sc$chip_tf_queried[i]
      r$screen_available <- sc$chip_available[i]
      r$screen_n_exp     <- sc$chip_n_exp[i]
      r$screen_n_bound   <- sc$chip_n_bound[i]
      r$screen_frac_exp  <- sc$chip_frac_exp[i]
    }
  }
  bound <- r$shortlisted & !is.na(r$screen_n_bound) & r$screen_n_bound >= min_bound
  r$selected <- FALSE; r$route <- NA_character_
  for (p in pids) {
    ib <- which(r$promoter_id == p & bound)
    ib <- ib[order(-r$screen_n_bound[ib], r$p_rank[ib])]
    ## TFs screened with the same ChIP-Atlas data (JUN, JUN::JUNB, FOS::JUN ...)
    ## count once: the one with the best motif at this promoter
    key <- r$screen_tf_queried[ib]
    ib <- ib[order(key, r$p_rank[ib])][!duplicated(key[order(key, r$p_rank[ib])])]
    ib <- ib[order(-r$screen_n_bound[ib], r$p_rank[ib])]
    ib <- utils::head(ib, n_chip)
    im <- which(r$promoter_id == p & r$shortlisted & !bound)
    im <- utils::head(im[order(r$p_rank[im], r$pvalue[im])], n_motif)
    r$selected[c(ib, im)] <- TRUE
    r$route[ib] <- "ChIP-bound"; r$route[im] <- "motif"
  }
  if (length(extra_tfs)) {
    ex <- r$tf %in% extra_tfs & !r$selected
    r$selected[ex] <- TRUE; r$route[ex] <- "requested"
  }
  r$status <- ifelse(r$selected, "scored",
              ifelse(bound, "ChIP-bound, not scored",
              ifelse(r$shortlisted, "shortlisted", "weak motif")))
  r$rank_in_promoter <- stats::ave(r$p_rank, r$promoter_id, FUN = function(x) rank(x, ties.method = "first"))
  attr(r, "rule") <- list(mode = "regulators", shortlist_bg_p = shortlist_bg_p,
                          shortlist_p = shortlist_p, n_chip = n_chip, n_motif = n_motif,
                          min_bound = min_bound, screen = screen, screen_kb = screen_kb,
                          use_bg = use_bg, multi = multi,
                          n_per_promoter = n_chip + n_motif, max_bg_p = shortlist_bg_p,
                          min_pref = NA)
  rownames(r) <- NULL
  r
}

## Main discovery plot: the scored regulators of each promoter, ranked by their
## best proMotif score there. Open circle = the same site's motif-only score,
## filled = proMotif score with ChIP-Atlas evidence (the line between them is
## what ChIP-seq adds). Label: ChIP-Atlas experiments with a peak over the site
## and the motif's background p; "(other promoter)" = binds better elsewhere.
#' Plot the scored regulators of each promoter
#'
#' One panel per promoter; the scored TFs ranked by their best proMotif score
#' there. Open circle = motif-only score of the promoter's best locus, filled =
#' proMotif score with ChIP-Atlas evidence. Labels: ChIP-Atlas experiments with a
#' peak over that site, experiments bound near the gene, and the motif's
#' background p.
#'
#' @param cand Result of \code{\link{select_regulators}} after
#'   \code{\link{summarise_candidates}}.
#' @param disc Result of \code{\link{discover_regulators}}.
#' @param top_n TFs shown per promoter (10).
#' @param base_size,point_size Plot font and marker size.
#' @return A ggplot, or \code{NULL} when nothing was scored.
#' @export
plot_discovery_regulators <- function(cand, disc, top_n = 10, base_size = 12, point_size = 3) {
  d <- cand[cand$selected & !is.na(cand$score_best), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  m <- disc$meta; rule <- attr(cand, "rule")
  d <- do.call(rbind, lapply(split(d, d$promoter_id), function(x)
    utils::head(x[order(-x$score_best, -x$chip_exp_best %||% 0), , drop = FALSE], top_n)))
  pids <- unique(d$promoter_id); pids <- pids[order(as.integer(sub("^P", "", pids)))]
  pl <- .disc_panel_labels(disc, pids)
  d$panel <- factor(pl[d$promoter_id], levels = pl[pids])
  d$ykey <- paste(d$tf, d$promoter_id, sep = "___")
  d$yf <- factor(d$ykey, levels = unique(d$ykey[order(d$panel, d$score_best)]))
  d$motif0 <- ifelse(is.na(d$score_motif_only), d$score_best, pmin(d$score_motif_only, d$score_best))
  d$evid <- ifelse(!is.na(d$chip_exp_best) & d$chip_exp_best > 0, "ChIP-Atlas peak at the site",
                   "motif evidence only")
  d$lab <- sprintf("%s | %s | bg p %s%s",
                   ifelse(d$chip_exp_best %in% NA | d$chip_exp_best == 0, "no peak at site",
                          sprintf("site: %d exp", as.integer(d$chip_exp_best))),
                   ifelse(is.na(d$screen_n_bound), "near gene: no ChIP data",
                          sprintf("near gene: %d/%d exp", as.integer(d$screen_n_bound),
                                  as.integer(d$screen_n_exp))),
                   formatC(d$bg_pvalue, format = "g", digits = 2),
                   ifelse(d$binding_check %in% "binding favours other promoter", " | other promoter", ""))
  sub <- sprintf(paste0("per promoter: TFs with a motif match (background p <= %g, site p <= %g);\n",
                        "scored: the %d with most ChIP-Atlas experiments bound near %s (+/- %d kb) and the %d best motif-only TFs\n",
                        "open circle = motif-only score, filled = proMotif score (motif + ChIP-Atlas) of the TF's best site at that promoter;\n",
                        "label: ChIP-Atlas experiments with a peak over that site | bound near the gene / all experiments of the TF | motif background p"),
                 rule$shortlist_bg_p, rule$shortlist_p, rule$n_chip, m$gene,
                 as.integer(rule$screen_kb), rule$n_motif)
  base_theme <- theme_bw(base_size = .tsz(base_size)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          plot.title = element_text(size = rel(0.95)),
          plot.subtitle = element_text(size = rel(0.72)),
          strip.background = element_rect(fill = "grey92", color = NA),
          legend.position = "bottom")
  ggplot(d, aes(y = yf)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60") +
    geom_segment(aes(x = motif0, xend = score_best, yend = yf), color = "grey60") +
    geom_point(aes(x = motif0), shape = 21, fill = "white", color = "grey35", size = point_size * 0.9) +
    geom_point(aes(x = score_best, color = evid), size = point_size) +
    geom_text(aes(x = 1.02, label = lab), hjust = 0, size = .tsz(base_size) / 4.4, color = "grey20") +
    scale_color_manual(values = c(`ChIP-Atlas peak at the site` = "#C8553D",
                                  `motif evidence only` = "grey50"), name = NULL) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 2),
                       expand = expansion(mult = c(0.01, 0))) +
    facet_wrap(~ panel, nrow = 1, scales = "free_y") +
    labs(x = "proMotif score (best site at this promoter)", y = NULL,
         title = sprintf("Most likely regulators of %s, by promoter", m$gene), subtitle = sub) +
    base_theme
}

## ---- Discovery candidates: one rule per promoter -------------------------------
## For each promoter of the gene, a TF is a candidate when
##   (1) its best site there is unusual for a promoter: background p <= max_bg_p
##       (0.01; the share of 20,094 human promoters with an equal or better site);
##   (2) with more than one promoter, it is specific to this promoter: its
##       -log10 background p here is at least min_pref (1 = 10x) above its best
##       at every other promoter of the gene.
## The n_per_promoter best candidates of each promoter (background p, then FIMO
## p) are then scored with the targeted scan and ChIP-Atlas (proMotif score).
## One row per TF per promoter (the TF's best site there). status:
##   "scored"                 - a top candidate, scored with the proMotif score
##   "candidate, not scored"  - passes (1) and (2), beyond the top n_per_promoter
##   "shared with P<k>"       - passes (1), but as good (< 10x) at promoter P<k>
##   "background p > 0.01"    - the motif match is ordinary for a promoter
## extra_tfs (expand_tfs) are scored at every promoter as well ("requested").
#' Select promoter-specific candidate regulators (discovery mode, rule "specific")
#'
#' For each promoter, a TF is a candidate when its best site there is unusual
#' for a promoter (background \eqn{p \le} \code{max_bg_p}) and, with more than
#' one promoter, at least \code{min_pref} log10 units better than at every other
#' promoter; the \code{n_per_promoter} best per promoter are marked for scoring.
#'
#' @param disc Result of \code{\link{discover_regulators}}.
#' @param n_per_promoter Candidates scored per promoter (3).
#' @param max_bg_p Promoter-background p cut-off (0.01).
#' @param min_pref Specificity: minimum -log10 p difference to every other
#'   promoter (1 = 10x).
#' @param extra_tfs TFs to score in addition ("requested").
#' @return One row per TF per promoter with the rule and \code{status}.
#' @export
select_candidates <- function(disc, n_per_promoter = 3, max_bg_p = 0.01, min_pref = 1,
                              extra_tfs = NULL) {
  r <- disc$result
  if (is.null(r) || !nrow(r)) return(NULL)
  if (is.null(r$promoter_id)) r$promoter_id <- "P1"
  use_bg <- !is.null(r$bg_pvalue) && any(is.finite(r$bg_pvalue))
  r$p_rank <- if (use_bg) ifelse(is.na(r$bg_pvalue), 1, r$bg_pvalue) else r$pvalue
  r$lp <- -log10(pmax(r$p_rank, 1e-300))
  ## best -log10 p of each TF at each promoter (all rows, representative or not)
  best <- tapply(r$lp, list(r$tf, r$promoter_id), max)
  pids <- colnames(best); multi <- length(pids) > 1
  r$lp_other <- NA_real_; r$other_promoter <- NA_character_
  if (multi) for (i in seq_len(nrow(r))) {
    oth <- setdiff(pids, r$promoter_id[i])
    o <- best[r$tf[i], oth]; names(o) <- oth
    o[is.na(o)] <- 0
    r$lp_other[i] <- max(o); r$other_promoter[i] <- names(o)[which.max(o)]
  }
  r$pref     <- if (multi) r$lp - r$lp_other else NA_real_
  r$specific <- if (multi) r$pref >= min_pref else TRUE
  ## one row per TF per promoter: its best site there (representative rows first,
  ## so a TF that leads a shared site keeps the "(+N)" annotation)
  r <- r[order(r$promoter_id, r$p_rank, r$pvalue, !r$is_representative), , drop = FALSE]
  r <- r[!duplicated(r[c("promoter_id", "tf")]), , drop = FALSE]
  r$eligible <- r$p_rank <= max_bg_p & r$specific
  ## top n per promoter among the eligible TFs, one TF per shared site
  r$selected <- FALSE
  for (p in pids) {
    i <- which(r$promoter_id == p & r$eligible)
    i <- i[!duplicated(r$site_group[i])]
    r$selected[utils::head(i, n_per_promoter)] <- TRUE
  }
  r$status <- ifelse(r$selected, "scored",
              ifelse(r$eligible, "candidate, not scored",
              ifelse(r$p_rank <= max_bg_p, paste("shared with", r$other_promoter),
                     sprintf("background p > %g", max_bg_p))))
  if (length(extra_tfs)) {
    ex <- r$tf %in% extra_tfs & !r$selected
    r$status[ex] <- "requested"; r$selected[ex] <- TRUE
  }
  r$rank_in_promoter <- stats::ave(r$p_rank, r$promoter_id, FUN = function(x) rank(x, ties.method = "first"))
  attr(r, "rule") <- list(n_per_promoter = n_per_promoter, max_bg_p = max_bg_p,
                          min_pref = min_pref, use_bg = use_bg, multi = multi)
  rownames(r) <- NULL
  r
}

## Add the targeted results (proMotif score, ChIP-Atlas) of the scored candidates
## to the candidate table: for each scored TF at each promoter, its best binding
## locus there, the motif-only score of that promoter's best locus, the number of
## significant loci, the ChIP-Atlas experiments with a peak over the best locus,
## and a binding check against the TF's other promoters (best score, then
## significant loci; a score difference < 0.05 is a tie).
#' Add proMotif scores to discovery candidates
#'
#' For each scored TF at each promoter: its best binding locus there
#' (\code{score_best}), the motif-only score, the number of significant loci,
#' the ChIP-Atlas experiments with a peak over the best locus, and a binding
#' check against the TF's other promoters.
#'
#' @param cand Result of \code{\link{select_regulators}} or
#'   \code{\link{select_candidates}}.
#' @param expanded Named list of scans (\code{\link{expand_regulators}}).
#' @return \code{cand} with the score columns added.
#' @export
summarise_candidates <- function(cand, expanded) {
  if (is.null(cand) || !nrow(cand)) return(cand)
  new <- c("score_best", "score_motif_only", "n_significant", "chip_exp_best",
           "best_locus", "best_locus_motif", "best_locus_pos", "score_other_best",
           "binding_check")
  for (cc in new) cand[[cc]] <- NA
  for (tf in intersect(names(expanded), unique(cand$tf))) {
    s <- expanded[[tf]]
    st <- if (!is.null(s)) .compat_scores(s)$sites else NULL
    if (is.null(st) || !nrow(st)) next
    inw <- strsplit(as.character(st$in_promoter_windows %||% st$promoter_id), ";", fixed = TRUE)
    at  <- function(p) vapply(inw, function(v) p %in% trimws(v), logical(1))
    pids <- unique(cand$promoter_id)
    per <- lapply(pids, function(p) {
      k <- which(at(p))
      if (!length(k)) return(NULL)
      b <- k[which.max(st$promotif_score[k])]
      list(best = st$promotif_score[b],
           motif = if (!is.null(st$motif_only_score)) max(st$motif_only_score[k], na.rm = TRUE) else NA,
           nsig = as.numeric(sum(st$significant[k] %in% TRUE)),
           chip = if (!is.null(st$chip_n_exp)) st$chip_n_exp[b] else NA,
           id = st$site_id[b], mot = st$site_motif[b], pos = st$position_rel_tss[b])
    })
    names(per) <- pids
    for (i in which(cand$tf == tf & cand$selected)) {
      p <- cand$promoter_id[i]; x <- per[[p]]
      oth <- per[setdiff(pids, p)]; oth <- oth[!vapply(oth, is.null, logical(1))]
      if (!is.null(x)) {
        cand$score_best[i] <- x$best; cand$score_motif_only[i] <- x$motif
        cand$n_significant[i] <- x$nsig; cand$chip_exp_best[i] <- x$chip
        cand$best_locus[i] <- x$id; cand$best_locus_motif[i] <- x$mot
        cand$best_locus_pos[i] <- x$pos
      }
      if (length(oth)) {
        ob <- max(vapply(oth, `[[`, numeric(1), "best"))
        os <- max(vapply(oth, `[[`, numeric(1), "nsig"))
        cand$score_other_best[i] <- ob
        own_b <- x$best %||% 0; own_s <- x$nsig %||% 0
        v <- (if (abs(own_b - ob) < 0.05) 0 else sign(own_b - ob)) + sign(own_s - os)
        cand$binding_check[i] <- if (v > 0) "binding agrees" else
                                 if (v < 0) "binding favours other promoter" else "binding mixed"
      } else if (!is.null(x)) cand$binding_check[i] <- "single promoter"
    }
  }
  cand
}

## Panel label of a promoter in discovery plots: "P1 (canonical) | 1 TSS, 6 tx | 1:3,652,516"
.disc_panel_labels <- function(disc, pids) {
  m <- disc$meta; pr <- disc$promoters
  if (is.null(pr)) return(stats::setNames(pids, pids))
  pr$n_ensembl_tx <- .ensembl_tx_count(pr, m$gene, m$species %||% "homo_sapiens")
  lab <- vapply(pids, function(pid) {
    i <- match(pid, pr$promoter_id)
    if (is.na(i)) return(pid)
    sprintf("%s%s | %s | %s:%s", pid, if (isTRUE(pr$has_canonical[i])) " (canonical)" else "",
            .tss_tx_label(pr[i, , drop = FALSE], pr$n_ensembl_tx[i], short = TRUE),
            m$chr, .fmt_bp(pr$tss[i]))
  }, character(1))
  stats::setNames(lab, pids)
}

## Discovery plot: the candidates of each promoter (one panel per promoter).
## x = -log10 promoter-background p of the TF's best site at that promoter;
## colour = status of the selection rule (see select_candidates); scored
## candidates are labelled with their best proMotif score at that promoter and
## the number of ChIP-Atlas experiments with a peak over that site.
#' Plot the motif shortlist / candidates of each promoter
#'
#' One panel per promoter; x = -log10 promoter-background p of each TF's best
#' site, colour = selection status; scored TFs labelled with their best proMotif
#' score and ChIP-Atlas experiments.
#'
#' @inheritParams plot_discovery_regulators
#' @param top_n TFs shown per promoter (plus every scored TF; 15).
#' @return A ggplot.
#' @export
plot_discovery_candidates <- function(cand, disc, top_n = 15, base_size = 12, point_size = 3) {
  if (is.null(cand) || !nrow(cand)) return(NULL)
  m <- disc$meta; rule <- attr(cand, "rule")
  ## per promoter: the top_n TFs by background p, plus every scored TF
  d <- do.call(rbind, lapply(split(cand, cand$promoter_id), function(x) {
    ## one row per shared site (its representative TF, "(+N)" = the others)
    x <- x[order(x$p_rank, x$pvalue), , drop = FALSE]
    x <- x[x$is_representative %in% TRUE | x$selected, , drop = FALSE]
    keep <- seq_len(nrow(x)) <= top_n | x$selected
    x[keep, , drop = FALSE]
  }))
  pids <- unique(d$promoter_id); pids <- pids[order(as.integer(sub("^P", "", pids)))]
  pl <- .disc_panel_labels(disc, pids)
  d$panel <- factor(pl[d$promoter_id], levels = pl[pids])
  nshare <- d$n_at_site %||% rep(1L, nrow(d))
  d$disp <- ifelse(!is.na(nshare) & nshare > 1, sprintf("%s (+%d)", d$tf, nshare - 1L), d$tf)
  d$ykey <- paste(d$disp, d$promoter_id, sep = "___")
  o <- order(d$panel, d$lp)
  d$yf <- factor(d$ykey, levels = unique(d$ykey[o]))
  ## status colours: scored = red, candidate = orange, shared = grey, ordinary = light
  d$status_lab <- ifelse(grepl("^shared", d$status), "as good at another promoter", d$status)
  regs <- identical(rule$mode, "regulators")
  lev <- if (regs) c("scored", "ChIP-bound, not scored", "shortlisted", "weak motif")
         else c("scored", "requested", "candidate, not scored", "as good at another promoter",
                sprintf("background p > %g", rule$max_bg_p))
  cols <- stats::setNames(if (regs) c("#C8553D", "#E9A03B", "grey55", "grey82")
                          else c("#C8553D", "#8E5572", "#E9A03B", "grey55", "grey82"), lev)
  lev <- lev[lev %in% d$status_lab]
  d$status_lab <- factor(d$status_lab, levels = lev)
  ## proMotif labels of the scored TFs
  sc <- d[d$selected, , drop = FALSE]
  if (nrow(sc)) {
    sc$lab <- ifelse(is.na(sc$score_best), "no site in the targeted scan",
              sprintf("score %.2f%s%s", sc$score_best,
                      ifelse(is.na(sc$chip_exp_best) | sc$chip_exp_best == 0, " | no ChIP peak",
                             sprintf(" | ChIP %d", as.integer(sc$chip_exp_best))),
                      ifelse(sc$binding_check %in% "binding favours other promoter", " | other promoter", "")))
  }
  xmax <- max(d$lp, na.rm = TRUE)
  thr  <- -log10(rule$max_bg_p)
  plab <- if (rule$use_bg) "-log10 promoter-background p (best site at this promoter)" else
                           "-log10 FIMO per-site p (best site)"
  sub <- if (regs) sprintf(paste0("shortlist = background p <= %g and site p <= %g | ChIP-bound = >= %d ChIP-Atlas experiments with a peak near %s | ",
                                  "scored = top %d ChIP-bound + %d motif-only per promoter\n",
                                  "label = best proMotif score at this promoter and ChIP-Atlas experiments with a peak there | (+N) = more TFs share that site | %d motifs, %s background promoters"),
                           rule$shortlist_bg_p, rule$shortlist_p, rule$min_bound, m$gene,
                           rule$n_chip, rule$n_motif, m$n_tested %||% NA,
                           format(m$bg_n %||% 20094, big.mark = ",")) else
         sprintf(paste0("candidate = background p <= %g at this promoter%s | top %d per promoter scored with the proMotif score (motif + ChIP-Atlas)\n",
                        "label = best proMotif score at this promoter and ChIP-Atlas experiments with a peak there | (+N) = more TFs share that site | %d motifs, %s background promoters"),
                 rule$max_bg_p,
                 if (rule$multi) sprintf(" and >= %gx better than at every other promoter", 10^rule$min_pref) else "",
                 rule$n_per_promoter, m$n_tested %||% NA,
                 format(m$bg_n %||% 20094, big.mark = ","))
  base_theme <- theme_bw(base_size = .tsz(base_size)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          plot.title = element_text(size = rel(0.95)),
          plot.subtitle = element_text(size = rel(0.72)),
          strip.background = element_rect(fill = "grey92", color = NA),
          legend.position = "bottom")
  p <- ggplot(d, aes(lp, yf)) +
    geom_vline(xintercept = thr, linetype = "dashed", color = "grey50") +
    geom_segment(aes(x = 0, xend = lp, yend = yf), color = "grey80") +
    geom_point(aes(color = status_lab), size = point_size) +
    scale_color_manual(values = cols[lev], name = NULL, drop = TRUE) +
    scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(0, xmax * 1.9), clip = "off") +
    facet_wrap(~ panel, nrow = 1, scales = "free_y") +
    labs(x = plab, y = NULL,
         title = sprintf(if (regs) "Motif shortlist of %s, by promoter" else "Candidate regulators of %s, by promoter", m$gene),
         subtitle = sub) + base_theme
  if (nrow(sc))
    p <- p + geom_text(data = sc, aes(x = lp, y = yf, label = lab), hjust = -0.08,
                       size = .tsz(base_size) / 4.2, color = "grey15")
  p
}

## Convenience wrapper: discover + select + score + plot + write files.
##   rule = "regulators" (default): the most likely regulators of each promoter,
##          common or not - motif shortlist, ChIP-Atlas Target Genes screen, then
##          proMotif score (see select_regulators); n_chip + n_motif per promoter
##   rule = "specific": TFs whose motif is unusual for promoters and specific to
##          one promoter of the gene (see select_candidates); n_per_promoter
##   chip           : ChIP-Atlas evidence (screen + promoter-resolved peaks in the
##                    targeted scans); chip = FALSE ranks by motif evidence only
##   expand_tfs     : TFs to score in addition to the selected ones
## Writes to <out_root>/<GENE>_regulators/:
##   <GENE>_discovery_regulators.pdf  the scored regulators of each promoter,
##                                    ranked by proMotif score (rule "regulators")
##   <GENE>_discovery_candidates.pdf  the motif shortlist / candidates per promoter
##   <GENE>_discovery_candidates.csv  one row per TF per promoter: motif, screen,
##                                    status, and for scored TFs the proMotif score
##   <GENE>_candidate_regulators.csv  every motif x promoter (best site, p-values)
##   top_TFs/<TF>/                    top-sites plot (proMotif score) and table
#' Discover the regulators of a gene, score them, plot, and write files
#'
#' Runs \code{\link{discover_regulators}}, selects TFs per promoter
#' (\code{rule = "regulators"}: \code{\link{select_regulators}};
#' \code{rule = "specific"}: \code{\link{select_candidates}}), scores them with
#' the targeted scan and ChIP-Atlas (\code{\link{expand_regulators}},
#' \code{\link{summarise_candidates}}), and writes to
#' \code{<out_root>/<GENE>_regulators/}: \code{_discovery_regulators.pdf},
#' \code{_discovery_candidates.pdf/.csv}, \code{_candidate_regulators.csv} and
#' \code{top_TFs/<TF>/} (top-sites plot and table per scored TF).
#'
#' @inheritParams discover_regulators
#' @param top_n Rows per promoter in the candidate plot (15).
#' @param out_root Output folder (default \code{"results"}).
#' @param write_files Write the CSVs and PDFs (default TRUE).
#' @param base_size,point_size Plot font and marker size.
#' @param chip ChIP-Atlas evidence: Target Genes screen and promoter-resolved
#'   peaks in the targeted scans (default FALSE = motif evidence only).
#' @param chip_cell_pattern Optional regular expression restricting ChIP
#'   experiments to cell types.
#' @param chip_max_mb ChIP-Atlas download budget per TF (MB).
#' @param rule \code{"regulators"} (default) or \code{"specific"}.
#' @param n_chip,n_motif,shortlist_bg_p,shortlist_p,min_bound,screen_kb See
#'   \code{\link{select_regulators}}.
#' @param n_per_promoter,max_bg_p,min_pref See \code{\link{select_candidates}}.
#' @param expand_threshold PWM threshold of the targeted scans (0.75).
#' @param expand_tfs TFs to score in addition to the selected ones.
#' @param label_top Labelled sites in the top-sites plots (10).
#' @param color_seed Seed for the database colours of the top-sites plots.
#' @param chip_distance,expand_top,expand_by No longer used (kept for old calls).
#' @param ... Passed to \code{\link{discover_regulators}}.
#' @return Invisibly, a list with \code{result}, \code{candidates}, \code{top},
#'   \code{plots}, \code{meta}, \code{promoters}, \code{expanded}.
#' @examples
#' \dontrun{
#' run_regulator_discovery("TFF1", tss_source = "epd", chip = TRUE)
#' }
#' @export
run_regulator_discovery <- function(gene, species = "homo_sapiens",
                                    upstream = 2500, downstream = 500,
                                    top_n = 15, out_root = "results",
                                    write_files = TRUE,
                                    base_size = 12, point_size = 3,
                                    chip = FALSE, chip_cell_pattern = NULL, chip_max_mb = 300,
                                    rule = c("regulators", "specific"),
                                    n_chip = 5, n_motif = 2, shortlist_bg_p = 0.2,
                                    shortlist_p = 1e-4, min_bound = 2, screen_kb = 5,
                                    n_per_promoter = 3, max_bg_p = 0.01, min_pref = 1,
                                    expand_threshold = 0.75, expand_tfs = NULL,
                                    label_top = 10, color_seed = NULL,
                                    chip_distance = NULL, expand_top = NULL, expand_by = NULL, ...) {
  rule <- match.arg(rule)
  if (!is.null(expand_top))
    message("  expand_top is no longer used: see n_chip / n_motif (rule = 'regulators') or n_per_promoter.")
  disc <- discover_regulators(gene, species = species, upstream = upstream,
                              downstream = downstream, top_n = top_n, ...)
  cand <- if (rule == "regulators")
    select_regulators(disc, shortlist_bg_p = shortlist_bg_p, shortlist_p = shortlist_p,
                      n_chip = if (isTRUE(chip)) n_chip else 0,
                      n_motif = if (isTRUE(chip)) n_motif else n_chip + n_motif,
                      min_bound = min_bound, screen = isTRUE(chip), screen_kb = screen_kb,
                      chip_cell_pattern = chip_cell_pattern, screen_max_mb = chip_max_mb,
                      extra_tfs = expand_tfs)
  else select_candidates(disc, n_per_promoter = n_per_promoter, max_bg_p = max_bg_p,
                         min_pref = min_pref, extra_tfs = expand_tfs)
  sel <- cand[cand$selected, , drop = FALSE]
  for (p in unique(cand$promoter_id))
    message(sprintf("  %s candidates: %s", p,
                    if (any(sel$promoter_id == p)) paste(sel$tf[sel$promoter_id == p], collapse = ", ")
                    else "none pass the rule"))
  tfs <- unique(c(sel$tf, expand_tfs))
  expanded <- if (length(tfs))
    expand_regulators(disc, tfs = tfs, chip = chip, threshold_frac = expand_threshold,
                      chip_cell_pattern = chip_cell_pattern, chip_max_mb = chip_max_mb)
    else list()
  cand  <- summarise_candidates(cand, expanded)
  plots <- list(regulators = plot_discovery_regulators(cand, disc, base_size = base_size,
                                                       point_size = point_size),
                candidates = plot_discovery_candidates(cand, disc, top_n = top_n,
                                                       base_size = base_size,
                                                       point_size = point_size))
  if (write_files) {
    out_dir <- file.path(out_root, paste0(gene, "_regulators"))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ## outputs of the old discovery plots
    unlink(file.path(out_dir, sprintf("%s_regulators_%s.pdf", gene,
                                      c("ranking", "map", "needle", "chip"))))
    write.csv(disc$result, file.path(out_dir, sprintf("%s_candidate_regulators.csv", gene)),
              row.names = FALSE)
    ct <- cand; ct$p_rank <- NULL; ct$lp <- NULL; ct$eligible <- NULL
    write.csv(ct, file.path(out_dir, sprintf("%s_discovery_candidates.csv", gene)), row.names = FALSE)
    message("Wrote: ", file.path(out_dir, sprintf("%s_discovery_candidates.csv", gene)),
            sprintf("  (%d TF x promoter rows, %d scored)", nrow(cand), sum(cand$selected)))
    npan <- length(unique(cand$promoter_id))
    if (!is.null(plots$regulators)) {
      nreg <- max(table(cand$promoter_id[cand$selected]))
      .ggsave_safe(file.path(out_dir, sprintf("%s_discovery_regulators.pdf", gene)), plots$regulators,
                   width = if (npan > 1) min(7.5 * npan, 30) else 12,
                   height = max(4.5, 0.38 * nreg + 3))
    }
    if (!is.null(plots$candidates)) {
      nrow_max <- max(table(cand$promoter_id[cand$rank_in_promoter <= top_n | cand$selected]))
      .ggsave_safe(file.path(out_dir, sprintf("%s_discovery_candidates.pdf", gene)), plots$candidates,
                   width = if (npan > 1) min(7.5 * npan, 30) else 12,
                   height = max(5, 0.3 * nrow_max + 3))
    }
    ## the scored candidates: top-sites plot (proMotif score) and table each
    if (is.null(color_seed)) color_seed <- sample.int(1e6, 1)
    unlink(file.path(out_dir, "top_TFs"), recursive = TRUE)
    for (tf in names(expanded)) {
      s <- expanded[[tf]]
      if (is.null(s) || !nrow(s$result)) {
        message("  top_TFs: ", tf, " skipped (no sites at threshold ", expand_threshold, ").")
        next
      }
      d <- file.path(out_dir, "top_TFs", gsub("[^A-Za-z0-9_.-]+", "-", tf))
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      stem <- file.path(d, sprintf("%s_%s", gene, gsub("[^A-Za-z0-9_.-]+", "-", tf)))
      write.csv(top_binding_sites(s), paste0(stem, "_top_sites.csv"), row.names = FALSE)
      with_chip <- !is.null(s$chip) && nrow(s$chip$peaks) > 0
      tp <- plot_top_sites(s, top_n = label_top, base_size = base_size, color_seed = color_seed)
      if (!is.null(tp))
        .ggsave_safe(paste0(stem, "_top_sites.pdf"), tp, width = 12, height = if (with_chip) 8 else 6.5)
    }
    message("Wrote plots to: ", out_dir,
            if (length(expanded)) sprintf(" (+ top_TFs/ for %d TF(s))", length(expanded)) else "")
  }
  message("Done (discovery).")
  invisible(list(result = disc$result, candidates = cand, top = disc$top, plots = plots,
                 meta = disc$meta, promoters = disc$promoters, expanded = expanded))
}
