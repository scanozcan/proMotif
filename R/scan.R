#' Scan a promoter for TF binding sites (data only)
#'
#' Enumerates the gene's transcripts from Ensembl, groups their TSSs into
#' promoters, fetches each promoter window (Ensembl, with a UCSC fallback), and
#' scans both strands for JASPAR and HOCOMOCO matches to \code{tf}. This is the
#' slow (network) step; pass the result to \code{\link{plot_binding_sites}} or
#' \code{\link{plot_promoter_map}} to build figures.
#'
#' Each hit gets a FIMO-style per-site p-value (PWM score under a
#' promoter-derived, GC-aware background, via \pkg{TFMPvalue}) and a
#' Benjamini-Hochberg q-value over all positions tested. Both are \code{NA} if
#' \pkg{TFMPvalue} is not installed.
#'
#' Promoters: TSSs whose scan windows overlap by at least \code{max_overlap}
#' (i.e. lie within \code{(1 - max_overlap) * window} bp) form one promoter,
#' numbered P1, P2, ... 5' to 3'.
#'
#' @param gene HGNC gene symbol (e.g. \code{"ATP7B"}).
#' @param tf Transcription-factor symbol to scan for (e.g. \code{"MTF1"}).
#' @param species Ensembl species token (default \code{"homo_sapiens"}).
#' @param upstream,downstream Basepairs up/downstream of each TSS (default 2500/500).
#' @param threshold_frac PWM cutoff as a fraction of the matrix min-max score (default 0.80).
#' @param hocomoco_version,hocomoco_species HOCOMOCO bundle (default \code{"v12"}/\code{"human"}).
#' @param transcripts \code{"protein_coding"} (default) or \code{"all"}.
#' @param max_transcripts Cap on TSS windows scanned (default 8). The canonical
#'   TSS is kept first, then one TSS per promoter.
#' @param merge_dist TSSs closer than this (bp) form one promoter. Default
#'   \code{NULL} = derived from \code{max_overlap} and the window.
#' @param scan_unit \code{"tss"} (default; one window per distinct TSS) or
#'   \code{"promoter"} (one window per promoter).
#' @param max_overlap Window overlap fraction at or above which TSSs belong to
#'   the same promoter (default 0.2).
#' @return A list with \code{result} (data.frame of predicted sites, including
#'   \code{pvalue}, \code{qvalue}, \code{promoter_id}, \code{in_promoter_windows}),
#'   \code{reps} (scanned representative transcripts), \code{tss_all} (every
#'   TSS), \code{promoters}, \code{promoter_overlap} and \code{meta}.
#'   \code{result} has zero rows if no sites pass the threshold.
#' @examples
#' \dontrun{
#' s <- scan_binding_sites("ATP7B", "MTF1")
#' head(s$result)
#' s$promoters
#' }
#' @export
scan_binding_sites <- function(gene, tf,
                               species          = "homo_sapiens",
                               upstream         = 2500,
                               downstream       = 500,
                               threshold_frac   = 0.80,
                               hocomoco_version = "v12",
                               hocomoco_species = "human",
                               transcripts      = c("protein_coding", "all"),
                               max_transcripts  = 8,
                               merge_dist       = NULL,
                               scan_unit        = c("tss", "promoter"),
                               max_overlap      = 0.2) {
  transcripts <- match.arg(transcripts)
  scan_unit   <- match.arg(scan_unit)
  message("== ", gene, " / ", tf, " ==")

  js <- .ensembl_lookup_expand(gene, species)
  gene_chr <- js$seq_region_name; gene_strand <- js$strand
  txs <- js$Transcript %||% list()
  if (length(txs) == 0) stop("No transcripts returned for ", gene)
  txdf <- do.call(rbind, lapply(txs, function(tx) data.frame(
    transcript_id = tx$id %||% NA_character_,
    biotype       = tx$biotype %||% NA_character_,
    is_canonical  = isTRUE((tx$is_canonical %||% 0) == 1),
    start         = as.numeric(tx$start), end = as.numeric(tx$end),
    stringsAsFactors = FALSE)))
  txdf$tss <- if (gene_strand == 1) txdf$start else txdf$end   # per-transcript 5' end
  n_tx_all <- nrow(txdf); n_tss_all <- length(unique(txdf$tss))
  if (transcripts == "protein_coding" &&
      any(txdf$biotype == "protein_coding", na.rm = TRUE))
    txdf <- txdf[txdf$biotype == "protein_coding", , drop = FALSE]

  ## group nearby TSSs into promoters P1, P2, ... (section 5). scan_unit = "tss"
  ## scans every distinct TSS (original behaviour); "promoter" scans one window
  ## per promoter (canonical TSS, else the TSS used by most transcripts).
  merge_dist <- merge_dist %||% round((1 - max_overlap) * (upstream + downstream + 1))
  txdf <- .cluster_tss(txdf, gene_strand, merge_dist)
  unit <- if (scan_unit == "promoter") txdf$promoter_id else txdf$tss
  reps <- do.call(rbind, lapply(split(txdf, unit), function(d) {
    r <- if (any(d$is_canonical)) d[d$is_canonical, ][1, ] else {
      tt <- table(d$tss); d[d$tss == as.numeric(names(tt)[which.max(tt)]), ][1, ] }
    r$n_sharing <- nrow(d); r
  }))
  ## max_transcripts cap. Priority: canonical TSS, then one TSS per promoter (its
  ## most-used TSS), then the remaining TSSs by number of transcripts - so no
  ## promoter is dropped just because of its coordinate (the old rule kept the
  ## lowest coordinates, which on the - strand are the most downstream TSSs).
  reps <- reps[order(!reps$is_canonical, -reps$n_sharing, reps$tss), ]
  pri  <- ifelse(reps$is_canonical, 0L, ifelse(!duplicated(reps$promoter_id), 1L, 2L))
  reps <- reps[order(pri, -reps$n_sharing, reps$tss), ]
  n_dropped <- max(0L, nrow(reps) - max_transcripts)
  reps <- head(reps, max_transcripts)
  reps <- reps[order(!reps$is_canonical, reps$tss), ]
  if (n_dropped)
    message(sprintf("  max_transcripts = %d: %d TSS window(s) not scanned (each promoter keeps >= 1 scanned TSS while max_transcripts allows; raise max_transcripts to scan all).",
                    as.integer(max_transcripts), n_dropped))
  reps$label <- sprintf("%s: %s%s | %s:%s%s", reps$promoter_id, reps$transcript_id,
                        ifelse(reps$is_canonical, " (canonical)", ""),
                        gene_chr, .fmt_bp(reps$tss),
                        ifelse(reps$n_sharing > 1,
                               sprintf("  +%d", reps$n_sharing - 1), ""))
  message(sprintf("  strand %s | chr %s", gene_strand, gene_chr))
  message(sprintf("  %d transcript(s) kept [%s] -> %d distinct TSS scanned",
                  nrow(txdf), transcripts, nrow(reps)))
  if (transcripts == "protein_coding" && n_tss_all > nrow(reps))
    message(sprintf("  note: all %d transcripts span %d distinct TSS; use transcripts=\"all\" to include alternative-promoter isoforms.",
                    n_tx_all, n_tss_all))
  prom <- .promoter_table(txdf, gene_strand, upstream, downstream, max_overlap)
  .report_promoters(prom$promoters, gene_chr, merge_dist, max_overlap, scan_unit)

  pwms <- list(
    JASPAR   = tryCatch(GLproxScape::fetch_jaspar_pwm(tf), error = function(e) NULL),
    HOCOMOCO = tryCatch(GLproxScape::fetch_hocomoco_pwm(tf, version = hocomoco_version,
                                           species = hocomoco_species),
                        error = function(e) NULL))
  for (nm in names(pwms))
    message("  ", nm, ": ",
            if (is.null(pwms[[nm]])) paste0("no ", tf, " matrix found")
            else sprintf("matrix %s (%s)", pwms[[nm]]$id %||% "NA",
                         pwms[[nm]]$name %||% "NA"))
  if (all(vapply(pwms, is.null, logical(1))))
    stop("No ", tf, " matrix in either database - nothing to scan.")

  rows <- list()
  for (i in seq_len(nrow(reps))) {
    rp <- reps[i, ]
    gene_info <- list(name = gene, chr = gene_chr, strand = gene_strand,
                      tss = rp$tss, start = rp$start, end = rp$end,
                      species = species, transcript_id = rp$transcript_id)
    pinfo <- .with_retry(
      function() GLproxScape::fetch_promoter_seq(gene_info, upstream = upstream,
                                    downstream = downstream),
      tries = 2, waits = c(3, 8),
      what = paste("Ensembl seq", rp$transcript_id))
    if (is.null(pinfo)) {
      pinfo <- .fetch_seq_ucsc(gene_chr, rp$tss, gene_strand,
                               upstream, downstream, species)
      if (!is.null(pinfo)) message("  used UCSC fallback for ", rp$transcript_id)
    }
    if (is.null(pinfo)) {
      message("  could not fetch sequence for ", rp$transcript_id,
              " (Ensembl and UCSC both failed) - skipping"); next
    }
    bg <- .seq_bg(pinfo$seq)
    for (db in names(pwms)) {
      if (is.null(pwms[[db]])) next
      df <- .scan_promoter(pinfo, pwms[[db]], threshold_frac)
      if (nrow(df) == 0) next
      n_pos <- 2L * max(1L, nchar(pinfo$seq) - pwms[[db]]$len + 1L)
      df$pvalue <- .site_pvalues(df$score, pwms[[db]]$pwm, bg)
      df$qvalue <- .bh_qvalue(df$pvalue, n_pos)  # FIMO-style BH/FDR over positions tested
      df$genomic_position <- if (gene_strand == 1) rp$tss + df$position
                             else rp$tss - df$position
      rows[[paste(rp$transcript_id, db)]] <- cbind(
        transcript_id = rp$transcript_id, promoter_id = rp$promoter_id,
        is_canonical = rp$is_canonical, transcript_label = rp$label, chr = gene_chr, gene_strand = gene_strand,
        tss_genomic = rp$tss, database = db, tf = tf,
        matrix_id = pwms[[db]]$id %||% NA_character_,
        matrix_name = pwms[[db]]$name %||% NA_character_,
        df, stringsAsFactors = FALSE)
    }
  }
  result <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(result)) result <- result[order(result$transcript_id, result$position), ]
  ## which promoters' windows contain each site ("P1", "P1;P2", ...)
  if (nrow(result))
    result$in_promoter_windows <- .covering_promoters(result$genomic_position,
                                                      prom$promoters)

  ## every TSS of the kept transcripts (for drawing all TSS lines)
  tss_all <- do.call(rbind, lapply(split(txdf, txdf$tss), function(d) data.frame(
    tss = d$tss[1], promoter_id = d$promoter_id[1], is_canonical = any(d$is_canonical),
    n_transcripts = nrow(d), stringsAsFactors = FALSE)))
  tss_all$scanned <- tss_all$tss %in% reps$tss
  rownames(tss_all) <- NULL

  list(result = result, reps = reps, tss_all = tss_all,
       promoters = prom$promoters, promoter_overlap = prom$overlap,
       meta = list(gene = gene, tf = tf, species = species, chr = gene_chr,
                   strand = gene_strand, upstream = upstream, downstream = downstream,
                   threshold_frac = threshold_frac, n_sites = nrow(result),
                   n_facets = nrow(reps), merge_dist = merge_dist,
                   scan_unit = scan_unit, max_overlap = max_overlap,
                   n_promoters = nrow(prom$promoters),
                   n_regions = length(unique(prom$promoters$region))))
}
