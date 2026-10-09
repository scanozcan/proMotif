#' Scan a promoter for TF binding sites (data only)
#'
#' Enumerates the gene's transcripts from Ensembl, fetches each distinct
#' promoter window (Ensembl, with a UCSC fallback), and scans both strands for
#' JASPAR and HOCOMOCO matches to \code{tf}. This is the slow (network) step;
#' pass the result to \code{\link{plot_binding_sites}} to build figures.
#'
#' Each hit also gets an analytic per-site p-value (from the PWM score under a
#' promoter-derived, GC-aware background, via \pkg{TFMPvalue}) and a Bonferroni
#' adjusted p-value over the positions tested. These are \code{NA} if
#' \pkg{TFMPvalue} is not installed.
#'
#' @param gene HGNC gene symbol (e.g. \code{"ATP7B"}).
#' @param tf Transcription-factor symbol to scan for (e.g. \code{"MTF1"}).
#' @param species Ensembl species token (default \code{"homo_sapiens"}).
#' @param upstream,downstream Basepairs up/downstream of each TSS (default 2500/500).
#' @param threshold_frac PWM cutoff as a fraction of the matrix min-max score (default 0.80).
#' @param hocomoco_version,hocomoco_species HOCOMOCO bundle (default \code{"v12"}/\code{"human"}).
#' @param transcripts \code{"protein_coding"} (default) or \code{"all"}.
#' @param max_transcripts Cap on distinct TSS windows scanned (canonical always kept; default 8).
#' @return A list with \code{result} (data.frame of predicted sites, including
#'   \code{pvalue}/\code{p_adj}), \code{reps} (the per-TSS representative
#'   transcripts), and \code{meta} (gene/TF/coordinate metadata). \code{result}
#'   has zero rows if no sites pass the threshold.
#' @examples
#' \dontrun{
#' s <- scan_binding_sites("ATP7B", "MTF1")
#' head(s$result)
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
                               max_transcripts  = 8) {
  transcripts <- match.arg(transcripts)
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

  reps <- do.call(rbind, lapply(split(txdf, txdf$tss), function(d) {
    r <- if (any(d$is_canonical)) d[d$is_canonical, ][1, ] else d[1, ]
    r$n_sharing <- nrow(d); r
  }))
  reps <- reps[order(!reps$is_canonical, reps$tss), ]
  if (nrow(reps) > max_transcripts) {
    keep <- union(which(reps$is_canonical), seq_len(max_transcripts))
    reps <- reps[sort(unique(keep))[seq_len(max_transcripts)], ]
  }
  reps$label <- sprintf("%s%s | %s:%s%s", reps$transcript_id,
                        ifelse(reps$is_canonical, " (canonical)", ""),
                        gene_chr, .fmt_bp(reps$tss),
                        ifelse(reps$n_sharing > 1,
                               sprintf("  +%d", reps$n_sharing - 1), ""))
  message(sprintf("  strand %s | chr %s", gene_strand, gene_chr))
  message(sprintf("  %d transcript(s) kept [%s] -> %d distinct TSS scanned",
                  nrow(txdf), transcripts, nrow(reps)))
  if (n_tss_all > nrow(reps))
    message(sprintf("  note: all %d transcripts span %d distinct TSS; use transcripts=\"all\" to include alternative-promoter isoforms.",
                    n_tx_all, n_tss_all))

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
      df$p_adj  <- pmin(1, df$pvalue * n_pos)   # Bonferroni over positions tested
      df$genomic_position <- if (gene_strand == 1) rp$tss + df$position
                             else rp$tss - df$position
      rows[[paste(rp$transcript_id, db)]] <- cbind(
        transcript_id = rp$transcript_id, is_canonical = rp$is_canonical,
        transcript_label = rp$label, chr = gene_chr, gene_strand = gene_strand,
        tss_genomic = rp$tss, database = db, tf = tf,
        matrix_id = pwms[[db]]$id %||% NA_character_,
        matrix_name = pwms[[db]]$name %||% NA_character_,
        df, stringsAsFactors = FALSE)
    }
  }
  result <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(result)) result <- result[order(result$transcript_id, result$position), ]

  list(result = result, reps = reps,
       meta = list(gene = gene, tf = tf, species = species, chr = gene_chr,
                   strand = gene_strand, upstream = upstream, downstream = downstream,
                   threshold_frac = threshold_frac, n_sites = nrow(result),
                   n_facets = nrow(reps)))
}
