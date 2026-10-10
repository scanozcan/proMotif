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
#' Scans the canonical promoter of \code{gene} against the whole JASPAR CORE
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
#' @param jaspar_db Optional pre-opened JASPAR database/connection to use
#'   instead of auto-detecting an installed JASPAR package.
#'
#' @return A list with \code{result} (the full per-TF table: best-site position,
#'   strand, score, \code{score_frac}, \code{pvalue}, \code{qvalue},
#'   \code{genomic_position}, \code{family}, \code{site_group},
#'   \code{is_representative}, \code{n_at_site}, \code{shares_site_with}),
#'   \code{top} (the \code{top_n} representatives) and \code{meta}.
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
                                use_motifmatchr = TRUE, jaspar_db = NULL) {
  message("== discovery: candidate regulators of ", gene, " ==")
  if (!requireNamespace("TFMPvalue", quietly = TRUE))
    stop("Discovery mode ranks by the FIMO p-value and needs TFMPvalue: ",
         "install.packages('TFMPvalue').")
  cp   <- .canonical_promoter(gene, species, upstream, downstream)
  seq  <- cp$pinfo$seq; tss_i <- cp$pinfo$tss_offset + 1
  bg   <- .seq_bg(seq)
  pwms <- .jaspar_core_pwms(collection, tax_group, jaspar_db)

  ## Optional fast shortlist: which motifs have any match at p < 5e-4.
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
      dna <- Biostrings::DNAStringSet(seq)
      mm  <- motifmatchr::matchMotifs(pwmlist, dna, bg = bg,
                                      out = "matches", p.cutoff = 5e-4)
      which(as.logical(motifmatchr::motifMatches(mm)[1, ]))
    }, error = function(e) {
      message("  motifmatchr shortlist skipped: ", conditionMessage(e)); NULL })
    if (length(sl)) { keep_idx <- sl
      message("  motifmatchr shortlisted ", length(sl), " motifs.") }
  }

  n_tested <- length(keep_idx)
  rows <- lapply(keep_idx, function(k) {
    p  <- pwms[[k]]
    bs <- .best_site(seq, tss_i, p)
    if (is.null(bs)) return(NULL)
    pv <- .fimo_pvalue(p$pwm, bs$score, bg)
    data.frame(tf = p$name %||% p$id, matrix_id = p$id %||% NA_character_,
               family = p$family, strand = bs$strand,
               position = bs$position, score = bs$score, score_frac = bs$score_frac,
               pvalue = pv, site = bs$site,
               genomic_position = if (cp$strand == 1) cp$tss + bs$position
                                  else cp$tss - bs$position,
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)
  if (is.null(res) || !nrow(res)) stop("No motif matches found for ", gene, ".")

  ## One row per TF name (best matrix), best-first.
  res <- res[order(res$pvalue, -res$score_frac), ]
  res <- res[!duplicated(res$tf), ]
  res$n_at_site <- 1L
  res$shares_site_with <- ""
  res$is_representative <- TRUE
  res$site_group <- seq_len(nrow(res))

  ## Group TFs whose best site is the SAME locus (within site_window bp). Nothing
  ## is dropped: the full table keeps a site_group id and an is_representative
  ## flag; the representative (smallest p) also lists the folded members.
  if (collapse_site && nrow(res) > 1) {
    grp <- integer(nrow(res)); reps_pos <- numeric(0)
    for (i in seq_len(nrow(res))) {
      gp  <- res$genomic_position[i]
      hit <- which(abs(reps_pos - gp) <= site_window)
      if (length(hit)) grp[i] <- hit[1]
      else { reps_pos <- c(reps_pos, gp); grp[i] <- length(reps_pos) }
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
  res <- res[order(res$pvalue, -res$score_frac), ]
  rownames(res) <- NULL
  top <- head(res[res$is_representative, , drop = FALSE], top_n)

  list(result = res, top = top,
       meta = list(gene = gene, species = species, chr = cp$chr, strand = cp$strand,
                   tss = cp$tss, transcript_id = cp$transcript_id,
                   upstream = upstream, downstream = downstream,
                   n_tested = n_tested, top_n = top_n, collection = collection,
                   tax_group = tax_group))
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
  top$neglogp <- -log10(pmax(top$pvalue, 1e-300))
  nshare <- if (!is.null(top$n_at_site)) top$n_at_site else rep(1L, nrow(top))
  disp   <- ifelse(nshare > 1, sprintf("%s (+%d)", top$tf, nshare - 1L), as.character(top$tf))
  top$tf_disp <- factor(disp, levels = disp[order(top$neglogp)])
  famlab <- ifelse(is.na(top$family) | top$family == "", "unassigned", top$family)
  top$family_lab <- factor(famlab)
  base_theme <- theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = rel(0.95)),
          plot.subtitle = element_text(size = rel(0.78)))
  sub <- sprintf("canonical promoter -%d/+%d bp | %d motifs tested | FIMO per-site p | match-priority, not binding | (+N) = more TFs share that site",
                 m$upstream, m$downstream, m$n_tested)

  ranking <- ggplot(top, aes(neglogp, tf_disp)) +
    geom_segment(aes(x = 0, xend = neglogp, yend = tf_disp), color = "grey70") +
    geom_point(aes(color = family_lab), size = point_size) +
    labs(x = "-log10 FIMO per-site p (best site)", y = NULL, color = "TF family",
         title = sprintf("Top candidate regulators of %s", m$gene), subtitle = sub) +
    base_theme

  map <- ggplot(top, aes(position, tf_disp)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    geom_point(aes(color = family_lab, shape = strand), size = point_size) +
    scale_shape_manual(values = c("+" = 16, "-" = 17)) +
    labs(x = "Best-site position relative to TSS (bp)", y = NULL,
         color = "TF family", shape = "Motif strand",
         title = sprintf("Where the top %s regulators' best sites sit", m$gene),
         subtitle = sub) + base_theme

  needle <- ggplot(top, aes(position, neglogp)) +
    geom_segment(aes(xend = position, yend = 0), color = "grey75", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(aes(color = family_lab, shape = strand), size = point_size) +
    scale_shape_manual(values = c("+" = 16, "-" = 17)) +
    labs(x = "Position relative to TSS (bp)", y = "-log10 FIMO per-site p",
         color = "TF family", shape = "Motif strand",
         title = sprintf("Where %s's candidate regulators sit", m$gene),
         subtitle = sub) + base_theme
  lab_df <- head(top[order(-top$neglogp), , drop = FALSE], 10)
  needle <- if (requireNamespace("ggrepel", quietly = TRUE))
    needle + ggrepel::geom_text_repel(
      data = lab_df, aes(label = as.character(tf_disp)), size = base_size / 4,
      min.segment.length = 0, max.overlaps = Inf, seed = 1)
  else
    needle + geom_text(
      data = lab_df, aes(label = as.character(tf_disp)),
      vjust = -0.6, size = base_size / 4, check_overlap = TRUE)

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
                                    base_size = 12, point_size = 3, ...) {
  disc  <- discover_regulators(gene, species = species, upstream = upstream,
                               downstream = downstream, top_n = top_n, ...)
  plots <- plot_regulators(disc, base_size = base_size, point_size = point_size)
  if (write_files) {
    out_dir <- file.path(out_root, paste0(gene, "_regulators"))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    csv_path <- file.path(out_dir, sprintf("%s_candidate_regulators.csv", gene))
    utils::write.csv(disc$result, csv_path, row.names = FALSE)
    message("Wrote: ", csv_path, "  (", nrow(disc$result), " TFs)")
    h <- max(4, 0.32 * nrow(disc$top) + 1.5)
    if (!is.null(plots$ranking))
      ggsave(file.path(out_dir, sprintf("%s_regulators_ranking.pdf", gene)),
             plots$ranking, width = 9, height = h, limitsize = FALSE)
    if (!is.null(plots$map))
      ggsave(file.path(out_dir, sprintf("%s_regulators_map.pdf", gene)),
             plots$map, width = 9, height = h, limitsize = FALSE)
    if (!is.null(plots$needle))
      ggsave(file.path(out_dir, sprintf("%s_regulators_needle.pdf", gene)),
             plots$needle, width = 10, height = 5)
    message("Wrote plots to: ", out_dir)
  }
  message("Done (discovery).")
  invisible(list(result = disc$result, top = disc$top, plots = plots, meta = disc$meta))
}
