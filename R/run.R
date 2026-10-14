#' Scan, plot, and write binding-site results for a gene + TF
#'
#' Convenience wrapper: runs \code{\link{scan_binding_sites}} then
#' \code{\link{plot_binding_sites}} and \code{\link{plot_promoter_map}}, and (by
#' default) writes the site CSV, a per-promoter summary CSV and the PDF figures
#' to \code{out_root/<GENE>_<TF>/}. With \code{chip = TRUE}, sites are annotated
#' with ChIP-Atlas peaks of the same TF (\code{\link{add_chip_peaks}}).
#'
#' @inheritParams scan_binding_sites
#' @param out_root Output directory root (default \code{"results"}).
#' @param write_files Logical; write CSV + PDFs (default \code{TRUE}). Set
#'   \code{FALSE} to only return objects (as the Shiny app does).
#' @param color_seed Optional integer to fix the database colours reproducibly.
#' @param base_size,point_size Plot font and marker size.
#' @param top_n Number of candidate regulators to return in discovery mode
#'   (used only when \code{tf} is \code{NULL}; see Details).
#' @param chip Add ChIP-Atlas evidence (default \code{FALSE}; needs internet on
#'   first use, files are cached).
#' @param chip_threshold MACS2 q-value class of the peaks: \code{"05"},
#'   \code{"10"} or \code{"20"}.
#' @param chip_cell_pattern Optional regular expression restricting ChIP
#'   experiments to cell types (e.g. \code{"HepG2|liver"}).
#' @param chip_max_mb Skip ChIP-Atlas files larger than this (MB; default 300).
#' @param chip_ring_p,chip_ring_min_exp Which in-peak sites get a black ring on
#'   the promoter map (see \code{\link{plot_promoter_map}}).
#' @details If \code{tf} is \code{NULL} or empty, the function switches to
#'   discovery mode and calls \code{\link{run_regulator_discovery}} to rank the
#'   gene's candidate regulators instead of scanning one named TF.
#' @return Invisibly, a list with \code{result} (data.frame), \code{plots}
#'   (the ggplots, or \code{NULL} if no sites), \code{meta}, \code{promoters},
#'   \code{promoter_summary} and \code{chip}.
#' @examples
#' \dontrun{
#' run_binding_site_analysis("ATP7B", "MTF1", threshold_frac = 0.75)
#' run_binding_site_analysis("TP73", "TP53", transcripts = "all",
#'                           scan_unit = "promoter", chip = TRUE, chip_max_mb = 600)
#' run_binding_site_analysis("SNCA")   # no TF -> discovery mode
#' }
#' @export
run_binding_site_analysis <- function(gene, tf = NULL,
                                      species          = "homo_sapiens",
                                      upstream         = 2500,
                                      downstream       = 500,
                                      threshold_frac   = 0.80,
                                      hocomoco_version = "v12",
                                      hocomoco_species = "human",
                                      transcripts      = c("protein_coding", "all"),
                                      max_transcripts  = 8,
                                      out_root         = "results",
                                      write_files      = TRUE,
                                      color_seed       = NULL,
                                      base_size        = 12,
                                      point_size       = 2.6,
                                      top_n            = 20,
                                      merge_dist       = NULL,
                                      scan_unit        = c("tss", "promoter"),
                                      max_overlap      = 0.2,
                                      chip             = FALSE,
                                      chip_threshold   = "05",
                                      chip_cell_pattern = NULL,
                                      chip_max_mb      = 300,
                                      chip_ring_p      = 1e-3,
                                      chip_ring_min_exp = 1) {
  transcripts <- match.arg(transcripts)
  scan_unit   <- match.arg(scan_unit)

  ## ---- no TF supplied -> DISCOVERY MODE (top candidate regulators) ---------
  if (is.null(tf) || !nzchar(tf)) {
    message("No TF supplied -> discovery mode: ranking candidate regulators of ", gene, ".")
    return(run_regulator_discovery(gene, species = species,
             upstream = upstream, downstream = downstream,
             top_n = top_n, out_root = out_root, write_files = write_files,
             base_size = base_size, point_size = point_size,
             chip = chip, chip_cell_pattern = chip_cell_pattern,
             chip_max_mb = chip_max_mb))
  }

  scan <- scan_binding_sites(gene, tf, species, upstream, downstream, threshold_frac,
                             hocomoco_version, hocomoco_species, transcripts, max_transcripts,
                             merge_dist = merge_dist, scan_unit = scan_unit,
                             max_overlap = max_overlap)
  out_dir  <- file.path(out_root, paste0(gene, "_", tf))
  csv_path <- file.path(out_dir, sprintf("%s_%s_binding_sites.csv", gene, tf))

  if (nrow(scan$result) == 0) {
    if (write_files) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      write.csv(data.frame(), csv_path, row.names = FALSE)
    }
    message("No ", tf, " sites at threshold ", threshold_frac,
            ". Lower threshold_frac or widen the window. No plots drawn.")
    return(invisible(list(result = scan$result, plots = NULL, meta = scan$meta)))
  }

  if (isTRUE(chip))
    scan <- add_chip_peaks(scan, threshold = chip_threshold,
                           cell_pattern = chip_cell_pattern, max_mb = chip_max_mb)

  plots <- plot_binding_sites(scan, base_size = base_size, point_size = point_size,
                              color_seed = color_seed)
  plots$promoter_map <- plot_promoter_map(scan, base_size = base_size,
                                          point_size = point_size, color_seed = color_seed,
                                          chip_ring_p = chip_ring_p,
                                          chip_ring_min_exp = chip_ring_min_exp)
  psum <- summarise_promoters(scan)

  if (write_files) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(scan$result, csv_path, row.names = FALSE)
    message("Wrote: ", csv_path, "  (", nrow(scan$result), " site rows)")
    if (!is.null(plots$canonical))
      ggsave(file.path(out_dir, sprintf("%s_%s_canonical.pdf", gene, tf)),
             plots$canonical, width = 9, height = 4.5)
    ggsave(file.path(out_dir, sprintf("%s_%s_by_transcript.pdf", gene, tf)),
           plots$by_transcript, width = 9,
           height = max(4, 1.4 * scan$meta$n_facets + 1.2), limitsize = FALSE)
    ggsave(file.path(out_dir, sprintf("%s_%s_genomic.pdf", gene, tf)),
           plots$genomic, width = 10, height = 5)
    if (!is.null(plots$canonical_sig))
      ggsave(file.path(out_dir, sprintf("%s_%s_canonical_significance.pdf", gene, tf)),
             plots$canonical_sig, width = 9, height = 4.5)
    if (!is.null(plots$by_transcript_sig))
      ggsave(file.path(out_dir, sprintf("%s_%s_by_transcript_significance.pdf", gene, tf)),
             plots$by_transcript_sig, width = 9,
             height = max(4, 1.4 * scan$meta$n_facets + 1.2), limitsize = FALSE)
    if (!is.null(plots$genomic_sig))
      ggsave(file.path(out_dir, sprintf("%s_%s_genomic_significance.pdf", gene, tf)),
             plots$genomic_sig, width = 10, height = 5)
    if (!is.null(psum))
      write.csv(psum, file.path(out_dir, sprintf("%s_%s_promoters.csv", gene, tf)),
                row.names = FALSE)
    if (!is.null(plots$promoter_map)) {
      with_chip <- !is.null(scan$chip) && nrow(scan$chip$peaks) > 0
      ggsave(file.path(out_dir, sprintf("%s_%s_promoter_map%s.pdf", gene, tf,
                                        if (with_chip) "_chip" else "")),
             plots$promoter_map, width = 11, height = if (with_chip) 7.5 else 6)
    }
    message("Wrote plots to: ", out_dir)
  }
  message("Done.")
  invisible(list(result = scan$result, plots = plots, meta = scan$meta,
                 promoters = scan$promoters, promoter_summary = psum, chip = scan$chip))
}
