#' Scan, plot, and write binding-site results for a gene + TF
#'
#' Convenience wrapper: runs \code{\link{scan_binding_sites}} then
#' \code{\link{plot_binding_sites}}, and (by default) writes a CSV and the
#' score-based and significance PDF figures to \code{out_root/<GENE>_<TF>/}.
#'
#' @inheritParams scan_binding_sites
#' @param out_root Output directory root (default \code{"results"}).
#' @param write_files Logical; write CSV + PDFs (default \code{TRUE}). Set
#'   \code{FALSE} to only return objects (as the Shiny app does).
#' @param color_seed Optional integer to fix the database colours reproducibly.
#' @param base_size,point_size Plot font and marker size.
#' @param top_n Number of candidate regulators to return in discovery mode
#'   (used only when \code{tf} is \code{NULL}; see Details).
#' @details If \code{tf} is \code{NULL} or empty, the function switches to
#'   discovery mode and calls \code{\link{run_regulator_discovery}} to rank the
#'   gene's candidate regulators instead of scanning one named TF.
#' @return Invisibly, a list with \code{result} (data.frame), \code{plots}
#'   (the ggplots, or \code{NULL} if no sites), and \code{meta}.
#' @examples
#' \dontrun{
#' run_binding_site_analysis("ATP7B", "MTF1", threshold_frac = 0.75)
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
                                      top_n            = 20) {
  transcripts <- match.arg(transcripts)

  ## No TF supplied -> discovery mode (rank candidate regulators).
  if (is.null(tf) || !nzchar(tf)) {
    message("No TF supplied -> discovery mode: ranking candidate regulators of ", gene, ".")
    return(run_regulator_discovery(gene, species = species,
             upstream = upstream, downstream = downstream,
             top_n = top_n, out_root = out_root, write_files = write_files,
             base_size = base_size, point_size = point_size))
  }

  scan <- scan_binding_sites(gene, tf, species, upstream, downstream, threshold_frac,
                             hocomoco_version, hocomoco_species, transcripts, max_transcripts)
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

  plots <- plot_binding_sites(scan, base_size = base_size, point_size = point_size,
                              color_seed = color_seed)

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
    message("Wrote plots to: ", out_dir)
  }
  message("Done.")
  invisible(list(result = scan$result, plots = plots, meta = scan$meta))
}
