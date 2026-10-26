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
#' @param usage Add FANTOM5 CAGE promoter usage across cell types and the
#'   normal vs cancer comparison (default \code{FALSE}; needs the FANTOM5
#'   folder, see \code{\link{add_promoter_usage}}).
#' @param usage_flank Window (bp) around each promoter's TSSs for CAGE peaks.
#' @details If \code{tf} is \code{NULL} or empty, the function switches to
#'   discovery mode and calls \code{\link{run_regulator_discovery}} to rank the
#'   gene's candidate regulators instead of scanning one named TF.
#' @param chip_mode `"auto"` (default): the merged ChIP-Atlas file if it is within `chip_max_mb`, else the newest per-experiment files up to `chip_max_mb`; `"assembled"`: merged file only; `"experiments"`: per-experiment files only.
#' @param chip_max_experiments Maximum number of ChIP-Atlas experiments (default `Inf`).
#' @param chip_sigma Summit-proximity scale in bp for ChIP support (default 50).
#' @param chip_support Normalisation of ChIP support: `"locus"` (experiments binding the locus, default) or `"all"` (all queried experiments).
#' @param tss_source `"ensembl"` (default) or `"epd"`; see [scan_binding_sites()].
#' @param flank Bases on each side of a site in the sequence context (default 10).
#' @param top_sites_min_prob Binding probability for the top-sites table (default 0.5).
#' @param top_sites_min_n Minimum number of sites in the top-sites table (default 10).
#' @param label_top Number of sites labelled in the top-sites plot (default 10).
#' @return Invisibly, a list with \code{result} (data.frame), \code{plots}
#'   (the ggplots, or \code{NULL} if no sites), \code{meta}, \code{promoters},
#'   \code{promoter_summary}, \code{chip} and \code{usage}.
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
                                      hocomoco_version = "v14",
                                      jaspar_release   = "2026",
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
                                      chip_mode        = c("auto", "assembled", "experiments"),
                                      chip_max_experiments = Inf,
                                      chip_sigma       = 50,
                                      chip_support     = c("locus", "all"),
                                      chip_ring_p      = 1e-3,
                                      chip_ring_min_exp = 1,
                                      usage            = FALSE,
                                      usage_flank      = 500,
                                      tss_source       = c("ensembl", "epd"),
                                      flank            = 10,
                                      top_sites_min_prob = 0.5,
                                      top_sites_min_n  = 10,
                                      label_top        = 10) {
  chip_mode <- match.arg(chip_mode); tss_source <- match.arg(tss_source)
  ## one colour seed per run, so every plot of this run uses the same two
  ## database colours (color_seed = a number makes them reproducible across runs)
  if (is.null(color_seed)) color_seed <- sample.int(1e6, 1)
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
             chip_max_mb = chip_max_mb, jaspar_release = jaspar_release,
             transcripts = transcripts, merge_dist = merge_dist, max_overlap = max_overlap))
  }

  scan <- scan_binding_sites(gene, tf, species, upstream, downstream, threshold_frac,
                             hocomoco_version, hocomoco_species, transcripts, max_transcripts,
                             merge_dist = merge_dist, scan_unit = scan_unit,
                             max_overlap = max_overlap, jaspar_release = jaspar_release,
                             tss_source = tss_source, flank = flank)
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
                           cell_pattern = chip_cell_pattern, max_mb = chip_max_mb,
                           mode = chip_mode, max_experiments = chip_max_experiments,
                           sigma = chip_sigma, support = match.arg(chip_support))

  if (isTRUE(usage)) scan <- add_promoter_usage(scan, flank = usage_flank)
  top_sites <- top_binding_sites(scan, min_prob = top_sites_min_prob, min_n = top_sites_min_n)

  plots <- plot_binding_sites(scan, base_size = base_size, point_size = point_size,
                              color_seed = color_seed)
  if (!is.null(scan$usage)) {
    plots$usage        <- plot_promoter_usage(scan, base_size = base_size)
    plots$usage_cancer <- plot_promoter_usage_cancer(scan, base_size = base_size)
  }
  plots$promoter_map <- plot_promoter_map(scan, base_size = base_size,
                                          point_size = point_size, color_seed = color_seed,
                                          chip_ring_p = chip_ring_p,
                                          chip_ring_min_exp = chip_ring_min_exp)
  plots$top_sites <- plot_top_sites(scan, top_n = label_top, base_size = base_size,
                                    point_size = point_size, color_seed = color_seed)
  plots$top_sites_nolabels <- plot_top_sites(scan, top_n = 0, base_size = base_size,
                                             point_size = point_size, color_seed = color_seed)
  plots$motif_logo <- plot_motif_logo(scan, base_size = base_size)
  psum <- summarise_promoters(scan)

  if (write_files) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(scan$result, csv_path, row.names = FALSE)
    message("Wrote: ", csv_path, "  (", nrow(scan$result), " site rows)")
    ggsave(file.path(out_dir, sprintf("%s_%s_genomic.pdf", gene, tf)),
           plots$genomic, width = 10, height = 5)
    if (!is.null(psum))
      write.csv(psum, file.path(out_dir, sprintf("%s_%s_promoters.csv", gene, tf)),
                row.names = FALSE)
    if (!is.null(scan$usage)) {
      write.csv(scan$usage$table, file.path(out_dir, sprintf("%s_%s_promoter_usage.csv", gene, tf)),
                row.names = FALSE)
      if (!is.null(plots$usage))
        ggsave(file.path(out_dir, sprintf("%s_%s_promoter_usage.pdf", gene, tf)), plots$usage,
               width = 11, height = max(5, 0.27 * length(unique(plots$usage$data$label)) + 2.2),
               limitsize = FALSE)
      if (!is.null(plots$usage_cancer))
        ggsave(file.path(out_dir, sprintf("%s_%s_promoter_usage_normal_vs_cancer.pdf", gene, tf)),
               plots$usage_cancer, width = 12,
               height = 2.6 * length(unique(plots$usage_cancer$data$panel)) + 1.8)
    }
    with_chip <- !is.null(scan$chip) && nrow(scan$chip$peaks) > 0
    if (!is.null(plots$promoter_map))
      ggsave(file.path(out_dir, sprintf("%s_%s_promoter_map%s.pdf", gene, tf,
                                        if (with_chip) "_chip" else "")),
             plots$promoter_map, width = 11, height = if (with_chip) 7.5 else 6)
    ## ranked likely binding sites: exact coordinates, sequences, probability
    if (!is.null(top_sites) && nrow(top_sites)) {
      write.csv(top_sites, file.path(out_dir, sprintf("%s_%s_top_sites.csv", gene, tf)),
                row.names = FALSE)
      write.csv(scan$sites, file.path(out_dir, sprintf("%s_%s_all_sites_ranked.csv", gene, tf)),
                row.names = FALSE)
      message(sprintf("Wrote: %s_%s_top_sites.csv (%d site(s) with binding probability >= %.2f or top %d; %d high, %d medium)",
                      gene, tf, nrow(top_sites), top_sites_min_prob, top_sites_min_n,
                      sum(top_sites$confidence == "high", na.rm = TRUE),
                      sum(top_sites$confidence == "medium", na.rm = TRUE)))
    }
    if (!is.null(plots$top_sites))
      ggsave(file.path(out_dir, sprintf("%s_%s_top_sites%s.pdf", gene, tf,
                                        if (with_chip) "_chip" else "")),
             plots$top_sites, width = 12, height = if (with_chip) 8 else 6.5)
    if (!is.null(plots$top_sites_nolabels))
      ggsave(file.path(out_dir, sprintf("%s_%s_top_sites%s_nolabels.pdf", gene, tf,
                                        if (with_chip) "_chip" else "")),
             plots$top_sites_nolabels, width = 12, height = if (with_chip) 8 else 6.5)
    if (!is.null(plots$motif_logo)) {
      ggsave(file.path(out_dir, sprintf("%s_%s_motif_logo.pdf", gene, tf)), plots$motif_logo,
             width = max(5, 0.45 * max(vapply(scan$pwms, `[[`, numeric(1), "len")) + 2),
             height = 1.6 * 2 * length(scan$pwms) + 1)
      write.csv(motif_matrices(scan), file.path(out_dir, sprintf("%s_%s_motif_matrices.csv", gene, tf)),
                row.names = FALSE)
    }
    message("Wrote plots to: ", out_dir)
  }
  message("Done.")
  invisible(list(result = scan$result, sites = scan$sites, top_sites = top_sites,
                 plots = plots, meta = scan$meta,
                 promoters = scan$promoters, promoter_summary = psum, chip = scan$chip,
                 usage = scan$usage, scan = scan))
}
