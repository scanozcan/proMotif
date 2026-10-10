#' Build binding-site plots from a scan result
#'
#' Fast, local plotting step: turns a \code{\link{scan_binding_sites}} result
#' into ggplot objects. Separating this from the scan lets a UI re-style
#' (text/point size, colours) instantly without re-running the scan.
#'
#' @param scan A list returned by \code{\link{scan_binding_sites}}.
#' @param base_size Base font size for the plots (default 12).
#' @param point_size Size of the site markers (default 2.6).
#' @param color_seed Optional integer; fixes the two database colours
#'   reproducibly. \code{NULL} (default) draws fresh colours each call.
#' @return A list of ggplot objects: \code{canonical}, \code{by_transcript},
#'   \code{genomic} (score-based), plus \code{canonical_sig},
#'   \code{by_transcript_sig}, \code{genomic_sig} (-log10 adjusted p; \code{NULL}
#'   when no p-values are available). Elements are \code{NULL} if no sites.
#' @examples
#' \dontrun{
#' s <- scan_binding_sites("ATP7B", "MTF1")
#' p <- plot_binding_sites(s, base_size = 16)
#' p$genomic
#' }
#' @export
plot_binding_sites <- function(scan, base_size = 12, point_size = 2.6,
                               color_seed = NULL) {
  result <- scan$result; reps <- scan$reps; m <- scan$meta
  if (is.null(result) || nrow(result) == 0)
    return(list(canonical = NULL, by_transcript = NULL, genomic = NULL,
                canonical_sig = NULL, by_transcript_sig = NULL, genomic_sig = NULL))

  gene <- m$gene; tf <- m$tf; gene_chr <- m$chr; gene_strand <- m$strand
  upstream <- m$upstream; downstream <- m$downstream; threshold_frac <- m$threshold_frac
  ## wrap the long transcript label onto 2-3 lines for the facet strips
  .wrap_lab <- function(x) gsub(" \\| ", "\n", sub(" \\(canonical\\)", "\n(canonical)", x))
  wlev <- .wrap_lab(reps$label)
  result$transcript_label <- factor(.wrap_lab(result$transcript_label), levels = wlev)

  shape_vals <- c("+" = 16, "-" = 17)
  base_theme <- theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(size = rel(0.95)),
          plot.subtitle = element_text(size = rel(0.78)))
  if (!is.null(color_seed)) set.seed(as.integer(color_seed))
  db_cols <- setNames(sample(.palette_pool, 2), c("JASPAR", "HOCOMOCO"))

  landscape <- function(d, title, subtitle) {
    ggplot(d, aes(position, score_frac, color = database)) +
      geom_segment(aes(xend = position, yend = threshold_frac), alpha = 0.4) +
      geom_point(aes(shape = strand), size = point_size) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
      scale_shape_manual(values = shape_vals) +
      scale_color_manual(values = db_cols, name = "Database") +
      scale_x_continuous(n.breaks = 6) +
      coord_cartesian(ylim = c(threshold_frac, 1.03)) +
      labs(x = "Position relative to TSS (bp)", y = "PWM score (fraction of max)",
           shape = "Motif strand", title = title, subtitle = subtitle) + base_theme
  }

  can_lab <- (wlev[reps$is_canonical])[1] %||% wlev[1]
  can <- result[result$transcript_label == can_lab, ]
  p1 <- NULL
  if (nrow(can) > 0)
    p1 <- landscape(can,
      sprintf("Predicted %s binding sites in the %s promoter (canonical)", tf, gene),
      sprintf("%s | window -%d/+%d bp | threshold %.2f | %d site(s)",
              can_lab, upstream, downstream, threshold_frac, nrow(can)))

  p2 <- landscape(result,
      sprintf("Predicted %s binding sites across %s transcripts", tf, gene),
      sprintf("window -%d/+%d bp around each TSS | threshold %.2f",
              upstream, downstream, threshold_frac)) +
    facet_wrap(~ transcript_label, ncol = 1, drop = FALSE) +
    theme(strip.text = element_text(size = base_size - 4),
          strip.background = element_rect(fill = "grey92", color = NA))

  tss_lines <- data.frame(tss = reps$tss, is_canonical = reps$is_canonical)

  ## Genomic view: collapse the same physical site re-detected across overlapping
  ## promoter windows to a single marker. A site inside several transcripts' windows
  ## is one genomic locus, so without this it would stack one point per transcript.
  ## Keep per-database and per-strand distinctions.
  gkeys   <- c("genomic_position", "strand", "database", "matrix_id")
  gresult <- result[!duplicated(result[gkeys]), , drop = FALSE]

  grng  <- range(c(result$genomic_position, reps$tss))
  dir   <- if (gene_strand == 1) 1 else -1          # transcription direction on genomic axis
  y_arr <- max(result$score_frac) + 0.015
  can_tss <- reps$tss[reps$is_canonical]
  if (length(can_tss) == 0) can_tss <- reps$tss[1]   # fallback if none flagged
  tss_arrows <- data.frame(x = can_tss, xend = can_tss + dir * 0.05 * diff(grng),
                           y = y_arr, yend = y_arr)
  p3 <- ggplot(gresult, aes(genomic_position, score_frac)) +
    geom_vline(data = tss_lines, aes(xintercept = tss, linetype = is_canonical),
               color = "grey45") +
    geom_segment(data = tss_arrows, aes(x = x, xend = xend, y = y, yend = yend),
                 arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"),
                 color = "grey30", linewidth = 0.6, inherit.aes = FALSE) +
    geom_point(aes(color = database, shape = strand), size = point_size) +
    scale_shape_manual(values = shape_vals) +
    scale_color_manual(values = db_cols, name = "Database") +
    scale_linetype_manual(values = c(`TRUE` = "dashed", `FALSE` = "dotted"),
                          labels = c(`TRUE` = "canonical TSS", `FALSE` = "alt TSS"),
                          name = "TSS") +
    scale_x_continuous(labels = .fmt_bp, n.breaks = 5) +
    labs(x = sprintf("Absolute position on chromosome %s (bp)", gene_chr),
         y = "PWM score (fraction of max)", shape = "Motif strand",
         title = sprintf("Predicted %s binding sites at the %s locus (genomic)", tf, gene),
         subtitle = sprintf("strand %s | lines = TSSs | threshold %.2f | %d site(s), shared promoters shown once",
                            ifelse(gene_strand == 1, "+", "-"), threshold_frac, nrow(gresult))) +
    base_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ## ---- significance plots (the three plots above are unchanged) ------------
  sig <- list(canonical_sig = NULL, by_transcript_sig = NULL, genomic_sig = NULL)
  if ("pvalue" %in% names(result) && any(is.finite(result$pvalue))) {
    result$neglogq <- -log10(pmax(result$pvalue, 1e-300))
    sig_line <- -log10(0.05)
    sig_sub  <- "y = -log10 per-site p (TFMPvalue, GC-aware background); higher = stronger match; dotted line = 0.05"
    sig_landscape <- function(d, title) {
      ggplot(d, aes(position, neglogq, color = database)) +
        geom_hline(yintercept = sig_line, linetype = "dotted", color = "grey50") +
        geom_point(aes(shape = strand), size = point_size) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
        scale_shape_manual(values = shape_vals) +
        scale_color_manual(values = db_cols, name = "Database") +
        scale_x_continuous(n.breaks = 6) +
        labs(x = "Position relative to TSS (bp)", y = "-log10 per-site p",
             shape = "Motif strand", title = title, subtitle = sig_sub) + base_theme
    }
    if (nrow(can) > 0)
      sig$canonical_sig <- sig_landscape(result[result$transcript_label == can_lab, ],
        sprintf("%s binding-site significance in the %s promoter (canonical)", tf, gene))
    sig$by_transcript_sig <- sig_landscape(result,
        sprintf("%s binding-site significance across %s transcripts", tf, gene)) +
      facet_wrap(~ transcript_label, ncol = 1, drop = FALSE) +
      theme(strip.text = element_text(size = base_size - 4),
            strip.background = element_rect(fill = "grey92", color = NA))
    gresult$neglogq <- -log10(pmax(gresult$pvalue, 1e-300))
    sig_arrows <- data.frame(x = can_tss, xend = can_tss + dir * 0.05 * diff(grng),
                             y = max(result$neglogq), yend = max(result$neglogq))
    sig$genomic_sig <- ggplot(gresult, aes(genomic_position, neglogq)) +
      geom_hline(yintercept = sig_line, linetype = "dotted", color = "grey50") +
      geom_vline(data = tss_lines, aes(xintercept = tss, linetype = is_canonical),
                 color = "grey45") +
      geom_segment(data = sig_arrows, aes(x = x, xend = xend, y = y, yend = yend),
                   arrow = grid::arrow(length = grid::unit(0.18, "cm"), type = "closed"),
                   color = "grey30", linewidth = 0.6, inherit.aes = FALSE) +
      geom_point(aes(color = database, shape = strand), size = point_size) +
      scale_shape_manual(values = shape_vals) +
      scale_color_manual(values = db_cols, name = "Database") +
      scale_linetype_manual(values = c(`TRUE` = "dashed", `FALSE` = "dotted"),
                            labels = c(`TRUE` = "canonical TSS", `FALSE` = "alt TSS"),
                            name = "TSS") +
      scale_x_continuous(labels = .fmt_bp, n.breaks = 5) +
      labs(x = sprintf("Absolute position on chromosome %s (bp)", gene_chr),
           y = "-log10 per-site p", shape = "Motif strand",
           title = sprintf("%s binding-site significance at the %s locus (genomic)", tf, gene),
           subtitle = sig_sub) +
      base_theme + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  }

  c(list(canonical = p1, by_transcript = p2, genomic = p3), sig)
}
