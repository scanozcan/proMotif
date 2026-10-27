#' Promoter track of a gene
#'
#' Draws one gene on the genome: FANTOM5 CAGE signal in two groups of samples
#' (normal vs cancer by default, or any two sets given in `groups`), the gene's
#' promoters (arrows labelled with expression and usage share in the two
#' groups) and its Ensembl transcripts, with exons coloured by the promoter each
#' transcript starts from.
#'
#' @param gene HGNC gene symbol.
#' @param species Ensembl species name (default `"homo_sapiens"`).
#' @param tss_source `"epd"` (EPDnew mapped TSSs, default) or `"ensembl"`.
#' @param transcripts `"all"` (default) or `"protein_coding"`.
#' @param normal FANTOM5 category used as normal (default `"primary cells"`).
#' @param lineage Optional lineage to restrict normal and cancer samples to.
#' @param cage Draw the CAGE panel (default `TRUE`; needs the FANTOM5 folder).
#' @param max_transcripts Maximum transcripts drawn (default 30).
#' @param tx_match Distance (bp) for assigning a transcript to a promoter (default 500).
#' @param flank Extra bp shown on each side (default 2000).
#' @param max_overlap Window overlap above which TSSs form one promoter (default 0.2).
#' @param base_size Base font size (default 11).
#' @param dir FANTOM5 folder (default: `options(promotif.fantom5_dir)`).
#' @param groups Optional named list of two sets of FANTOM5 sample types that
#'   replaces the normal / cancer comparison, e.g.
#'   `list("cell types using P2" = c(...), "other primary cells" = c(...))`.
#'   Sample types may be given as `"sample type"` or `"sample type [category]"`.
#'   Default `NULL` (normal vs cancer).
#' @return A patchwork plot; `attr(p, "height")` gives a suggested page height
#'   and `attr(p, "data")` the drawn promoters, transcripts, exons and peaks.
#' @export
#' @examples
#' \dontrun{
#' p <- plot_promoter_track("TP73")
#' ggplot2::ggsave("TP73_track.pdf", p, width = 11, height = attr(p, "height"))
#'
#' # any two groups of FANTOM5 sample types, e.g. the primary cells that use
#' # GATA2 P1 vs all other primary cells
#' s <- add_promoter_usage(scan_binding_sites("GATA2", "GATA2", tss_source = "epd"))
#' prim <- s$usage$table[s$usage$table$category == "primary cells", ]
#' p1 <- prim$sample_type[prim$dominant == "P1"]
#' other <- setdiff(prim$sample_type, p1)
#' p2 <- plot_promoter_track("GATA2",
#'                           groups = list("cells using P1" = p1, "other primary cells" = other))
#' }
plot_promoter_track <- function(gene, species = "homo_sapiens", tss_source = c("epd", "ensembl"),
                                transcripts = c("all", "protein_coding"),
                                normal = "primary cells", lineage = NULL, cage = TRUE,
                                max_transcripts = 30, tx_match = 500, flank = 2000,
                                max_overlap = 0.2, base_size = 11, dir = NULL, groups = NULL) {
  tss_source <- match.arg(tss_source); transcripts <- match.arg(transcripts)
  if (!is.null(groups) && (length(groups) != 2 || is.null(names(groups))))
    stop("groups must be a named list of two sets of sample types")
  gnames <- if (is.null(groups)) c("normal", "cancer") else names(groups)
  js <- .ensembl_lookup_expand(gene, species)
  strand <- js$strand; chr <- js$seq_region_name
  gt <- .gene_transcripts(gene, species, transcripts, tss_source)
  merge_dist <- round((1 - max_overlap) * 3001)
  txp  <- .cluster_tss(gt$txdf, strand, merge_dist)
  prom <- .promoter_table(txp, strand, 2500, 500, max_overlap)$promoters

  ## Ensembl transcripts with exons
  ens <- Filter(function(t) transcripts == "all" || identical(t$biotype, "protein_coding"),
                js$Transcript %||% list())
  tx <- do.call(rbind, lapply(ens, function(t) data.frame(
    transcript_id = t$id, name = t$display_name %||% t$id, biotype = t$biotype %||% "",
    start = t$start, end = t$end, canonical = isTRUE((t$is_canonical %||% 0) == 1),
    stringsAsFactors = FALSE)))
  tx$tss <- if (strand == 1) tx$start else tx$end
  tx$promoter_id <- vapply(tx$tss, function(s) {
    hit <- which(s >= prom$tss_min - tx_match & s <= prom$tss_max + tx_match)
    if (length(hit)) prom$promoter_id[hit[1]] else "none"
  }, character(1))
  tx <- tx[order(match(tx$promoter_id, c(prom$promoter_id, "none")), !tx$canonical,
                 tx$biotype != "protein_coding"), , drop = FALSE]
  if (nrow(tx) > max_transcripts) {
    message("  showing ", max_transcripts, " of ", nrow(tx), " transcripts (canonical and coding first per promoter)")
    keep <- unlist(lapply(split(seq_len(nrow(tx)), tx$promoter_id), function(i) i[seq_len(min(length(i), 3))]))
    rest <- setdiff(seq_len(nrow(tx)), keep)
    tx <- tx[sort(c(keep, rest[seq_len(max(0, max_transcripts - length(keep)))])), , drop = FALSE]
  }
  tx$row <- rev(seq_len(nrow(tx)))
  ex <- do.call(rbind, lapply(ens[match(tx$transcript_id, vapply(ens, `[[`, "", "id"))], function(t)
    do.call(rbind, lapply(t$Exon, function(e) data.frame(transcript_id = t$id, start = e$start,
                                                           end = e$end)))))
  ex$row <- tx$row[match(ex$transcript_id, tx$transcript_id)]
  ex$promoter_id <- tx$promoter_id[match(ex$transcript_id, tx$transcript_id)]

  lo <- min(tx$start, prom$tss_min) - flank; hi <- max(tx$end, prom$tss_max) + flank
  kb <- function(x) x / 1000
  pcol <- setNames(c(grDevices::hcl.colors(max(3, nrow(prom)), "Dark 3")[seq_len(nrow(prom))], "grey60"),
                   c(prom$promoter_id, "none"))

  ## CAGE peaks and promoter shares by group
  pk <- NULL; share_lab <- setNames(rep("", nrow(prom)), prom$promoter_id)
  if (cage && identical(species, "homo_sapiens")) {
    ch <- .fantom5_chunk(paste0("chr", sub("^chr", "", chr)), if (strand == 1) "+" else "-", dir)
    if (!is.null(ch)) {
      if (is.null(groups)) {
        cls <- .classify_samples(ch$samples$sample_type, ch$samples$category)
        grp <- ifelse(cls$condition == "cancer" & ch$samples$category == "cell lines", "cancer",
               ifelse(cls$condition != "cancer" & ch$samples$category %in% normal, "normal", NA))
        if (!is.null(lineage)) grp[is.na(cls$lineage) | !(cls$lineage %in% lineage)] <- NA
      } else {
        lab2 <- sprintf("%s [%s]", ch$samples$sample_type, ch$samples$category)
        grp <- rep(NA_character_, nrow(ch$samples))
        for (k in 2:1) grp[ch$samples$sample_type %in% groups[[k]] | lab2 %in% groups[[k]]] <- gnames[k]
      }
      i <- which(ch$peaks$summit >= lo & ch$peaks$summit <= hi)
      if (length(i) && any(grp == gnames[1], na.rm = TRUE) && any(grp == gnames[2], na.rm = TRUE)) {
        M <- ch$tpm100[i, , drop = FALSE] / 100
        mean_g <- function(g) rowMeans(M[, which(grp == g), drop = FALSE])
        pk <- data.frame(summit = ch$peaks$summit[i], A = mean_g(gnames[1]), B = mean_g(gnames[2]))
        pk$promoter_id <- vapply(pk$summit, function(s) {
          hit <- which(s >= prom$tss_min - tx_match & s <= prom$tss_max + tx_match)
          if (length(hit)) prom$promoter_id[hit[1]] else "none"
        }, character(1))
        tot <- c(A = sum(pk$A[pk$promoter_id != "none"]), B = sum(pk$B[pk$promoter_id != "none"]))
        for (p in prom$promoter_id)
          share_lab[p] <- sprintf("  %.1f -> %.1f TPM, share %.0f%% -> %.0f%%",
                                  sum(pk$A[pk$promoter_id == p]), sum(pk$B[pk$promoter_id == p]),
                                  100 * sum(pk$A[pk$promoter_id == p]) / max(tot["A"], 1e-9),
                                  100 * sum(pk$B[pk$promoter_id == p]) / max(tot["B"], 1e-9))
        attr(pk, "n_types") <- c(A = sum(grp == gnames[1], na.rm = TRUE), B = sum(grp == gnames[2], na.rm = TRUE))
      }
    }
  }

  th <- theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = "grey92", color = NA),
          strip.text.y = element_text(angle = 0), legend.position = "none")
  span <- hi - lo
  ## room for transcript IDs on the 5' side
  xr <- kb(c(lo - 0.22 * span * (strand == 1), hi + 0.22 * span * (strand == -1)))
  xs <- scale_x_continuous(limits = xr, expand = c(0, 0))
  plots <- list()
  if (!is.null(pk)) {
    glab <- if (is.null(groups)) c("normal", "cancer lines") else gnames
    pl <- rbind(data.frame(pk[, c("summit", "promoter_id")], tpm = pk$A,
                           group = sprintf("%s (%d types)", glab[1], attr(pk, "n_types")["A"])),
                data.frame(pk[, c("summit", "promoter_id")], tpm = pk$B,
                           group = sprintf("%s (%d types)", glab[2], attr(pk, "n_types")["B"])))
    pl$group <- factor(pl$group, unique(pl$group))
    plots$cage <- ggplot(pl) +
      geom_segment(aes(x = kb(summit), xend = kb(summit), y = 0, yend = tpm, color = promoter_id),
                   linewidth = 0.9) +
      facet_grid(group ~ .) + scale_color_manual(values = pcol) + xs +
      scale_y_sqrt() + labs(x = NULL, y = "CAGE TPM\n(mean, sqrt)",
        title = sprintf("%s  |  chr%s:%s-%s (%s strand)  |  promoters from %s", gene, chr,
                        .fmt_bp(lo), .fmt_bp(hi), if (strand == 1) "+" else "-",
                        if (gt$source == "epd") "EPDnew (mapped TSSs)" else "Ensembl TSSs"),
        subtitle = if (is.null(groups))
          sprintf("FANTOM5 CAGE peaks; normal = %s%s; cancer = cancer cell lines | promoter labels: summed mean TPM and share, normal -> cancer",
                  paste(normal, collapse = " + "),
                  if (!is.null(lineage)) paste0("; lineage: ", paste(lineage, collapse = ", ")) else "")
          else "FANTOM5 CAGE peaks, mean over sample types | promoter labels: summed mean TPM and share, top -> bottom panel") +
      th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  arr <- 0.025 * span * strand
  np <- nrow(prom)
  pa <- data.frame(promoter_id = prom$promoter_id, x = prom$tss, xend = prom$tss + arr,
                   ly = 1.5 + 0.55 * ((seq_len(np) - 1) %% 3),          # stagger labels
                   hj = ifelse(prom$tss < mean(xr) * 1000, 0, 1),        # keep labels inside
                   lab = sprintf("%s%s (%d TSS)%s", prom$promoter_id, ifelse(prom$has_canonical, "*", ""),
                                 prom$n_tss, share_lab[prom$promoter_id]))
  plots$prom <- ggplot(pa) +
    geom_rect(data = prom, aes(xmin = kb(tss_min - tx_match), xmax = kb(tss_max + tx_match),
                               ymin = -0.4, ymax = 1.2, fill = promoter_id), alpha = 0.15, inherit.aes = FALSE) +
    geom_segment(aes(x = kb(x), xend = kb(x), y = 0, yend = 1, color = promoter_id), linewidth = 0.9) +
    geom_segment(aes(x = kb(x), xend = kb(xend), y = 1, yend = 1, color = promoter_id), linewidth = 0.9,
                 arrow = grid::arrow(length = grid::unit(0.12, "inches"), type = "closed")) +
    geom_segment(aes(x = kb(x), xend = kb(x), y = 1, yend = ly - 0.15, color = promoter_id),
                 linewidth = 0.3, linetype = 3) +
    geom_text(aes(x = kb(x), y = ly, label = lab, color = promoter_id, hjust = hj),
              size = base_size / 3.4, fontface = "bold") +
    scale_color_manual(values = pcol) + scale_fill_manual(values = pcol) + xs +
    scale_y_continuous(limits = c(-0.5, 1.5 + 0.55 * min(np, 3))) + labs(x = NULL, y = "promoters") + th +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
  plots$tx <- ggplot() +
    geom_segment(data = tx, aes(x = kb(start), xend = kb(end), y = row, yend = row, color = promoter_id),
                 linewidth = 0.4) +
    geom_rect(data = ex, aes(xmin = kb(start), xmax = kb(end), ymin = row - 0.3, ymax = row + 0.3,
                             fill = promoter_id), color = NA) +
    geom_text(data = tx, aes(x = kb(tss) - kb(strand * 0.004 * span), y = row,
                             label = paste0(transcript_id, ifelse(canonical, "*", ""), .bt_abbrev(biotype))),
              hjust = if (strand == 1) 1 else 0, size = base_size / 4.2, color = "grey20") +
    scale_color_manual(values = pcol) + scale_fill_manual(values = pcol) + xs +
    labs(x = sprintf("chr%s position (kb)", chr), y = NULL,
         caption = "* canonical (Ensembl) | NMD nonsense-mediated decay, RI retained intron, nc non-coding | colours: promoter a transcript starts from (grey: none)") +
    th + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  h <- c(cage = 2.2, prom = 0.6 + 0.3 * min(np, 3), tx = max(1.5, 0.22 * nrow(tx)))[names(plots)]
  p <- patchwork::wrap_plots(plots, ncol = 1, heights = h)
  attr(p, "data") <- list(promoters = prom, transcripts = tx, exons = ex, peaks = pk)
  attr(p, "height") <- 1.2 + sum(h) * 1.1
  p
}

.bt_abbrev <- function(b) {
  a <- c(protein_coding = "", nonsense_mediated_decay = " (NMD)", retained_intron = " (RI)",
         processed_transcript = " (nc)", lncRNA = " (nc)", protein_coding_CDS_not_defined = " (nc)")
  out <- a[b]; out[is.na(out)] <- paste0(" (", b[is.na(out)], ")"); unname(out)
}
