## ===========================================================================
## 6. ChIP-Atlas EVIDENCE (optional; internet on first use, then cached)
## ---------------------------------------------------------------------------
## Confirm mode : add_chip_peaks(scan)   - which predicted sites lie inside real
##                ChIP-seq peaks for the same TF (all cell types or a subset).
## Discovery    : add_chip_targets(disc) - for each top candidate TF, the share
##                of its ChIP-seq experiments with a peak within +/- 1/5/10 kb
##                of this gene's TSS (ChIP-Atlas "Target Genes").
## Caveats: a peak (~200-500 bp) shows the TF binds that REGION in some cell
## type - not necessarily at the predicted motif, and not in your cells unless
## you restrict cell types (cell_pattern). HOCOMOCO matrices were partly trained
## on public ChIP-seq, so HOCOMOCO-vs-ChIP agreement is partly circular.
## Human only for peak coordinates (ChIP-Atlas mouse = mm10, Ensembl = GRCm39).
## File formats: https://github.com/inutano/chip-atlas/wiki
## ===========================================================================

.chipatlas_base <- "https://chip-atlas.dbcls.jp/data"

.chipatlas_genome <- function(species, coordinates = TRUE) {
  if (species == "homo_sapiens") return("hg38")
  if (species == "mus_musculus") {
    if (coordinates) {
      message("  ChIP-Atlas mouse peaks are on mm10 but Ensembl mouse is GRCm39 - coordinates do not match; ChIP peak layer skipped.")
      return(NULL)
    }
    return("mm10")
  }
  message("  ChIP-Atlas layer not available for ", species, " in this script.")
  NULL
}

.chip_cache_dir <- function() {
  d <- file.path(tools::R_user_dir("proMotif", "cache"), "chipatlas")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

## Download a ChIP-Atlas file once into the cache. Returns the local path,
## NA if the file does not exist (HTTP 404, e.g. no data for that TF), or NULL
## if ChIP-Atlas could not be reached or the file is larger than max_mb.
.chipatlas_get <- function(url, max_mb = 300) {
  dest <- file.path(.chip_cache_dir(),
                    gsub("/", "__", sub(paste0("^", .chipatlas_base, "/"), "", url)))
  if (file.exists(dest) && file.size(dest) > 0) return(dest)
  h <- tryCatch(httr::HEAD(url, httr::timeout(30)), error = function(e) NULL)
  if (is.null(h)) { message("  ChIP-Atlas not reachable - ChIP layer skipped."); return(NULL) }
  sc <- httr::status_code(h)
  if (sc == 404) return(NA_character_)
  if (sc != 200) {
    message("  ChIP-Atlas returned HTTP ", sc, " for ", basename(url), " - skipped.")
    return(NULL)
  }
  mb <- suppressWarnings(as.numeric(httr::headers(h)[["content-length"]])) / 1e6
  known <- length(mb) == 1 && is.finite(mb)
  if (known && mb > max_mb) {
    message(sprintf("  %s is %.0f MB (> chip_max_mb = %s) - skipped. Raise chip_max_mb to fetch it once (it is then cached).",
                    basename(url), mb, format(max_mb)))
    return(NULL)
  }
  message(sprintf("  downloading %s%s (cached for next time) ...", basename(url),
                  if (known) sprintf(" (%.1f MB)", mb) else ""))
  tmp <- paste0(dest, ".part")
  ok <- .with_retry(function() {
    r <- httr::GET(url, httr::write_disk(tmp, overwrite = TRUE), httr::timeout(3600))
    if (httr::status_code(r) != 200) stop("HTTP ", httr::status_code(r))
    TRUE
  }, tries = 3, waits = c(5, 20), what = paste("ChIP-Atlas", basename(url)))
  if (is.null(ok)) { unlink(tmp); return(NULL) }
  file.rename(tmp, dest)
  dest
}

## Cell-type matching that ignores spaces, "-" and "_" (ChIP-Atlas writes e.g.
## "Hep_G2", "K-562"), so cell_pattern = "HepG2|K562" works as expected.
.cell_match <- function(x, pattern)
  grepl(gsub("[-_ ]", "", pattern), gsub("[-_ ]", "", x), ignore.case = TRUE)

.empty_peaks <- function()
  data.frame(chr = character(0), start = numeric(0), end = numeric(0),
             name = character(0), score = numeric(0), srx = character(0),
             cell_type = character(0), cell_class = character(0),
             stringsAsFactors = FALSE)

## Stream a (large) BED file and keep only peaks overlapping chrom:from-to.
.read_bed_region <- function(path, chrom, from, to, chunk = 500000L) {
  con <- file(path, "r"); on.exit(close(con))
  key <- paste0(chrom, "\t"); keep <- character(0)
  repeat {
    ln <- readLines(con, n = chunk, warn = FALSE)
    if (!length(ln)) break
    ln <- ln[startsWith(ln, key)]
    if (!length(ln)) next
    f <- strsplit(ln, "\t", fixed = TRUE)
    s <- as.numeric(vapply(f, `[`, "", 2)); e <- as.numeric(vapply(f, `[`, "", 3))
    keep <- c(keep, ln[e >= from & s <= to])
  }
  if (!length(keep)) return(.empty_peaks())
  f  <- strsplit(keep, "\t", fixed = TRUE)
  nm <- vapply(f, function(x) if (length(x) >= 4 && !is.na(x[4]))
                 utils::URLdecode(x[4]) else "", "")
  df <- data.frame(chr = chrom,
                   start = as.numeric(vapply(f, `[`, "", 2)),
                   end   = as.numeric(vapply(f, `[`, "", 3)),
                   name  = nm,
                   score = suppressWarnings(as.numeric(vapply(f, `[`, "", 5))),
                   stringsAsFactors = FALSE)
  ## name = GFF3-style attributes, e.g. "ID=SRX...;Name=MTF1 (@ HepG2);..."
  df$srx <- ifelse(grepl("ID=", nm), sub("^.*ID=([^;]+).*$", "\\1", nm),
                   paste0("peak", seq_len(nrow(df))))
  df$cell_type  <- ifelse(grepl("@ ", nm), sub("^.*@ ([^);]+)\\).*$", "\\1", nm),
                          NA_character_)
  df$cell_class <- ifelse(grepl("Cell group=", nm),
                          sub("^.*Cell group=([^;]+).*$", "\\1", nm), NA_character_)
  df
}

## ChIP-Atlas peaks for one TF in a genomic interval (all experiments merged).
## threshold: MACS2 q-value class "05", "10" or "20" (q < 1e-5, 1e-10, 1e-20).
## cell_pattern: optional regex matched against each peak's sample annotation
## (cell line / tissue / cell class), e.g. "HepG2|liver".
## Returns list(peaks, available) or NULL if ChIP-Atlas was not reachable.
#' ChIP-Atlas peaks of a TF in a genomic interval
#'
#' Downloads (once; cached under \code{tools::R_user_dir("proMotif", "cache")})
#' the ChIP-Atlas file merging all ChIP-seq experiments for \code{tf} and
#' returns the peaks overlapping \code{chr:from-to}. Human (hg38) only.
#'
#' @param tf TF symbol as used by ChIP-Atlas (HGNC for human).
#' @param chr,from,to Genomic interval (Ensembl chromosome name, 1-based).
#' @param species Ensembl species token; only \code{"homo_sapiens"} is supported.
#' @param threshold MACS2 q-value class: \code{"05"}, \code{"10"} or \code{"20"}
#'   (q < 1e-5, 1e-10, 1e-20).
#' @param cell_pattern Optional regular expression on the sample annotation
#'   (e.g. \code{"HepG2|liver"}); spaces, "-" and "_" are ignored.
#' @param max_mb Skip files larger than this many MB (default 300).
#' @return \code{list(peaks, available)}, or \code{NULL} if ChIP-Atlas was not
#'   reachable or the file was too large.
#' @export
fetch_chip_peaks <- function(tf, chr, from, to, species = "homo_sapiens",
                             threshold = "05", cell_pattern = NULL, max_mb = 300) {
  genome <- .chipatlas_genome(species, coordinates = TRUE)
  if (is.null(genome)) return(NULL)
  url  <- sprintf("%s/%s/assembled/Oth.ALL.%s.%s.AllCell.bed",
                  .chipatlas_base, genome, threshold, tf)
  path <- .chipatlas_get(url, max_mb)
  if (is.null(path)) return(NULL)
  if (is.na(path)) {
    message("  ChIP-Atlas has no ", tf, " ChIP-seq peak file (", basename(url), ").")
    return(list(peaks = .empty_peaks(), available = FALSE))
  }
  chrom <- if (grepl("^chr", chr)) chr else paste0("chr", chr)
  pk <- .read_bed_region(path, chrom, from, to)
  if (!is.null(cell_pattern) && nrow(pk))
    pk <- pk[.cell_match(pk$name, cell_pattern), , drop = FALSE]
  list(peaks = pk, available = TRUE)
}

## Annotate a confirm-mode scan with ChIP-Atlas peaks of the same TF.
## Adds to scan$result: chip_n_exp (experiments with a peak covering the site),
## chip_max_score (-10*log10 MACS2 q), chip_cell_types, chip_supported.
#' Annotate a confirm-mode scan with ChIP-Atlas peaks
#'
#' Adds \code{chip_n_exp}, \code{chip_max_score}, \code{chip_cell_types} and
#' \code{chip_supported} to \code{scan$result} and stores the peaks in
#' \code{scan$chip}. A peak shows the TF binds that region in some cell type,
#' not necessarily at the predicted motif or in your cells.
#'
#' @inheritParams summarise_promoters
#' @inheritParams fetch_chip_peaks
#' @return The scan, annotated (unchanged if ChIP-Atlas was unavailable).
#' @export
add_chip_peaks <- function(scan, threshold = "05", cell_pattern = NULL, max_mb = 300) {
  m <- scan$meta; prom <- scan$promoters
  from <- min(prom$win_start); to <- max(prom$win_end)
  message("  ChIP-Atlas: ", m$tf, " peaks (MACS2 q < 1e-", threshold, ")",
          if (!is.null(cell_pattern)) paste0(", cells matching '", cell_pattern, "'") else "",
          " ...")
  ck <- fetch_chip_peaks(m$tf, m$chr, from, to, m$species, threshold, cell_pattern, max_mb)
  if (is.null(ck)) return(scan)
  pk <- ck$peaks; res <- scan$result
  if (nrow(res)) {
    hit <- lapply(res$genomic_position, function(g) which(pk$start < g & pk$end >= g))
    res$chip_n_exp      <- vapply(hit, function(h) length(unique(pk$srx[h])), integer(1))
    res$chip_max_score  <- vapply(hit, function(h) if (length(h)) max(pk$score[h])
                                  else NA_real_, numeric(1))
    res$chip_cell_types <- vapply(hit, function(h) {
      ct <- sort(unique(stats::na.omit(pk$cell_type[h])))
      if (length(ct) > 10) paste0(paste(ct[1:10], collapse = ", "), ", +",
                                  length(ct) - 10, " more")
      else paste(ct, collapse = ", ")
    }, character(1))
    res$chip_supported <- res$chip_n_exp > 0
    scan$result <- res
  }
  n_exp <- length(unique(pk$srx))
  if (ck$available)
    message(sprintf("  ChIP-Atlas: %d peak(s) from %d experiment(s) overlap the promoter windows; %d of %d site row(s) lie inside a peak.",
                    nrow(pk), n_exp, sum(res$chip_supported %in% TRUE), nrow(res)))
  scan$chip <- list(peaks = pk, threshold = threshold, cell_pattern = cell_pattern,
                    n_experiments = n_exp, available = ck$available)
  scan
}

## Peak coverage along a region: number of experiments with a peak at each step.
.peak_coverage <- function(peaks, from, to, step = NULL) {
  step <- step %||% max(1, round((to - from) / 600))
  x <- unique(c(seq(from, to, by = step), to))
  n <- vapply(x, function(p) length(unique(peaks$srx[peaks$start < p & peaks$end >= p])),
              integer(1))
  data.frame(x = x, n = n)
}

## ChIP-Atlas "Target Genes" support for one TF at one gene.
.chip_target_row <- function(tf, gene, genome, distance = 1, cell_pattern = NULL,
                             max_mb = 300) {
  out <- data.frame(chip_tf_queried = tf, chip_available = NA, chip_n_exp = NA_integer_,
                    chip_n_bound = NA_integer_, chip_frac_exp = NA_real_,
                    chip_mean_score = NA_real_, chip_max_score = NA_real_,
                    stringsAsFactors = FALSE)
  url  <- sprintf("%s/%s/target/%s.%d.tsv", .chipatlas_base, genome, tf, as.integer(distance))
  path <- .chipatlas_get(url, max_mb)
  if (is.null(path)) return(out)
  if (is.na(path)) { out$chip_available <- FALSE; return(out) }
  con <- file(path, "r"); on.exit(close(con))
  hdr <- strsplit(readLines(con, n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
  exp_cols <- grep("RX[0-9]+", hdr)             # SRX/ERX/DRX experiment columns
  if (!length(exp_cols)) {
    message("  unexpected Target Genes header for ", tf, " - skipped."); return(out)
  }
  if (!is.null(cell_pattern))
    exp_cols <- exp_cols[.cell_match(hdr[exp_cols], cell_pattern)]
  key <- paste0(gene, "\t"); row <- NULL
  repeat {
    ln <- readLines(con, n = 20000L, warn = FALSE)
    if (!length(ln)) break
    hit <- ln[startsWith(ln, key)]
    if (length(hit)) { row <- strsplit(hit[1], "\t", fixed = TRUE)[[1]]; break }
  }
  vals <- if (is.null(row)) rep(0, length(exp_cols))
          else suppressWarnings(as.numeric(row[exp_cols]))
  vals[is.na(vals)] <- 0
  out$chip_available <- length(exp_cols) > 0
  out$chip_n_exp     <- length(exp_cols)
  if (length(exp_cols)) {
    out$chip_n_bound    <- sum(vals > 0)
    out$chip_frac_exp   <- mean(vals > 0)
    out$chip_mean_score <- mean(vals)
    out$chip_max_score  <- max(vals)
  }
  out
}

## Annotate discovery-mode candidates (disc$top, and the matching rows of
## disc$result) with ChIP-Atlas Target Genes support. Heterodimer names such as
## "MAX::MYC" are split; the better-supported partner is reported.
#' Annotate discovery-mode candidates with ChIP-Atlas Target Genes
#'
#' For each top candidate TF, the share of its ChIP-Atlas experiments with a
#' peak within +/- \code{distance} kb of the gene's TSS. Heterodimer names
#' (\code{"MAX::MYC"}) are split and the better-supported partner is reported.
#'
#' @param disc A list returned by \code{\link{discover_regulators}}.
#' @param distance 1, 5 or 10 (kb from the TSS).
#' @inheritParams fetch_chip_peaks
#' @return \code{disc} with \code{chip_*} columns added to \code{top} and
#'   \code{result}.
#' @export
add_chip_targets <- function(disc, distance = 1, cell_pattern = NULL, max_mb = 300) {
  m <- disc$meta
  genome <- .chipatlas_genome(m$species, coordinates = FALSE)
  if (is.null(genome)) return(disc)
  if (!distance %in% c(1, 5, 10)) stop("chip_distance must be 1, 5 or 10 (kb).")
  top <- disc$top
  message(sprintf("  ChIP-Atlas Target Genes (+/- %d kb of TSS)%s for the top %d TF(s) ...",
                  as.integer(distance),
                  if (!is.null(cell_pattern)) paste0(", cells matching '", cell_pattern, "'") else "",
                  nrow(top)))
  fix_case <- if (m$species == "homo_sapiens") toupper else
    function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))
  ann <- do.call(rbind, lapply(unique(top$tf), function(tf) {
    parts <- fix_case(unique(trimws(strsplit(tf, "::", fixed = TRUE)[[1]])))
    rows  <- do.call(rbind, lapply(parts, .chip_target_row, gene = m$gene, genome = genome,
                                   distance = distance, cell_pattern = cell_pattern,
                                   max_mb = max_mb))
    ok   <- which(rows$chip_available %in% TRUE)
    best <- if (length(ok)) ok[which.max(rows$chip_frac_exp[ok])] else 1L
    cbind(tf = tf, rows[best, , drop = FALSE], stringsAsFactors = FALSE)
  }))
  rownames(ann) <- NULL
  cols <- setdiff(names(ann), "tf")
  for (cc in cols) {
    disc$top[[cc]]    <- ann[[cc]][match(disc$top$tf, ann$tf)]
    disc$result[[cc]] <- ann[[cc]][match(disc$result$tf, ann$tf)]
  }
  disc$meta$chip <- list(distance = distance, cell_pattern = cell_pattern)
  message(sprintf("  ChIP-Atlas: %d of %d top TF(s) have ChIP-seq data; %d have a peak near %s in >= 1 experiment.",
                  sum(disc$top$chip_available %in% TRUE), nrow(disc$top),
                  sum(disc$top$chip_n_bound > 0, na.rm = TRUE), m$gene))
  disc
}

## NEW discovery plot: motif evidence (x) vs ChIP-Atlas binding evidence (y).
## Top-right = supported by both. Crosses at the bottom = TF has no ChIP data.
#' Motif vs ChIP-Atlas evidence plot for discovery mode
#'
#' @param disc A list returned by \code{\link{add_chip_targets}}.
#' @param base_size,point_size Plot font and marker size.
#' @return A ggplot, or \code{NULL} if there is no ChIP annotation.
#' @export
plot_regulators_chip <- function(disc, base_size = 12, point_size = 3) {
  top <- disc$top; m <- disc$meta
  if (is.null(top) || !nrow(top) || !"chip_frac_exp" %in% names(top)) return(NULL)
  top$neglogp    <- -log10(pmax(top$pvalue, 1e-300))
  top$family_lab <- ifelse(is.na(top$family) | top$family == "", "unassigned", top$family)
  has <- top$chip_available %in% TRUE
  av <- top[has, , drop = FALSE]; na <- top[!has, , drop = FALSE]
  na$y <- rep(-0.08, nrow(na))
  base_theme <- theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = rel(0.95)),
          plot.subtitle = element_text(size = rel(0.78)))
  cp  <- m$chip$cell_pattern
  sub <- sprintf("x = motif evidence (best-site FIMO p, canonical promoter)\ny = share of the TF's ChIP-Atlas experiments with a peak within +/- %d kb of the TSS%s | crosses = no ChIP-Atlas data",
                 as.integer(m$chip$distance),
                 if (!is.null(cp)) paste0(" (cells ~ '", cp, "')") else "")
  p <- ggplot() + geom_hline(yintercept = 0, color = "grey80")
  if (nrow(av)) p <- p + geom_point(data = av, aes(neglogp, chip_frac_exp,
                                                   color = family_lab, size = chip_n_exp),
                                    alpha = 0.85)
  if (nrow(na)) p <- p + geom_point(data = na, aes(neglogp, y), shape = 4,
                                    size = point_size * 0.8, color = "grey45")
  lab <- rbind(if (nrow(av)) data.frame(x = av$neglogp, y = av$chip_frac_exp, l = av$tf),
               if (nrow(na)) data.frame(x = na$neglogp, y = na$y, l = na$tf))
  p <- p + (if (requireNamespace("ggrepel", quietly = TRUE))
    ggrepel::geom_text_repel(data = lab, aes(x, y, label = l), size = base_size / 4,
                             min.segment.length = 0, max.overlaps = Inf, seed = 1)
  else geom_text(data = lab, aes(x, y, label = l), vjust = -0.7, size = base_size / 4,
                 check_overlap = TRUE))
  p + scale_y_continuous(breaks = c(0, .25, .5, .75, 1),
                         labels = function(x) paste0(round(100 * x), "%")) +
    coord_cartesian(ylim = c(-0.12, 1.02)) +
    scale_size_continuous(range = c(1.5, point_size * 2.5), name = "ChIP experiments") +
    labs(x = "-log10 FIMO per-site p (best site)",
         y = "Experiments with a peak at this promoter", color = "TF family",
         title = sprintf("Motif vs ChIP-Atlas evidence for %s's candidate regulators", m$gene),
         subtitle = sub) + base_theme
}

