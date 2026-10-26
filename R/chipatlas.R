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
             summit = numeric(0), summit_exact = logical(0),
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
  ## merged ("assembled") files carry no summit: the peak midpoint stands in
  df$summit <- (df$start + df$end) / 2; df$summit_exact <- FALSE
  df
}

## ChIP-Atlas peaks for one TF in a genomic interval (all experiments merged).
## threshold: MACS2 q-value class "05", "10" or "20" (q < 1e-5, 1e-10, 1e-20).
## cell_pattern: optional regex matched against each peak's sample annotation
## (cell line / tissue / cell class), e.g. "HepG2|liver".
## Returns list(peaks, available) or NULL if ChIP-Atlas was not reachable.
#' ChIP-Atlas peaks of a TF in a genomic interval
#'
#' Returns the ChIP-Atlas peaks of \code{tf} overlapping \code{chr:from-to}.
#' Files are downloaded once and cached under
#' \code{tools::R_user_dir("proMotif", "cache")}. Human (hg38) only.
#'
#' @details \code{max_mb} is a data budget. With \code{mode = "auto"}, the
#'   ChIP-Atlas file merging all experiments of the TF is used if it is cached
#'   or not larger than \code{max_mb}; otherwise the TF's experiments are
#'   listed (\code{\link{list_chip_experiments}}), filtered by
#'   \code{cell_pattern}, ordered newest first, and per-experiment peak files
#'   are used until \code{max_mb} is reached
#'   (\code{\link{fetch_chip_peaks_budget}}). Per-experiment files also give the
#'   MACS2 summit; merged files give the peak interval only.
#'   \code{mode = "assembled"} uses only the merged file (skipped if larger than
#'   \code{max_mb}); \code{mode = "experiments"} always uses per-experiment files.
#'
#' @param tf TF symbol as used by ChIP-Atlas (HGNC for human).
#' @param chr,from,to Genomic interval (Ensembl chromosome name, 1-based).
#' @param species Ensembl species token; only \code{"homo_sapiens"} is supported.
#' @param threshold MACS2 q-value class: \code{"05"}, \code{"10"} or \code{"20"}
#'   (q < 1e-5, 1e-10, 1e-20).
#' @param cell_pattern Optional regular expression on the sample annotation
#'   (e.g. \code{"HepG2|liver"}); spaces, "-" and "_" are ignored.
#' @param max_mb Data budget in MB (default 300); see Details.
#' @param mode `"auto"` (default), `"assembled"` or `"experiments"`; see Details.
#' @param max_experiments Maximum number of per-experiment files (default `Inf`;
#'   a finite value switches to per-experiment files).
#' @return A list with \code{peaks} (data frame with \code{summit}),
#'   \code{available}, \code{mode}, \code{n_total} and \code{n_used}
#'   (experiments), or \code{NULL} if ChIP-Atlas was not reachable.
#' @export
fetch_chip_peaks <- function(tf, chr, from, to, species = "homo_sapiens",
                             threshold = "05", cell_pattern = NULL, max_mb = 300,
                             mode = c("auto", "assembled", "experiments"),
                             max_experiments = Inf) {
  mode   <- match.arg(mode)
  genome <- .chipatlas_genome(species, coordinates = TRUE)
  if (is.null(genome)) return(NULL)
  url  <- sprintf("%s/%s/assembled/Oth.ALL.%s.%s.AllCell.bed",
                  .chipatlas_base, genome, threshold, tf)
  if (mode == "experiments" || is.finite(max_experiments))
    return(fetch_chip_peaks_budget(tf, chr, from, to, species, threshold, cell_pattern,
                                   budget_mb = max_mb, max_experiments = max_experiments))
  if (mode == "auto") {
    dest <- file.path(.chip_cache_dir(),
                      gsub("/", "__", sub(paste0("^", .chipatlas_base, "/"), "", url)))
    if (!(file.exists(dest) && file.size(dest) > 0)) {
      h  <- tryCatch(httr::HEAD(url, httr::timeout(30)), error = function(e) NULL)
      if (is.null(h)) { message("  ChIP-Atlas not reachable - ChIP layer skipped."); return(NULL) }
      if (httr::status_code(h) == 404) {
        message("  ChIP-Atlas has no ", tf, " ChIP-seq peak file (", basename(url), ").")
        return(list(peaks = .empty_peaks(), available = FALSE, mode = "assembled",
                    n_total = 0L, n_used = 0L))
      }
      mb <- suppressWarnings(as.numeric(httr::headers(h)[["content-length"]])) / 1e6
      if (length(mb) == 1 && is.finite(mb) && mb > max_mb) {
        message(sprintf("  %s is %.0f MB (> chip_max_mb = %s): switching to per-experiment files, newest first.",
                        basename(url), mb, format(max_mb)))
        return(fetch_chip_peaks_budget(tf, chr, from, to, species, threshold, cell_pattern,
                                       budget_mb = max_mb))
      }
    }
  }
  path <- .chipatlas_get(url, max_mb)
  if (is.null(path)) return(NULL)
  if (is.na(path)) {
    message("  ChIP-Atlas has no ", tf, " ChIP-seq peak file (", basename(url), ").")
    return(list(peaks = .empty_peaks(), available = FALSE, mode = "assembled",
                n_total = 0L, n_used = 0L))
  }
  chrom <- if (grepl("^chr", chr)) chr else paste0("chr", chr)
  pk <- .read_bed_region(path, chrom, from, to)
  if (!is.null(cell_pattern) && nrow(pk))
    pk <- pk[.cell_match(pk$name, cell_pattern), , drop = FALSE]
  ## number of experiments queried (for the support fraction): all of the TF's
  ## experiments (matching cell_pattern), from its experiment list
  ex <- tryCatch(list_chip_experiments(tf, species, cell_pattern), error = function(e) NULL)
  n_tot <- max(if (is.null(ex)) 0L else nrow(ex), length(unique(pk$srx)))
  list(peaks = pk, available = TRUE, mode = "assembled", n_total = n_tot, n_used = n_tot,
       mb_used = file.size(path) / 1e6)
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
#' @param mode ChIP-Atlas data source; see [fetch_chip_peaks()] (default `"auto"`).
#' @param max_experiments Maximum number of experiments (default `Inf`).
#' @param sigma Summit-proximity scale in bp (default 50; doubled for merged files, which give only the peak midpoint).
#' @param support `"locus"` (default): normalise ChIP support by the experiments with a peak in the scanned windows (+ `prior`); `"all"`: by all queried experiments.
#' @param prior Pseudo-count added to the number of locus experiments (default 1, so one experiment gives at most 0.5).
#' @param rep_half Number of experiments at which the replication term reaches half its maximum (default 20).
#' @return The scan, annotated (unchanged if ChIP-Atlas was unavailable).
#' @export
add_chip_peaks <- function(scan, threshold = "05", cell_pattern = NULL, max_mb = 300,
                           mode = c("auto", "assembled", "experiments"),
                           max_experiments = Inf, sigma = 50,
                           support = c("locus", "all"), prior = 1, rep_half = 20) {
  mode <- match.arg(mode); support <- match.arg(support)
  m <- scan$meta; prom <- scan$promoters
  from <- min(prom$win_start); to <- max(prom$win_end)
  message("  ChIP-Atlas: ", m$tf, " peaks (MACS2 q < 1e-", threshold, ")",
          if (!is.null(cell_pattern)) paste0(", cells matching '", cell_pattern, "'") else "",
          " ...")
  ## heterodimers ("MAX::MYC"): peaks of both partners are combined
  parts <- unique(trimws(strsplit(m$tf, "::", fixed = TRUE)[[1]]))
  if (identical(m$species, "homo_sapiens")) parts <- toupper(parts)
  if (length(parts) > 1)
    message("  heterodimer: combining ChIP-Atlas peaks of ", paste(parts, collapse = " and "))
  cks <- Filter(Negate(is.null), lapply(parts, function(tf1)
    fetch_chip_peaks(tf1, m$chr, from, to, m$species, threshold, cell_pattern, max_mb,
                     mode = mode, max_experiments = max_experiments)))
  if (!length(cks)) return(scan)
  ck <- list(peaks = do.call(rbind, lapply(cks, `[[`, "peaks")),
             available = any(vapply(cks, `[[`, logical(1), "available")))
  n_queried <- sum(vapply(cks, function(x) as.numeric(x$n_used %||% 0), numeric(1)))
  modes <- unique(vapply(cks, function(x) x$mode %||% "assembled", ""))
  pk <- ck$peaks; res <- scan$result
  if (nrow(res)) {
    ctr <- res$site_center %||% res$genomic_position
    hit <- lapply(ctr, function(g) which(pk$start <= g & pk$end >= g))
    res$chip_n_exp      <- vapply(hit, function(h) length(unique(pk$srx[h])), integer(1))
    ## summit proximity: per experiment, the best Gaussian weight of its peaks
    ## covering the site (sigma bp; doubled when only the peak midpoint is known)
    res$chip_summit_dist <- vapply(seq_along(hit), function(i) {
      h <- hit[[i]]; if (!length(h)) NA_real_ else min(abs(pk$summit[h] - ctr[i])) }, numeric(1))
    ## credit of experiment e to site i: summit proximity w (Gaussian, sigma bp;
    ## doubled when only the peak midpoint is known) x motif strength relative to
    ## the best motif (same database) inside e's peaks here, so weak matches next
    ## to a strong one share less of the peak (CentriMo-like assignment)
    long <- do.call(rbind, lapply(seq_along(hit), function(i) {
      h <- hit[[i]]; if (!length(h)) return(NULL)
      sg <- ifelse(pk$summit_exact[h] %in% TRUE, sigma, 2 * sigma)
      w  <- exp(-(pk$summit[h] - ctr[i])^2 / (2 * sg^2))
      w  <- tapply(w, pk$srx[h], max)
      data.frame(i = i, srx = names(w), w = as.numeric(w), stringsAsFactors = FALSE) }))
    sw <- numeric(nrow(res)); rel_mean <- numeric(nrow(res))
    if (!is.null(long) && nrow(long)) {
      mstr <- if ("pvalue" %in% names(res) && any(is.finite(res$pvalue)))
                -log10(pmax(res$pvalue, 1e-300)) else res$score_frac
      long$m  <- mstr[long$i]; long$db <- res$database[long$i]
      long$mx <- ave(long$m, long$srx, long$db, FUN = function(v) max(v, na.rm = TRUE))
      long$rel    <- ifelse(long$mx > 0, long$m / long$mx, 1)
      long$credit <- long$w * long$rel
      s  <- tapply(long$credit, long$i, sum)
      sw[as.integer(names(s))] <- s
      rr <- tapply(long$rel, long$i, mean)
      rel_mean[as.integer(names(rr))] <- rr
    }
    ## replication: how many experiments have a peak over the site - Hill curve,
    ## negligible for 1-3 experiments, 0.5 at `rep_half` (20), ~1 from ~60 on;
    ## scaled by the site's motif strength relative to the best motif in those
    ## peaks, so weak matches inside a broad peak are not all lifted equally
    k <- res$chip_n_exp
    res$chip_replication <- (k^2 / (k^2 + rep_half^2)) * rel_mean
    ## two normalisations: among all queried experiments, and among the
    ## experiments that bind this locus at all (+ pseudo-count `prior`, so one
    ## experiment gives at most 1/(1 + prior)). The second is the default
    ## (support = "locus"): ChIP-Atlas mixes many cell types, most of which do
    ## not use a given promoter, which would cap the first near zero.
    n_locus <- length(unique(pk$srx))
    res$chip_support_all   <- if (n_queried > 0) pmin(1, sw / n_queried) else NA_real_
    res$chip_support_locus <- pmin(1, sw / (n_locus + prior))
    res$chip_support <- if (support == "all") res$chip_support_all else res$chip_support_locus
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
    message(sprintf("  ChIP-Atlas: %d peak(s) from %d of %d queried experiment(s) overlap the promoter windows; %d of %d site row(s) lie inside a peak.",
                    nrow(pk), n_exp, as.integer(n_queried), sum(res$chip_supported %in% TRUE), nrow(res)))
  scan$chip <- list(peaks = pk, threshold = threshold, cell_pattern = cell_pattern,
                    n_experiments = n_exp, n_queried = n_queried, available = ck$available,
                    mode = paste(modes, collapse = "+"), max_mb = max_mb, sigma = sigma,
                    support = support, prior = prior,
                    n_total = sum(vapply(cks, function(x) as.numeric(x$n_total %||% 0), numeric(1))),
                    mb_used = sum(vapply(cks, function(x) as.numeric(x$mb_used %||% 0), numeric(1))))
  score_binding_sites(scan)
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
#' @param max_mb Skip Target Genes files larger than this many MB (default 300).
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
  use_bg <- identical(m$rank_by, "background") && !is.null(top$bg_pvalue)
  top$neglogp    <- -log10(pmax(if (use_bg) top$bg_pvalue else top$pvalue, 1e-300))
  plab <- if (use_bg) "-log10 promoter-background p (best site)" else
                      "-log10 FIMO per-site p (best site)"
  top <- top[order(-top$neglogp), , drop = FALSE]       # one point per TF
  top <- top[!duplicated(top$tf), , drop = FALSE]        # (best over promoters)
  top$family_lab <- ifelse(is.na(top$family) | top$family == "", "unassigned", top$family)
  has <- top$chip_available %in% TRUE
  av <- top[has, , drop = FALSE]; na <- top[!has, , drop = FALSE]
  na$y <- rep(-0.08, nrow(na))
  base_theme <- theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = rel(0.95)),
          plot.subtitle = element_text(size = rel(0.78)))
  cp  <- m$chip$cell_pattern
  sub <- sprintf("x = motif evidence (best site over the scanned promoters)\ny = share of the TF's ChIP-Atlas experiments with a peak within +/- %d kb of the TSS%s | crosses = no ChIP-Atlas data",
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
    labs(x = plab,
         y = "Experiments with a peak at this promoter", color = "TF family",
         title = sprintf("Motif vs ChIP-Atlas evidence for %s's candidate regulators", m$gene),
         subtitle = sub) + base_theme
}


#' List the ChIP-Atlas experiments of a TF
#'
#' Reads the experiment list of a TF (accession and cell type) from the header
#' line of its ChIP-Atlas "Target Genes" table; only the first ~256 kB of the
#' file are requested. Experiments are ordered newest first (highest accession
#' number, a proxy for submission date). The list is cached.
#'
#' @param tf TF symbol (e.g. `"ESR1"`).
#' @param species `"homo_sapiens"` (hg38) or `"mus_musculus"` (mm10).
#' @param cell_pattern Optional regular expression matched against the cell
#'   type (spaces, `-` and `_` are ignored), e.g. `"MCF-7|T-47D"`.
#' @return A data frame with `srx`, `cell_type` and `accession_number`, or
#'   `NULL` if ChIP-Atlas could not be reached.
#' @export
#' @examples
#' \dontrun{
#' head(list_chip_experiments("ESR1", cell_pattern = "MCF-7"))
#' }
list_chip_experiments <- function(tf, species = "homo_sapiens", cell_pattern = NULL) {
  genome <- .chipatlas_genome(species, coordinates = FALSE)
  if (is.null(genome)) return(NULL)
  dest <- file.path(.chip_cache_dir(), sprintf("experiments_%s_%s.tsv", genome, tf))
  if (file.exists(dest) && file.size(dest) > 0) {
    ex <- utils::read.delim(dest, stringsAsFactors = FALSE, colClasses = "character")
  } else {
    url <- sprintf("%s/%s/target/%s.1.tsv", .chipatlas_base, genome, tf)
    hdr <- NULL
    for (kb in c(256, 2048)) {
      r <- tryCatch(httr::GET(url, httr::add_headers(Range = sprintf("bytes=0-%d", kb * 1024 - 1)),
                              httr::timeout(60)), error = function(e) NULL)
      if (is.null(r)) { message("  ChIP-Atlas not reachable."); return(NULL) }
      if (httr::status_code(r) == 404) {
        message("  ChIP-Atlas has no experiment list for ", tf, ".")
        return(data.frame(srx = character(0), cell_type = character(0)))
      }
      txt <- rawToChar(httr::content(r, as = "raw"))
      if (grepl("\n", txt, fixed = TRUE)) { hdr <- sub("\n.*$", "", txt); break }
    }
    if (is.null(hdr)) { message("  could not read the experiment list of ", tf, "."); return(NULL) }
    f  <- strsplit(hdr, "\t", fixed = TRUE)[[1]]
    f  <- f[grepl("^[SDE]RX[0-9]+", f)]
    ex <- data.frame(srx = sub("\\|.*$", "", f),
                     cell_type = gsub("_", " ", sub("^[^|]*\\|?", "", f)),
                     stringsAsFactors = FALSE)
    ex <- ex[!duplicated(ex$srx), , drop = FALSE]
    utils::write.table(ex, dest, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  ex$accession_number <- as.numeric(sub("^[SDE]RX", "", ex$srx))
  ex <- ex[order(-ex$accession_number), , drop = FALSE]      # newest first
  if (!is.null(cell_pattern) && nrow(ex))
    ex <- ex[.cell_match(ex$cell_type, cell_pattern) | .cell_match(ex$srx, cell_pattern), , drop = FALSE]
  rownames(ex) <- NULL
  ex
}

## One experiment's peaks overlapping chrom:from-to (narrowPeak: col 10 = summit offset).
.read_each_bed_region <- function(path, chrom, from, to, srx, cell_type) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[startsWith(ln, paste0(chrom, "\t"))]
  if (!length(ln)) return(.empty_peaks())
  f <- strsplit(ln, "\t", fixed = TRUE)
  s <- as.numeric(vapply(f, `[`, "", 2)); e <- as.numeric(vapply(f, `[`, "", 3))
  k <- which(e >= from & s <= to)
  if (!length(k)) return(.empty_peaks())
  f <- f[k]
  sm <- suppressWarnings(as.numeric(vapply(f, function(x) if (length(x) >= 10) x[10] else NA_character_, "")))
  data.frame(chr = chrom, start = s[k], end = e[k],
             name = sprintf("ID=%s;Name=(@ %s)", srx, cell_type),
             score = suppressWarnings(as.numeric(vapply(f, `[`, "", 5))),
             srx = srx, cell_type = cell_type, cell_class = NA_character_,
             summit = ifelse(is.na(sm), (s[k] + e[k]) / 2, s[k] + sm),
             summit_exact = !is.na(sm), stringsAsFactors = FALSE)
}

#' ChIP-Atlas peaks of a TF within a data budget
#'
#' Downloads per-experiment ChIP-Atlas peak files of a TF, newest experiments
#' first, until `budget_mb` of data (or `max_experiments` files) is reached,
#' and returns the peaks overlapping a region. Per-experiment files are small
#' (about 0.5-2 MB), carry the MACS2 summit and are cached, so later genes with
#' the same TF need no download. [fetch_chip_peaks()] uses this function
#' automatically when the merged peak file of a TF is larger than its `max_mb`.
#'
#' @inheritParams fetch_chip_peaks
#' @param budget_mb Maximum amount of peak data to use, in MB (default 1000).
#' @param max_experiments Maximum number of experiments (default `Inf`).
#' @param batch Number of files downloaded in parallel per round (default 24;
#'   parallel downloads need the curl package >= 5.0).
#' @return A list with `peaks` (data frame incl. `summit`), `available`,
#'   `mode = "experiments"`, `n_total`, `n_used`, `mb_used` and `experiments`
#'   (one row per experiment, with `used`).
#' @export
#' @examples
#' \dontrun{
#' pk <- fetch_chip_peaks_budget("ESR1", "21", 42366000, 42369000,
#'                               cell_pattern = "MCF-7", budget_mb = 300)
#' }
fetch_chip_peaks_budget <- function(tf, chr, from, to, species = "homo_sapiens",
                                    threshold = "05", cell_pattern = NULL,
                                    budget_mb = 1000, max_experiments = Inf,
                                    batch = 24L) {
  genome <- .chipatlas_genome(species, coordinates = TRUE)
  if (is.null(genome)) return(NULL)
  ex <- list_chip_experiments(tf, species, cell_pattern)
  if (is.null(ex)) return(NULL)
  if (!nrow(ex)) {
    message("  no ", tf, " experiments", if (!is.null(cell_pattern))
            paste0(" matching '", cell_pattern, "'") else "", ".")
    return(list(peaks = .empty_peaks(), available = FALSE, mode = "experiments",
                n_total = 0L, n_used = 0L, mb_used = 0, experiments = ex))
  }
  dir <- file.path(.chip_cache_dir(), "eachData", genome, paste0("bed", threshold))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ex$path <- file.path(dir, sprintf("%s.%s.bed", ex$srx, threshold))
  ex$url  <- sprintf("%s/%s/eachData/bed%s/%s.%s.bed", .chipatlas_base, genome, threshold,
                     ex$srx, threshold)
  ex$mb <- NA_real_; ex$used <- FALSE
  cached <- file.exists(ex$path) & file.size(ex$path) > 0
  ex$mb[cached] <- file.size(ex$path[cached]) / 1e6
  use_multi <- requireNamespace("curl", quietly = TRUE) &&
               "multi_download" %in% getNamespaceExports("curl")
  message(sprintf("  ChIP-Atlas %s: %d experiment(s)%s; using newest first up to %s MB%s ...",
                  tf, nrow(ex), if (!is.null(cell_pattern)) paste0(" matching '", cell_pattern, "'") else "",
                  format(budget_mb), if (is.finite(max_experiments))
                    sprintf(" / %d experiments", as.integer(max_experiments)) else ""))
  used_mb <- 0; n_dl <- 0L; i <- 1L
  while (i <= nrow(ex) && used_mb < budget_mb && sum(ex$used) < max_experiments) {
    idx  <- i:min(nrow(ex), i + batch - 1L)
    need <- idx[!(file.exists(ex$path[idx]) & file.size(ex$path[idx]) > 0)]
    if (length(need)) {
      tmp <- paste0(ex$path[need], ".part")
      if (use_multi) {
        r <- tryCatch(curl::multi_download(ex$url[need], tmp, progress = FALSE, resume = FALSE),
                      error = function(e) NULL)
        okk <- if (is.null(r)) rep(FALSE, length(need)) else r$success %in% TRUE & r$status_code == 200
      } else {
        okk <- vapply(seq_along(need), function(j) {
          rr <- tryCatch(httr::GET(ex$url[need[j]], httr::write_disk(tmp[j], overwrite = TRUE),
                                   httr::timeout(600)), error = function(e) NULL)
          !is.null(rr) && httr::status_code(rr) == 200 }, logical(1))
      }
      for (j in seq_along(need)) {
        if (okk[j]) { file.rename(tmp[j], ex$path[need[j]]); n_dl <- n_dl + 1L } else unlink(tmp[j])
      }
      got <- need[okk]; ex$mb[got] <- file.size(ex$path[got]) / 1e6
    }
    for (k in idx) {                          # keep newest-first order within the budget
      if (is.na(ex$mb[k])) next               # missing / failed file
      if (used_mb + ex$mb[k] > budget_mb || sum(ex$used) >= max_experiments) { used_mb <- budget_mb; break }
      ex$used[k] <- TRUE; used_mb <- used_mb + ex$mb[k]
    }
    i <- max(idx) + 1L
  }
  chrom <- if (grepl("^chr", chr)) chr else paste0("chr", chr)
  u  <- ex[ex$used, , drop = FALSE]
  pk <- if (nrow(u)) do.call(rbind, lapply(seq_len(nrow(u)), function(k)
          .read_each_bed_region(u$path[k], chrom, from, to, u$srx[k], u$cell_type[k])))
        else .empty_peaks()
  avg <- if (nrow(u)) mean(u$mb) else NA
  message(sprintf("  ChIP-Atlas %s: used %d of %d experiment(s) (%.0f MB; %d newly downloaded)%s.",
                  tf, nrow(u), nrow(ex), sum(u$mb), n_dl,
                  if (nrow(u) < nrow(ex) && is.finite(avg))
                    sprintf("; all %d would be ~%.0f MB - raise chip_max_mb to use more", nrow(ex), avg * nrow(ex))
                  else ""))
  list(peaks = pk, available = nrow(u) > 0, mode = "experiments",
       n_total = nrow(ex), n_used = nrow(u), mb_used = sum(u$mb),
       experiments = ex[, c("srx", "cell_type", "accession_number", "mb", "used")])
}
