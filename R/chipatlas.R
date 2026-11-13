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
## Pulldown type (chip_pulldown): TF-antibody ("native") and epitope-tag pulldowns
## ("tagged": FLAG, HA, V5, GFP, HALO, biotin ...) are told apart from ChIP-Atlas
## metadata; "all" (default) uses both, including tag pulldowns that ChIP-Atlas
## files under "Epitope tags" / "GFP" / "Biotin" (see build_chip_pulldown_table).
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
             pulldown = character(0), stringsAsFactors = FALSE)

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
  df$pulldown <- NA_character_
  df
}

## ChIP-Atlas peaks for one TF in a genomic interval (all experiments merged).
## threshold: MACS2 q-value class "05", "10" or "20" (q < 1e-5, 1e-10, 1e-20).
## cell_pattern: optional regex matched against the ChIP-Atlas cell type name
## (or accession), e.g. "Hep G2|Hepatocytes", for peaks and experiment lists alike.
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
#' @param cell_pattern Optional regular expression matched against the
#'   ChIP-Atlas cell type name (e.g. \code{"Hep G2|Hepatocytes"}) or experiment
#'   accession; spaces, "-" and "_" are ignored. ChIP-Atlas names some lines
#'   differently from common usage (501mel = "501A"); see
#'   \code{table(list_chip_experiments(tf)$cell_type)}.
#' @param max_mb Data budget in MB (default 300); see Details.
#' @param mode `"auto"` (default), `"assembled"` or `"experiments"`; see Details.
#' @param max_experiments Maximum number of per-experiment files (default `Inf`;
#'   a finite value switches to per-experiment files).
#' @param pulldown Which pulldowns count as the TF's ChIP-seq: `"all"` (default;
#'   TF-antibody and epitope-tag pulldowns, including tag pulldowns that
#'   ChIP-Atlas files under "Epitope tags", "GFP" or "Biotin"), `"native"`
#'   (TF-antibody pulldowns only) or `"tagged"` (tag pulldowns only). Set for a
#'   whole session with `options(promotif.chip_pulldown = ...)`. See
#'   [build_chip_pulldown_table()].
#' @return A list with \code{peaks} (data frame with \code{summit} and
#'   \code{pulldown}),
#'   \code{available}, \code{mode}, \code{n_total} and \code{n_used}
#'   (experiments), or \code{NULL} if ChIP-Atlas was not reachable.
#' @export
fetch_chip_peaks <- function(tf, chr, from, to, species = "homo_sapiens",
                             threshold = "05", cell_pattern = NULL, max_mb = 300,
                             mode = c("auto", "assembled", "experiments"),
                             max_experiments = Inf,
                             pulldown = getOption("promotif.chip_pulldown", "all")) {
  mode   <- match.arg(mode); pulldown <- .check_pulldown(pulldown)
  genome <- .chipatlas_genome(species, coordinates = TRUE)
  if (is.null(genome)) return(NULL)
  url  <- sprintf("%s/%s/assembled/Oth.ALL.%s.%s.AllCell.bed",
                  .chipatlas_base, genome, threshold, tf)
  if (mode == "experiments" || is.finite(max_experiments))
    return(fetch_chip_peaks_budget(tf, chr, from, to, species, threshold, cell_pattern,
                                   budget_mb = max_mb, max_experiments = max_experiments,
                                   pulldown = pulldown))
  if (mode == "auto") {
    dest <- file.path(.chip_cache_dir(),
                      gsub("/", "__", sub(paste0("^", .chipatlas_base, "/"), "", url)))
    if (!(file.exists(dest) && file.size(dest) > 0)) {
      h  <- tryCatch(httr::HEAD(url, httr::timeout(30)), error = function(e) NULL)
      if (is.null(h)) { message("  ChIP-Atlas not reachable - ChIP layer skipped."); return(NULL) }
      if (httr::status_code(h) == 404) {
        ## no experiment filed under the TF: tag pulldowns of it may still exist
        if (pulldown != "native" && nrow(.chip_pulldown_of(tf, genome)$extra))
          return(fetch_chip_peaks_budget(tf, chr, from, to, species, threshold, cell_pattern,
                                         budget_mb = max_mb, pulldown = pulldown))
        message("  ChIP-Atlas has no ", tf, " ChIP-seq peak file (", basename(url), ").")
        return(list(peaks = .empty_peaks(), available = FALSE, mode = "assembled",
                    n_total = 0L, n_used = 0L))
      }
      mb <- suppressWarnings(as.numeric(httr::headers(h)[["content-length"]])) / 1e6
      if (length(mb) == 1 && is.finite(mb) && mb > max_mb) {
        message(sprintf("  %s is %.0f MB (> chip_max_mb = %s): switching to per-experiment files, newest first.",
                        basename(url), mb, format(max_mb)))
        return(fetch_chip_peaks_budget(tf, chr, from, to, species, threshold, cell_pattern,
                                       budget_mb = max_mb, pulldown = pulldown))
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
    ## matched on the ChIP-Atlas cell type (and accession), as in the experiment
    ## list, so that "experiments with a peak" and "experiments queried" use the
    ## same experiments
    pk <- pk[.cell_match(ifelse(is.na(pk$cell_type), "", pk$cell_type), cell_pattern) |
             .cell_match(pk$srx, cell_pattern), , drop = FALSE]
  ## pulldown type of each experiment; keep the requested ones
  pd <- .chip_pulldown_of(tf, .chipatlas_genome(species, coordinates = FALSE))
  pk$pulldown <- ifelse(pk$srx %in% pd$tagged, "tagged", "native")
  if (pulldown != "all") pk <- pk[pk$pulldown == pulldown, , drop = FALSE]
  ## number of experiments queried (for the support fraction): all of the TF's
  ## experiments (matching cell_pattern and pulldown), from its experiment list
  ex <- tryCatch(list_chip_experiments(tf, species, cell_pattern, pulldown = pulldown),
                 error = function(e) NULL)
  ## tag pulldowns filed under "Epitope tags" / "GFP" / "Biotin" are not in the
  ## TF's merged file: their per-experiment files are added (few, small, cached)
  mb_extra <- 0
  if (!is.null(ex) && pulldown != "native") {
    xe <- ex[ex$srx %in% pd$extra$srx, , drop = FALSE]
    if (nrow(xe)) {
      fx <- .each_files(xe, genome, threshold, budget_mb = Inf)
      u  <- fx$ex[fx$ex$used, , drop = FALSE]
      if (nrow(u)) {
        px <- do.call(rbind, lapply(seq_len(nrow(u)), function(k)
                .read_each_bed_region(u$path[k], chrom, from, to, u$srx[k], u$cell_type[k])))
        px$pulldown <- rep("tagged", nrow(px))
        pk <- rbind(pk, px); mb_extra <- sum(u$mb)
      }
    }
  }
  n_tot <- max(if (is.null(ex)) 0L else nrow(ex), length(unique(pk$srx)))
  list(peaks = pk, available = TRUE, mode = "assembled", n_total = n_tot, n_used = n_tot,
       n_tagged = if (is.null(ex)) NA_integer_ else sum(ex$pulldown == "tagged"),
       mb_used = file.size(path) / 1e6 + mb_extra)
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
#' @param pulldown `"all"` (default), `"native"` or `"tagged"`; see [fetch_chip_peaks()].
#'   `chip_n_exp_tagged` counts the tag pulldowns among `chip_n_exp`.
#' @return The scan, annotated (unchanged if ChIP-Atlas was unavailable).
#' @export
add_chip_peaks <- function(scan, threshold = "05", cell_pattern = NULL, max_mb = 300,
                           mode = c("auto", "assembled", "experiments"),
                           max_experiments = Inf, sigma = 50,
                           support = c("locus", "all"), prior = 1, rep_half = 20,
                           pulldown = getOption("promotif.chip_pulldown", "all")) {
  mode <- match.arg(mode); support <- match.arg(support); pulldown <- .check_pulldown(pulldown)
  m <- scan$meta; prom <- scan$promoters
  from <- min(prom$win_start); to <- max(prom$win_end)
  message("  ChIP-Atlas: ", m$tf, " peaks (MACS2 q < 1e-", threshold, ")",
          if (!is.null(cell_pattern)) paste0(", cells matching '", cell_pattern, "'") else "",
          c(all = ", TF-antibody and tag pulldowns", native = ", TF-antibody pulldowns only",
            tagged = ", tag pulldowns only")[[pulldown]], " ...")
  ## heterodimers ("MAX::MYC"): peaks of both partners are combined
  parts <- unique(trimws(strsplit(m$tf, "::", fixed = TRUE)[[1]]))
  if (identical(m$species, "homo_sapiens")) parts <- toupper(parts)
  if (length(parts) > 1)
    message("  heterodimer: combining ChIP-Atlas peaks of ", paste(parts, collapse = " and "))
  cks <- Filter(Negate(is.null), lapply(parts, function(tf1)
    fetch_chip_peaks(tf1, m$chr, from, to, m$species, threshold, cell_pattern, max_mb,
                     mode = mode, max_experiments = max_experiments, pulldown = pulldown)))
  if (!length(cks)) return(scan)
  ck <- list(peaks = do.call(rbind, lapply(cks, `[[`, "peaks")),
             available = any(vapply(cks, `[[`, logical(1), "available")))
  n_queried <- sum(vapply(cks, function(x) as.numeric(x$n_used %||% 0), numeric(1)))
  modes <- unique(vapply(cks, function(x) x$mode %||% "assembled", ""))
  pk <- ck$peaks; res <- scan$result
  ## peaks are fetched for the whole span from the first to the last promoter
  ## window; keep only those overlapping a promoter window, so that peaks between
  ## distant promoters (e.g. enhancers in a long intron) neither count as
  ## "experiments with a peak here" nor enlarge n_locus
  if (nrow(pk) && nrow(prom)) {
    inw <- vapply(seq_len(nrow(pk)), function(j)
      any(pk$start[j] <= prom$win_end & pk$end[j] >= prom$win_start), logical(1))
    pk <- pk[inw, , drop = FALSE]
  }
  if (nrow(res)) {
    ctr <- res$site_center %||% res$genomic_position
    hit <- lapply(ctr, function(g) which(pk$start <= g & pk$end >= g))
    res$chip_n_exp      <- vapply(hit, function(h) length(unique(pk$srx[h])), integer(1))
    res$chip_n_exp_tagged <- vapply(hit, function(h)
      length(unique(pk$srx[h][pk$pulldown[h] %in% "tagged"])), integer(1))
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
  n_tag <- length(unique(pk$srx[pk$pulldown %in% "tagged"]))
  if (ck$available)
    message(sprintf("  ChIP-Atlas: %d peak(s) from %d of %d queried experiment(s)%s overlap the promoter windows; %d of %d site row(s) lie inside a peak.",
                    nrow(pk), n_exp, as.integer(n_queried),
                    if (n_tag) sprintf(" (%d tag pulldown(s))", n_tag) else "",
                    sum(res$chip_supported %in% TRUE), nrow(res)))
  scan$chip <- list(peaks = pk, threshold = threshold, cell_pattern = cell_pattern,
                    n_experiments = n_exp, n_queried = n_queried, available = ck$available,
                    mode = paste(modes, collapse = "+"), max_mb = max_mb, sigma = sigma,
                    support = support, prior = prior, pulldown = pulldown, n_experiments_tagged = n_tag,
                    n_queried_tagged = sum(vapply(cks, function(x) as.numeric(x$n_tagged %||% 0), numeric(1)), na.rm = TRUE),
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
## Genes in the screen's universe (as the promoter background): the share of
## all genes bound by a TF in at least as many experiments as this gene is
## chip_target_q (rank / universe).
.screen_universe <- 20094L
## Is a ChIP-Atlas file already in the cache? / Its size in MB (HEAD; NA if unknown)
.chipatlas_cached <- function(url) {
  dest <- file.path(.chip_cache_dir(),
                    gsub("/", "__", sub(paste0("^", .chipatlas_base, "/"), "", url)))
  file.exists(dest) && file.size(dest) > 0
}
.chipatlas_size_mb <- function(url) {
  h <- tryCatch(httr::HEAD(url, httr::timeout(30)), error = function(e) NULL)
  if (is.null(h) || httr::status_code(h) != 200) return(NA_real_)
  mb <- suppressWarnings(as.numeric(httr::headers(h)[["content-length"]])) / 1e6
  if (length(mb) == 1 && is.finite(mb)) mb else NA_real_
}

## ChIP screen within a budget: when a TF's Target Genes file is larger than the
## budget (e.g. CTCF, 113 MB at +/- 5 kb), its experiments are used newest first,
## as for the ChIP-seq scans: per-experiment peak files (cached) are read until
## budget_mb is used, and an experiment counts as bound when a peak lies within
## +/- distance kb of one of the gene's TSSs. The genome-wide target rank needs
## the whole file and is left NA.
.chip_target_row_budget <- function(tf, region, distance, cell_pattern, budget_mb,
                                    pulldown, file_mb, out) {
  ex <- list_chip_experiments(tf, region$species, cell_pattern, pulldown = pulldown)
  if (is.null(ex)) return(out)
  if (!nrow(ex)) { out$chip_available <- FALSE; out$chip_n_exp <- 0L; return(out) }
  genome <- .chipatlas_genome(region$species, coordinates = TRUE)
  if (is.null(genome)) return(out)
  message(sprintf("  %s Target Genes file is %.0f MB (> %s MB): screening its newest experiments within the budget instead ...",
                  tf, file_mb, format(budget_mb)))
  fx <- .each_files(ex, genome, "05", budget_mb)
  u  <- fx$ex[fx$ex$used, , drop = FALSE]
  if (!nrow(u)) return(out)
  d     <- as.numeric(distance) * 1000
  lo    <- region$tss_min - d; hi <- region$tss_max + d
  chrom <- paste0("chr", sub("^chr", "", region$chr))
  bound <- vapply(seq_len(nrow(u)), function(k) {
    p <- .read_each_bed_region(u$path[k], chrom, min(lo), max(hi), u$srx[k], u$cell_type[k])
    nrow(p) > 0 && any(vapply(seq_along(lo), function(j) any(p$end >= lo[j] & p$start <= hi[j]), logical(1)))
  }, logical(1))
  out$chip_available <- TRUE
  out$chip_n_exp     <- nrow(u)
  out$chip_n_bound   <- sum(bound)
  out$chip_frac_exp  <- mean(bound)
  out$chip_source    <- sprintf("newest %d of %d experiments (%s MB budget)", nrow(u), nrow(ex), format(budget_mb))
  out
}

## Where the gene's TSSs are, for the budgeted ChIP screen
.screen_region <- function(disc) {
  m <- disc$meta; pr <- disc$promoters
  if (is.null(m$chr) || is.null(pr) || !nrow(pr)) return(NULL)
  list(species = m$species %||% "homo_sapiens", chr = m$chr,
       tss_min = pr$tss_min %||% pr$tss, tss_max = pr$tss_max %||% pr$tss)
}

.chip_target_row <- function(tf, gene, genome, distance = 1, cell_pattern = NULL,
                             max_mb = 300, pulldown = getOption("promotif.chip_pulldown", "all"),
                             region = NULL) {
  pulldown <- .check_pulldown(pulldown)
  out <- data.frame(chip_tf_queried = tf, chip_available = NA, chip_n_exp = NA_integer_,
                    chip_n_bound = NA_integer_, chip_frac_exp = NA_real_,
                    chip_mean_score = NA_real_, chip_max_score = NA_real_,
                    chip_target_rank = NA_integer_, chip_target_q = NA_real_,
                    chip_source = NA_character_, stringsAsFactors = FALSE)
  url  <- sprintf("%s/%s/target/%s.%d.tsv", .chipatlas_base, genome, tf, as.integer(distance))
  ## Target Genes file larger than the budget and not cached: use the TF's newest
  ## experiments within the budget instead of skipping the TF (needs the gene's TSSs)
  if (!is.null(region) && !.chipatlas_cached(url)) {
    mb <- .chipatlas_size_mb(url)
    if (isTRUE(mb > max_mb))
      return(.chip_target_row_budget(tf, region, distance, cell_pattern, max_mb, pulldown, mb, out))
  }
  path <- .chipatlas_get(url, max_mb)
  if (is.null(path)) return(out)
  if (is.na(path)) { out$chip_available <- FALSE; return(out) }
  out$chip_source <- "Target Genes file"
  con <- file(path, "r"); on.exit(close(con))
  hdr <- strsplit(readLines(con, n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]]
  exp_cols <- grep("RX[0-9]+", hdr)             # SRX/ERX/DRX experiment columns
  if (!length(exp_cols)) {
    message("  unexpected Target Genes header for ", tf, " - skipped."); return(out)
  }
  if (!is.null(cell_pattern))
    exp_cols <- exp_cols[.cell_match(hdr[exp_cols], cell_pattern)]
  ## pulldown type: the Target Genes table holds the experiments filed under the
  ## TF; tag pulldowns filed under "Epitope tags" / "GFP" / "Biotin" are not in it
  if (pulldown != "all") {
    tg <- sub("\\|.*$", "", hdr[exp_cols]) %in% .chip_pulldown_of(tf, genome)$tagged
    exp_cols <- exp_cols[if (pulldown == "tagged") tg else !tg]
  }
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
    ## where this gene ranks among the TF's targets genome-wide: experiments
    ## bound at every listed gene (counted once per file and cell filter, cached)
    cpath <- sub("\\.tsv$", paste0(if (is.null(cell_pattern)) ".counts" else
                 sprintf(".counts_%s", substr(gsub("[^A-Za-z0-9]+", "_", cell_pattern), 1, 40)),
                 if (pulldown == "all") "" else paste0("_", pulldown), ".rds"), path)
    cnt <- if (file.exists(cpath)) tryCatch(readRDS(cpath), error = function(e) NULL) else NULL
    if (is.null(cnt)) {
      ln <- readLines(path, warn = FALSE)[-1]
      sp <- strsplit(ln, "\t", fixed = TRUE)
      cnt <- vapply(sp, function(v) sum(suppressWarnings(as.numeric(v[exp_cols])) > 0, na.rm = TRUE), numeric(1))
      names(cnt) <- vapply(sp, `[`, "", 1)
      tryCatch(saveRDS(cnt, cpath), error = function(e) NULL)
    }
    g <- out$chip_n_bound
    if (g > 0) {
      out$chip_target_rank <- as.integer(sum(cnt >= g))
      out$chip_target_q    <- min(1, out$chip_target_rank / .screen_universe)
    } else out$chip_target_q <- 1
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
#' @details Target Genes tables contain the experiments ChIP-Atlas files under
#'   the TF: with `pulldown = "all"` these are TF-antibody and tag pulldowns
#'   filed under the TF; tag pulldowns filed under "Epitope tags" are not included.
#' @return \code{disc} with \code{chip_*} columns added to \code{top} and
#'   \code{result}.
#' @export
add_chip_targets <- function(disc, distance = 1, cell_pattern = NULL, max_mb = 300,
                             pulldown = getOption("promotif.chip_pulldown", "all")) {
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
                                   max_mb = max_mb, pulldown = pulldown,
                                   region = .screen_region(disc)))
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
  base_theme <- theme_bw(base_size = .tsz(base_size)) +
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
    ggrepel::geom_text_repel(data = lab, aes(x, y, label = l), size = .tsz(base_size) / 4,
                             min.segment.length = 0, max.overlaps = Inf, seed = 1)
  else geom_text(data = lab, aes(x, y, label = l), vjust = -0.7, size = .tsz(base_size) / 4,
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
#' @param pulldown `"all"` (default), `"native"` or `"tagged"`; see
#'   [fetch_chip_peaks()]. With `"all"` or `"tagged"`, tag pulldowns of the TF
#'   that ChIP-Atlas files under "Epitope tags", "GFP" or "Biotin" are added.
#' @return A data frame with `srx`, `cell_type`, `pulldown` (`"native"` or
#'   `"tagged"`) and `accession_number`, or `NULL` if ChIP-Atlas could not be reached.
#' @export
#' @examples
#' \dontrun{
#' head(list_chip_experiments("ESR1", cell_pattern = "MCF-7"))
#' }
list_chip_experiments <- function(tf, species = "homo_sapiens", cell_pattern = NULL,
                                  pulldown = getOption("promotif.chip_pulldown", "all")) {
  pulldown <- .check_pulldown(pulldown)
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
      if (httr::status_code(r) == 404) { hdr <- ""; break }   # nothing filed under the TF
      txt <- rawToChar(httr::content(r, as = "raw"))
      if (grepl("\n", txt, fixed = TRUE)) { hdr <- sub("\n.*$", "", txt); break }
    }
    if (is.null(hdr)) { message("  could not read the experiment list of ", tf, "."); return(NULL) }
    f  <- if (nzchar(hdr)) strsplit(hdr, "\t", fixed = TRUE)[[1]] else character(0)
    f  <- f[grepl("^[SDE]RX[0-9]+", f)]
    ex <- data.frame(srx = sub("\\|.*$", "", f),
                     cell_type = gsub("_", " ", sub("^[^|]*\\|?", "", f)),
                     stringsAsFactors = FALSE)
    ex <- ex[!duplicated(ex$srx), , drop = FALSE]
    utils::write.table(ex, dest, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  ## pulldown type; tag pulldowns filed under "Epitope tags" / "GFP" / "Biotin"
  pd <- .chip_pulldown_of(tf, genome)
  ex$pulldown <- ifelse(ex$srx %in% pd$tagged, "tagged", "native")
  if (pulldown != "native" && nrow(pd$extra))
    ex <- rbind(ex, data.frame(srx = pd$extra$srx, cell_type = pd$extra$cell_type,
                               pulldown = rep("tagged", nrow(pd$extra)), stringsAsFactors = FALSE))
  ex <- ex[!duplicated(ex$srx), , drop = FALSE]
  if (pulldown != "all") ex <- ex[ex$pulldown == pulldown, , drop = FALSE]
  if (!nrow(ex)) message("  ChIP-Atlas has no ", if (pulldown == "all") "" else paste0(pulldown, " "),
                         "experiments for ", tf, ".")
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
             summit_exact = !is.na(sm), pulldown = NA_character_, stringsAsFactors = FALSE)
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
#' @param pulldown `"all"` (default), `"native"` or `"tagged"`; see [fetch_chip_peaks()].
#' @return A list with `peaks` (data frame incl. `summit` and `pulldown`), `available`,
#'   `mode = "experiments"`, `n_total`, `n_used`, `n_tagged`, `mb_used` and `experiments`
#'   (one row per experiment, with `pulldown` and `used`).
#' @export
#' @examples
#' \dontrun{
#' pk <- fetch_chip_peaks_budget("ESR1", "21", 42366000, 42369000,
#'                               cell_pattern = "MCF-7", budget_mb = 300)
#' }
fetch_chip_peaks_budget <- function(tf, chr, from, to, species = "homo_sapiens",
                                    threshold = "05", cell_pattern = NULL,
                                    budget_mb = 1000, max_experiments = Inf,
                                    batch = 24L,
                                    pulldown = getOption("promotif.chip_pulldown", "all")) {
  pulldown <- .check_pulldown(pulldown)
  genome <- .chipatlas_genome(species, coordinates = TRUE)
  if (is.null(genome)) return(NULL)
  ex <- list_chip_experiments(tf, species, cell_pattern, pulldown = pulldown)
  if (is.null(ex)) return(NULL)
  if (!nrow(ex)) {
    message("  no ", tf, " experiments", if (!is.null(cell_pattern))
            paste0(" matching '", cell_pattern, "'") else "", ".")
    return(list(peaks = .empty_peaks(), available = FALSE, mode = "experiments",
                n_total = 0L, n_used = 0L, n_tagged = 0L, mb_used = 0, experiments = ex))
  }
  message(sprintf("  ChIP-Atlas %s: %d experiment(s)%s%s; using newest first up to %s MB%s ...",
                  tf, nrow(ex), if (any(ex$pulldown == "tagged")) sprintf(" (%d tag pulldown(s))", sum(ex$pulldown == "tagged")) else "",
                  if (!is.null(cell_pattern)) paste0(" matching '", cell_pattern, "'") else "",
                  format(budget_mb), if (is.finite(max_experiments))
                    sprintf(" / %d experiments", as.integer(max_experiments)) else ""))
  fx <- .each_files(ex, genome, threshold, budget_mb, max_experiments, batch)
  ex <- fx$ex
  chrom <- if (grepl("^chr", chr)) chr else paste0("chr", chr)
  u  <- ex[ex$used, , drop = FALSE]
  pk <- if (nrow(u)) do.call(rbind, lapply(seq_len(nrow(u)), function(k) {
          p <- .read_each_bed_region(u$path[k], chrom, from, to, u$srx[k], u$cell_type[k])
          p$pulldown <- rep(u$pulldown[k], nrow(p)); p }))
        else .empty_peaks()
  avg <- if (nrow(u)) mean(u$mb) else NA
  message(sprintf("  ChIP-Atlas %s: used %d of %d experiment(s) (%.0f MB; %d newly downloaded)%s.",
                  tf, nrow(u), nrow(ex), sum(u$mb), fx$n_dl,
                  if (nrow(u) < nrow(ex) && is.finite(avg))
                    sprintf("; all %d would be ~%.0f MB - raise chip_max_mb to use more", nrow(ex), avg * nrow(ex))
                  else ""))
  list(peaks = pk, available = nrow(u) > 0, mode = "experiments",
       n_total = nrow(ex), n_used = nrow(u), n_tagged = sum(u$pulldown == "tagged"), mb_used = sum(u$mb),
       experiments = ex[, c("srx", "cell_type", "pulldown", "accession_number", "mb", "used")])
}

## Download per-experiment peak files (newest first, cached) until budget_mb or
## max_experiments is reached; ex gains path, mb and used.
.each_files <- function(ex, genome, threshold = "05", budget_mb = Inf, max_experiments = Inf,
                        batch = 24L) {
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
  list(ex = ex, n_dl = n_dl)
}

## ---- pulldown type of ChIP-Atlas experiments: TF antibody or epitope tag ------
## ChIP-Atlas files an experiment under the protein it curated as the antigen.
## Pulldowns of a tagged protein (FLAG-, HA-, V5-, GFP-, HALO-, biotin-tagged)
## are mostly filed under "Epitope tags", "GFP" or "Biotin", and some under the
## TF itself. The pulldown table lists every tagged experiment of a TF:
##   tagged pulldowns filed under the TF (antibody metadata names a tag), and
##   tag-filed pulldowns whose title names the TF next to a tag word (or names
##   only that TF). Experiments of a TF not in the table are TF-antibody
##   ("native") pulldowns.
## chip_pulldown = "all" (default): native + tagged experiments;
##   "native": TF-antibody pulldowns only; "tagged": tagged pulldowns only.
## The table ships with the tool (built from ChIP-Atlas experimentList.tab with
## build_chip_pulldown_table()); update_chip_pulldown() rebuilds it in the cache.

.tag_words <- c(flag = "FLAG", "3xflag" = "FLAG", "3x flag" = "FLAG", ha = "HA", "3xha" = "HA",
                hemagglutinin = "HA", v5 = "V5", "myc-tag" = "MYC-tag", "myc tag" = "MYC-tag",
                gfp = "GFP", egfp = "GFP", yfp = "GFP", venus = "GFP", mcherry = "mCherry",
                mneongreen = "GFP", halo = "HALO", halotag = "HALO", bira = "BIOTIN",
                biotin = "BIOTIN", biotinylated = "BIOTIN", avitag = "BIOTIN", "avi-tag" = "BIOTIN",
                t7 = "T7", "t7-tag" = "T7", strep = "STREP", "strep-tag" = "STREP",
                "twin-strep" = "STREP", tagged = "tag")
.tag_rx <- paste0("(?i)(?<![A-Za-z0-9])(", paste(gsub("([-+.])", "\\\\\\1", names(.tag_words)), collapse = "|"),
                  ")(?![A-Za-z0-9])")
## knock-down / control wording that mentions a tag protein without pulling it down;
## input / IgG / mock samples filed under a tag are not pulldowns of a TF
.tag_noise_rx <- "(?i)(sh|si)[-_ ]?(gfp|egfp)|gfp[-_ ](knock|control|kd)|epitope[- ]mapp"
.not_pulldown_rx <- "(?i)(?<![a-z])(input|igg|mock|empty[- ]?vector|untransfected|no[- ]antibody|luciferase|luc)(?![a-z])"

.first_tag <- function(x) {
  m <- regmatches(x, regexpr(.tag_rx, x, perl = TRUE))
  if (!length(m)) return(NA_character_)
  unname(.tag_words[tolower(m)])
}

#' Build the ChIP-Atlas pulldown table (tagged experiments of each TF)
#'
#' Reads ChIP-Atlas `experimentList.tab` and lists, for one genome, every
#' experiment that pulled down an epitope-tagged TF: experiments filed under
#' the TF whose antibody is a tag antibody, and experiments filed under
#' "Epitope tags", "GFP" or "Biotin" whose title names the TF. Used by the
#' ChIP-Atlas layer through `chip_pulldown`.
#'
#' @param experiment_list Path or URL of `experimentList.tab` (default: the
#'   ChIP-Atlas URL; about 400 MB, read line by line).
#' @param genome `"hg38"` or `"mm10"`.
#' @param out Output file (`.tsv` or `.tsv.gz`); `NULL` = return the table only.
#' @return A data frame: `srx`, `filed_as` (ChIP-Atlas antigen), `target`
#'   (TF), `tag`, `match` (`"antibody"`, `"title-adjacent"` or `"title-only"`),
#'   `cell_type`, `title`.
#' @export
build_chip_pulldown_table <- function(experiment_list = "https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab",
                                      genome = "hg38", out = NULL) {
  con <- if (grepl("^https?://", experiment_list)) url(experiment_list, "r") else file(experiment_list, "r")
  on.exit(close(con))
  keep <- character(0)
  repeat {
    ln <- readLines(con, n = 200000L, warn = FALSE, encoding = "UTF-8")
    if (!length(ln)) break
    ln <- ln[startsWith(ln, "") & grepl(paste0("^[^\t]+\t", genome, "\tTFs and others\t"), ln, useBytes = TRUE)]
    keep <- c(keep, ln)
  }
  f <- strsplit(keep, "\t", fixed = TRUE)
  srx  <- vapply(f, `[`, "", 1); ag <- vapply(f, `[`, "", 4); ct <- vapply(f, `[`, "", 6)
  ttl  <- vapply(f, function(x) if (length(x) >= 9) x[9] else "", "")
  meta <- vapply(f, function(x) if (length(x) >= 10) paste(x[10:length(x)], collapse = " | ") else "", "")
  ttl  <- iconv(ttl, "UTF-8", "ASCII", sub = " "); meta <- iconv(meta, "UTF-8", "ASCII", sub = " ")
  tagfiled <- ag %in% c("Epitope tags", "GFP", "Biotin")
  ## 1. tagged pulldowns filed under a TF: the antibody field names a tag
  ##    (or, without an antibody field, the title / metadata do)
  ab <- vapply(strsplit(meta, " | ", fixed = TRUE), function(v) {
    a <- v[grepl("antibod|^ab=|chip_ab", v, ignore.case = TRUE)]
    if (length(a)) paste(a, collapse = " ") else NA_character_ }, "")
  txt1 <- ifelse(is.na(ab), paste(ttl, meta), ab)
  txt1 <- gsub(.tag_noise_rx, " ", txt1, perl = TRUE)
  ## "Myc" is a tag word except in MYC / MYCN pulldowns (9E10 recognises human MYC)
  own_tag <- vapply(seq_along(txt1), function(i) {
    if (tagfiled[i]) return(NA_character_)
    t <- txt1[i]
    if (toupper(ag[i]) %in% c("MYC", "MYCN")) t <- gsub("(?i)myc[- ]tag", " ", t, perl = TRUE)
    .first_tag(t) }, "")
  own <- which(!tagfiled & !is.na(own_tag) & own_tag != "tag")
  ## 2. tag-filed pulldowns: the TF named in the title (or the source name)
  syms <- setdiff(unique(ag[!tagfiled]), "")
  syms <- syms[grepl("^[A-Za-z0-9][A-Za-z0-9.-]*$", syms) & nchar(syms) >= 2]
  sym_uc <- toupper(syms)
  tok_split <- function(s) { s <- gsub("[^A-Za-z0-9.]+", " ", s); s <- trimws(gsub("\\.+( |$)", " ", s))
    strsplit(s, " +")[[1]] }
  find_target <- function(s, own_tag_word) {
    tk <- tok_split(s); if (!length(tk)) return(NULL)
    ## short symbols (<= 3 characters) must match in case, longer ones in any case
    is_tf <- vapply(tk, function(x) {
      u <- toupper(x); j <- match(u, sym_uc)
      !is.na(j) && (nchar(u) > 3 || identical(x, syms[j]) || identical(x, u)) }, logical(1))
    is_tag <- grepl(.tag_rx, tk, perl = TRUE) | tolower(tk) %in% c("tag", "tagged", "myc")
    is_tf[is_tag & toupper(tk) %in% c("MYC", "HA", "T7")] <- FALSE
    if (!any(is_tf)) return(NULL)
    cand <- unique(toupper(tk[is_tf]))
    if (any(is_tag)) {
      d <- vapply(which(is_tf), function(i) min(abs(i - which(is_tag))), numeric(1))
      near <- unique(toupper(tk[which(is_tf)[d <= 2]]))
      if (length(near) == 1) return(c(syms[match(near, sym_uc)], "title-adjacent"))
      if (length(near) > 1) return(NULL)
    }
    if (length(cand) == 1) return(c(syms[match(cand, sym_uc)], "title-only"))
    NULL
  }
  ext <- lapply(which(tagfiled & !grepl(.not_pulldown_rx, ttl, perl = TRUE)), function(i) {
    r <- find_target(ttl[i])
    if (is.null(r)) {
      sn <- regmatches(meta[i], regexpr("source_name=[^|]*", meta[i]))
      if (length(sn)) r <- find_target(sn)
    }
    if (is.null(r)) return(NULL)
    tg <- .first_tag(paste(ttl[i], ab[i] %||% "", meta[i]))
    if (is.na(tg) || tg == "tag") tg <- c("Epitope tags" = "tag", GFP = "GFP", Biotin = "BIOTIN")[[ag[i]]]
    data.frame(srx = srx[i], filed_as = ag[i], target = r[1], tag = tg, match = r[2],
               cell_type = ct[i], title = substr(ttl[i], 1, 120), stringsAsFactors = FALSE)
  })
  tab <- rbind(
    data.frame(srx = srx[own], filed_as = ag[own], target = ag[own], tag = own_tag[own],
               match = "antibody", cell_type = ct[own], title = substr(ttl[own], 1, 120),
               stringsAsFactors = FALSE),
    do.call(rbind, ext))
  tab <- tab[order(tab$target, tab$srx), , drop = FALSE]; rownames(tab) <- NULL
  attr(tab, "genome") <- genome
  if (!is.null(out)) {
    con2 <- if (grepl("\\.gz$", out)) gzfile(out, "w") else file(out, "w")
    utils::write.table(tab, con2, sep = "\t", quote = FALSE, row.names = FALSE)
    close(con2)
  }
  tab
}

## The pulldown table of a genome: a rebuilt copy in the cache wins over the
## one shipped with the tool. Read once per session.
.pulldown_env <- new.env(parent = emptyenv())
.chip_pulldown_table <- function(genome) {
  if (!is.null(.pulldown_env[[genome]])) return(.pulldown_env[[genome]])
  fn   <- sprintf("chipatlas_pulldown_%s.tsv.gz", genome)
  cand <- c(file.path(.chip_cache_dir(), fn),
            system.file("extdata", fn, package = "proMotif"),
            file.path(getOption("promotif.bg_dir", "."), fn))
  f <- cand[nzchar(cand) & file.exists(cand)][1]
  tab <- if (is.na(f)) {
    message("  ChIP-Atlas pulldown table ", fn, " not found - tagged pulldowns not recognised.")
    data.frame(srx = character(0), filed_as = character(0), target = character(0), tag = character(0),
               match = character(0), cell_type = character(0), title = character(0))
  } else utils::read.delim(f, stringsAsFactors = FALSE, quote = "", comment.char = "",
                           colClasses = "character")
  .pulldown_env[[genome]] <- tab
  tab
}

#' Rebuild the ChIP-Atlas pulldown table in the cache
#'
#' Downloads ChIP-Atlas `experimentList.tab` (about 400 MB) and rebuilds the
#' table of tagged pulldowns (see [build_chip_pulldown_table()]) for `genome`,
#' so experiments added to ChIP-Atlas after the tool's release are labelled.
#' @inheritParams build_chip_pulldown_table
#' @return The table, invisibly.
#' @export
update_chip_pulldown <- function(genome = c("hg38", "mm10")) {
  genome <- match.arg(genome)
  out <- file.path(.chip_cache_dir(), sprintf("chipatlas_pulldown_%s.tsv.gz", genome))
  tab <- build_chip_pulldown_table(genome = genome, out = out)
  .pulldown_env[[genome]] <- NULL
  message(sprintf("  %d tagged pulldown(s) of %d TF(s) written to %s", nrow(tab), length(unique(tab$target)), out))
  invisible(tab)
}

## Pulldown type of a TF's experiments and the extra tag-filed experiments.
## Returns list(tagged = SRXs of tagged pulldowns filed under the TF,
##              extra  = data frame srx / cell_type of tag-filed pulldowns of the TF)
.chip_pulldown_of <- function(tf, genome) {
  tab <- .chip_pulldown_table(genome)
  hit <- tab[toupper(tab$target) == toupper(tf), , drop = FALSE]
  own <- hit$filed_as == hit$target
  list(tagged = hit$srx[own],
       extra  = data.frame(srx = hit$srx[!own], cell_type = hit$cell_type[!own],
                           tag = hit$tag[!own], stringsAsFactors = FALSE))
}

.check_pulldown <- function(x) {
  x <- x %||% "all"
  if (!x %in% c("all", "native", "tagged"))
    stop("chip_pulldown must be \"all\", \"native\" or \"tagged\".", call. = FALSE)
  x
}
