## ===========================================================================
## 8. MOTIF DATABASES (JASPAR 2026, HOCOMOCO v14)
## ---------------------------------------------------------------------------
## JASPAR: the release's non-redundant collection is downloaded once as one
## TRANSFAC-format file (counts + TF family/class) and cached. It serves both
## confirm mode (.jaspar_pwm) and discovery mode (.jaspar_core_pwms). There is no
## JASPAR2026 Bioconductor package, so the default path needs neither TFBSTools
## nor a JASPAR data package; jaspar_release = "2024" / "2022" also work, and a
## Bioconductor JASPAR object can still be passed as jaspar_db (discovery).
## HOCOMOCO v14 (H14CORE: human and mouse TFs in one collection) is downloaded
## once in MEME format; v11 / v12 still come from GLproxScape.
## Cache: tools::R_user_dir("proMotif", "cache").
## ===========================================================================

.motif_db_env <- new.env(parent = emptyenv())

.motif_cache_dir <- function(sub) {
  d <- file.path(tools::R_user_dir("proMotif", "cache"), sub)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

.download_once <- function(url, dest, what) {
  if (file.exists(dest) && file.size(dest) > 0) return(dest)
  message("  downloading ", what, " (once; cached) ...")
  tmp <- paste0(dest, ".part")
  ok <- .with_retry(function() {
    r <- httr::GET(url, httr::write_disk(tmp, overwrite = TRUE), httr::timeout(600))
    if (httr::status_code(r) != 200) stop("HTTP ", httr::status_code(r))
    TRUE
  }, tries = 3, waits = c(5, 20), what = what)
  if (is.null(ok)) { unlink(tmp); return(NULL) }
  file.rename(tmp, dest)
  dest
}

## JASPAR TRANSFAC export -> list of (id, name, counts 4 x L, family, class).
.parse_jaspar_transfac <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ends <- which(ln == "//")
  if (!length(ends)) return(list())
  starts <- c(1L, utils::head(ends, -1) + 1L)
  out <- lapply(seq_along(ends), function(k) {
    b <- ln[starts[k]:ends[k]]
    id   <- sub("^AC\\s+", "", b[startsWith(b, "AC ")][1])
    name <- sub("^ID\\s+", "", b[startsWith(b, "ID ")][1])
    rows <- b[grepl("^[0-9]+\\s", b)]
    if (is.na(id) || !length(rows)) return(NULL)
    m <- do.call(rbind, lapply(strsplit(rows, "\\s+"), function(x) as.numeric(x[2:5])))
    counts <- t(m); rownames(counts) <- c("A", "C", "G", "T")
    tag <- function(key) {
      v <- sub(paste0("^CC ", key, ":"), "", b[startsWith(b, paste0("CC ", key, ":"))])
      v <- v[nzchar(v)]
      if (length(v)) paste(unique(v), collapse = "/") else NA_character_
    }
    list(id = id, name = name, counts = counts,
         family = tag("tf_family"), class = tag("tf_class"))
  })
  Filter(Negate(is.null), out)
}

## One JASPAR release / collection / taxon, downloaded once and kept in memory.
.jaspar_collection <- function(release = "2026", collection = "CORE",
                               tax_group = "vertebrates") {
  key <- paste("jaspar", release, collection, tax_group, sep = "_")
  if (!is.null(.motif_db_env[[key]])) return(.motif_db_env[[key]])
  fn   <- sprintf("JASPAR%s_%s_%s_non-redundant_pfms_transfac.txt",
                  release, collection, tax_group)
  url  <- sprintf("https://jaspar.elixir.no/download/data/%s/%s/%s", release, collection, fn)
  path <- .download_once(url, file.path(.motif_cache_dir("jaspar"), fn),
                         paste("JASPAR", release, collection, tax_group))
  if (is.null(path)) return(NULL)
  mot <- .parse_jaspar_transfac(path)
  if (!length(mot)) return(NULL)
  .motif_db_env[[key]] <- mot
  mot
}

## Confirm mode: the TF's JASPAR matrix from a given release (latest version of
## an exact, case-insensitive name match). PFM -> log-odds with the same
## convention as GLproxScape::fetch_jaspar_pwm (+0.25 per cell, flat background),
## so scores are unchanged. Falls back to the JASPAR REST API (latest release)
## if the release file cannot be downloaded or does not contain the TF.
.jaspar_pwm <- function(tf, release = "2026", tax_group = "vertebrates") {
  mot <- .jaspar_collection(release, "CORE", tax_group)
  hit <- if (is.null(mot)) integer(0)
         else which(toupper(vapply(mot, `[[`, "", "name")) == toupper(tf))
  if (!length(hit)) {
    message("  ", tf, " not found in JASPAR ", release, " ", tax_group,
            " CORE file - trying the JASPAR REST API (latest release).")
    return(tryCatch(GLproxScape::fetch_jaspar_pwm(tf), error = function(e) NULL))
  }
  ids <- vapply(mot[hit], `[[`, "", "id")
  ver <- suppressWarnings(as.numeric(sub("^.*\\.", "", ids)))
  k   <- hit[order(-ver, -xtfrm(ids))][1]
  m   <- mot[[k]]$counts + 0.25
  m   <- sweep(m, 2, colSums(m), "/")
  list(id = mot[[k]]$id, name = tf, pwm = log2(m / 0.25), len = ncol(m),
       family = mot[[k]]$family)
}

## MEME motif file -> list of (id, ppm 4 x w, w).
.parse_meme <- function(path) {
  ln <- readLines(path, warn = FALSE)
  out <- list()
  for (i in grep("^MOTIF\\s", ln)) {
    id <- strsplit(trimws(sub("^MOTIF\\s+", "", ln[i])), "\\s+")[[1]][1]
    j <- i + 1L
    while (j <= length(ln) && !startsWith(ln[j], "letter-probability matrix")) j <- j + 1L
    if (j > length(ln)) next
    w <- suppressWarnings(as.integer(sub("^.*w=\\s*([0-9]+).*$", "\\1", ln[j])))
    if (is.na(w) || j + w > length(ln)) next
    ppm <- vapply(strsplit(trimws(ln[(j + 1L):(j + w)]), "\\s+"),
                  function(x) as.numeric(x[1:4]), numeric(4))
    out[[id]] <- list(id = id, ppm = ppm, w = w)
  }
  out
}

## HOCOMOCO v14 H14CORE, one matrix per TF: best quality letter (A > B > ...),
## then the lowest subtype index - the same rule GLproxScape uses for v11/v12.
.hocomoco14 <- function() {
  key <- "hocomoco_v14"
  if (!is.null(.motif_db_env[[key]])) return(.motif_db_env[[key]])
  url  <- "https://hocomoco14.autosome.org/final_bundle/hocomoco14/H14CORE/formatted_motifs/H14CORE_meme_format.meme"
  path <- .download_once(url, file.path(.motif_cache_dir("hocomoco"),
                                        "H14CORE_meme_format.meme"),
                         "HOCOMOCO v14 (H14CORE)")
  if (is.null(path)) return(NULL)
  mot <- .parse_meme(path)
  if (!length(mot)) return(NULL)
  parts <- strsplit(names(mot), ".", fixed = TRUE)
  gene  <- toupper(vapply(parts, `[`, "", 1))
  qual  <- vapply(parts, function(p) {
    q <- p[length(p)]
    if (nchar(q) == 1 && q %in% LETTERS) match(q, LETTERS) else 99L }, integer(1))
  subi  <- suppressWarnings(as.integer(vapply(parts, function(p)
    if (length(p) >= 3) p[3] else NA_character_, "")))
  o    <- order(gene, qual, subi)
  keep <- o[!duplicated(gene[o])]
  res  <- stats::setNames(mot[keep], gene[keep])
  message("  HOCOMOCO v14: ", length(res), " TFs from ", length(mot), " matrices")
  .motif_db_env[[key]] <- res
  res
}

## Confirm mode: the TF's HOCOMOCO matrix. PPM -> log-odds vs a flat background
## (same convention as GLproxScape::fetch_hocomoco_pwm).
.hocomoco_pwm <- function(tf, version = "v14", species = "human") {
  if (version %in% c("v11", "v12"))
    return(GLproxScape::fetch_hocomoco_pwm(tf, version = version, species = species))
  if (version != "v14") stop("hocomoco_version must be \"v14\", \"v12\" or \"v11\".")
  db <- .hocomoco14()
  if (is.null(db)) return(NULL)
  m <- db[[toupper(tf)]]
  if (is.null(m)) return(NULL)
  pwm <- log2(pmax(m$ppm, 1e-6) / 0.25)
  rownames(pwm) <- c("A", "C", "G", "T")
  list(id = m$id, name = tf, pwm = pwm, len = m$w)
}
