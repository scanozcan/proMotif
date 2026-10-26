#' Clear the sequence cache
#'
#' Promoter sequences and Ensembl gene lookups are cached on disk under
#' `tools::R_user_dir("proMotif", "cache")` the first time they are fetched,
#' so repeated analyses of the same gene do not query Ensembl or UCSC again.
#' Set `options(promotif.refresh_cache = TRUE)` to re-fetch and overwrite
#' cached entries, or call this function to remove them. Motif databases and
#' ChIP-Atlas files are cached separately and are not removed.
#'
#' @return The number of removed files (invisibly).
#' @export
#' @examples
#' \dontrun{
#' clear_sequence_cache()
#' }
clear_sequence_cache <- function() {
  d <- file.path(tools::R_user_dir("proMotif", "cache"), c("promoter_seqs", "ensembl_lookup"))
  n <- length(list.files(d, full.names = TRUE))
  unlink(d, recursive = TRUE)
  message("removed ", n, " cached sequence / lookup file(s)")
  invisible(n)
}
