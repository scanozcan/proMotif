#' Launch the proMotif Shiny app
#'
#' Opens an interactive front-end with three modes: a gene and a TF (binding
#' sites), discovery of the likely regulators of each promoter, and comparison
#' of a gene's promoters, each with a gene track (optionally with FANTOM5 CAGE).
#' The analysis runs once per "Run"; text size, colours, the significance rings,
#' the ChIP-seq track, the promoter map, logo, usage and gene-track options
#' re-style the plots instantly.
#'
#' @param launch.browser Logical; open in the system browser (default TRUE).
#' @return Called for its side effect (runs the app); does not return.
#' @examples
#' \dontrun{
#' launch_app()
#' }
#' @export
launch_app <- function(launch.browser = TRUE) {
  if (!requireNamespace("shiny", quietly = TRUE))
    stop("The 'shiny' package is required: install.packages('shiny')")
  app_dir <- system.file("shiny", package = "proMotif")
  if (!nzchar(app_dir))
    stop("Shiny app not found - reinstall proMotif.")
  shiny::runApp(app_dir, launch.browser = launch.browser)
}
