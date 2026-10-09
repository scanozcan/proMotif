#' Launch the proMotif Shiny app
#'
#' Opens an interactive front-end: type a gene and TF, adjust the promoter
#' window / threshold, and view the plots. The network scan runs once per
#' "Run scan"; text/point size and colours re-style instantly.
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
