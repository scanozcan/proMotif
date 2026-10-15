## proMotif Shiny app (launched via proMotif::launch_app()).
## Modes: gene + TF, discover regulators, compare promoters. The analysis runs
## once per "Run"; size, colours and the ChIP ring threshold re-style instantly.
library(shiny)
library(proMotif)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
APP_TITLE <- "proMotif - promoter binding-site scanner (JASPAR + HOCOMOCO)"

## ---- helpers ----------------------------------------------------------------
## run an expression, collecting its console messages for the "Log" tab
.with_log <- function(expr, log) {
  withCallingHandlers(expr, message = function(m) {
    log(c(log(), sub("\n$", "", conditionMessage(m))))
    message(conditionMessage(m), appendLF = FALSE)   # still echo to the console
    invokeRestart("muffleMessage")
  })
}
.need_plot <- function(p, msg = "Nothing to show for this run.") {
  validate(need(!is.null(p), msg)); p
}

ui <- fluidPage(
  titlePanel(APP_TITLE),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      radioButtons("mode", NULL, inline = FALSE,
                   c("Gene + TF (binding sites)" = "confirm",
                     "Discover regulators (gene only)" = "discover",
                     "Compare promoters (gene only)" = "compare")),
      textInput("gene", "Gene (HGNC symbol)", "ATP7B"),
      conditionalPanel("input.mode == 'confirm'",
        textInput("tf", "Transcription factor", "MTF1")),
      fluidRow(
        column(6, numericInput("upstream",   "Upstream (bp)",   2500, min = 0, step = 100)),
        column(6, numericInput("downstream", "Downstream (bp)",  500, min = 0, step = 100))),
      conditionalPanel("input.mode != 'discover'",
        selectInput("transcripts", "Transcripts", c("protein_coding", "all"))),
      conditionalPanel("input.mode == 'confirm'",
        sliderInput("threshold", "PWM threshold (fraction of max)",
                    min = 0.60, max = 0.95, value = 0.80, step = 0.01),
        fluidRow(
          column(6, numericInput("max_transcripts", "Max TSS windows", 8, min = 1, step = 1)),
          column(6, selectInput("scan_unit", "Scan per", c("tss", "promoter"))))),
      conditionalPanel("input.mode != 'confirm'",
        numericInput("top_n", "TFs to show", 20, min = 5, step = 5)),
      conditionalPanel("input.mode == 'discover'",
        selectInput("rank_by", "Rank TFs by",
                    c("promoter background" = "background", "FIMO per-site p" = "fimo"))),
      tags$details(tags$summary("Motif databases"),
        fluidRow(
          column(6, selectInput("jaspar_release", "JASPAR", c("2026", "2024", "2022"))),
          column(6, selectInput("hocomoco_version", "HOCOMOCO", c("v14", "v12", "v11"))))),
      conditionalPanel("input.mode != 'compare'",
        checkboxInput("chip", "Add ChIP-Atlas evidence (human)", FALSE),
        conditionalPanel("input.chip",
          textInput("chip_cell", "Cell types (regex, optional)", "", placeholder = "HepG2|liver"),
          fluidRow(
            column(6, numericInput("chip_max_mb", "Max file (MB)", 300, min = 10, step = 50)),
            column(6, conditionalPanel("input.mode == 'discover'",
              selectInput("chip_distance", "TSS +/- kb", c(1, 5, 10))))))),
      actionButton("run", "Run", class = "btn-primary", width = "100%"),
      tags$hr(),
      tags$b("Display (applies instantly, no re-run)"),
      sliderInput("base_size", "Plot text / point size", min = 10, max = 26, value = 16, step = 1),
      conditionalPanel("input.mode == 'confirm' && input.chip",
        selectInput("ring_p", "Black ring: in ChIP peak and p <=",
                    c("1e-4", "1e-3", "1e-2", "1"), selected = "1e-3")),
      actionButton("recolour", "Shuffle colours", width = "100%"),
      numericInput("color_seed", "Colour seed (blank = random per run)", value = NA),
      helpText("Gene + TF takes seconds; discovery and promoter comparison scan",
               "~1,000 JASPAR motifs and take a minute or more. The first run",
               "downloads JASPAR / HOCOMOCO (and ChIP-Atlas files) once.")
    ),
    mainPanel(
      width = 9,
      div(style = "margin-bottom:8px; font-weight:600;", textOutput("status")),
      uiOutput("tabs")
    )
  )
)

server <- function(input, output, session) {
  log_rv <- reactiveVal(character(0))

  ## ---- slow step: run the chosen analysis (only on Run) ---------------------
  run_rv <- eventReactive(input$run, {
    g <- trimws(input$gene); mode <- input$mode
    validate(need(nzchar(g), "Enter a gene symbol."))
    if (mode == "confirm") validate(need(nzchar(trimws(input$tf)), "Enter a TF symbol."))
    log_rv(character(0))
    cell <- if (nzchar(trimws(input$chip_cell))) trimws(input$chip_cell) else NULL
    msg <- switch(mode, confirm = sprintf("Scanning %s / %s ...", g, trimws(input$tf)),
                  discover = sprintf("Ranking candidate regulators of %s ...", g),
                  compare  = sprintf("Comparing the promoters of %s ...", g))
    withProgress(message = msg, value = 0.4, {
      res <- tryCatch(.with_log({
        if (mode == "confirm") {
          s <- scan_binding_sites(g, trimws(input$tf),
                 upstream = input$upstream, downstream = input$downstream,
                 threshold_frac = input$threshold, transcripts = input$transcripts,
                 max_transcripts = input$max_transcripts, scan_unit = input$scan_unit,
                 hocomoco_version = input$hocomoco_version,
                 jaspar_release = input$jaspar_release)
          if (isTRUE(input$chip) && nrow(s$result))
            s <- add_chip_peaks(s, cell_pattern = cell, max_mb = input$chip_max_mb)
          list(mode = mode, obj = s)
        } else if (mode == "discover") {
          d <- discover_regulators(g, upstream = input$upstream,
                 downstream = input$downstream, top_n = input$top_n,
                 jaspar_release = input$jaspar_release, rank_by = input$rank_by)
          if (isTRUE(input$chip))
            d <- add_chip_targets(d, distance = as.numeric(input$chip_distance),
                                  cell_pattern = cell, max_mb = input$chip_max_mb)
          list(mode = mode, obj = d)
        } else {
          list(mode = mode, obj = compare_promoter_regulators(g,
                 upstream = input$upstream, downstream = input$downstream,
                 transcripts = input$transcripts, jaspar_release = input$jaspar_release))
        }
      }, log_rv), error = function(e) list(mode = mode, error = conditionMessage(e)))
    })
    res
  }, ignoreNULL = TRUE)

  ok <- reactive({
    r <- run_rv()
    if (!is.null(r$error)) return(FALSE)
    switch(r$mode, confirm = nrow(r$obj$result) > 0,
           discover = !is.null(r$obj$top) && nrow(r$obj$top) > 0,
           compare = !is.null(r$obj$comparison))
  })

  ## ---- colour seed: fixed per run so size changes don't reshuffle -----------
  seed_rv <- reactiveVal(NULL)
  observeEvent(input$run, {
    seed_rv(if (is.na(input$color_seed)) sample.int(1e6, 1) else input$color_seed)
  })
  observeEvent(input$recolour, { seed_rv(sample.int(1e6, 1)) })
  observeEvent(input$color_seed, {
    if (!is.na(input$color_seed)) seed_rv(input$color_seed)
  })

  ## ---- fast step: (re)build plots on display changes, no re-run -------------
  plots_rv <- reactive({
    req(ok()); r <- run_rv(); bs <- input$base_size
    switch(r$mode,
      confirm = c(plot_binding_sites(r$obj, base_size = bs, point_size = bs / 4,
                                     color_seed = seed_rv()),
                  list(promoter_map = plot_promoter_map(r$obj, base_size = bs,
                         point_size = bs / 4, color_seed = seed_rv(),
                         chip_ring_p = as.numeric(input$ring_p)))),
      discover = c(plot_regulators(r$obj, base_size = bs, point_size = bs / 4),
                   list(chip = plot_regulators_chip(r$obj, base_size = bs,
                                                    point_size = bs / 4))),
      compare = list(cmp = plot_promoter_regulators(r$obj, top_n = input$top_n,
                                                    base_size = bs)))
  })

  output$status <- renderText({
    r <- run_rv()
    if (!is.null(r$error)) return(paste("Error:", r$error))
    m <- r$obj$meta
    switch(r$mode,
      confirm = if (!nrow(r$obj$result))
        sprintf("No %s sites found for %s at threshold %.2f - try lowering the threshold or widening the window.",
                m$tf, m$gene, input$threshold)
      else sprintf("%s / %s: %d site row(s), %d promoter(s), %d TSS window(s) on chr%s (strand %s) | JASPAR %s, HOCOMOCO %s%s.",
                   m$gene, m$tf, m$n_sites, m$n_promoters, m$n_facets, m$chr,
                   ifelse(m$strand == 1, "+", "-"), m$jaspar_release %||% "?",
                   m$hocomoco_version %||% "?",
                   if (!is.null(r$obj$chip)) sprintf(" | ChIP-Atlas: %d experiment(s) in region",
                                                    r$obj$chip$n_experiments) else ""),
      discover = sprintf("%s: %d motifs tested, top %d TFs shown, ranked by %s.",
                         m$gene, m$n_tested, nrow(r$obj$top),
                         if (identical(m$rank_by, "background")) "promoter background" else "FIMO p"),
      compare = if (is.null(r$obj$comparison))
        sprintf("%s has a single promoter at this window - nothing to compare (try Transcripts = all).", m$gene)
      else sprintf("%s: %d promoters compared; %d TFs with a call.", m$gene, m$n_regions,
                   sum(r$obj$comparison$class != "none")))
  })

  ## ---- tabs depend on the mode of the last run --------------------------------
  output$tabs <- renderUI({
    r <- run_rv(); mode <- r$mode
    log_tab <- tabPanel("Log", verbatimTextOutput("log"))
    if (!isTRUE(ok())) return(tabsetPanel(log_tab))
    nf <- if (mode == "confirm") r$obj$meta$n_facets else 1
    h  <- paste0(max(600, 240 * nf), "px")
    switch(mode,
      confirm = tabsetPanel(
        tabPanel("Promoter map", plotOutput("p_promoter_map", height = "680px")),
        tabPanel("Canonical",    plotOutput("p_canonical", height = "600px")),
        tabPanel("By transcript", plotOutput("p_by_transcript", height = h)),
        tabPanel("Genomic",      plotOutput("p_genomic", height = "640px")),
        tabPanel("Canonical (sig)",     plotOutput("p_canonical_sig", height = "600px")),
        tabPanel("By transcript (sig)", plotOutput("p_by_transcript_sig", height = h)),
        tabPanel("Genomic (sig)",       plotOutput("p_genomic_sig", height = "640px")),
        tabPanel("Promoters", br(), tableOutput("tbl_prom")),
        tabPanel("Table", br(), downloadButton("dl_csv", "Download CSV"), br(), br(),
                 tableOutput("tbl")),
        log_tab),
      discover = tabsetPanel(
        tabPanel("Ranking", plotOutput("p_ranking", height = "700px")),
        tabPanel("Needle",  plotOutput("p_needle", height = "600px")),
        tabPanel("Map",     plotOutput("p_map", height = "700px")),
        tabPanel("Motif vs ChIP", plotOutput("p_chip", height = "650px")),
        tabPanel("Table", br(), downloadButton("dl_csv", "Download CSV"), br(), br(),
                 tableOutput("tbl")),
        log_tab),
      compare = tabsetPanel(
        tabPanel("Promoter-specific TFs", plotOutput("p_cmp", height = "750px")),
        tabPanel("Promoters", br(), tableOutput("tbl_prom")),
        tabPanel("Table", br(), downloadButton("dl_csv", "Download CSV"), br(), br(),
                 tableOutput("tbl")),
        log_tab))
  })

  .sig_need <- "Significance plots need the 'TFMPvalue' package: install.packages('TFMPvalue')."
  output$p_promoter_map  <- renderPlot(.need_plot(plots_rv()$promoter_map))
  output$p_canonical     <- renderPlot(.need_plot(plots_rv()$canonical))
  output$p_by_transcript <- renderPlot(.need_plot(plots_rv()$by_transcript))
  output$p_genomic       <- renderPlot(.need_plot(plots_rv()$genomic))
  output$p_canonical_sig     <- renderPlot(.need_plot(plots_rv()$canonical_sig, .sig_need))
  output$p_by_transcript_sig <- renderPlot(.need_plot(plots_rv()$by_transcript_sig, .sig_need))
  output$p_genomic_sig       <- renderPlot(.need_plot(plots_rv()$genomic_sig, .sig_need))
  output$p_ranking <- renderPlot(.need_plot(plots_rv()$ranking))
  output$p_needle  <- renderPlot(.need_plot(plots_rv()$needle))
  output$p_map     <- renderPlot(.need_plot(plots_rv()$map))
  output$p_chip    <- renderPlot(.need_plot(plots_rv()$chip,
                        "Tick 'Add ChIP-Atlas evidence' and run again to see this plot."))
  output$p_cmp     <- renderPlot(.need_plot(plots_rv()$cmp,
                        "No TF passes the present/absent cut-offs for these promoters."))
  output$log <- renderText(paste(log_rv(), collapse = "\n"))

  table_rv <- reactive({
    req(ok()); r <- run_rv()
    switch(r$mode, confirm = r$obj$result, discover = r$obj$result,
           compare = r$obj$comparison)
  })
  output$tbl <- renderTable({
    t <- table_rv(); r <- run_rv()
    keep <- switch(r$mode,
      confirm = c("promoter_id", "transcript_id", "is_canonical", "database", "matrix_id",
                  "position", "strand", "score_frac", "pvalue", "qvalue",
                  "genomic_position", "site", "chip_n_exp", "chip_cell_types"),
      discover = c("tf", "family", "matrix_id", "position", "strand", "score_frac",
                   "pvalue", "bg_pvalue", "bg_qvalue", "n_at_site", "shares_site_with",
                   "chip_n_exp", "chip_frac_exp"),
      compare = names(t))
    utils::head(t[, intersect(keep, names(t)), drop = FALSE], 500)
  }, digits = 4)
  output$tbl_prom <- renderTable({
    req(ok()); r <- run_rv()
    if (r$mode == "confirm") summarise_promoters(r$obj) else
      r$obj$promoters[, c("promoter_id", "tss", "n_transcripts", "n_tss", "span_bp",
                          "has_canonical", "dist_to_canonical_bp")]
  })
  output$dl_csv <- downloadHandler(
    filename = function() {
      r <- run_rv(); m <- r$obj$meta
      sprintf("%s_%s.csv", m$gene, switch(r$mode, confirm = paste0(m$tf, "_binding_sites"),
                                          discover = "candidate_regulators",
                                          compare = "promoter_regulators"))
    },
    content = function(file) utils::write.csv(table_rv(), file, row.names = FALSE)
  )
}

shinyApp(ui, server)
