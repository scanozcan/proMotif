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
      selectInput("transcripts", "Transcripts", c("protein_coding", "all")),
      conditionalPanel("input.mode != 'compare'",
        sliderInput("threshold", "PWM threshold (fraction of max)",
                    min = 0.60, max = 0.95, value = 0.80, step = 0.01)),
      conditionalPanel("input.mode == 'confirm'",
        fluidRow(
          column(6, numericInput("max_transcripts", "Max TSS windows", 8, min = 1, step = 1)),
          column(6, selectInput("scan_unit", "Scan per", c("tss", "promoter"))))),
      conditionalPanel("input.mode != 'confirm'",
        numericInput("top_n", "TFs to show", 20, min = 5, step = 5)),
      conditionalPanel("input.mode == 'discover'",
        selectInput("rank_by", "Rank TFs by",
                    c("promoter background" = "background", "FIMO per-site p" = "fimo")),
        numericInput("expand_n", "Targeted scans of top TFs (0 = off)", 3, min = 0, step = 1),
        selectInput("expand_by", "Pick them by",
                    c("ChIP support first (if ChIP on)" = "auto", "motif rank" = "rank")),
        textInput("expand_tfs", "...or these TFs (comma-separated)", "")),
      tags$details(tags$summary("Motif databases"),
        fluidRow(
          column(6, selectInput("jaspar_release", "JASPAR", c("2026", "2024", "2022"))),
          column(6, selectInput("hocomoco_version", "HOCOMOCO", c("v14", "v12", "v11"))))),
      conditionalPanel("input.mode == 'confirm'",
        checkboxInput("usage", "Promoter usage across cell types (FANTOM5 CAGE)", FALSE)),
      conditionalPanel("input.mode != 'compare'",
        checkboxInput("chip", "Add ChIP-Atlas evidence (human)", FALSE),
        conditionalPanel("input.chip",
          textInput("chip_cell", "Cell types (regex, optional)", "", placeholder = "HepG2|liver"),
          fluidRow(
            column(6, numericInput("chip_max_mb", "ChIP data budget (MB)", 300, min = 10, step = 50)),
            column(6, conditionalPanel("input.mode == 'discover'",
              selectInput("chip_distance", "TSS +/- kb", c(1, 5, 10))))),
          helpText("If a TF's merged ChIP-Atlas file is larger than the budget, the newest",
                   "experiments are used until the budget is reached (files are cached)."))),
      actionButton("run", "Run", class = "btn-primary", width = "100%"),
      tags$hr(),
      tags$b("Display (applies instantly, no re-run)"),
      sliderInput("base_size", "Plot text / point size", min = 10, max = 26, value = 16, step = 1),
      conditionalPanel("input.mode != 'compare' && input.chip",
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
          if (isTRUE(input$usage)) s <- add_promoter_usage(s)
          list(mode = mode, obj = s)
        } else if (mode == "discover") {
          d <- discover_regulators(g, upstream = input$upstream,
                 downstream = input$downstream, top_n = input$top_n,
                 transcripts = input$transcripts,
                 jaspar_release = input$jaspar_release, rank_by = input$rank_by)
          if (isTRUE(input$chip))
            d <- add_chip_targets(d, distance = as.numeric(input$chip_distance),
                                  cell_pattern = cell, max_mb = input$chip_max_mb)
          etf <- trimws(strsplit(input$expand_tfs, ",")[[1]]); etf <- etf[nzchar(etf)]
          ex <- if (isTRUE(input$expand_n > 0) || length(etf))
            expand_regulators(d, tfs = if (length(etf)) etf else NULL, n = input$expand_n,
                              by = input$expand_by, chip = isTRUE(input$chip),
                              threshold_frac = input$threshold,
                              hocomoco_version = input$hocomoco_version,
                              chip_cell_pattern = cell, chip_max_mb = input$chip_max_mb)
            else list()
          list(mode = mode, obj = d, expanded = ex)
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
      discover = sprintf("%s: %d promoter(s) scanned (%s), %d motifs, top %d TFs per promoter, ranked by %s.",
                         m$gene, m$n_promoters %||% 1L,
                         paste(m$promoters_scanned %||% "P1", collapse = ", "),
                         m$n_tested, m$top_n,
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
        tabPanel("Top sites",
                 helpText("Hover a point to see the site, its exact genomic coordinates and sequence."),
                 div(style = "position:relative",
                     plotOutput("p_top", height = "620px",
                                hover = hoverOpts("top_hover", delay = 80, delayType = "debounce")),
                     uiOutput("top_hover_info")),
                 fluidRow(column(4, numericInput("label_top", "Label the top N sites", 10, min = 0, step = 1)),
                          column(4, br(), downloadButton("dl_top", "Download all ranked sites"))),
                 tableOutput("tbl_top")),
        tabPanel("Motif", plotOutput("p_logo", height = "560px")),
        tabPanel("Promoter map", plotOutput("p_promoter_map", height = "680px")),
        tabPanel("Promoter usage", plotOutput("p_usage", height = "800px"),
                 br(), tableOutput("tbl_usage")),
        tabPanel("Normal vs cancer", plotOutput("p_usage_cancer", height = "750px")),
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
        tabPanel("Top TF maps", br(), uiOutput("exp_sel"),
                 plotOutput("p_exp", height = "720px")),
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

  ## ---- top sites: probability map without the ChIP track (one ggplot, so
  ## hover coordinates map to the data); tooltip = site, coordinates, sequence
  top_plot_rv <- reactive({
    req(ok()); r <- run_rv(); req(r$mode == "confirm")
    bs <- input$base_size
    plot_top_sites(r$obj, top_n = input$label_top %||% 10, base_size = bs, point_size = bs / 4,
                   color_seed = seed_rv(), chip_ring_p = as.numeric(input$ring_p),
                   chip_track = FALSE)
  })
  output$p_top <- renderPlot(.need_plot(top_plot_rv()))
  output$p_logo <- renderPlot({
    r <- run_rv(); req(ok(), r$mode == "confirm")
    .need_plot(plot_motif_logo(r$obj, base_size = input$base_size))
  })
  output$top_hover_info <- renderUI({
    h <- input$top_hover; p <- top_plot_rv(); req(h, p)
    pts <- attr(p, "points"); req(pts)
    pv <- if (!is.null(h$mapping$panelvar1)) "panel" else NULL
    near <- nearPoints(pts, h, xvar = "genomic_position", yvar = "y", panelvar1 = pv,
                       threshold = 8, maxpoints = 1)
    if (!nrow(near)) return(NULL)
    s <- near[1, ]
    num <- function(x, f = "%.3g") if (is.null(x) || is.na(x)) "NA" else sprintf(f, x)
    sid <- ""
    st <- run_rv()$obj$sites
    if (!is.null(st)) {
      k <- which(st$genomic_start == s$genomic_start & st$genomic_end == s$genomic_end &
                 st$motif_strand_genomic == s$motif_strand_genomic)
      if (length(k)) sid <- paste0(st$site_id[k[1]], " &middot; ")
    }
    html <- paste0(
      "<b>", sid, s$database, " ", s$matrix_id, "</b><br>",
      sprintf("chr%s:%s-%s (%s strand)", sub("^chr", "", run_rv()$obj$meta$chr),
              format(s$genomic_start, big.mark = ","), format(s$genomic_end, big.mark = ","),
              s$motif_strand_genomic), "<br>",
      "motif: <span style='font-family:monospace'>", s$site_motif, "</span><br>",
      "genome (+): <span style='font-family:monospace'>", s$context_genomic, "</span><br>",
      sprintf("%s bp from the %s TSS (%s) | p = %s, q = %s", s$position, s$promoter_id,
              s$transcript_id, num(s$pvalue), num(s$qvalue)), "<br>",
      if (!is.null(s$chip_n_exp)) sprintf("ChIP: %s experiment(s), support %s, summit %s bp away<br>",
                                          s$chip_n_exp, num(s$chip_support, "%.2f"),
                                          num(s$chip_summit_dist, "%.0f")) else "",
      sprintf("<b>binding probability %s</b> (%s)", num(s$binding_prob, "%.2f"), s$confidence))
    left <- h$coords_css$x + 14; top <- h$coords_css$y + 14
    div(style = sprintf(paste0("position:absolute; left:%dpx; top:%dpx; z-index:100;",
                               "background:rgba(255,255,255,0.96); border:1px solid #bbb;",
                               "border-radius:4px; padding:6px 9px; font-size:12px;",
                               "pointer-events:none; max-width:460px;"), round(left), round(top)),
        HTML(html))
  })
  output$tbl_top <- renderTable({
    r <- run_rv(); req(ok(), r$mode == "confirm")
    t <- top_binding_sites(r$obj, min_n = max(10, input$label_top %||% 10))
    req(t)
    t[, intersect(c("site_id", "chrom", "genomic_start", "genomic_end", "motif_strand_genomic",
                    "promoter_id", "position_rel_tss", "databases", "site_motif", "context_genomic",
                    "pvalue", "qvalue", "chip_n_exp", "chip_support", "binding_prob", "confidence"),
                  names(t)), drop = FALSE]
  }, digits = 3)
  output$dl_top <- downloadHandler(
    filename = function() { m <- run_rv()$obj$meta; sprintf("%s_%s_top_sites.csv", m$gene, m$tf) },
    content = function(file) utils::write.csv(run_rv()$obj$sites, file, row.names = FALSE))

  .sig_need <- "Significance plots need the 'TFMPvalue' package: install.packages('TFMPvalue')."
  output$p_promoter_map  <- renderPlot(.need_plot(plots_rv()$promoter_map))
  output$p_usage <- renderPlot({
    r <- run_rv(); req(ok(), r$mode == "confirm")
    validate(need(!is.null(r$obj$usage),
      "Tick 'Promoter usage across cell types' and run again (needs the FANTOM5 folder, see ?add_promoter_usage)."))
    plot_promoter_usage(r$obj, base_size = input$base_size)
  })
  output$p_usage_cancer <- renderPlot({
    r <- run_rv(); req(ok(), r$mode == "confirm")
    validate(need(!is.null(r$obj$usage),
      "Tick 'Promoter usage across cell types' and run again."))
    p <- plot_promoter_usage_cancer(r$obj, base_size = input$base_size)
    validate(need(!is.null(p),
      "Not enough normal and cancer samples of the same lineage express this gene."))
    p
  })
  output$tbl_usage <- renderTable({
    r <- run_rv(); req(ok(), r$mode == "confirm", !is.null(r$obj$usage))
    utils::head(r$obj$usage$table, 100)
  }, digits = 2)
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

  ## targeted scans of the top discovery TFs (promoter map, ChIP track, rings)
  output$exp_sel <- renderUI({
    ex <- run_rv()$expanded
    if (!length(ex)) return(helpText("Set 'Targeted scans of top TFs' above 0 and run again."))
    selectInput("exp_tf", "TF", names(ex), width = "300px")
  })
  output$p_exp <- renderPlot({
    ex <- run_rv()$expanded; req(length(ex), input$exp_tf)
    s <- ex[[input$exp_tf]]
    validate(need(!is.null(s) && nrow(s$result) > 0,
                  "No sites for this TF at the chosen threshold."))
    bs <- input$base_size
    plot_promoter_map(s, base_size = bs, point_size = bs / 4, color_seed = seed_rv(),
                      chip_ring_p = as.numeric(input$ring_p))
  })

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
                  "genomic_start", "genomic_end", "site_motif", "chip_n_exp", "chip_support",
                  "binding_prob", "confidence", "chip_cell_types"),
      discover = c("promoter_id", "tf", "family", "matrix_id", "position", "strand", "score_frac",
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
