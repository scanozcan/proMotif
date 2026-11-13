## proMotif Shiny app (launched via proMotif::launch_app()).
## Modes: gene + TF, discover regulators, compare promoters. Promoters from
## EPDnew (default) or Ensembl. Discovery selects the likely regulators of each
## promoter and scores them with the proMotif score (run_regulator_discovery).
## The analysis runs once per "Run"; text size, colours, the ChIP rings, the
## promoter-map, logo, usage and gene-track options re-style instantly.
library(shiny)
library(proMotif)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
## FANTOM5 lineages of the normal vs cancer comparison (gene track)
LINEAGES <- c("blood", "brain", "liver", "colon", "stomach", "esophagus", "lung", "breast", "kidney",
              "prostate", "pancreas", "skin", "female_reproductive", "testis", "bone_soft_tissue",
              "thyroid", "bladder", "head_neck", "eye", "endothelium")
APP_TITLE <- "proMotif - promoter binding-site scanner (JASPAR + HOCOMOCO)"

## ---- helpers ----------------------------------------------------------------
## run an expression, collecting its console messages for the "Log" tab
.with_log <- function(expr, log) {
  withCallingHandlers(expr, message = function(m) {
    log(c(isolate(log()), sub("\n$", "", conditionMessage(m))))   # isolate: no reactive loop
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
      selectInput("tss_source", "Promoters (TSSs) from",
                  c("EPDnew (curated promoters)" = "epd", "Ensembl transcripts" = "ensembl")),
      selectInput("transcripts", "Transcripts", c("protein_coding", "all")),
      conditionalPanel("input.mode == 'confirm'",
        sliderInput("threshold", "PWM threshold (fraction of max)",
                    min = 0.60, max = 0.95, value = 0.80, step = 0.01)),
      conditionalPanel("input.mode == 'confirm'",
        fluidRow(
          column(6, numericInput("max_transcripts", "Max TSS windows", 8, min = 1, step = 1)),
          column(6, selectInput("scan_unit", "Scan per", c("tss", "promoter"))))),
      conditionalPanel("input.mode != 'confirm'",
        numericInput("top_n", "Motifs to show per promoter", 15, min = 5, step = 5)),
      conditionalPanel("input.mode == 'discover'",
        selectInput("rule", "Select TFs",
                    c("Likely regulators (motif + ChIP-Atlas)" = "regulators",
                      "Promoter-specific motifs" = "specific")),
        conditionalPanel("input.rule == 'regulators'",
          checkboxInput("auto_n", "Number of TFs: automatic (decided by the scores)", TRUE),
          conditionalPanel("!input.auto_n",
            fluidRow(
              column(6, numericInput("n_chip", "ChIP-bound TFs / promoter", 5, min = 0, step = 1)),
              column(6, numericInput("n_motif", "Motif-only TFs / promoter", 2, min = 0, step = 1))))),
        conditionalPanel("input.rule == 'specific'",
          numericInput("n_per_promoter", "TFs per promoter", 3, min = 1, step = 1)),
        textInput("expand_tfs", "Also score these TFs (comma-separated)", ""),
        numericInput("expand_threshold", "PWM threshold of the TF scans", 0.75,
                     min = 0.60, max = 0.95, step = 0.01),
        helpText("Likely regulators: TFs whose motif is at least as good here as at most",
                 "promoters; with ChIP-Atlas on, the TFs with the largest share of their",
                 "experiments bound near the gene are scored too. Each scored TF gets the",
                 "proMotif score; automatic shows the TFs with score >= 0.8 plus 2 more",
                 "(5-10 per promoter).")),
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
          selectInput("chip_pulldown", "Pulldowns",
                      c("TF antibody + epitope tag" = "all", "TF antibody only" = "native",
                        "Epitope tag only (FLAG, HA, GFP ...)" = "tagged"), "all"),
          numericInput("chip_max_mb", "ChIP data budget per TF (MB)", 300, min = 10, step = 50),
          helpText("If a TF's ChIP-Atlas data are larger than the budget, the newest",
                   "experiments are used until the budget is reached (files are cached).",
                   "In discovery, the first genes download the screen files of many TFs."))),
      checkboxInput("track", "Gene track (promoters and transcripts)", TRUE),
      conditionalPanel("input.track",
        checkboxInput("track_cage", "...with FANTOM5 CAGE, normal vs cancer (needs the FANTOM5 folder)", FALSE)),
      actionButton("run", "Run", class = "btn-primary", width = "100%"),
      tags$hr(),
      tags$b("Display (applies instantly, no re-run)"),
      sliderInput("base_size", "Plot text / point size", min = 10, max = 26, value = 16, step = 1),
      conditionalPanel("input.mode != 'compare' && input.chip",
        fluidRow(
          column(7, selectInput("ring_p", "Black ring: in ChIP peak and p <=",
                                c("1e-4", "1e-3", "1e-2", "1"), selected = "1e-3")),
          column(5, numericInput("ring_min_exp", "...and >= experiments", 1, min = 1, step = 1))),
        checkboxInput("chip_track", "ChIP-seq coverage track under the maps", TRUE)),
      conditionalPanel("input.mode == 'confirm'",
        tags$details(tags$summary("Promoter map, logo and usage"),
          selectInput("map_metric", "Promoter map: y axis",
                      c("motif -log10 p" = "p", "proMotif score" = "score")),
          numericInput("map_label", "Promoter map: label the top N sites (score axis)", 0, min = 0, step = 1),
          checkboxInput("map_facet", "Promoter map: one panel per distant promoter", TRUE),
          checkboxInput("logo_revcomp", "Motif logo: also reverse complement", TRUE),
          fluidRow(
            column(6, numericInput("usage_top_n", "Usage: sample types per promoter", 10, min = 3, step = 1)),
            column(6, numericInput("usage_min_n", "Normal vs cancer: min samples per group", 3, min = 2, step = 1))))),
      conditionalPanel("input.track",
        tags$details(tags$summary("Gene track"),
          fluidRow(
            column(6, selectInput("track_tx", "Transcripts", c("all", "protein_coding"))),
            column(6, numericInput("track_max_tx", "Max transcripts", 30, min = 1, step = 1))),
          conditionalPanel("input.track_cage",
            selectInput("track_lineage", "CAGE: normal vs cancer in lineage",
                        c("all lineages" = "", LINEAGES)),
            uiOutput("track_group_ui")))),
      actionButton("recolour", "Shuffle colours", width = "100%"),
      numericInput("color_seed", "Colour seed (blank = random per run)", value = NA),
      helpText("Gene + TF takes seconds; discovery and promoter comparison scan",
               "~1,000 JASPAR motifs and take a few minutes (more with ChIP-Atlas).",
               "The first run downloads JASPAR / HOCOMOCO (and ChIP-Atlas files) once.")
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
                  discover = sprintf("Finding the likely regulators of %s ...", g),
                  compare  = sprintf("Comparing the promoters of %s ...", g))
    withProgress(message = msg, value = 0.4, {
      res <- tryCatch(.with_log({
        if (mode == "confirm") {
          s <- scan_binding_sites(g, trimws(input$tf),
                 upstream = input$upstream, downstream = input$downstream,
                 threshold_frac = input$threshold, transcripts = input$transcripts,
                 max_transcripts = input$max_transcripts, scan_unit = input$scan_unit,
                 hocomoco_version = input$hocomoco_version,
                 jaspar_release = input$jaspar_release, tss_source = input$tss_source)
          if (isTRUE(input$chip) && nrow(s$result))
            s <- add_chip_peaks(s, cell_pattern = cell, max_mb = input$chip_max_mb,
                                pulldown = input$chip_pulldown %||% "all")
          if (isTRUE(input$usage)) s <- add_promoter_usage(s)
          list(mode = mode, obj = s)
        } else if (mode == "discover") {
          etf <- trimws(strsplit(input$expand_tfs, ",")[[1]]); etf <- etf[nzchar(etf)]
          d <- run_regulator_discovery(g, upstream = input$upstream,
                 downstream = input$downstream, top_n = input$top_n, write_files = FALSE,
                 chip = isTRUE(input$chip), chip_cell_pattern = cell,
                 chip_max_mb = input$chip_max_mb, chip_pulldown = input$chip_pulldown %||% "all",
                 rule = input$rule,
                 n_chip = if (isTRUE(input$auto_n)) "auto" else input$n_chip,
                 n_motif = if (isTRUE(input$auto_n)) "auto" else input$n_motif,
                 n_per_promoter = input$n_per_promoter,
                 expand_threshold = input$expand_threshold,
                 expand_tfs = if (length(etf)) etf else NULL,
                 transcripts = input$transcripts, jaspar_release = input$jaspar_release,
                 tss_source = input$tss_source)
          list(mode = mode, obj = d, expanded = d$expanded, chip = isTRUE(input$chip))
        } else {
          list(mode = mode, obj = compare_promoter_regulators(g,
                 upstream = input$upstream, downstream = input$downstream,
                 transcripts = input$transcripts, jaspar_release = input$jaspar_release,
                 tss_source = input$tss_source))
        }
      }, log_rv), error = function(e) list(mode = mode, error = conditionMessage(e)))
      res$gene <- g; res$tss_source <- input$tss_source
    })
    res
  }, ignoreNULL = TRUE)

  ok <- reactive({
    r <- run_rv()
    if (!is.null(r$error)) return(FALSE)
    switch(r$mode, confirm = nrow(r$obj$result) > 0,
           discover = !is.null(r$obj$candidates) && nrow(r$obj$candidates) > 0,
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
                         chip_ring_p = as.numeric(input$ring_p),
                         chip_ring_min_exp = input$ring_min_exp %||% 1,
                         facet_regions = isTRUE(input$map_facet %||% TRUE),
                         y_metric = input$map_metric %||% "p", label_top = input$map_label %||% 0,
                         chip_track = isTRUE(input$chip_track %||% TRUE)))),
      discover = list(
        regulators = plot_discovery_regulators(r$obj$candidates, r$obj, base_size = bs,
                                               point_size = bs / 4),
        candidates = plot_discovery_candidates(r$obj$candidates, r$obj, top_n = input$top_n,
                                               base_size = bs, point_size = bs / 4),
        specific = plot_discovery_promoter_specific(r$obj$specific, r$obj, base_size = bs)),
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
      discover = {
        cd <- r$obj$candidates; sel <- cd[cd$selected, , drop = FALSE]
        sprintf("%s: %d promoter(s) (%s, %s), %s motifs scanned; %d TF(s) scored with the proMotif score, %d shown (%s)%s.",
                m$gene, length(unique(cd$promoter_id)), paste(sort(unique(cd$promoter_id)), collapse = ", "),
                if (identical(m$tss_source, "epd")) "EPDnew" else "Ensembl",
                format(m$n_tested %||% NA), length(unique(sel$tf)),
                length(unique(cd$tf[cd$reported %in% TRUE])),
                if (identical(attr(cd, "rule")$mode, "regulators")) "likely-regulator rule" else "promoter-specific rule",
                if (isTRUE(r$chip)) " | ChIP-Atlas screen and peaks" else " | motif only (tick ChIP-Atlas for the screen)")
      },
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
    track_tab <- tabPanel("Gene track", plotOutput("p_track", height = "auto"))
    if (mode == "discover") {
      cd <- r$obj$candidates
      nreg <- max(c(1, table(cd$promoter_id[if (!is.null(cd$reported)) cd$reported else cd$selected])))
      h_reg  <- paste0(max(450, 42 * nreg + 220), "px")
      h_cand <- paste0(max(500, 30 * (input$top_n + nreg) + 220), "px")
      h_spec <- paste0(max(420, 38 * length(unique(r$obj$specific$tf)) + 220), "px")
    }
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
        track_tab,
        tabPanel("Table", br(), downloadButton("dl_csv", "Download CSV"), br(), br(),
                 tableOutput("tbl")),
        log_tab),
      discover = tabsetPanel(
        tabPanel("Regulators",
                 helpText("Scored TFs of each promoter. Open circle = motif-only score, filled",
                          "= proMotif score (motif + ChIP-Atlas); labels: experiments with a peak",
                          "over the site, experiments bound near the gene, motif background p."),
                 plotOutput("p_regulators", height = h_reg)),
        tabPanel("Promoter-specific",
                 helpText("Best proMotif score of each scored TF at every promoter (targeted scans).",
                          "Framed = binding specific to that promoter (score >= 0.5 and >= 0.3 above",
                          "every other promoter). Genes with one promoter have nothing to compare."),
                 plotOutput("p_specific", height = h_spec)),
        tabPanel("Candidates",
                 helpText("Motif shortlist of each promoter (best motifs vs a 20,094-promoter",
                          "background), coloured by what happened to each TF."),
                 plotOutput("p_candidates", height = h_cand)),
        tabPanel("Top sites per TF", br(), uiOutput("exp_sel"),
                 fluidRow(column(4, radioButtons("exp_view", NULL, inline = TRUE,
                                                 c("Top sites (proMotif score)" = "top",
                                                   "Promoter map (motif -log10 p)" = "map"))),
                          column(4, numericInput("label_top_d", "Label the top N sites", 10, min = 0, step = 1))),
                 plotOutput("p_exp", height = "720px")),
        track_tab,
        tabPanel("Table", br(), downloadButton("dl_csv", "Download candidates (CSV)"), " ",
                 downloadButton("dl_motifs", "Download all motifs (CSV)"), br(), br(),
                 tableOutput("tbl")),
        log_tab),
      compare = tabsetPanel(
        tabPanel("Promoter-specific TFs", plotOutput("p_cmp", height = "750px")),
        tabPanel("Promoters", br(), tableOutput("tbl_prom")),
        track_tab,
        tabPanel("Table", br(), downloadButton("dl_csv", "Download CSV"), br(), br(),
                 tableOutput("tbl")),
        log_tab))
  })

  ## ---- top sites: score map without the ChIP track (one ggplot, so
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
    .need_plot(plot_motif_logo(r$obj, base_size = input$base_size,
                               revcomp = isTRUE(input$logo_revcomp %||% TRUE)))
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
      if (!is.null(s$chip_n_exp)) sprintf("ChIP: %s experiment(s)%s, support %s, summit %s bp away<br>",
                                          s$chip_n_exp,
                                          if (isTRUE(s$chip_n_exp_tagged > 0)) sprintf(" (%d tag pulldown)", s$chip_n_exp_tagged) else "",
                                          num(s$chip_support, "%.2f"),
                                          num(s$chip_summit_dist, "%.0f")) else "",
      sprintf("<b>proMotif score %s</b> (%s)", num(s$promotif_score, "%.2f"), s$confidence))
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
                    "pvalue", "qvalue", "chip_n_exp", "chip_n_exp_tagged", "chip_support", "promotif_score", "confidence"),
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
    plot_promoter_usage(r$obj, top_n = input$usage_top_n %||% 10, base_size = input$base_size)
  })
  output$p_usage_cancer <- renderPlot({
    r <- run_rv(); req(ok(), r$mode == "confirm")
    validate(need(!is.null(r$obj$usage),
      "Tick 'Promoter usage across cell types' and run again."))
    p <- plot_promoter_usage_cancer(r$obj, min_n = input$usage_min_n %||% 3, base_size = input$base_size)
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
  output$p_regulators <- renderPlot(.need_plot(plots_rv()$regulators,
                           "No TF was scored for this gene (no motif passed the rule)."))
  output$p_candidates <- renderPlot(.need_plot(plots_rv()$candidates))
  output$p_specific <- renderPlot(.need_plot(plots_rv()$specific,
                         "Only one promoter (or no scored TF): nothing to compare between promoters."))
  ## ---- gene track: promoters, transcripts and optional CAGE, rebuilt when its
  ## options or the text size change (Ensembl and FANTOM5 data are cached)
  output$track_group_ui <- renderUI({
    r <- run_rv(); u <- r$obj$usage
    if (is.null(u)) return(helpText("Tick 'Promoter usage' (Gene + TF) to compare the cells",
                                     "that use one promoter with the other primary cells."))
    pids <- rownames(u$tpm)
    selectInput("track_groups", "CAGE: compare",
                c("normal vs cancer" = "",
                  stats::setNames(pids, sprintf("primary cells using %s vs other primary cells", pids))))
  })
  track_rv <- reactive({
    r <- run_rv(); req(is.null(r$error), isTRUE(input$track), r$gene)
    grp <- NULL
    tg <- input$track_groups %||% ""
    if (isTRUE(input$track_cage) && nzchar(tg) && !is.null(r$obj$usage)) {
      tab <- r$obj$usage$table; prim <- tab[tab$category == "primary cells", , drop = FALSE]
      p1 <- prim$sample_type[prim$dominant == tg]
      if (length(p1))
        grp <- stats::setNames(list(p1, setdiff(prim$sample_type, p1)),
                               c(sprintf("primary cells using %s", tg), "other primary cells"))
    }
    lin <- input$track_lineage %||% ""
    tryCatch(.with_log(plot_promoter_track(r$gene, tss_source = r$tss_source,
                                           transcripts = input$track_tx %||% "all",
                                           cage = isTRUE(input$track_cage),
                                           lineage = if (nzchar(lin)) lin else NULL, groups = grp,
                                           max_transcripts = input$track_max_tx %||% 30,
                                           base_size = round(input$base_size * 0.8)), log_rv),
             error = function(e) { log_rv(c(isolate(log_rv()), paste("Gene track:", conditionMessage(e)))); NULL })
  })
  output$p_track <- renderPlot({
    validate(need(isTRUE(input$track), "Tick 'Gene track' to draw it."))
    p <- track_rv()
    validate(need(!is.null(p), "The gene track could not be drawn (see the Log tab)."))
    p
  }, height = function() {
    p <- tryCatch(track_rv(), error = function(e) NULL)
    round(80 * (attr(p, "height") %||% 8))
  })
  output$p_cmp     <- renderPlot(.need_plot(plots_rv()$cmp,
                        "No TF passes the present/absent cut-offs for these promoters."))
  output$log <- renderText(paste(log_rv(), collapse = "\n"))

  ## the scored discovery TFs: top-sites plot (proMotif score) of each,
  ## listed by their best score (with the promoter where it is reached)
  output$exp_sel <- renderUI({
    r <- run_rv(); ex <- r$expanded
    if (!length(ex)) return(helpText("No TF was scored for this gene."))
    cd <- r$obj$candidates; cd <- cd[cd$selected & cd$tf %in% names(ex) & !is.na(cd$score_best), , drop = FALSE]
    cd <- cd[order(-cd$score_best), , drop = FALSE]; cd <- cd[!duplicated(cd$tf), , drop = FALSE]
    tfs <- c(cd$tf, setdiff(names(ex), cd$tf))
    lab <- ifelse(tfs %in% cd$tf,
                  sprintf("%s  (%s, score %.2f)", tfs, cd$promoter_id[match(tfs, cd$tf)],
                          cd$score_best[match(tfs, cd$tf)]), tfs)
    selectInput("exp_tf", "TF", stats::setNames(tfs, lab), width = "360px")
  })
  output$p_exp <- renderPlot({
    ex <- run_rv()$expanded; req(length(ex), input$exp_tf)
    s <- ex[[input$exp_tf]]
    validate(need(!is.null(s) && nrow(s$result) > 0,
                  "No sites for this TF at the chosen threshold."))
    bs <- input$base_size
    if (identical(input$exp_view, "map"))
      plot_promoter_map(s, base_size = bs, point_size = bs / 4, color_seed = seed_rv(),
                        chip_ring_p = as.numeric(input$ring_p),
                        chip_ring_min_exp = input$ring_min_exp %||% 1,
                        label_top = input$label_top_d %||% 10,
                        chip_track = isTRUE(input$chip_track %||% TRUE))
    else
      plot_top_sites(s, top_n = input$label_top_d %||% 10, base_size = bs, point_size = bs / 4,
                     color_seed = seed_rv(), chip_ring_p = as.numeric(input$ring_p),
                     chip_track = isTRUE(input$chip_track %||% TRUE))
  })

  table_rv <- reactive({
    req(ok()); r <- run_rv()
    switch(r$mode, confirm = r$obj$result,
           discover = { cd <- r$obj$candidates
                        cd[order(!cd$selected, cd$promoter_id, -ifelse(is.na(cd$score_best), -1, cd$score_best),
                                 cd$p_rank %||% seq_len(nrow(cd))), , drop = FALSE] },
           compare = r$obj$comparison)
  })
  output$tbl <- renderTable({
    t <- table_rv(); r <- run_rv()
    keep <- switch(r$mode,
      confirm = c("promoter_id", "transcript_id", "is_canonical", "database", "matrix_id",
                  "position", "strand", "score_frac", "pvalue", "qvalue",
                  "genomic_start", "genomic_end", "site_motif", "chip_n_exp", "chip_n_exp_tagged", "chip_support",
                  "promotif_score", "confidence", "chip_cell_types"),
      discover = c("promoter_id", "tf", "status", "reported", "route", "screen_n_bound", "screen_n_exp", "screen_target_q",
                   "screen_source",
                   "bg_pvalue", "pvalue", "position", "site", "pref", "score_best", "score_motif_only",
                   "n_significant", "chip_exp_best", "best_locus", "binding_check"),
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
                                          discover = "discovery_candidates",
                                          compare = "promoter_regulators"))
    },
    content = function(file) {
      t <- table_rv(); t$p_rank <- NULL; t$lp <- NULL; t$eligible <- NULL
      utils::write.csv(t, file, row.names = FALSE)
    }
  )
  output$dl_motifs <- downloadHandler(
    filename = function() sprintf("%s_candidate_regulators.csv", run_rv()$obj$meta$gene),
    content = function(file) utils::write.csv(run_rv()$obj$result, file, row.names = FALSE))
}

shinyApp(ui, server)
