## proMotif Shiny app (launched via proMotif::launch_app()).
## Scan runs once per "Run scan"; text/point size and colours re-style instantly.
library(shiny)
library(proMotif)

ui <- fluidPage(
  titlePanel("proMotif - promoter binding-site scanner (JASPAR + HOCOMOCO)"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      textInput("gene", "Gene (HGNC symbol)", "ATP7B"),
      textInput("tf",   "Transcription factor", "MTF1"),
      fluidRow(
        column(6, numericInput("upstream",   "Upstream (bp)",   2500, min = 0, step = 100)),
        column(6, numericInput("downstream", "Downstream (bp)",  500, min = 0, step = 100))
      ),
      sliderInput("threshold", "PWM threshold (fraction of max)",
                  min = 0.60, max = 0.95, value = 0.80, step = 0.01),
      selectInput("transcripts", "Transcripts",
                  c("protein_coding", "all"), selected = "protein_coding"),
      numericInput("max_transcripts", "Max transcripts (facets)", 8, min = 1, step = 1),
      actionButton("run", "Run scan", class = "btn-primary", width = "100%"),
      tags$hr(),
      tags$b("Display (applies instantly, no re-scan)"),
      sliderInput("base_size", "Plot text / point size", min = 10, max = 26, value = 16, step = 1),
      actionButton("recolour", "Shuffle colours", width = "100%"),
      numericInput("color_seed", "Colour seed (blank = random per run)", value = NA),
      helpText("Scan queries Ensembl / JASPAR / HOCOMOCO live and takes a few",
               "seconds (longer on the first HOCOMOCO download). Size + colour",
               "changes redraw the existing results immediately.")
    ),
    mainPanel(
      width = 9,
      div(style = "margin-bottom:8px; font-weight:600;", textOutput("status")),
      tabsetPanel(
        tabPanel("Canonical",     plotOutput("p_canonical", height = "600px")),
        tabPanel("By transcript", uiOutput("p_by_transcript_ui")),
        tabPanel("Genomic",       plotOutput("p_genomic", height = "640px")),
        tabPanel("Canonical (sig)",     plotOutput("p_canonical_sig", height = "600px")),
        tabPanel("By transcript (sig)", uiOutput("p_by_transcript_sig_ui")),
        tabPanel("Genomic (sig)",       plotOutput("p_genomic_sig", height = "640px")),
        tabPanel("Table",
                 br(), downloadButton("dl_csv", "Download CSV"),
                 br(), br(), tableOutput("tbl"))
      )
    )
  )
)

server <- function(input, output, session) {

  scan_rv <- eventReactive(input$run, {
    g <- trimws(input$gene); tf <- trimws(input$tf)
    validate(need(nzchar(g), "Enter a gene symbol."),
             need(nzchar(tf), "Enter a TF symbol."))
    withProgress(message = sprintf("Scanning %s / %s ...", g, tf), value = 0.4, {
      tryCatch(
        scan_binding_sites(
          gene = g, tf = tf,
          upstream = input$upstream, downstream = input$downstream,
          threshold_frac = input$threshold, transcripts = input$transcripts,
          max_transcripts = input$max_transcripts),
        error = function(e) structure(list(error = conditionMessage(e)),
                                      class = "scan_error"))
    })
  }, ignoreNULL = TRUE)

  scan_ok <- reactive({
    s <- scan_rv(); !inherits(s, "scan_error") && nrow(s$result) > 0
  })

  seed_rv <- reactiveVal(NULL)
  observeEvent(input$run, {
    seed_rv(if (is.na(input$color_seed)) sample.int(1e6, 1) else input$color_seed)
  })
  observeEvent(input$recolour, { seed_rv(sample.int(1e6, 1)) })
  observeEvent(input$color_seed, {
    if (!is.na(input$color_seed)) seed_rv(input$color_seed)
  })

  plots_rv <- reactive({
    req(scan_ok())
    plot_binding_sites(scan_rv(), base_size = input$base_size,
                       point_size = input$base_size / 4, color_seed = seed_rv())
  })

  output$status <- renderText({
    s <- scan_rv()
    if (inherits(s, "scan_error")) return(paste("Error:", s$error))
    if (nrow(s$result) == 0)
      return(sprintf("No %s sites found for %s at threshold %.2f - try lowering the threshold or widening the window.",
                     s$meta$tf, s$meta$gene, input$threshold))
    sprintf("%s / %s: %d site(s) across %d transcript window(s) on chr%s (strand %s).",
            s$meta$gene, s$meta$tf, s$meta$n_sites, s$meta$n_facets, s$meta$chr,
            ifelse(s$meta$strand == 1, "+", "-"))
  })

  output$p_canonical <- renderPlot({ p <- plots_rv()$canonical; req(!is.null(p)); p })
  output$p_genomic   <- renderPlot({ plots_rv()$genomic })
  output$p_by_transcript_ui <- renderUI({
    req(scan_ok())
    plotOutput("p_by_transcript",
               height = paste0(max(600, 240 * scan_rv()$meta$n_facets), "px"))
  })
  output$p_by_transcript <- renderPlot({ plots_rv()$by_transcript })

  .sig_need <- "Significance plots need the 'TFMPvalue' package: install.packages('TFMPvalue')."
  output$p_canonical_sig <- renderPlot({
    p <- plots_rv()$canonical_sig; validate(need(!is.null(p), .sig_need)); p
  })
  output$p_genomic_sig <- renderPlot({
    p <- plots_rv()$genomic_sig; validate(need(!is.null(p), .sig_need)); p
  })
  output$p_by_transcript_sig_ui <- renderUI({
    req(scan_ok())
    plotOutput("p_by_transcript_sig",
               height = paste0(max(600, 240 * scan_rv()$meta$n_facets), "px"))
  })
  output$p_by_transcript_sig <- renderPlot({
    p <- plots_rv()$by_transcript_sig; validate(need(!is.null(p), .sig_need)); p
  })

  output$tbl <- renderTable({
    req(scan_ok())
    r <- scan_rv()$result
    cols <- c("transcript_id", "is_canonical", "database", "matrix_id",
              "position", "strand", "score", "score_frac",
              "pvalue", "qvalue", "genomic_position", "site")
    r[, intersect(cols, names(r))]
  })
  output$dl_csv <- downloadHandler(
    filename = function() sprintf("%s_%s_binding_sites.csv",
                                  scan_rv()$meta$gene, scan_rv()$meta$tf),
    content  = function(file) write.csv(scan_rv()$result, file, row.names = FALSE)
  )
}

shinyApp(ui, server)
