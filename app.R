# app.R — DE + Expression Matrix Viewer
# ------------------------------------------------------------------------------
# Upload:
#   - Metadata (CSV/TXT/TSV/XLSX) with a 'sample' column (IDs)
#   - Expression matrix (CSV/XLSX): MUST HAVE a 'gene' column (rows=genes; other columns = samples)
#     (can be TPM, counts, or normalized counts)
#   - DE results:
#       * RDS (DESeqResults or data.frame with rownames=genes), OR
#       * CSV/XLSX table that MUST HAVE a 'gene' column (gene symbols or IDs)
# Plot:
#   - Choose ANY metadata column as X; optional Facet; optional "X × Facet" as X
#   - Jitter aligned to boxes (shape 21, black outline)
#   - Custom main title and axis titles
# DE table:
#   - Auto-maps common column aliases; requires 'gene' for non-RDS
# ------------------------------------------------------------------------------
options(shiny.maxRequestSize = 1024^3)  # ~1 GB
suppressPackageStartupMessages({
  library(shiny); library(DT); library(dplyr); library(tidyr)
  library(tibble); library(ggplot2); library(readr); library(tools); library(readxl)
})

# ---------------------- Helpers ----------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
nz_or <- function(x, default) if (is.null(x) || (is.character(x) && !nzchar(x))) default else x

# Read CSV/XLSX/TSV/TXT into data.frame 
read_any_table <- function(path, sheet = 1) {
  ext <- tolower(file_ext(path))
  if (ext == "csv") {
    read.csv(path, check.names = FALSE)
  } else if (ext %in% c("tsv","txt")) {
    read.delim(path, check.names = FALSE)
  } else if (ext %in% c("xlsx","xls")) {
    as.data.frame(readxl::read_excel(path, sheet = sheet))
  } else {
    validate(need(FALSE, paste("Unsupported format:", ext, "| Allowed: .csv, .tsv, .txt, .xlsx/.xls")))
  }
}

# Read DE results from .rds; standardize columns
read_de_rds <- function(path, design_name) {
  x <- readRDS(path)
  if (inherits(x, "DESeqResults") || inherits(x, "DFrame") || inherits(x, "DataFrame")) x <- as.data.frame(x)
  if (!("gene" %in% names(x))) x <- tibble::rownames_to_column(x, "gene")
  ren <- c("logFC"="log2FoldChange", "lfc"="log2FoldChange",
           "FDR"="padj", "qvalue"="padj", "PValue"="pvalue", "p"="pvalue")
  for (nm in names(ren)) if (nm %in% names(x) && !(ren[[nm]] %in% names(x))) names(x)[names(x)==nm] <- ren[[nm]]
  tibble(
    design   = design_name,
    gene     = x$gene,
    log2FC   = suppressWarnings(as.numeric(x$log2FoldChange)),
    padj     = if ("padj" %in% names(x)) suppressWarnings(as.numeric(x$padj)) else NA_real_,
    baseMean = if ("baseMean" %in% names(x)) suppressWarnings(as.numeric(x$baseMean)) else NA_real_
  )
}

# Read DE results from CSV/XLSX; MUST HAVE a 'gene' column (no fallback)
read_de_table <- function(path, design_name, sheet = 1) {
  x <- read_any_table(path, sheet = sheet)
  validate(need("gene" %in% names(x),
                "DE table must include a 'gene' column (symbols or IDs)."))
  ren <- c("logFC"="log2FoldChange", "lfc"="log2FoldChange",
           "FDR"="padj", "qvalue"="padj", "PValue"="pvalue", "p"="pvalue")
  for (nm in names(ren)) if (nm %in% names(x) && !(ren[[nm]] %in% names(x))) names(x)[names(x)==nm] <- ren[[nm]]
  validate(need("log2FoldChange" %in% names(x),
                "DE table must have 'log2FoldChange' (or alias: logFC/lfc)."))
  tibble(
    design   = design_name,
    gene     = x$gene,
    log2FC   = suppressWarnings(as.numeric(x$log2FoldChange)),
    padj     = if ("padj" %in% names(x)) suppressWarnings(as.numeric(x$padj)) else NA_real_,
    baseMean = if ("baseMean" %in% names(x)) suppressWarnings(as.numeric(x$baseMean)) else NA_real_
  )
}

# Plot builder
plot_base <- function(df, xvar, facet_var = NULL, main_title = NULL,
                      xlab = NULL, ylab = "Expression", jitter_on = TRUE,
                      facet_scales = "free_y") {
  stopifnot(xvar %in% names(df))
  df[[xvar]] <- as.factor(df[[xvar]])
  p <- ggplot(df, aes(x = .data[[xvar]], y = Expression, fill = gene)) +
    geom_boxplot(
      aes(group = interaction(.data[[xvar]], gene)),
      width = 0.55, outlier.shape = NA, alpha = 0.9,
      position = position_dodge(width = 0.70)
    )
  if (isTRUE(jitter_on)) {
    p <- p + geom_jitter(
      aes(group = interaction(.data[[xvar]], gene), fill = gene),
      shape = 21, color = "black", stroke = 0.25,
      size = 2, alpha = 0.7,
      position = position_jitterdodge(jitter.width = 0.05, jitter.height = 0, dodge.width = 0.70)
    )
  }
  p <- p +
    labs(title = main_title, x = xlab %||% xvar, y = ylab) +
    scale_x_discrete(drop = TRUE) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  if (!is.null(facet_var) && facet_var != "None" && facet_var %in% names(df)) {
    df[[facet_var]] <- as.factor(df[[facet_var]])
    p <- p + facet_wrap(stats::as.formula(paste0("~", facet_var)), scales = facet_scales)
  }
  p
}

# ---------------------------- UI -----------------------------------
ui <- navbarPage(
  title = "DE + Expression Matrix Viewer",
  id = "topnav",
  
  # ---- Main App tab ----
  tabPanel(
    "App",
    sidebarLayout(
      sidebarPanel(
        h4("Upload data"),
        fileInput("meta_file", "Metadata", accept = c(".csv",".tsv",".txt",".xlsx",".xls")),
        helpText("CSV/TSV/TXT/XLSX; must include a 'sample' column (IDs)."),
        tags$hr(),
        fileInput("tpm_file",  "Expression matrix", accept = c(".csv",".xlsx",".xls")),
        checkboxInput("tpm_is_log", "Values already log2-transformed? (e.g., log2(x+1))", TRUE),
        helpText("CSV/XLSX; MUST include a 'gene' column (rows=genes; columns=sample IDs)."),
        tags$hr(),
        h5("DESeq2 results"),
        fluidRow(
          column(7, fileInput("de1", NULL, accept = c(".rds",".csv",".xlsx",".xls"),
                              buttonLabel = "Add DE file")),
          column(5, textInput("de1_name", "Design label", value = "Design1"))
        ),
        fluidRow(
          column(7, fileInput("de2", NULL, accept = c(".rds",".csv",".xlsx",".xls"),
                              buttonLabel = "Add DE file")),
          column(5, textInput("de2_name", "Design label", value = "Design2"))
        ),
        fluidRow(
          column(7, fileInput("de3", NULL, accept = c(".rds",".csv",".xlsx",".xls"),
                              buttonLabel = "Add DE file")),
          column(5, textInput("de3_name", "Design label", value = "Design3"))
        ),
        helpText("If not RDS, DE table must include 'gene' and 'log2FoldChange' (or alias logFC/lfc)."),
        tags$hr(),
        h4("Filters & options"),
        selectizeInput("genes", "Gene(s)", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "Type gene IDs/symbols")),
        numericInput("alpha",  "FDR threshold (padj ≤)", value = 0.05, step = 0.01, min = 0),
        numericInput("lfc_thr","|log2FC| ≥", value = 0, step = 0.1, min = 0),
        tags$hr(),
        h4("Plot mappings"),
        uiOutput("meta_map_ui"),
        checkboxInput("use_interaction", "Use (X × Facet) as X-axis", FALSE),
        tags$hr(),
        h4("Plot titles"),
        textInput("plot_title", "Main title", ""),
        textInput("plot_xlabel", "X-axis label (blank = use column name)", ""),
        textInput("plot_ylabel", "Y-axis label", "Expression"),
        checkboxInput("jitter", "Add jitter on boxes", TRUE),
        downloadButton("dl_tbl",  "Download DE table (CSV)"),
        downloadButton("dl_plot", "Download Plot (PDF)")
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("DE table", DTOutput("de_table")),
          tabPanel("Plot", plotOutput("main_plot", height = "500px"))
        )
      )
    )
  ),
  
  # ---- Info tab ----
  tabPanel(
    "Info",
    fluidPage(
      h2("About this app"),
      p("This viewer loads a metadata table, an expression matrix (TPM / counts / normalized counts), ",
        "and optional DE results to let you:"),
      tags$ul(
        tags$li("Search genes and view their DE across designs."),
        tags$li("Plot per-gene expression across any metadata factor (with optional faceting)."),
        tags$li("Download DE tables and plots.")
      ),
      h3("File formats"),
      tags$ul(
        tags$li(strong("Metadata:"), " CSV/TSV/TXT/XLSX; must include a 'sample' column (IDs)."),
        tags$li(strong("Expression matrix:"), " CSV/XLSX; rows=genes with a 'gene' column; columns=sample IDs."),
        tags$li(strong("DE results:"), " RDS (DESeqResults or data.frame) or CSV/XLSX with 'gene' and ",
                "'log2FoldChange' (or 'logFC'/'lfc'); optional 'padj'.")
      ),
      h3("Notes"),
      tags$ul(
        tags$li("If expression values are not log2(x+1), uncheck the log box to auto-transform."),
        tags$li("When using 'X × Facet as X', facets use ", code("free_x"), " and drop unused x-levels.")
      )
    )
  ),
  
  # ---- Contact tab ----
  tabPanel(
    "Contact",
    fluidPage(
      h2("Contact"),
      p("Questions, bugs, or feature requests?"),
      tags$ul(
        tags$li(HTML("Email: <a href='mailto:aysevilpektas@gmail.com'>aysevilpektas@gmail.com</a>")),
        tags$li(HTML("GitHub: <a href='https://github.com/aysevllpkts' target='_blank' rel='noopener noreferrer'>@aysevllpkts</a>"))
      ),
    )
  )
)

# --------------------------- Server --------------------------------
server <- function(input, output, session){
  # ---- Metadata ----
  meta_raw <- reactive({
    req(input$meta_file)
    df <- read_any_table(input$meta_file$datapath)
    validate(need("sample" %in% names(df), "Metadata must include a 'sample' column (IDs)."))
    df
  })
  
  output$meta_map_ui <- renderUI({
    req(meta_raw())
    cn <- names(meta_raw())
    tagList(
      selectInput("x_var",     "X-axis (categorical)", choices = cn),
      selectInput("facet_var", "Facet by (optional)",  choices = c("None", cn), selected = "None")
    )
  })
  
  coldata <- reactive({
    df <- meta_raw()
    names(df) <- tolower(names(df))
    df$sample <- as.character(df$sample)
    df
  })
  
  # ---- Expression matrix — MUST HAVE 'gene' column ----
  expr_mat <- reactive({
    req(input$tpm_file)
    df <- read_any_table(input$tpm_file$datapath)
    validate(need("gene" %in% names(df),
                  "Expression matrix must include a 'gene' column (rows=genes)."))
    rownames(df) <- df$gene
    df <- df[, setdiff(names(df), "gene"), drop = FALSE]
    validate(need("sample" %in% names(coldata()), "Metadata not loaded yet."))
    common <- intersect(colnames(df), coldata()$sample)
    validate(need(length(common) > 0,
                  "No overlapping sample IDs between expression matrix and metadata."))
    df[, common, drop = FALSE]
  })
  
  observe({
    if (!is.null(expr_mat())) updateSelectizeInput(session, "genes", choices = rownames(expr_mat()), server = TRUE)
  })
  
  expr_long <- reactive({
    df <- expr_mat() %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "sample", values_to = "Expression") %>%
      mutate(Expression = as.numeric(Expression)) %>%
      left_join(coldata(), by = "sample")
    if (!isTRUE(input$tpm_is_log)) df <- df %>% mutate(Expression = log2(Expression + 1))
    df
  })
  
  # ---- DE results (RDS or CSV/XLSX) ----
  de_all <- reactive({
    files <- list(
      list(inp = input$de1, name = nz_or(input$de1_name, "Design1")),
      list(inp = input$de2, name = nz_or(input$de2_name, "Design2")),
      list(inp = input$de3, name = nz_or(input$de3_name, "Design3"))
    )
    pieces <- list()
    for (f in files) {
      if (is.null(f$inp$datapath)) next
      ext <- tolower(file_ext(f$inp$datapath))
      if (ext == "rds") {
        pieces <- c(pieces, list(read_de_rds(f$inp$datapath, f$name)))
      } else if (ext %in% c("csv","xlsx","xls")) {
        pieces <- c(pieces, list(read_de_table(f$inp$datapath, f$name)))
      } else {
        validate(need(FALSE, paste("Unsupported DE file type:", ext, "| Allowed: .rds, .csv, .xlsx/.xls")))
      }
    }
    if (length(pieces) == 0) return(NULL)
    bind_rows(pieces)
  })
  
  de_filtered <- reactive({
    req(de_all())
    df <- de_all()
    if (length(input$genes)) df <- df %>% filter(gene %in% input$genes)
    df %>%
      mutate(is_sig = !is.na(padj) & padj <= input$alpha & abs(log2FC) >= input$lfc_thr) %>%
      arrange(desc(is_sig), padj, desc(abs(log2FC)))
  })
  
  output$de_table <- renderDT({
    validate(need(!is.null(de_all()), "Upload one or more DE results (.rds/.csv/.xlsx) to see the DE table."))
    df <- de_filtered()
    datatable(df, rownames = FALSE, filter = "top",
              options = list(pageLength = 25, dom = "Bfrtip",
                             buttons = c("copy","csv","excel"),
                             scrollX = TRUE),
              extensions = c("Buttons")) %>%
      formatRound("log2FC", 3) %>%
      formatSignif("padj", 3) %>%
      formatStyle("is_sig", target = "row",
                  backgroundColor = styleEqual(c(TRUE, FALSE), c("#e8ffe8", "transparent"))) %>%
      formatStyle("log2FC", color = styleInterval(0, c("blue","red")), fontWeight = "bold")
  })
  
  output$dl_tbl <- downloadHandler(
    filename = function() paste0("DE_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content  = function(f) readr::write_csv(de_filtered(), f)
  )
  
  # ---- Plot ----
  make_plot <- reactive({
    req(input$genes, expr_long(), input$x_var)
    df <- expr_long() %>% filter(gene %in% input$genes)
    xv <- input$x_var
    fv <- if (isTruthy(input$facet_var) && input$facet_var != "None") input$facet_var else NULL
    xv <- tolower(xv); if (!is.null(fv)) fv <- tolower(fv)
    validate(need(xv %in% names(df), "Selected X-axis column not found in metadata."))
    ylab_txt <- if (nzchar(input$plot_ylabel)) input$plot_ylabel else "Expression"
    main_t   <- if (nzchar(input$plot_title)) input$plot_title else NULL
    xlab_txt <- if (nzchar(input$plot_xlabel)) input$plot_xlabel else NULL
    if (isTRUE(input$use_interaction) && !is.null(fv) && fv %in% names(df)) {
      df <- df %>% mutate(.interaction = paste0(.data[[xv]], "_", .data[[fv]]))
      plot_base(df, xvar = ".interaction", facet_var = fv,
                main_title = main_t, xlab = xlab_txt, ylab = ylab_txt,
                jitter_on = input$jitter, facet_scales = "free_x")
    } else {
      plot_base(df, xvar = xv, facet_var = fv,
                main_title = main_t, xlab = xlab_txt, ylab = ylab_txt,
                jitter_on = input$jitter, facet_scales = "free_y")
    }
  })
  
  output$main_plot <- renderPlot({ make_plot() })
  
  # Plot download (PDF)
  output$dl_plot <- downloadHandler(
    filename = function() paste0("Expression_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"),
    content = function(file) {
      g <- make_plot()
      ggsave(file, plot = g, width = 9, height = 5, dpi = 300, device = "pdf")
    }
  )
}

shinyApp(ui, server)