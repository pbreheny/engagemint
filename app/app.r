library(dplyr)
library(DT)
library(magrittr)
library(shiny)
library(shinyjs)
library(shinyWidgets)

# Setup
hits_clip <- readRDS('rds/hits-clip.rds')
e_clip <- readRDS('rds/eclip.rds')
mircount <- readRDS('rds/abundance.rds')
hits_genes <- lapply(hits_clip, "[", "Gene Symbol") %>%
  unlist(use.names = FALSE) %>%
  unique()
eclip_genes <- lapply(e_clip, "[", "Gene Symbol") %>%
  unlist(use.names = FALSE) %>%
  unique()
all_genes <- c(hits_genes, eclip_genes) %>% 
  na.omit() %>%
  unique()
algorithms <- c("TargetScan 1", "PITA 1")
all_mirs <- lapply(hits_clip, "[", algorithms) %>%
  unlist(use.names = FALSE) %>%
  na.omit() %>%
  unique()

# create data frame for tracking datasets - needed for hierarchical selectInputs
metadat <- data.frame(generation = c(rep("HITS-CLIP", 2), 
                                     rep("eCLIP", 6)),
                      dataset = c("Heart", "Brain", "Aorta", "Cerebellum", "Frontal Cortex",
                                  "Gastrointestinal", "Hippocampus", "Left Ventricle"))

# Tooltips JS string
tooltips <- "header = table.columns().header();
                $(header[3]).attr('title', 
                'This is a -log10 transformation of the p-value. Often preferred for ' + 
                'visualization in genetics studies, a larger value of -log10(p) ' +
                'corresponds to a greater level of statistical significance. For ' +
                'example, the common GWAS threshold p-value of 5 C 10^-8 would equal ' +
                '7.30 in -log10 form.')"

# ui -----
ui <- fluidPage(
  useShinyjs(),
  div(style = "height:25px"), 
  
  sidebarLayout(
    sidebarPanel(
      width = 2,
      # Select generation and dataset
      selectInput("gen", "Select Generation:",
                  choices = metadat$generation),
      selectInput("dataset", "Select Dataset:", choices = NULL)
    ),
    
    mainPanel(
      width = 10,
      tabsetPanel(
        tabPanel("Gene Selection",
                 
                 # Drop down list for gene selection - server side, the list is huge
                 selectizeInput("symbol", label = "",
                                choices = NULL, selected = NULL, multiple = TRUE),
                 
                 fluidRow(
                   column(8, align = "center",
                          imageOutput("tracks", height = "40%"))),
                 
                 # Dynamic data table of microRNAs
                 DTOutput("geneList"),
                 
                 # Download .csv file
                 downloadButton("geneDownload", "Download .csv"),
                 
                 div(style = "height:25px")
        ),
        
        tabPanel("miRNA Selection",
                 # Similar drop down list for miRNAs
                 selectizeInput("miRNA", label = "",
                                choices = NULL, selected = NULL, multiple = TRUE),
                 
                 # Dynamic data table of gene hits
                 DTOutput("mirList", height = "50%"),
                 
                 # Download .csv file
                 downloadButton("mirDownload", "Download .csv"),
                 
                 div(style = "height:25px")
        ),
        
        tabPanel("miRNA Abundance",
                 
                 div(style = "height:25px"),
                 
                 DTOutput("mirAbun"),
                 
                 div(style = "height:25px")
        )
      ) 
    )
  )
)

# server -----

server <- function(input, output, session) {
  
  # Retrieve the currently selected generation and dataset
  meta <- reactive(
    filter(metadat, generation == input$gen)
  )

  observeEvent(meta(), {
    updateSelectInput(session, "dataset", choices = unique(meta()$dataset)) 
  })
  
  dataset <- reactive(switch(input$dataset,
                             Heart = hits_clip$`Heart 1st gen`,
                             Brain = hits_clip$`Brain 1st gen`,
                             Aorta = e_clip$Aorta_cluster.bed,
                             Cerebellum = e_clip$CB_cluster.bed,
                             `Frontal Cortex` = e_clip$FC_cluster.bed,
                             Gastrointestinal = e_clip$Gastroc_cluster.bed,
                             Hippocampus = e_clip$CA_cluster.bed,
                             `Left Ventricle` = e_clip$LV_cluster.bed
  ))
  
  ##### Gene Selection Tab #####
  
  # Update selectize options as the user types
  updateSelectizeInput(session, "symbol", server = TRUE,
                       choices = all_genes,
                       selected = NULL,
                       options = list(
                         placeholder = "Select a gene symbol"))
  observe(input$dataset)
  
  # Pull selected gene symbol, turn gene symbol into hyperlink
  selected_sym <- reactive(
    dataset() %>%
      filter(`Gene Symbol` %in% input$symbol) %>%
      arrange(desc(`-Log10 P Value`)))
  
  observe(input$symbol)

  # Render track image
  output$tracks <- renderImage({
    if(length(input$symbol) == 1 & input$dataset == "Heart") {
      srcstring <- paste0("Images/", input$symbol, ".png")
    } else{
      srcstring <- ""
    }
    list(src = srcstring)
  }, deleteFile = FALSE)
  
  # Table of miRNAs
  output$geneList <- renderDT(
    datatable(selected_sym() %>%
                select_if(!(names(.) %in% "Genomic Coordinates (Hg19)")),
              rownames = FALSE,
              escape = FALSE,
              callback = JS(tooltips),
              filter='top'
    )
  )
  
  # Download button
  output$geneDownload <- downloadHandler(
    filename = function() {
      paste0(paste0(input$dataset, "_"),
             paste(input$symbol, collapse = "_"),
             "_CLIP.csv")
    },
    content = function(file) {
      write.csv(selected_sym() %>%
                  select_if(!(names(.) %in% "Genomic Coordinates Link (Hg19)")),
                file,
                quote = FALSE, row.names = FALSE)
    }
  )
  
  # Toggle button / table
  observe({
    toggle(id = "geneDownload", condition = !is.null(input$symbol))
  })
  
  observe({
    toggle(id = "geneList", condition = !is.null(input$symbol))
  })
  
  # ##### miRNA Selection Tab #####
  
  updateSelectizeInput(session, "miRNA", server = TRUE, choices = all_mirs,
                       options = list(placeholder = "Select a miRNA"))
  
  # Pull selected miRNA, turn gene symbol into hyperlink
  selected_mir <- reactive(dataset() %>%
                             filter(if_any(matches(c("PITA 1", "TargetScan 1")),
                                           any_vars(. %in% input$miRNA))) %>%
                             arrange(desc(`-Log10 P Value`)))
  
  # Table of gene hits
  output$mirList <- renderDT(
    if("TargetScan 1" %in% colnames(dataset())) {
      datatable(selected_mir() %>%
                  select_if(!(names(.) %in% "Genomic Coordinates (Hg19)")),
                rownames = FALSE,
                escape = FALSE,
                callback = JS(tooltips),
                filter='top')
    } else {
      datatable(data.frame("miRNA target calling is not available for the currently selected data."),
                rownames = FALSE, colnames = c(""))
    }
  )
  
  # Download button
  output$mirDownload <- downloadHandler(
    filename = function() {
      paste0(paste0(input$dataset, "_"),
             paste(input$miRNA, collapse = "_"),
             "_CLIP.csv")
    },
    content = function(file) {
      write.csv(selected_mir() %>%
                  select_if(!(names(.) %in% "Genomic Coordinates Link (Hg19)")),
                file,
                quote = FALSE, row.names = FALSE)
    }
  )
  
  # Toggle button / table
  observe({
    toggle(id = "mirDownload", condition = !is.null(input$miRNA))
  })
  
  observe({
    toggle(id = "mirList", condition = !is.null(input$miRNA))
  })
  
  
  ##### miRNA Abundance Tab #####
  
  output$mirAbun <- renderDT(
    if(!("TargetScan 1" %in% colnames(dataset()))) {
      datatable(mircount,
                rownames = FALSE,
                class = "stripe nowrap",
                escape = FALSE,
                options = list(autoWidth = TRUE,
                               scrollX = TRUE),
                filter='top')
    } else {
      datatable(data.frame("miRNA abundance is not available for the currently selected data."),
                rownames = FALSE, colnames = c(""))        
    }
  )
  
  observe({
    toggle(id = "mirAbun", condition = !is.null(dataset()))
  })
}

shinyApp(ui = ui, server = server)
