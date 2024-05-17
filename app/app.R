library(dplyr)
library(DT)
library(magrittr)
library(shiny)
library(shinyjs)
library(shinyWidgets)

# create RDS objects
# library(readxl)
# path <- "R:/rdss_pbreheny/grade/CLIP_bg/miR-CLIP database build example table_Boudreau lab.xlsx"
# miRdata <- lapply(excel_sheets(path), read_excel, 
#                    path = path, na = "NA")
# names(miRdata) <- excel_sheets(path)
# saveRDS(miRdata, "miRdata.rds")
# 
# library(vroom)
# mircount <- vroom("R:/rdss_pbreheny/grade/CLIP_bg/miR quant/master_miRcount.csv",
#                   col_select = 1:26) %>% na.omit
# saveRDS(mircount, "miRcount.rds")
# 
# bedpaths <- list.files("R:/rdss_pbreheny/grade/CLIP_bg/bed12",
#                        full.names = TRUE)
# bedfiles <- list()
# for(i in 1:length(bedpaths)) {
#   bedfiles[[i]] <- read.delim(bedpaths[i], sep = "\t", header = FALSE)
# }
# names(bedfiles) <- gsub(".*\\/", "", bedpaths)
# saveRDS(bedfiles, "bedfiles.rds")

##### Preliminaries #####

# create data frame for tracking datasets - needed for hierarchical selectInputs
metadat <- data.frame(generation = c(rep("HITS-CLIP", 2), 
                                     rep("eCLIP", 6)),
                      dataset = c("Heart", "Brain", "Aorta", "Cerebellum", "Frontal Cortex",
                                  "Gastrointestinal", "Hippocampus", "Left Ventricle"))

# Hyperlink function
hyperlink <- function(link, text) {
  paste0("<a href='",
         link,
         "'target='_blank'>",
         text)
}

# read in original HITS-CLIP peak data
mir <- readRDS("miRdata.rds") %>%
  lapply(function(x) {
    x %>%
      mutate(`P-value` = ifelse(`P-value` == 0 | `P-value` < 1e-100,
                                1e-100,
                                `P-value`)) %>%
      mutate(`Avg Peak Height` = round(`Avg Peak Height`, 2),
             `-Log10 P Value` = round(-log10(`P-value`), 2)) %>%
      relocate(`-Log10 P Value`, .after = `Avg Peak Height`) %>%
      select(-`P-value`) %>%
      mutate(`Genomic Coordinates Link (Hg19)` =
               hyperlink(`UCSC Link (can just be a button)`,
                         `Genomic Coordinates (Hg19)`)) %>%
      relocate(`Genomic Coordinates Link (Hg19)`) %>%
      select(-`UCSC Link (can just be a button)`)
  })

# read in new eCLIP initial bed files
bedfiles <- readRDS("bedfiles.rds") %>%
  lapply(function(x) {
    select(x, c(Chromosome = V1, Strand = V6, `Gene Symbol` = V4, 
                `-Log10 P Value` = V5, `Window Start` = V2, `Window End` = V3)) %>%
    mutate(`-Log10 P Value` = ifelse(`-Log10 P Value` > 100, 
                                     100.00, round(`-Log10 P Value`, 2)))
  })

# read in miRNA abundance
mircount <- readRDS("miRcount.rds") %>%
  select(-c(Species_ID,
            grep("fam_sum", colnames(.)),
            grep("geomean", colnames(.))))
colnames(mircount) <- c("miR Family 6", "miR Family 8", "Seed+m8",
                        "MiRBase ID", "Mature Sequence Family",
                        "Conservation", "MiRBase Accession", "Left Ventricle<br>Rank",
                        "Aorta<br>Rank", "Gastrointestinal<br>Rank", "Hippocampus<br>Rank",
                        "Cerebellum<br>Rank", "Frontal Cortex<br>Rank")
mircount <- mircount[,c("MiRBase ID", "Aorta<br>Rank", "Cerebellum<br>Rank", 
                        "Frontal Cortex<br>Rank", "Gastrointestinal<br>Rank",
                        "Hippocampus<br>Rank", "Left Ventricle<br>Rank",
                        "Conservation", "Seed+m8", "miR Family 6", "miR Family 8",
                        "Mature Sequence Family")]

# Generate list of all possible genes
hits_genes <- lapply(mir, "[", "Gene Symbol") %>%
  unlist(use.names = FALSE)
eclip_genes <- lapply(bedfiles, "[", "Gene Symbol") %>%
  unlist(use.names = FALSE)

all_genes <- c(hits_genes, eclip_genes) %>% 
  na.omit() %>%
  unique()

# Generate list of all possible miRNAs
algorithms <- c("TargetScan 1", "PITA 1")
all_miRNAs <- lapply(mir, "[", algorithms) %>%
  unlist(use.names = FALSE) %>%
  na.omit() %>%
  unique()

# Tooltips JS string
tooltips <- "header = table.columns().header();
                $(header[3]).attr('title', 
                'This is a -log10 transformation of the p-value. Often preferred for ' + 
                'visualization in genetics studies, a larger value of -log10(p) ' +
                'corresponds to a greater level of statistical significance. For ' +
                'example, the common GWAS threshold p-value of 5 C 10^-8 would equal ' +
                '7.30 in -log10 form.')"

##### App #####

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

server <- function(input, output, session) {
  
  # Retrieve the currently selected generation and dataset
  meta <- reactive(filter(metadat, generation == input$gen))
  
  observeEvent(meta(), {
    updateSelectInput(session, "dataset", choices = unique(meta()$dataset)) 
  })

  dataset <- reactive(switch(input$dataset,
                             Heart = mir$`Heart 1st gen`,
                             Brain = mir$`Brain 1st gen`,
                             Aorta = bedfiles$AOpval.bed,
                             Cerebellum = bedfiles$CBpval.bed,
                             `Frontal Cortex` = bedfiles$FCpval.bed,
                             Gastrointestinal = bedfiles$SKMpval.bed,
                             Hippocampus = bedfiles$CApval.bed,
                             `Left Ventricle` = bedfiles$LVpval.bed
                            ))

##### Gene Selection Tab #####

# Update selectize options as the user types
updateSelectizeInput(session, "symbol", server = TRUE,
                     choices = all_genes,
                     selected = NULL,
                     options = list(
                       placeholder = "Select a gene symbol"))

observe(input$dataset,
        )
  # Pull selected gene symbol, turn gene symbol into hyperlink
  selected_sym <- reactive(dataset() %>%
                             filter(`Gene Symbol` %in% input$symbol) %>%
                             arrange(desc(`-Log10 P Value`)))
  
  observe(input$symbol,
          )
  
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
              callback = JS(tooltips)
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

  updateSelectizeInput(session, "miRNA", server = TRUE, choices = all_miRNAs,
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
                callback = JS(tooltips)) 
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
                               scrollX = TRUE)
        )
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
