library(DT)
library(shiny)
library(shinyjs)
library(shinyWidgets)

# Setup
hits_clip <- readRDS('rds/hits-clip.rds')
e_clip <- readRDS('rds/eclip.rds')

# ui -----

ui <- fluidPage(
  useShinyjs(),
  
  div(style = "height:25px"), 
  
  tabsetPanel(
    tabPanel("eCLIP",
             DTOutput("eCLIPtable"),
             downloadButton("eCLIPdownload", "Download .csv"),
    ),
    
    tabPanel("HITSCLIP",
             DTOutput("HITSCLIPtable"),
             downloadButton("HITSCLIPdownload", "Download .csv"),
    )
  ),
  
  div(style = "height:25px")
) 

# server -----

server <- function(input, output, session) {
  
  # Table of eCLIP miRNAs
  output$eCLIPtable <- renderDT(
    datatable(e_clip,
              rownames = FALSE,
              filter = "top",
              escape = FALSE
    )
  )
  
  # eCLIP download button
  output$eCLIPdownload <- downloadHandler(
    filename = "eCLIPtable.csv",
    content = function(file) {
      write.csv(e_clip,
                file,
                quote = FALSE, row.names = FALSE)
    }
  )
  
  # Table of HITS-CLIP miRNAs
  output$HITSCLIPtable <- renderDT(
    datatable(hits_clip[,!(colnames(hits_clip) %in% "Genomic Coordinates (Hg19)")],
              rownames = FALSE,
              filter = "top",
              escape = FALSE
    )
  )
  
  # HITS-CLIP download button
  output$HITSCLIPdownload <- downloadHandler(
    filename = "HITSCLIPtable.csv",
    content = function(file) {
      write.csv(hits_clip[,!(colnames(hits_clip) %in% "Genomic Coordinates Link (Hg19)")],
                file,
                quote = FALSE, row.names = FALSE)
    }
  )
  
}

shinyApp(ui = ui, server = server)
