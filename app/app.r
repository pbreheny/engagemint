library(dplyr)
library(DT)
library(magrittr)
library(shiny)
library(shinyjs)
library(shinyWidgets)

# Setup
hits_clip <- readRDS('rds/hits-clip.rds')
e_clip <- readRDS('rds/eclip.rds')

# function to render DT
DTrender <- function(x) {
  renderDT(
    datatable(x,
              rownames = FALSE,
              filter = "top",
              escape = TRUE
    )
  )
}

# function to create download button
DLbutton <- function(table, filepath) {
  downloadHandler(
    filename = filepath,
    content = function(file) {
      write.csv(table,
                file,
                quote = FALSE, row.names = FALSE)
    }
  )
}

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
  output$eCLIPtable <- DTrender(e_clip)
  
  # eCLIP download button
  output$eCLIPdownload <- DLbutton(eCLIPtable, "eCLIPtable.csv")
  
  # Table of HITS-CLIP miRNAs
  output$HITSCLIPtable <- DTrender(hits_clip)
  
  # HITS-CLIP download button
  output$HITSCLIPdownload <- DLbutton(HITSCLIPtable, "HITSCLIPtable.csv")
}

shinyApp(ui = ui, server = server)
