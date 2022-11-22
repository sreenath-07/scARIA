library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(shinythemes)
library(dplyr)
library(shinycssloaders)
library(scATAC.Explorer)
library(scRNAseq)
library(shinyjs)
library(DT)
library(rsconnect)


atacdata_table<-as.data.frame(queryATAC(metadata_only=TRUE))
accession_atac<-as.vector(unique(atacdata_table$Accession))
atacdata_table<-select(atacdata_table, Reference, Accession, Sequencing_Name, Organism, Disease, Tissue_Cell_Type, Data_Summary)
scrna_table<-as.data.frame(listDatasets())
scrna_table<-select(scrna_table, Reference, Taxonomy, Part)
accession_scrna<-as.vector(unique(scrna_table$Reference))

######UI
ui <- fluidPage(theme=shinytheme("yeti"),
                titlePanel(
                  h1(HTML("scATAC-Seq and scRNA-Seq App"),
                     style={'background-color:#158FAD; padding:20px;align: top;'}
                  )),
                h4("You can upload your own data or analyze publcily available data"),
                sidebarLayout (
                  sidebarPanel(
                    h3(HTML("Upload datasets"), style={'text-align: center;'}),
                    h5(HTML("Upload scATAC-Seq data"), style={'text-align: center; font-weight: bold;'}),
                    fileInput("atacdata", "Choose scATAC-Seq fragments.tsv File", accept = ".tsv.gz"),
                    fileInput("atacdata", "Choose scATAC-Seq H5 File", accept = ".tsv.gz"),
                    h5(HTML("Upload scRNA-Seq data"), style={'text-align: center; font-weight: bold;'}),
                    fileInput("scrnadata", "Choose scRNA-Seq RDS/H5AD/Matrix File"),
                    fileInput("scrnadata", "Choose scRNA-Seq Features File"),
                    fileInput("scrnadata", "Choose scRNA-Seq Barcode File"),
                    actionButton('integrateupload_btn', 'Analyze and Integrate',
                                 style="color: #fff; background-color: #4CAF50; text-align: center"),
                    h4(HTML("OR"), style={'text-align: center;'}),
                    h3(HTML("Select publicly available datasets"), style={'text-align: center;'}),
                    h4(HTML("Click the buttons below to see available datasets"), style={'text-align: center;'}),
                    actionButton('atac_public', 'Show scATAC-Seq Datasets',
                                 style="text-align: center"),
                    actionButton('scrna_public', 'Show scRNA-Seq Datasets',
                                 style="text-align: center"),
                    selectInput("atac_public_data", "Choose the scATAC-Seq dataset of interest",
                                choices=accession_atac),
                    selectInput("scrna_public_data", "Choose the scRNA-Seq dataset of interest",
                                choices=accession_scrna),
                    actionButton('integratepublic_btn', 'Analyze and Integrate',
                                 style="color: #fff; background-color: #4CAF50; text-align: center")
                  ),
                  mainPanel(
                    DT::dataTableOutput("mytable")
                    #plotOutput("kmplot", width = "100%") %>% withSpinner(color="#0dc5c1")
                  )
                )
)



######Server
server <- function(input, output) {
  #To upload user datasets
  
  
  #To show public scATAC-Seq datasets
  observeEvent(input$atac_public, {
    output$mytable = DT::renderDataTable({
      unique(atacdata_table) 
    }, options = list(scrollX = TRUE, scrollY = TRUE), 
    caption=htmltools::tags$caption("scATAC-Seq Datasets", style="color:black; font-size: 20px; font-weight:bold"))
  })
  
  
  #To show public scRNA-Seq datasets
  observeEvent(input$scrna_public, {
    output$mytable = DT::renderDataTable({
      unique(scrna_table)
    }, options = list(scrollX = TRUE, scrollY = TRUE),
    caption=htmltools::tags$caption("scRNA-Seq Datasets", style="color:black; font-size: 20px; font-weight:bold"))
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
