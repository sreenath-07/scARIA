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
library(EnsDb.Hsapiens.v75)
library(Signac)
library(SeuratData)
library(Seurat)
library(ggplot2)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)

atacdata_table<-as.data.frame(queryATAC(metadata_only=TRUE))
accession_atac<-as.vector(unique(atacdata_table$Accession))
atacdata_table<-dplyr::select(atacdata_table, Reference, Accession, Sequencing_Name, Organism, Disease, Tissue_Cell_Type, Data_Summary)
scrna_table<-as.data.frame(listDatasets())
scrna_table<-dplyr::select(scrna_table, Reference, Taxonomy, Part)
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
                    fileInput("atacfrag", "Choose scATAC-Seq fragments.tsv File*"),
                    fileInput("atach5", "Choose scATAC-Seq H5 File*"),
                    fileInput("atacmeta", "Choose scATAC-Seq Metadata File*"),
                    h5(HTML("Upload scRNA-Seq data"), style={'text-align: center; font-weight: bold;'}),
                    fileInput("scrnamatrix", "Choose scRNA-Seq RDS/H5AD/Matrix File*"),
                    fileInput("scrnafeatures", "Choose scRNA-Seq Features File"),
                    fileInput("scrnabarcode", "Choose scRNA-Seq Barcode File"),
                    actionButton('integrateupload_btn', 'Analyze and Integrate',
                                 style="color: #fff; background-color: #4CAF50; text-align: center"),
                    h6(HTML("*=required"), style={'text-align: center; font-weight: bold; color: red;'}),
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
                    tabsetPanel(type = "tabs",
                                tabPanel("Refined Clusters", plotOutput("refined_clust")),
                                tabPanel("scATAC-Seq", plotOutput("atacplot")),
                                tabPanel("scRNA-Seq", plotOutput("screeplot")))
                    #plotOutput("kmplot", width = "100%") %>% withSpinner(color="#0dc5c1")
                  )
                )
)



######Server
server <- function(input, output) {
  #To upload user datasets
  options(shiny.maxRequestSize=3000*1024^2)
  
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
  
  ##Analysis and Integration
  observeEvent(input$integrateupload_btn, {
    #Check if mandatory files have been uploaded
    
    if(is.null(input$atacfrag) | is.null(input$atach5) | is.null(input$scrnamatrix))
    {showModal(modalDialog(
      title = "Missing File(s)!",
      paste0("You have not uploaded one or more required files. Please ensure you have uploaded the mandatory input files"),
      easyClose = TRUE,
      footer = NULL
    ))
    }
   else{ 
   ##ATAC Pipeline
     print("hello")
     counts <- Read10X_h5(filename = input$atach5$datapath)
     metadata <- read.csv(
       file = input$atacmeta$datapath,
       header = TRUE,
       row.names = 1
     )
     
     chrom_assay <- CreateChromatinAssay(
       counts = counts,
       sep = c(":", "-"),
       genome = 'hg19',
       fragments = "C:/Users/varsh/Desktop/GT Research/BIOL 8803/atac_v1_pbmc_10k_fragments.tsv.gz",
       min.cells = 10,
       min.features = 200
     )
     
     pbmc <- CreateSeuratObject(
       counts = chrom_assay,
       assay = "peaks",
       meta.data = metadata
     ) 
     
     
     # extract gene annotations from EnsDb
     annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
     
     # change to UCSC style since the data was mapped to hg19
     seqlevelsStyle(annotations) <- 'UCSC'
     
     # add the gene information to the object
     Annotation(pbmc) <- annotations
     
     # compute nucleosome signal score per cell
     pbmc <- NucleosomeSignal(object = pbmc)
     
     # compute TSS enrichment score per cell
     pbmc <- TSSEnrichment(object = pbmc, fast = TRUE)
     
     # add blacklist ratio and fraction of reads in peaks
     pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
     pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
     
     #Remove cells that are outliers for these QC metrics
     pbmc <- subset(
       x = pbmc,
       subset = peak_region_fragments > 3000 &
         peak_region_fragments < 20000 &
         pct_reads_in_peaks > 15 &
         blacklist_ratio < 0.05 &
         nucleosome_signal < 4 &
         TSS.enrichment > 2
     )
     
     #Normalization and linear dimensional reduction
     pbmc <- RunTFIDF(pbmc)
     pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
     pbmc <- RunSVD(pbmc)
     
     #Nonlinear dimenison reduction and clustering
     pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
     pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
     pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
     
     gene.activities <- GeneActivity(pbmc)
     # add the gene activity matrix to the Seurat object as a new assay and normalize it
     pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
     pbmc <- NormalizeData(
       object = pbmc,
       assay = 'RNA',
       normalization.method = 'LogNormalize',
       scale.factor = median(pbmc$nCount_RNA)
     )
    
     
     DefaultAssay(pbmc) <- 'RNA'
     FeaturePlot(
       object = pbmc,
       features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
       pt.size = 0.1,
       max.cutoff = 'q95',
       ncol = 3
     )
     
     
    ##RNA-Seq
     pbmc_rna <- readRDS(input$scrnamatrix$datapath)
     
     transfer.anchors <- FindTransferAnchors(
       reference = pbmc_rna,
       query = pbmc,
       reduction = 'cca'
     )
    #DimPlot(object = pbmc, label = TRUE) + NoLegend()
    
   
   ##Integration
   
   }
  })
  ####Results: 3 plots and one table
  output$atacplot<-renderPlot({
  DimPlot(object = pbmc, label = TRUE) + NoLegend()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
