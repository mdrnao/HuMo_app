## --------------
##
## Purpose: Shiny app to process new HCC transcriptomics into HuMo clusters
##
## Author: Holly Hall
##
## 2024
##
## packages --------------

options(shiny.maxRequestSize = 1000 * 1024^2)

## Load in packages
### Shiny bits
require(shiny)
require(shinythemes)
require(DT)
require(shinyWidgets)
### Processing and vis
require(uwot)
require(dplyr)
require(DESeq2)
require(ggplot2)
require(fst)
require(cowplot)
require(colorspace)


## load in data --------------

## Load in data
mouse_human_lookup <- readRDS("data/mouse_human_lookup.rds")
mouse_human_lookup <- mouse_human_lookup[, c(1:4)]
source("./mapToSpace.R")

umap.model <- load_uwot("./data/umap_model")
human_svd <- readRDS("./data/LIHC_TCGA_SVD_100dim.rds")
rownames <- readRDS("./data/gene_names.rds")

human_data <- list(
  embedding = read.fst("./data/human_data_umap_embedding.fst"),
  sample_info = read.fst("./data/human_data_sample_info.fst")
)
colnames(human_data$sample_info) <- "membership"

mouse_data <- list(
  embedding = read.fst("./data/mouse_data_umap_embedding.fst"),
  sample_info = read.fst("./data/mouse_data_sample_info.fst")
)

counts = read.fst("./data/all_counts.fst")

umap <- data.frame(
  rbind(human_data$embedding, mouse_data$embedding),
  cluster = c(
    human_data$sample_info$membership,
    mouse_data$sample_info$membership
  ),
  species = c(rep("Human", 371), rep("Mouse", 187))
)
colnames(umap)[1:2] <- c("umap1", "umap2")


hallmark <- read.fst("./data/hallmark_enrichments.fst")

cols <- c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")

## UI setup --------------

import_data_tab <- {
  sidebarLayout(
    sidebarPanel(
      h3("Upload data"),
      p("First column must be the gene name",
        class = "text-success"
      ),
      p("File format can be CSV or other delimited file."),
      fileInput("filedata", "Upload a file", accept = c(".csv", ".tsv")),
      h3("Process data"),
      p("Select the input gene identifier from the dropdown box, then click process to transform."),
      pickerInput(
        inputId = "gene_column",
        label = "Input gene identifier",
        choices = colnames(mouse_human_lookup),
        options = list(style = "btn-primary")
      ),
      actionBttn(
        inputId = "process",
        label = "Process",
        style = "bordered",
        color = "success",
        block = T
      ),
      h4("Once processed (bottom table appears) proceed to visualise or download tabs.")
    ),
    mainPanel(
      DTOutput(outputId = "userdataf"),
      DTOutput(outputId = "checkprog")
    )
  )
}

visualise_umap_tab <- {
  sidebarLayout(
    sidebarPanel(
      h3("View data in UMAP plot"),
      p("New data appears as boxed crosses, existing data is mouse (crosses) and TCGA data (circles)."),
      p("Existing data is coloured by HuMo group."),
      p("New data HuMo classification can be downloaded on the next tab.")
    ),
    mainPanel(
      plotOutput("umap.plot")
    )
  )
}

download_result_tab <- {
  sidebarLayout(
    sidebarPanel(
      h3("Download HuMo membership and X/Y coordinates."),
      downloadButton("download", "Download .tsv"),
      h3("Summary"),
      plotOutput("summary")
    ),
    mainPanel(
      DTOutput(outputId = "preview")
    )
  )
}

explore_hallmarks_tab <- {
  sidebarLayout(
    sidebarPanel(
      h2("Select hallmark enrichment"),
      p("In this section you can view the enrichments either aggregated by cluster (Violin),
                                    or overlaid on the UMAP", class = "text-success"),
      p("Hallmark pathway enrichment of individual samples was performed by ssGSEA. The graph will automatically
                                    update as you choose a pathway. "),
      prettyRadioButtons(
        inputId = "umapORviolin",
        label = "View hallmark enrichments by:",
        choices = c("UMAP", "Violin"),
        icon = icon("check"),
        bigger = TRUE,
        status = "primary",
        animation = "jelly"
      ),
      pickerInput(
        inputId = "hallmark_pathway", label = "Hallmark pathway", choices = colnames(hallmark),
        options = list(`live-search` = TRUE, style = "btn-primary")
      ),
    ),
    mainPanel(
      plotOutput(outputId = "selected_hallmark", width = "100%")
    )
  )
}

explore_counts_tab <- {
  sidebarLayout(
    sidebarPanel(
      h2("Gene expression by cluster"),
      p("In this section you can view the gene expression of your favourite gene either aggregated by cluster (Violin),
                                    or overlaid on the UMAP", class = "text-success"),
      p("Normalised gene expression: VST and scaled data. The graph will automatically
                                    update as you choose a gene "),
      
      prettyRadioButtons(
        inputId = "umapORviolin2",
        label = "View gene expression by:",
        choices = c("UMAP", "Violin"),
        icon = icon("check"),
        bigger = TRUE,
        status = "primary",
        animation = "jelly"
      ),
      
      pickerInput(
        inputId = "gene", 
        label = "Gene symbol", choices = rownames,
        options = list(`live-search` = TRUE, style = "btn-primary")
      ),
    ),
    mainPanel(
      plotOutput(outputId = "gene_by_humo", width = "100%")
    )
  )
}


## UI define --------------

# Define UI for application
ui <- navbarPage("HuMo",
  tabPanel("Import data", import_data_tab),
  tabPanel("Visualise", visualise_umap_tab),
  tabPanel("Download", download_result_tab),
  tabPanel("Explore: Hallmarks", explore_hallmarks_tab),
  tabPanel("Explore: Gene expression", explore_counts_tab),
  
  theme = shinytheme("flatly")
)

## server define --------------

# Define server logic required
server <- function(input, output) {
  ## Import data
  userdata <- reactive({
    req(input$filedata)
    ext <- tools::file_ext(input$filedata$name)

    switch(ext,
      csv = read.csv(input$filedata$datapath),
      tsv = read.delim(input$filedata$datapath),
      validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })

  output$userdataf <- renderDT({
    req(input$filedata)
    data.frame(userdata())
  })

  ## Process

  prepro_userdata <- eventReactive(input$process, {
    ## Pre-processing
    # Load in user data
    userdata.raw <- userdata()
    userdata.raw <- data.frame(userdata.raw)
    showNotification("Userdata in")
    # Move gene IDs to the rowname
    rownames(userdata.raw) <- userdata.raw[, 1]
    userdata.raw2 <- userdata.raw[, -1]
    # Create a new column with the right geneID name
    userdata.raw2[[input$gene_column]] <- rownames(userdata.raw2)
    showNotification("Merging")
    genenamevar <- as.character(input$gene_column)
    # Add all the gene name information to the table
    temp.userdata <- merge(userdata.raw2, mouse_human_lookup, by = genenamevar)
    rownames(temp.userdata) <- make.unique(temp.userdata$hsapiens_homolog_associated_gene_name)
    # Make the rownames the gene names, remove the rest
    temp.userdata2 <- temp.userdata[rownames$gene, ]
    temp.userdata2 <- temp.userdata2[, colnames(userdata.raw)[-1]]
    # Replace any NAs with 0
    temp.userdata2[is.na(temp.userdata2)] <- 0

    ## Normalisation
    # Create dds object for VST
    # Check we have counts
    if(sum(temp.userdata2[,1]) == 0) {
      validate("No counts - check counts table or rowname selection",
               errorClass = "warning")
    }
    dds <- DESeqDataSetFromMatrix(
      countData = temp.userdata2,
      colData = data.frame(names = colnames(temp.userdata2)),
      design = ~1
    )
    showNotification("dds made, applying vst")
    # Apply VST normalisation
    countsTable_vst <- as.data.frame(assay(vst(dds, blind = T)))
    # Scale the data
    showNotification("Scaling")
    countsTable_vst_scaled <- data.frame(t(scale(t(countsTable_vst), center = T, scale = F)))
    countsTable_vst_scaled2 <- data.frame(countsTable_vst_scaled)
    ## Map each new sample into "human sample-type" space
    showNotification("Mapping")
    human_human_space <- mapToSpace(countsTable_vst_scaled2, human_svd, rank = 100)
    # Use umap transform to position
    newembedding <- umap_transform(t(human_human_space), umap.model,
      ret_extra = "nn"
    )
  })

  output$checkprog <- renderDT({
    req(prepro_userdata())
    df <- prepro_userdata()
    head(data.frame(umap1 = df$embedding[, 1], umap2 = df$embedding[, 2]))
  })

  output$umap.plot <- renderPlot({
    df <- prepro_userdata()
    df <- data.frame(umap1 = df$embedding[, 1], umap2 = df$embedding[, 2])
    ggplot(df, aes(umap1, umap2)) +
      geom_point(
        data = umap,
        aes(
          colour = paste0("HuMo", cluster),
          shape = species
        )
      ) +
      labs(colour = "HuMo") +
      scale_colour_manual(values = cols) +
      scale_shape_manual(values = c(
        Human = 16,
        Mouse = 4
      )) +
      geom_point(size = 4, shape = 7) +
      theme_minimal() +
      panel_border()
  }, height = 600)

  user_embedding <- reactive({
    req(prepro_userdata())
    # Return the NN index from the transform function
    tmp <- as.matrix(prepro_userdata()$nn$euclidean$idx)
    # Identify the HuMo classifications of the NN
    d <- data.frame(apply(tmp, 2, function(x) umap$cluster[x]))
    # Assign memembership based on the most frequent
    d$max <- apply(d, 1, function(x) names(which.max(table(x))))
    data.frame(
      membership = d$max,
      prepro_userdata()$embedding
    )
  })

  ## Download
  output$preview <- renderDT({
    req(user_embedding())
    data.frame(user_embedding())
  })
  output$download <- downloadHandler(
    filename = function() {
      paste0("user_embedding.tsv")
    },
    content = function(file) {
      write.table(user_embedding(), file, sep = "\t")
    }
  )

  output$summary <- renderPlot({
    df <- user_embedding()
    ggplot(df, aes(x = paste0("HuMo", membership))) +
      geom_bar() +
      theme_minimal() +
      panel_border() +
      labs(x = "HuMo class")
  })


  ## Hallmarks exploration

  output$selected_hallmark <- renderPlot(
    {
      if (input$umapORviolin == "UMAP") {
        cols <- cbind(umap, hallmark)

        ggplot(cols, aes(x = umap1, y = umap2)) +
          geom_point(size = 3, aes_string(colour = input$hallmark_pathway)) +
          labs(x = "UMAP 1", y = "UMAP 2", colour = "NES") +
          theme_minimal() +
          panel_border() +
          scale_colour_continuous_diverging(palette = "Blue-Red")
      } else {
        hallmark$cluster <- paste0("HuMo", umap$cluster)

        ggplot(hallmark, aes_string(y = input$hallmark_pathway)) +
          geom_violin(aes(x = cluster, colour = cluster)) +
          geom_boxplot(width = 0.2, aes(x = cluster, colour = cluster)) +
          scale_colour_manual(values = c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")) +
          facet_wrap(~ umap$species) +
          theme_minimal() +
          panel_border() +
          theme(legend.position = "none")
      }
    },
    height = 600
  )
  
  ## Gene counts exploration
  
  output$gene_by_humo <- renderPlot(
    {
      selected_gene = as.numeric(counts[which(input$gene == rownames$genes),])
      
      selected_gene.df <- data.frame(umap, 
                         gene = selected_gene)
      
  
      if (input$umapORviolin2 == "UMAP") {
        
        ggplot(selected_gene.df, aes(x = umap1, y = umap2)) +
          geom_point(size = 3, aes(colour = gene)) +
          labs(x = "UMAP 1", y = "UMAP 2", colour = "NES") +
          theme_minimal() +
          panel_border() +
          scale_colour_continuous_diverging(palette = "Blue-Red")
      } else {
        selected_gene.df$cluster <- paste0("HuMo", selected_gene.df$cluster)
        
        ggplot(selected_gene.df, aes(y = gene)) +
          geom_violin(aes(x = cluster, colour = cluster)) +
          geom_boxplot(width = 0.2, aes(x = cluster, colour = cluster)) +
          scale_colour_manual(values = c("#440154FF", "#31688EFF", "#35B779FF", "#FDE725FF")) +
          facet_wrap(~ umap$species) +
          theme_minimal() +
          panel_border() +
          theme(legend.position = "none") +
          labs(y= as.character(input$gene))
      }
    },
    height = 600
  )
}

# Run the application
shinyApp(ui = ui, server = server)
