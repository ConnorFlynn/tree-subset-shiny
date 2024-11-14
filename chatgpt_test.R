library(shiny)
library(ggtree)
library(ape)
library(dplyr)
library(tibble)
library(stringr)
library(shinyalert)
library(shinyjs)
# Define UI
ui <- fluidPage(
  navbarPage(
    "Phylogenetic Tree",
    tabPanel(
      "Explore Tree",
      fluidRow(
        class = "inputs",
        column(6, fileInput("upload_tree", "Select Tree File:")),
        column(6, selectInput("file_type", "Tree File Type", choices = c("tree", "beast", "mlc"), selected = "tree"))
      ),
      uiOutput("select_node_render"),
      uiOutput("subtree_render"),
      uiOutput("patristic_render")
    )
  )
)

# Define Server
server <- function(input, output, session) {
  session$onSessionEnded(stopApp)
  
  # Reactive values for resetting inputs
  rv <- reactiveValues(
    data = NULL,
    clear = FALSE
  )
  
  observeEvent(input$upload_tree, {
    rv$clear <- FALSE
  }, priority = 1000)
  
  observeEvent(input$file_type, {
    rv$data <- NULL
    rv$clear <- TRUE
    reset('upload_tree')
  }, priority = 1000)
  
  # Reactive function to read the tree file based on input file type
  tree <- reactive({
    req(input$upload_tree, input$file_type, !rv$clear)
    
    file <- input$upload_tree$datapath
    
    output <- switch(
      input$file_type,
      tree = possibly(read.tree, otherwise = NULL)(file),
      beast = possibly(read.beast, otherwise = NULL)(file),
      mlc = possibly(read.codeml_mlc, otherwise = NULL)(file),
      jplace = possibly(read.jplace, otherwise = NULL)(file),
      mrbayes = possibly(read.mrbayes, otherwise = NULL)(file),
      nhx = possibly(read.nhx, otherwise = NULL)(file),
      rst = possibly(read.paml_rst, otherwise = NULL)(file),
      phylip = possibly(read.phylip, otherwise = NULL)(file),
      r8s = possibly(read.r8s, otherwise = NULL)(file),
      raxml = possibly(read.raxml, otherwise = NULL)(file)
    )
    
    return(output)
  })
  
  # Data frame of the tree
  tree_df <- reactive({
    req(tree())
    as_data_frame(tree())
  })
  
  observe({
    req(input$upload_tree)
    
    if (is.null(tree())) {
      shinyalert("Tree import error", paste("There was an error when trying to read your tree!",
                                            "Did you select the correct file format?"),
                 type = "error")
    }
  })
  
  # UI for selecting node and plot settings
  output$select_node_render <- renderUI({
    req(input$upload_tree, tree())
    tagList(
      fluidRow(
        column(
          12,
          selectizeInput(
            inputId = "select_node",
            label = "Select Node:",
            choices = tree_df() %>% 
              select(label) %>% 
              arrange(label) %>% 
              pull(label),
            width = "100%"
          )
        )
      ),
      fluidRow(
        column(3, numericInput("subtree_levels_back", "Select Number of Ancestral Levels:", min = 1, value = 10)),
        column(3, numericInput("subtree_text_size", "Select label text size:", min = 2, value = 3)),
        column(3, numericInput("subtree_plot_height", "Select plot height", value = 1200)),
        column(3, numericInput("subtree_width_multiply", "Select plot width multiplier:", value = 1.4, min = 1, step = 0.1))
      )
    )
  })
  
  # Render subtree plot
  # Render subtree plot
  output$subtree <- renderPlot({
    req(input$select_node, tree(), input$subtree_width_multiply, input$subtree_text_size, input$subtree_plot_height)
    
    # Get the original tree and selected node
    full_tree <- tree()
    
    # Find the descendants of the selected node
    node_id <- match(input$select_node, full_tree$tip.label)
    descendants <- ape::Descendants(full_tree, node = node_id, type = "all")[[1]]
    
    # Drop tips that are NOT descendants of the selected node to create the subtree
    sub_tree <- ape::drop.tip(full_tree, setdiff(full_tree$tip.label, full_tree$tip.label[descendants]))
    
    # Labels and plotting logic
    labels <- sub_tree$tip.label
    
    labels_df <- tibble(
      label = labels,
      genus = str_extract(label, "[^;]+;[^;]+$") %>% str_replace(";[^;]+$", ""),
      species = str_extract(label, "[^;]+$")
    ) %>% 
      mutate(
        species = if_else(is.na(genus), "", str_replace(species, "s__", "")),
        genus = if_else(is.na(genus), label, str_replace(genus, "g__", ""))
      )
    
    p <- ggtree(sub_tree, aes(color = group)) %<+% labels_df +
      geom_tiplab(aes(label = paste(genus, species)), size = input$subtree_text_size) +
      theme_tree2() +
      scale_color_manual(values = c(`1` = "red", `0` = "black"))
    
    p + lims(x = c(0, max(p$data$x) * input$subtree_width_multiply))
  })
  
  
  # Patristic distance calculation
  patristic_distances <- reactive({
    req(tree())
    dist_matrix <- cophenetic(tree())  # Calculate patristic distances
    dist_matrix
  })
  
  # Render patristic distance table
  output$patristic_table <- renderTable({
    req(patristic_distances())
    as.data.frame(as.table(patristic_distances()))
  })
  
  # UI element for patristic distance table
  output$patristic_render <- renderUI({
    req(patristic_distances())
    tableOutput("patristic_table")
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
