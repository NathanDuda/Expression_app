





# Load libraries
library(shiny)
library(dplyr)
library(shinycssloaders)  # For loading animations
library(shinythemes)      # For theme support

# Define the UI
ui <- fluidPage(
  theme = shinytheme("flatly"),
  tags$head(
    tags$style(HTML("
    .shiny-input-container { margin-bottom: 20px; }
    h4 { color: #5c9fd7; }  /* Updated light blue color */
    .btn-primary { background-color: #5c9fd7; border-color: #5c9fd7; }
    .btn-primary:hover { background-color: #4a8cbb; border-color: #4a8cbb; }  /* Slightly darker for hover */
    table { margin-top: 20px; border-collapse: collapse; width: 100%; }
    th, td { padding: 10px; text-align: left; border: 1px solid #ddd; }
    th { background-color: #f2f2f2; }
    tr:nth-child(even) { background-color: #f9f9f9; }
    tr:hover { background-color: #f1f1f1; }
  "))
  ),
  
  titlePanel("Gene Expression Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Select up to 10 genes:"),
      selectizeInput(
        "selected_genes", "Genes:",
        choices = NULL,  # Initially empty, populated dynamically
        multiple = TRUE, 
        options = list(
          maxItems = 10,              # Limit selections to 10
          placeholder = 'Start typing to search for genes',
          create = FALSE              # Prevent users from entering new (non-existent) genes
        )
      ),
      actionButton("check_button", "Check Expression", class = "btn-primary")  # Use a primary button style
    ),
    
    mainPanel(
      h4("Gene Expression Results:"),
      withSpinner(tableOutput("expression_results"))  # Add loading spinner
    )
  )
)

# Define the Server logic
server <- function(input, output, session) {
  # Read all expression data files from 'preexisting_data' folder
  expression_data <- read.csv("C:/Users/17735/Downloads/Expression_app/preexisting_data/Expression_RPKM.tsv", sep = "")
  
  # Dynamically update selectize input with gene names as options
  observe({
    genes <- expression_data[[1]]  # Assume first column contains gene names
    updateSelectizeInput(session, "selected_genes", choices = genes, server = TRUE)
  })
  
  # Reactive expression to calculate expression summaries for selected genes
  expression_summary <- reactive({
    req(input$selected_genes)  # Ensure at least one gene is selected
    
    filtered_data <- expression_data %>%
      filter(ID %in% input$selected_genes) %>%
      mutate(across(where(is.numeric), ~ round(.x, 2)))
    
    filtered_data$Top_Donor <- apply(filtered_data[-1], 1, function(x) {
      max_value <- max(x, na.rm = TRUE)  
      names(x)[x == max_value]
    })
    
    filtered_data$Top_Donor <- sapply(filtered_data$Top_Donor, function(names) {
      if (length(names) == 0) NA else paste(names, collapse = ", ")
    })
    
    results <- filtered_data %>%
      rowwise() %>%
      mutate(
        Top_Expression_Value = max(c_across(where(is.numeric)), na.rm = TRUE),
        Total_Expression = sum(c_across(where(is.numeric))),
        Expressed = ifelse(Total_Expression > 0, "Yes", "No"),
        Top_Donor = ifelse(Expressed == 'No', '', Top_Donor)
      ) %>%
      select(Gene = ID, Expressed, Top_Donor, Top_Expression_Value)
    
    colnames(results) <- c('Gene', 'Expressed', 'Top Donor', 'Top Expression Value')
    
    return(results)
  })
  
  # Render the expression results table when the button is clicked
  output$expression_results <- renderTable({
    input$check_button  # Trigger table render on button click
    isolate(expression_summary())
  })
}

# Run the app
shinyApp(ui = ui, server = server)



