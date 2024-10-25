# Load libraries
library(shiny)
library(dplyr)
library(shinycssloaders)  # For loading animations
library(shinythemes)      # For theme support




aws_prefix <- '/mnt/efs/fs1/destination_folder/Expression_app/'
#aws_prefix <- 'C:/Users/17735/Downloads/Expression_app/'



# Define the UI
ui <- fluidPage(
  theme = shinytheme("flatly"),
  
  # Include Google Fonts and CSS styling for Montserrat
  tags$head(
    tags$link(
      href = "https://fonts.googleapis.com/css2?family=Montserrat:wght@300;400;700&display=swap", 
      rel = "stylesheet"
    ),
    tags$style(HTML("
      body { font-family: 'Montserrat', sans-serif; }
      h1, h2, h3, h4, h5, h6 { font-family: 'Montserrat', sans-serif; font-weight: 700; }
      .shiny-input-container { margin-bottom: 20px; }
      h4 { color: #5c9fd7; font-weight: 700; }
      .btn-primary { background-color: #5c9fd7; border-color: #5c9fd7; }
      .btn-primary:hover { background-color: #4a8cbb; border-color: #4a8cbb; }
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
      h4(""),
      selectizeInput(
        "selected_genes", 
        tags$span(style = "color: #5c9fd7;", "Select up to 10 genes:"),
        choices = NULL,  # Initially empty, populated dynamically
        multiple = TRUE, 
        options = list(
          maxItems = 10,
          placeholder = 'Start typing to search for genes',
          create = FALSE
        )
      ),
      actionButton("check_button", "Check Expression", class = "btn-primary")
    ),
    
    mainPanel(
      # Conditionally display the Results header
      conditionalPanel(
        condition = "output.showResults === true",
        h4("Results:")
      ),
      # Render the table and spinner only after the button press
      uiOutput("table_ui")
    )
  )
)

# Define the Server logic
server <- function(input, output, session) {
  # Read expression data
  expression_data <- read.csv(paste0(aws_prefix, "preexisting_data/Expression_RPKM.tsv"), sep = "")
  
  # Dynamically update selectize input with gene names as options
  observe({
    genes <- expression_data[[1]]  # Assume first column contains gene names
    updateSelectizeInput(session, "selected_genes", choices = genes, server = TRUE)
  })
  
  # Store the state of the button press
  button_pressed <- reactiveVal(FALSE)
  
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
  
  # Update the button press state and render the UI dynamically
  observeEvent(input$check_button, {
    button_pressed(TRUE)  # Set button as pressed
    
    output$table_ui <- renderUI({
      req(button_pressed())  # Ensure the button has been pressed
      withSpinner(tableOutput("expression_results"))
    })
    
    output$expression_results <- renderTable({
      expression_summary()
    })
  })
  
  # Check if the results should be shown
  output$showResults <- reactive({
    button_pressed() && nrow(expression_summary()) > 0
  })
  outputOptions(output, "showResults", suspendWhenHidden = FALSE)  # Ensure it updates immediately
}

# Run the app
shinyApp(ui = ui, server = server)

