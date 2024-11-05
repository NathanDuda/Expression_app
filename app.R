# Load libraries
library(shiny)
library(dplyr)
library(shinycssloaders)  # For loading animations
library(shinythemes)      # For theme support
library(DT)               # For paginated and searchable tables
library(tidyverse)
library(sysfonts)
library(showtext)
library(ggplot2)

# Define file path
#aws_prefix <- 'C:/Users/17735/Downloads/Expression_app/'
aws_prefix <- '/mnt/efs/fs1/destination_folder/Expression_app/'

font_add_google(name = "Montserrat", family = 'Montserrat') 
showtext_auto()


# Define the UI
ui <- fluidPage(
  theme = shinytheme("flatly"),
  
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
      #results-table { margin-top: 20px; }
      
      /* CSS for vertical column lines */
      table.dataTable {
        border-collapse: collapse; 
        width: 100%; 
      }
      table.dataTable th, table.dataTable td {
        border-left: 1px solid #ddd; 
      }
      table.dataTable th:first-child, table.dataTable td:first-child {
        border-left: none; 
      }
    "))
    
    
    
    
    
  ),
  
  titlePanel("Gene Expression Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Search Genes or Cell Lines"),
      
      selectizeInput(
        "selected_genes_temp", 
        "Select genes:",
        choices = NULL, 
        multiple = TRUE, 
        options = list(
          maxItems = 100,
          placeholder = 'Start typing to search for genes',
          create = FALSE
        )
      ),
      
      selectInput(
        "column_select_temp", 
        "Select a cell line:", 
        choices = NULL, 
        selected = "none"  # Default to no selection
      ),
      
      actionButton("check_button", "Check Expression", class = "btn-primary")
    ),
    
    # Inside Main Panel UI
    mainPanel(
      conditionalPanel(
        condition = "output.showResults === true",
        h4("Results:")
      ),
      
      # Use a fluidRow to position stats and plot side by side
      uiOutput("column_stats_ui"),  # UI for column stats and plot
      uiOutput("combined_ui"),
      # Add a results table below the column stats and plot
      uiOutput("table_ui")  # Main table UI
    )
    
  )
)



# Define the Server logic
server <- function(input, output, session) {
  # Read expression data
  expression_data <- read.csv(paste0(aws_prefix, "preexisting_data/Expression_RPKM.tsv"), sep = "")
  
  valid_genes_data <- expression_data %>%
    #filter(rowSums(select(., -ID)) > 0) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>% # Round values to 2 decimal places
    group_by(ID) %>% # remove all copies of genes with IDs that appear more than once (do not keep any duplicated gene id)
    filter(n() == 1) %>%
    ungroup()
  
  
  top_donor_stats <- valid_genes_data %>%
    pivot_longer(cols = colnames(expression_data)[-1], names_to = 'cell_line', values_to = 'exp') %>%
    group_by(cell_line) %>%
    summarise(
      median = median(exp, na.rm = TRUE),
      Q25 = round(quantile(exp, 0.25), 2),
      Q75 = round(quantile(exp, 0.75), 2),
      Q90 = round(quantile(exp, 0.90), 2),
      max = max(exp, na.rm = TRUE)
    )
  colnames(top_donor_stats) <- c('Top_Donor', 'Median', 'Q25', 'Q75', 'Q90', 'Max')
  selected_genes <- reactiveVal()
  selected_column <- reactiveVal()
  
  # Populate the gene selector with valid genes only
  observe({
    valid_genes <- valid_genes_data[[1]]  
    updateSelectizeInput(session, "selected_genes_temp", choices = valid_genes, server = TRUE)
  })
  
  # Populate the column selector
  observe({
    column_choices <- c("none", names(valid_genes_data)[-1])
    updateSelectInput(session, "column_select_temp", choices = column_choices, selected = "none")
  })
  
  
  
  
  # Handle the "Check Expression" button click
  observeEvent(input$check_button, {
    selected_genes(input$selected_genes_temp)
    selected_column(input$column_select_temp)
    
    output$table_ui <- renderUI({
      if (!is.null(selected_genes()) && length(selected_genes()) > 0 && !is.null(selected_column()) && selected_column() != "none") {
        tagList(
          h4(paste0("Expression of selected genes in ", selected_column())),  # Title for combined results table
          withSpinner(DTOutput("combined_results"))
        )
      } else if (!is.null(selected_genes()) && length(selected_genes()) > 0) {
        tagList(
          h4("Expression Results"),  # Title for main results table
          withSpinner(DTOutput("expression_results"))
        )
      } else {
        tagList(
          h4("All Genes"),  # Title for column-based table
          withSpinner(DTOutput("expression_results"))
        )
      }
    })
    
    
    output$combined_ui <- renderUI({
      if (!is.null(selected_column()) && selected_column() != "none") {
        fluidRow(
          column(
            width = 6,
            tagList(
              h4(paste0("Summary Statistics of ", selected_column())),  # Title for stats table
              withSpinner(DTOutput("column_stats_table"))
            )
          ),
          column(
            width = 6,
            tagList(
              h4(paste0("Expression Distribution of ", selected_column())),  # Title for distribution plot
              plotOutput("column_distribution_plot", height = "250px")
            )
          )
        )
      } else {
        NULL  # Return NULL to avoid rendering empty space when only genes are selected
      }
    })
    output$combined_results <- renderDT({
      if (!is.null(selected_genes()) && length(selected_genes()) > 0 && !is.null(selected_column()) && selected_column() != "none") {
        filtered_genes <- expression_summary_column() %>%
          filter(Gene %in% selected_genes())
        
        # Include the selected column in the results
        combined_data <- filtered_genes# %>%
        #  mutate(Selected_Column_Value = valid_genes_data[[selected_column()]][match(filtered_genes$Gene, valid_genes_data$ID)]) %>%
        # select(Gene, Expressed, Top_Donor, Top_Expression_Value, Selected_Column_Value)
        
        datatable(
          combined_data,
          options = list(
            pageLength = 10,
            dom = 't',  # Controls the elements in the table
            autoWidth = TRUE
          ),
          rownames = FALSE
        )
      }
    })
    
    output$expression_results <- renderDT({
      if (!is.null(selected_genes()) && length(selected_genes()) > 0) {
        filtered_data <- expression_summary_genes() %>%
          filter(Gene %in% selected_genes())
        
        datatable(
          filtered_data,
          options = list(
            pageLength = 10,
            dom = 't',  # Controls the elements in the table
            autoWidth = TRUE
          ),
          rownames = FALSE
        )
      } else {
        datatable(expression_summary_column(), options = list(pageLength = 10, searchable = TRUE),
                  rownames = F)
      }
    })
    # Conditionally render the column statistics table and plot
    output$column_stats_ui <- renderUI({
      if (!is.null(selected_genes()) && length(selected_genes()) > 0) {
        NULL  # Hide if genes are selected
      } else if (!is.null(selected_column()) && selected_column() != "") {
        NULL
      } else {
        NULL
      }
    })
    
    output$column_stats_table <- renderDT({
      if (!is.null(selected_column()) && selected_column() != "none") {
        req(selected_column())  # Ensure a column is selected
        datatable(column_statistics_summary(), options = list(pageLength = 5,
                                                              searching = F,
                                                              paging = F,
                                                              info = F),
                  rownames = F)
      }
    })
    
    
    
  })
  
  output$column_distribution_plot <- renderPlot({
    if (!is.null(selected_column()) && selected_column() != "none") {
      req(selected_column())  # Ensure a column is selected
      
      selected_col <- selected_column()
      values <- valid_genes_data[[selected_col]]  # Extract the column values
      
      # Create a data frame with sorted values and their ranks
      plot_data <- data.frame(
        Value = values  # Sorted expression values
      )
      
      # Calculate statistics
      median_val <- median(values, na.rm = TRUE)
      max_val <- max(values, na.rm = TRUE)
      lower_bound <- quantile(values, 0.25, na.rm = TRUE)  # 25th percentile
      upper_bound <- quantile(values, 0.75, na.rm = TRUE)  # 75th percentile
      
      # Create a line plot with a log-scaled x-axis using ggplot2
      ggplot(plot_data, aes(x = Value)) +
        geom_density(aes(y = ..density..), fill = "#8dbcE3", color = "#5c9fd7", size = 1, alpha = 0.6) +  # Fill with a lighter shade
        xlim(0, upper_bound*2) +
        geom_vline(aes(xintercept = median_val, linetype = "Median"), color = "gray20") +  # Median line
        geom_vline(aes(xintercept = lower_bound, linetype = "IQR"), color = "gray20") +  # 25th percentile line
        geom_vline(aes(xintercept = upper_bound, linetype = "IQR"), color = "gray20") +  # 75th percentile line
        scale_linetype_manual(
          name = "",
          values = c("Median" = "dashed", "IQR" = "dotted"),
          labels = c("IQR (25th and 75th Percentile)", "Median")
        ) +
        theme_minimal() +  # Apply a clean theme
        theme(
          legend.position = c(0.99, 0.99),  # Position legend in the upper right corner
          legend.justification = c("right", "top"),  # Align the legend to the top right
          text = element_text(family = "Montserrat")) +
        labs(
          title = '',
          x = "Expression Value", 
          y = "Density"
        ) 
    }
  })
  
  # Reactive summary logic (genes and column-based)
  expression_summary_genes <- reactive({
    if (!is.null(selected_genes()) && length(selected_genes()) > 0) {
      filtered_data <- valid_genes_data %>% 
        filter(ID %in% selected_genes())
      
      filtered_data$Top_Donor <- apply(filtered_data[-1], 1, function(x) names(x)[x == max(x, na.rm = TRUE)])
      
      out2 <- filtered_data %>%
        rowwise() %>%
        mutate(
          Top_Expression_Value = max(c_across(where(is.numeric)), na.rm = TRUE),
          Total_Expression = sum(c_across(where(is.numeric))),
          Expressed = ifelse(Total_Expression > 0, "Yes", "No"),
          Top_Donor = ifelse(Expressed == "No", "", Top_Donor)
        ) %>%
        merge(top_donor_stats, by = 'Top_Donor') %>%
        select(Gene = ID, Expressed, Top_Donor, Top_Expression_Value, 
               Q25, Median, Q75, Q90, Max)
      colnames(out2) <- c('Gene', 'Expressed', 'Top donor', 'Expression value', 
                          '25th Percentile', 'Median', '75th Percentile', '90th Percentile', 'Max')
      out2
    }
  })
  
  # Reactive expression for column-based summaries with percentile and statistics
  expression_summary_column <- reactive({
    if (!is.null(selected_column()) && selected_column() != "none") {
      selected_col <- selected_column()
      
      # Calculate column statistics
      filtered_data <- valid_genes_data %>%
        mutate(
          Value = get(selected_col),
          Expressed = ifelse(Value > 0, "Yes", "No"),
          Percentile = round(percent_rank(Value) * 100, 2)  # Calculate percentile
        ) %>%
        arrange(desc(Value))  # Sort descending by Value for clarity
      
      # get stats values for the selected column 
      filtered_data$Median <- top_donor_stats$Median[top_donor_stats$Top_Donor == selected_column()]
      filtered_data$Q25 <- top_donor_stats$Q25[top_donor_stats$Top_Donor == selected_column()]
      filtered_data$Q75 <- top_donor_stats$Q75[top_donor_stats$Top_Donor == selected_column()]
      filtered_data$Q90 <- top_donor_stats$Q90[top_donor_stats$Top_Donor == selected_column()]
      filtered_data$Max <- top_donor_stats$Max[top_donor_stats$Top_Donor == selected_column()]
      
      # Select columns to display
      out <- filtered_data %>%
        arrange(ID) %>%
        select(Gene = ID, Expressed, Value, Percentile)
      colnames(out)[3] <- 'Expression Value'
      out
    }
  })
  
  column_statistics_summary <- reactive({
    if (!is.null(selected_column()) && selected_column() != "none") {
      selected_col <- selected_column()
      
      # Filter and summarize statistics for the selected column
      stats <- top_donor_stats %>%
        filter(Top_Donor == selected_col) %>%
        select(Q25, Median, Q75, Q90, Max)
      
      # Create a simple summary table with labels
      data.frame(
        Statistic = c('25th Percentile', "Median", '75th Percentile', '90th Percentile', "Max"),
        Value = c(stats$Q25, stats$Median, stats$Q75, stats$Q90, stats$Max)
      )
    }
  })
  
  # Control when to show the results section
  output$showResults <- reactive({
    !is.null(selected_genes()) || !is.null(selected_column())
  })
  outputOptions(output, "showResults", suspendWhenHidden = FALSE)
}

# Run the app
shinyApp(ui = ui, server = server)

