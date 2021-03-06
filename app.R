# Maximal growth rate estimator - App version
# Author: Catharina Meyer, Quantitative Pharmacology
# Contact: c.meyer@lacdr.leidenuniv.nl
# Description: This shiny app fits time kill data in two steps. In step 1, 
#              the growth curves over time are fitted and the maximal growth rates 
#              are estimated. In step 2, the PD relation (maximal growth rates over 
#              drug concentration) is fitted.
# Date of last change: 22.04.2021

#--------------------------------------------------------------------------------------------------------------------------------------


# load R packages
library(shiny)
library(lattice)
library(deSolve)
library("growthrates") # functions from this package are used for the main fitting in step 1

#--------------------------------------------------------------------------------------------------------------------------------------
# ------------------- User interface object--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------

# This object defines the user input, side panels and output panels and is only necessary when run as app
ui <- fluidPage(
  
  titlePanel("Maximal growth rate estimator"),
  
  # Sidebar layout
  sidebarLayout(
    
    # Sidebar panel for inputs
    sidebarPanel(
      
      # Input: Select a file
      fileInput("file1", "Choose CSV File",
                multiple = FALSE,
                accept = c(".csv")),

      radioButtons(inputId = "all", label="Process the whole plate?", choices = c("yes","no"),selected="no"),
      conditionalPanel(
        condition = "input.all != 'yes'",
        selectInput('selected_strain', label = 'Select a strain', choices = 'No choices here yet'),
        selectInput('selected_drug', label = 'Select a drug', choices = 'No choices here yet'),
        selectInput('selected_media', label = 'Select a media', choices = 'No choices here yet')
      ),

      selectInput('selected_measurement', label = 'Select measurement data', choices = c('raw_measurement',"corrected_measurement"),selected="corrected_measurement"),
      selectInput('selected_time', label = 'Exclude time points before:', choices = 'No choices here yet'),
      
      helpText("This shiny app fits time kill data in two steps. In step 1, the growth curves over time are fitted and the maximal growth rates are estimated. In step 2, the relation of maximal growth rates over drug concentration is fitted."),
      
      # Input: actionButton() to defer the rendering of output
      # until the user explicitly clicks the button (rather than
      # doing it immediately when inputs change).
      fluidRow(
        column(3,
               actionButton("update_1", HTML("Step 1:<br/>Compute<br/>Growth Rates"))),
        column(3,
               actionButton("update_2", HTML("Step 2:<br/>Compute<br/>PD relation ")))),
      
      helpText("Note: Restart estimation when changing input parameters. If changing step 2 methods or parameters, it is only necessary to rerun step 2."),
      
      radioButtons(inputId = "pool", label="Compute growth rate", choices = c("per replicate","pooled"),selected="pooled"),

      # Input: Selector for choosing method
      selectInput(inputId = "method",
                  label = "Choose a fitting method for step 1:",
                  choices = c("smoothing", "baranyi", "huang")),
      
      # Numeric input of parameters step 1
      # show these parameters only when choosing a parametric method (huang or baranyi)
      conditionalPanel(
        condition = "input.method != 'smoothing'",
        radioButtons(inputId = "customrange", label="Customize lower and upper values:", choices = c("yes","no"),selected="no"),
        helpText("Choose parameter estimates for step 1:"),
        fluidRow(
          column(3,
                 numericInput("y0", "y0:",
                              min = -5, max = 5,
                              value = 0.1, step = 0.01)),
          column(3,
                 numericInput("mumax", "Maximal growth rate (mumax):",
                              min = 0, max = 5,
                              value = 0.1, step = 0.01)),
          column(3,
                 numericInput("K", "Capacity (K):",
                              min = 0, max = 5,
                              value = 0.2, step = 0.01))
          ),
        # only show range of initial parameters if the customrange is set to yes
        conditionalPanel(
          condition = "input.customrange == 'yes'",
          fluidRow(
            column(3,
                   numericInput("y0_lower", "Lower y0:",
                                min = -5, max = 5,
                                value = 0.01, step = 0.01)),
            column(3,
                   numericInput("mumax_lower", "Lower mumax:",
                                min = 0, max = 5,
                                value = 0.01, step = 0.01)),
            column(3,
                   numericInput("K_lower", "Lower K:",
                                min = 0, max = 5,
                                value = 0.01, step = 0.01))
            # )
          ),
          fluidRow(
            column(3,
                   numericInput("y0_upper", "Upper y0:",
                                min = -5, max = 5,
                                value = 0.6, step = 0.01)),
            column(3,
                   numericInput("mumax_upper", "Upper mumax:",
                                min = 0, max = 5,
                                value = 5, step = 0.01)),
            column(3,
                   numericInput("K_upper", "Upper K:",
                                min = 0, max = 5,
                                value = 2, step = 0.01)))
          )
        ),
      # if the method baranyi is chosen, additionally show these parameters and the ranges
      conditionalPanel(
        condition = "input.method == 'baranyi'",
        fluidRow(
          column(3,
                 numericInput("h0", "h0:",
                              min = 0, max = 10,
                              value = 1, step = 0.1))),
        conditionalPanel(
          condition = "input.customrange == 'yes'",
        fluidRow(
          column(3,
                 numericInput("h0_lower", "Lower h0:",
                              min = 0, max = 10,
                              value = 0, step = 0.1)),
          column(3,
                 numericInput("h0_upper", "Upper h0:",
                              min = 0, max = 10,
                              value = 10, step = 0.1)))
        )
      ),
      # if the method huang is chosen, additionally show these parameters and the ranges
      
      conditionalPanel(
        condition = "input.method == 'huang'",
        fluidRow(
          column(3,
                 numericInput("alpha", "Alpha:",
                              min = 0, max = 10,
                              value = 1, step = 0.1),
                 numericInput("lambda", "Lambda:",
                              min = 0, max = 5,
                              value = 0.1, step = 0.1))),
        conditionalPanel(
          condition = "input.customrange == 'yes'",
        fluidRow(
          column(3,
                 numericInput("alpha_lower", "Lower alpha:",
                              min = 0, max = 10,
                              value = 0, step = 0.1),
                 numericInput("lambda_lower", "Lower lambda:",
                              min = 0, max = 5,
                              value = 0, step = 0.1)),
          column(3,
                 numericInput("alpha_upper", "Upper alpha:",
                              min = 0, max = 10,
                              value = 7, step = 0.1),
                 numericInput("lambda_upper", "Upper lambda:",
                              min = 0, max = 5,
                              value = 4, step = 0.1)))
        )
      ),
      
      # Select input for method in step 2
      radioButtons(inputId = "model", label="Choose a model for step 2:", choices = c("exp","Emax","sigmoid_Emax", "capacity_Emax"), selected="sigmoid_Emax"),
      
      # for the capacity_Emax model, a shared fitting option is implemented
      conditionalPanel(
        condition = "input.model == 'capacity_Emax'",
        radioButtons(inputId = "fitting_type", label="Choose fitting type for parameters:", choices = c("individual","shared"), selected = "individual"),
        ),
      # Show input for parameters of exponential model if this is selected
      conditionalPanel(
        condition = "input.model == 'exp'",
        helpText("Exponential model equation: mumax[drug_concentration] = a*exp(b*drug_concentration)+c"),
        fluidRow(
          column(3,
                 numericInput("a", "a :",
                              min = -5, max = 5,
                              value = -0.4, step = 0.1)),
          column(3,
                 numericInput("b", "b :",
                              min = -5, max = 5,
                              value = -5, step = 0.1)),
          column(3,
                 numericInput("c", "c :",
                              min = -5, max = 5,
                              value = 0.6, step = 0.1)
          ))),
      # an automated guess of the inital parameters is available in step 2 for the variations of the Emax model
      conditionalPanel(
        condition = ("(input.model == 'Emax') | (input.model == 'sigmoid_Emax') | ((input.model == 'capacity_Emax') & (input.fitting_type=='individual'))"),
        radioButtons(inputId = "guess", label="Guess step 2 parameters:", choices = c("yes","no")),
        # show manual input option for inital parameters
        conditionalPanel(
          condition = "input.guess == 'no'",
          # Numeric input of parameters step 2
          helpText("Choose parameter estimates for step 2: "),
          helpText("(Sigmoid) E_max model formula: mumax ~ E0 + E_max * (((drug_concentration/EC50)^k)/(1 + ((drug_concentration/EC50)^k)))"),
          helpText("Note: for Emax, k is fixed to 1"),
        
          fluidRow(
            column(3,
                   numericInput("E0", "Effect without drug: E0",
                                min = -5, max = 5,
                                value = 0.4, step = 0.1)
            ),
            column(3,
                   numericInput("E_max", "Maximum effect: E_max",
                                min = -5, max = 5,
                                value = -0.4, step = 0.1)
            ),
            column(3,
                   numericInput("EC50", "Drug concentration: EC50",
                                min = -5, max = 5,
                                value = 0.2, step = 0.1)
            ),
            column(3,
                   numericInput("k", "Hill parameter: k",
                                min = -5, max = 10,
                                value = 5, step = 0.1))
          )
          
        )
      ),

      # download output file buttons
      downloadButton("downloadData_step1", "Download results step 1"),
      downloadButton("downloadData_step2", "Download results step 2")
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      
      # Output: Tabset
      tabsetPanel(type = "tabs",
                  tabPanel("Input Data", tableOutput("data")),
                  tabPanel(HTML("Step 1:<br/>Growth Rates<br/>Plots"), plotOutput("plots")),
                  tabPanel(HTML("Step 1:<br/>Growth Rates<br/>Results"), tableOutput("results")),
                  tabPanel(HTML("Step 2:<br/>PD relation<br/>Initial parameters"), plotOutput("plot1")),
                  tabPanel(HTML("Step 2:<br/>PD relation<br/>Plots"), plotOutput("Gconc")),
                  tabPanel(HTML("Step 2:<br/>PD relation<br/>Results"), tableOutput("results_PD")),
                  tabPanel(HTML("Step 2:<br/>PD relation<br/>Fit Summary"), verbatimTextOutput("fit_step2"))
                  )
      )
    )
  )


#--------------------------------------------------------------------------------------------------------------------------------------
#Server function ----------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------
# This function creates the output

server <- function(input, output, session) {
  
  source("Functions.R") # this imports the functions defined in the script Functions.R

#____________________event and input functions_____________________________________________________________________________________________
  
  # creates the input selection options from the input file
  observeEvent(input$file1, {
    df_all <- getData_all_R()
    updateSelectInput(session, "selected_strain", label = "Select a strain", choices = levels(df_all$strain_name))
    updateSelectInput(session, "selected_drug", label = "Select a drug", choices = levels(df_all$drug_name))
    updateSelectInput(session, "selected_media", label = "Select a media", choices = levels(df_all$media_name))
    updateSelectInput(session, "selected_time", label = 'Exclude time points before:', choices = levels(as.factor((df_all$time))))
    
  })
  
#_____________________________reactive wrapper functions_(everything input or app specific goes here)___________________________________________________________________________________
  # all the functions with "_R" are reactive functions to wrap the imported functions. 
  # This makes sure that the wrapped functions are only run when certain inputs are changed or when the action buttons are pressed.
  
  # reads in all the data from the input file
  getData_all_R <- reactive({
    req(input$file1)
    df_all <- getData_all(input$file1$datapath)
    return(df_all)
  })
  # selects a subset of the data according to user input
  getData_subset_R <- reactive({
    df_all <- getData_all_R()
    df_subset <- getData_subset(input, df_all)
  return(df_subset)
    })

  # fits the growthcurves to determine maximal growth rates (step 1)
  Fit_Growthrate_R <- eventReactive(input$update_1,{
    df <- getData_subset_R()
    many_fits <- Fit_Growthrate(input, df)
    return(many_fits)
  })
  
  # created all combintaions of drug/strain/media that is present in the data subset
  combinations_R <- reactive({
    df <- getData_subset_R()
    comb <- combinations(input, df)
    return(comb)
  }) 

  # fits the PD relation (maximal growth rates over drug concentration) for every combination of drug/strain/media individually or with shared parameters
  Fit_PD_R <- eventReactive(input$update_2,{
    comb <- combinations_R()
    res <- Results()
    m <- Fit_PD(input, comb, res)
    return(m)
  })

  # retrieves results from the PD fitting (step 2) in the correct data format
  results_PD_R <- reactive({
    comb <- combinations_R()
    many_fits <- Fit_Growthrate_R()
    # res <- results(many_fits) # replace with Results()???
    res <- Results()
    m <- Fit_PD_R()
    d <- results_PD(input, comb, res, m)
    return(d)
  })
  
  # retrieves results from the growthcurves fitting (step 1)
  Results <- function() {
    res <- results(Fit_Growthrate_R())
    return(res)
  }

#----Output--------------------------------------------------------------------------------------------------------------
  
  # Show input data in HTML table
  output$data <- renderTable({
    df_all <- getData_all_R()
    return(df_all)
    },digits=6)
  
  # show fitting results (step 1) in table
  output$results <- renderTable(Results(), digits=6)
 
  # plot with fitting
  output$plots <- renderPlot({
    many_fits <- Fit_Growthrate_R()
    comb <- combinations_R()
    plot_Growthrate(input, many_fits, comb)

  })
  
  # plot maximal growth rate estimates over drug concentration 
  # the model fit using the initial parameters before the fitting is shown in red
  output$plot1 <- renderPlot({
    res <- Results()
    comb <- combinations_R()
    plot_PD_init(input, res, comb)
  })  
  
  
  # plot maximal growth rate estimate over AB
  output$Gconc <- renderPlot({
    res <- Results()
    m <- Fit_PD_R()
    comb <- combinations_R()
    plot_PD(input, res, m, comb)
  })
  
  # show fitting results (step 2) in table
  output$results_PD <- renderTable({
    d <- results_PD_R()
    },digits=6)
    
  # print summary of fitting in step 2
  output$fit_step2 <- renderPrint({
    m <- Fit_PD_R()
    comb <- combinations_R()
    Print_PD(input, comb, m)
  })
  
  # defines output that can be downloaded
  output$downloadData_step1 <- downloadHandler(
    filename = function() {
      paste("Step1_results-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(Results(), file)
    }
  )
  
  # defines output that can be downloaded
  output$downloadData_step2 <- downloadHandler(
    filename = function() {
      paste("Step2_results-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(results_PD_R(), file)
    }
  )
}

# Create Shiny app
shinyApp(ui, server)
