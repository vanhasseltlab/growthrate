library(shiny)
library(lattice)
library(deSolve)
library("growthrates")
# library(ggplot2)

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
      
      helpText("This shiny app fits time kill data in two steps. In step 1, the growth curves over time are fitted and the maximal growth rates are estimated. In step 2, the relation of maximal growth rates over drug concentration is fitted."),
      
      # Input: actionButton() to defer the rendering of output ----
      # until the user explicitly clicks the button (rather than
      # doing it immediately when inputs change).
      fluidRow(
        column(3,
               actionButton("update_1", HTML("Step 1:<br/>Compute<br/>Growth Rates"))),
        column(3,
               actionButton("update_2", HTML("Step 2:<br/>Compute<br/>PD relation ")))),
      
      helpText("Note: Restart estimation when changing input parameters. If changing step 2 methods or parameters, it is only necessary to rerun step 2."),
      
      radioButtons(inputId = "pool", label="Compute growth rate", choices = c("per replicate","pooled"),selected="pooled"),
      #helpText("Note: Step 2 is only possible when data with different drug concentrations is not pooled together!"),

      # Input: Selector for choosing method
      selectInput(inputId = "method",
                  label = "Choose a fitting method for step 1:",
                  choices = c("smoothing", "baranyi", "huang")),
      
      # Numeric input of parameters step 1
      # radioButtons(inputId = "guess1", label="Guess initials step 1:", choices = c("yes","no"),selected="yes"),
      
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
        # #this doesn't work
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

      radioButtons(inputId = "model", label="Choose a model for step 2:", choices = c("exp","Emax","sigmoid_Emax", "capacity_Emax"), selected="sigmoid_Emax"),
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
      
      conditionalPanel(
        condition = "input.model != 'exp'",
        radioButtons(inputId = "guess", label="Guess step 2 parameters:", choices = c("yes","no")),
        conditionalPanel(
          condition = "input.guess == 'no'",
          # Numeric input of parameters step 2
          helpText("Choose parameter estimates for step 2: "),
          # helpText("(Sigmoid) E_max model equation: mumax[drug_concentration] = E0 + E_max *(((drug_concentration/EC50)**k)/(1+((drug_concentration/EC50)**k)))"),
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
                                min = -5, max = 5,
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
                  tabPanel(HTML("Step 2:<br/>PD relation<br/>Plots"), plotOutput("Gconc")),
                  tabPanel(HTML("Step 2:<br/>PD relation<br/>Results"), tableOutput("results_step2")),
                  tabPanel(HTML("Step 2:<br/>PD relation<br/>Fit Summary"), verbatimTextOutput("fit_step2"))
                  )
      )
    )
  )


#--------------------------------------------------------------------------------------------------------------------------------------
#Server function ----------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------
server <- function(input, output,session) {
  
  observeEvent(input$file1, {
    mytable <- read.csv(input$file1$datapath, stringsAsFactors=TRUE)
    updateSelectInput(session, "selected_strain", label = "Select a strain", choices = levels(mytable$strain_name))
    updateSelectInput(session, "selected_drug", label = "Select a drug", choices = levels(mytable$drug_name))
    updateSelectInput(session, "selected_media", label = "Select a media", choices = levels(mytable$media_name))
  })
  
  # input values step 1 fitting
  #Reactive expression to create data frame of step 1 input values, only updated when action button is pressed
  sliderValues <- eventReactive(input$update_1,{
    data.frame(
      Name = c("Y0", "Mumax","K","h0","Alpha","Lambda"),
      InitialEstimate = c(input$y0,
                          input$mumax,
                          input$K,
                          input$h0,
                          input$alpha,
                          input$lambda),
      LowerEst = c(input$y0_lower,
                   input$mumax_lower,
                   input$K_lower,
                   input$h0_lower,
                   input$alpha_lower,
                   input$lambda_lower),
      UpperEst = c(input$y0_upper,
                   input$mumax_upper,
                   input$K_upper,
                   input$h0_upper,
                   input$alpha_upper,
                   input$lambda_upper),
      stringsAsFactors = FALSE)
  })
  
  # input values for step 2 exp fitting
  sliderValues2 <- eventReactive(input$update_2,{
    data.frame(
      Name = c("a", "b","c"),
      Est_s = c(input$a,
              input$b,
              input$c),
      stringsAsFactors = FALSE)
  })
  
  # input values for step 2 Emax fitting
  sliderValues_Emax <- eventReactive(input$update_2,{
    data.frame(
      Name = c("E0", "E_max","EC50","k"),
      Est = c(input$E0,
                input$E_max,
                input$EC50,
                input$k),
      stringsAsFactors = FALSE)
  })

getData <- reactive({
  req(input$file1)
  df <- read.csv(input$file1$datapath, stringsAsFactors=TRUE)
  if (input$all=="no") {
    df_subset <- subset(df, (df$strain_name==input$selected_strain) & (df$drug_name==input$selected_drug) & (df$media_name==input$selected_media))
  }
  else if (input$all=="yes") {
    df_subset <- df
  }
  df_subset$corrected_measurement <- as.numeric(as.character(df_subset$corrected_measurement))
  return(df_subset)
})

  # step 1 fitting/ not modularized yet
  Fitting <- eventReactive(input$update_1,{
    df <- getData()
    
    # select subset of data according to user input
    # changes column name df, sorry about that
    if (input$selected_measurement=="raw_measurement"){
      colnames(df)[which((colnames(df)=="raw_measurement"))] <- "selected_measurement"
    } else if (input$selected_measurement=="corrected_measurement") {
      colnames(df)[which((colnames(df)=="corrected_measurement"))] <- "selected_measurement"
    }
    
    if (input$pool=="pooled") {
      if (input$method == "smoothing") {
        many_fits <- all_splines(selected_measurement ~ time | strain_name + drug_concentration + drug_name + media_name ,
                                 data = df, spar = 0.5)
      }
      
      else if (input$method == "baranyi") {
        p   <- c(y0 = sliderValues()$InitialEstimate[1], mumax = sliderValues()$InitialEstimate[2], K = sliderValues()$InitialEstimate[3], h0 = sliderValues()$InitialEstimate[4])
        lower   <- c(y0 = sliderValues()$LowerEst[1], mumax = sliderValues()$LowerEst[2], K = sliderValues()$LowerEst[3], h0 = sliderValues()$LowerEst[4])
        upper   <- c(y0 = sliderValues()$UpperEst[1], mumax = sliderValues()$UpperEst[2], K = sliderValues()$UpperEst[3], h0 = sliderValues()$UpperEst[4])
        
        many_fits <- all_growthmodels(
          selected_measurement ~ grow_baranyi(time, parms) | strain_name + drug_concentration + drug_name + media_name,
          data = df,
          p = p, lower = lower, upper = upper,
          transform = "log", ncores = 2)
      }
      
      else if (input$method == "huang") {
        p   <- c(y0 = sliderValues()$InitialEstimate[1], mumax = sliderValues()$InitialEstimate[2], K = sliderValues()$InitialEstimate[3], alpha = sliderValues()$InitialEstimate[5],lambda = sliderValues()$InitialEstimate[6] )
        lower   <- c(y0 = sliderValues()$LowerEst[1], mumax = sliderValues()$LowerEst[2], K = sliderValues()$LowerEst[3], alpha = sliderValues()$LowerEst[5],lambda = sliderValues()$LowerEst[6])
        upper   <- c(y0 = sliderValues()$UpperEst[1], mumax = sliderValues()$UpperEst[2], K = sliderValues()$UpperEst[3], alpha = sliderValues()$UpperEst[5],lambda = sliderValues()$UpperEst[6])
        
        many_fits <- all_growthmodels(
          selected_measurement ~ grow_huang(time, parms) |  strain_name + drug_concentration + drug_name + media_name,
          data = df,
          p = p, lower = lower, upper = upper,
          transform = "log", ncores = 2)
      }
    }
    
    else if (input$pool=="per replicate"){
      if (input$method == "smoothing") {
        # Smooting for mutliple data sets 
        many_fits <- all_splines(selected_measurement ~ time | strain_name + drug_concentration + replicate_ID + drug_name + media_name,
                                 data = df, spar = 0.5)
      }
      
      else if (input$method == "baranyi") {
        p   <- c(y0 = sliderValues()$InitialEstimate[1], mumax = sliderValues()$InitialEstimate[2], K = sliderValues()$InitialEstimate[3], h0 = sliderValues()$InitialEstimate[4])
        lower   <- c(y0 = sliderValues()$LowerEst[1], mumax = sliderValues()$LowerEst[2], K = sliderValues()$LowerEst[3], h0 = sliderValues()$LowerEst[4])
        upper   <- c(y0 = sliderValues()$UpperEst[1], mumax = sliderValues()$UpperEst[2], K = sliderValues()$UpperEst[3], h0 = sliderValues()$UpperEst[4])
        
        many_fits <- all_growthmodels(
          selected_measurement ~ grow_baranyi(time, parms) | strain_name + drug_concentration + replicate_ID + drug_name + media_name,
          data = df,
          p = p, lower = lower, upper = upper,
          transform = "log", ncores = 2)
      }
      
      else if (input$method == "huang") {
        p   <- c(y0 = sliderValues()$InitialEstimate[1], mumax = sliderValues()$InitialEstimate[2], K = sliderValues()$InitialEstimate[3], alpha = sliderValues()$InitialEstimate[5],lambda = sliderValues()$InitialEstimate[6] )
        lower   <- c(y0 = sliderValues()$LowerEst[1], mumax = sliderValues()$LowerEst[2], K = sliderValues()$LowerEst[3], alpha = sliderValues()$LowerEst[5],lambda = sliderValues()$LowerEst[6])
        upper   <- c(y0 = sliderValues()$UpperEst[1], mumax = sliderValues()$UpperEst[2], K = sliderValues()$UpperEst[3], alpha = sliderValues()$UpperEst[5],lambda = sliderValues()$UpperEst[6])
        
        many_fits <- all_growthmodels(
          selected_measurement ~ grow_huang(time, parms) | strain_name + drug_concentration + replicate_ID + drug_name + media_name,
          data = df,
          p = p, lower = lower, upper = upper,
          transform = "log", ncores = 2)
      }
      
    }
    
    # res <- results(many_fits) ### DO I NEED THIS???
    return(many_fits)
  })
  
  # this function determines initial parameter guesses for Emax and sigmoid Emax - step 2
  # guessInitials <- function() { 
  #   res <- results(Fitting())
  #   #set E0 guess as the highest value of conc lower than 0.5 or as the first one
  #   if (subset(res, res$mumax==max(res$mumax))$drug_concentration<0.5) {
  #     E0_x <- max(res$mumax)
  #   } else {E0_x<-res$mumax[1]}
  #   # E_max guess is the difference between the last value minus E0_x
  #   E_max_x <- res$mumax[length(res$mumax)]-E0_x
  #   # EC50 is the drug concentration for which half the maximal growth rate is reached
  #   EC50_yx <- E0_x-((E0_x-res$mumax[length(res$mumax)])/2)
  #   EC50_yx_dif <- abs(res$mumax-EC50_yx)
  #   which_x <- which(EC50_yx_dif==min(EC50_yx_dif))
  #   half_mumax <- res$mumax[which_x]
  #   EC50_x <-  res$drug_concentration[which_x]
  #   # for k the default slider value is used
  #   k_x<-sliderValues_Emax()$Est[4]
  #   # k_x<-5
  #   
  #   # return(c(E0_x,E_max_x,EC50_x,k_x))
  #   return(c(E0=E0_x,E_max=E_max_x,EC50=EC50_x,k=k_x))
  # }
  
  guessInitials <- function(res) { 
    #set E0 guess as the highest value of conc lower than 0.5 or as the first one
    # if (subset(res, res$mumax==max(res$mumax))$drug_concentration<0.5) {# CHANGE THI TO BE THE FIRST
    #   E0_x <- max(res$mumax)
    # } else {E0_x<-res$mumax[1]}
    E0_x<-res$mumax[1]
    
    if (input$model!="capacity_Emax") {
      # E_max guess is the difference between the last value minus E0_x
      E_max_x <- res$mumax[length(res$mumax)]-E0_x
    } else if (input$model=="capacity_Emax") {
      E_max_x <- (res$mumax[length(res$mumax)]/E0_x)-1
    }
    # EC50 is the drug concentration for which half the maximal growth rate is reached
    EC50_yx <- E0_x-((E0_x-res$mumax[length(res$mumax)])/2)
    EC50_yx_dif <- abs(res$mumax-EC50_yx)
    which_x <- which(EC50_yx_dif==min(EC50_yx_dif))
    half_mumax <- res$mumax[which_x]
    EC50_x <-  res$drug_concentration[which_x]
    # for k the default slider value is used
    k_x<-sliderValues_Emax()$Est[4]
    # k_x<-5
    
    return(c(E0=E0_x,E_max=E_max_x,EC50=EC50_x,k=k_x))
  }

  combinations <- function() {
    df <- getData()
    if (input$all == 'no') {
      comb <- data.frame(
        strain_names = input$selected_strain,
        drug_names= input$selected_drug,
        media_names=input$selected_media)
    } else if (input$all == 'yes') {
      comb <- expand.grid(levels(df$strain_name), levels(df$drug_name), levels(df$media_name))
    }
    return(comb)
  }  
  
  #step 2 fitting
  Fit_AB <- eventReactive(input$update_2,{
    comb <- combinations()
    m <- list()
    many_fits <- Fitting()
    for(row in 1:nrow(comb)) { # for every combinations of strain/drug/media
      part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      res <- results(part_many_fits)

      if (input$model=="exp"){
        part_m <- nls(mumax ~ a*exp(b*drug_concentration)+c,
                        data = res, start=c(a=sliderValues2()$Est[1], b= sliderValues2()$Est[2], c=sliderValues2()$Est[3]))
      }
      
      else if (input$model=="sigmoid_Emax") {
        if (input$guess=="yes"){
          Inits <- guessInitials(res)
        } else {Inits <-c(E0=sliderValues_Emax()$Est[1],E_max=sliderValues_Emax()$Est[2], EC50=sliderValues_Emax()$Est[3], k=sliderValues_Emax()$Est[4])}
        part_m <- nls(mumax ~ E0 + E_max *(((drug_concentration/EC50)**k)/(1+((drug_concentration/EC50)**k))),
                 data = res, start=c(Inits[1], Inits[2], Inits[3], Inits[4]))
      }
      
      else if (input$model=="Emax") {
        if (input$guess=="yes"){
          Inits <- guessInitials(res)
        } else {
          Inits <-c(E0=sliderValues_Emax()$Est[1],E_max=sliderValues_Emax()$Est[2], EC50=sliderValues_Emax()$Est[3], k=sliderValues_Emax()$Est[4])}
        part_m <- nls(mumax ~ E0 + E_max *(((drug_concentration/EC50))/(1+((drug_concentration/EC50)))),
                 data = res, start=c(Inits[1], Inits[2], Inits[3]))
      }
      
      else if (input$model=="capacity_Emax") {
        if (input$guess=="yes"){
          Inits <- guessInitials(res)
        } else {
          Inits <-c(E0=sliderValues_Emax()$Est[1],E_max=sliderValues_Emax()$Est[2], EC50=sliderValues_Emax()$Est[3], k=sliderValues_Emax()$Est[4])}
        part_m <- nls(mumax ~ E0*(1+ E_max * (((drug_concentration/EC50)**k)/(1 + ((drug_concentration/EC50)**k)))),
                      data = res, start=c(Inits[1], Inits[2], Inits[3], Inits[4]))
        # E0*(1- E_max * (((drug_concentration/EC50)^k)/(1 + ((drug_concentration/EC50)^k)))
      }
      
      m<-append(m,list(part_m))
      }
    return(m)
    
  })
  

  
#----Output--------------------------------------------------------------------------------------------------------------
  
  # # Show the values of the used input parameters in a table step 1----
  # output$parameters <- renderTable({
  #   paras <- sliderValues()
  #   if (input$method == "smoothing") {
  #     paras <- NA
  #   }
  #   if (input$method == "baranyi") {
  #     paras[5:6,] <- NA
  #   }
  #   if (input$method == "huang") {
  #     paras[4,] <- NA
  #   }
  #   return(paras)
  #   #sliderValues()
  #   },digits=2)

  # Show input data in HTML table
  output$data <- renderTable({
    req(input$file1)
    df <- read.csv(input$file1$datapath, stringsAsFactors=TRUE)
    return(df)
    },digits=6)
  
  # show fitting results in table
  output$results <- renderTable(results(Fitting()), digits=6)
 
  # plot with fitting
  output$plots <- renderPlot({
    #Plot results
    many_fits <- Fitting()
    
    l_plot <- ceiling((length(rownames(results(many_fits))))/5)
    par(mfrow = c(l_plot, 5))
    if (l_plot>6) { # for more than 6 rows
      par(mar = c(2, 1, 1, 1))
    } else if (l_plot<=6) {
      par(mar = c(2, 1, 1, 1))
    } else if (l_plot<=4) {
      par(mar = c(2.5, 4, 2, 1))
    } else if (l_plot<=2) {
      par(mar = c(4, 4, 2, 1))
    }
    
    comb <- combinations()
    for(row in 1:nrow(comb)) {
      part_plot <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      plot(part_plot,ylab="measurement", xlab="time")
    }
  })
  
  
  
  # plot maximal growth rate estimate over AB
  output$Gconc <- renderPlot({
    m <- Fit_AB()
    l_plot <- length(m)
    par(mfrow = c(l_plot, 2))
    par(mar = c(4, 4, 2, 1))
    comb <- combinations()
    many_fits <- Fitting()
    for(row in 1:nrow(comb)) { # for every combinations of strain/drug/media
      part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      part_res <- results(part_many_fits)
      ordered<-order(part_res$drug_concentration)
      
      plot(part_res$drug_concentration[ordered],part_res$mumax[ordered],xlab="drug concentration",ylab="Est. maximal growth rate", main = as.character(part_res$strain_name[[1]]))
      lines(part_res$drug_concentration[ordered],predict(m[[row]])[ordered],lty=2,col="red",lwd=2)
      
      plot(part_res$drug_concentration[ordered],part_res$mumax[ordered], log="x",xlab="Log-scale drug concentration",ylab="Est. maximal growth rate",  main = as.character(part_res$strain_name[[1]]))
      lines(part_res$drug_concentration[ordered],predict(m[[row]])[ordered],lty=2,col="red",lwd=2,log="x")
    }
  })
  
  results_step2 <- function() {
    comb <- combinations()
    many_fits <- Fitting()
    res <- results(many_fits)
    m <- Fit_AB()
    
    columns <- c( "drug_name", "strain_name", "media_name")
    
    if (input$model=="exp") {
      name = sliderValues_Emax()$Name
    } else if (input$model=="sigmoid_Emax") {
      name = c("E0","E_max","EC50","k")
    } else if (input$model=="Emax") {
      name = sliderValues_Emax()$Name[c(1,2,3)]
    } else if (input$model=="capacity_Emax") {
      name = c("E0","E_max","EC50","k")
    }
    
    columns <- c(columns, name)
    ncol <- length(columns)
    d <- setNames(data.frame(matrix(ncol = ncol, nrow = 0)), columns)
    
    for (row in 1:nrow(comb)) {
      part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      part_res <- results(part_many_fits)
      # columns_part_res <- colnames(part_res)
      
      vec1 <- list(drug_name = as.character(comb[row,2]),strain_name = as.character(comb[row,1]),media_name = as.character(comb[row,3]))
      vec2 <- c(coef(m[[row]]))
      vec <-c(vec1,vec2)
      
      #convert factor columns to character
      i <- sapply(d, is.factor)
      d[i] <- lapply(d[i], as.character)
      d <- rbind(d, vec)
    }
    return(d)
  }
  
  output$results_step2 <- renderTable({
    d <- results_step2()
    # return(d)
    },digits=6)
    
  
  output$fit_step2 <- renderPrint({
    m <- Fit_AB()
    comb <- combinations()
    for (i in 1:length(m)) {
      print(paste("Drug= ",comb[i,2], " , Strain= ", comb[i,1], " Media= ", comb[i,3] ))[[1]]
      print(summary(m[[i]]))
    }
  })
  
  output$downloadData_step1 <- downloadHandler(
    filename = function() {
      paste("Step1_results-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(results(Fitting()), file)
    }
  )
  
  output$downloadData_step2 <- downloadHandler(
    filename = function() {
      paste("Step2_results-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(results_step2(), file)
    }
  )
}

# Create Shiny app ----
shinyApp(ui, server)
