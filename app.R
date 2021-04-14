library(shiny)
library(lattice)
library(deSolve)
library("growthrates")

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
        condition = "input.model == 'capacity_Emax'",
        radioButtons(inputId = "fitting_type", label="Choose fitting type for parameters:", choices = c("individual","shared"), selected = "individual"),
        
        ),
      
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
        condition = ("(input.model == 'Emax') | (input.model == 'sigmoid_Emax') | ((input.model == 'capacity_Emax') & (input.fitting_type=='individual'))"),
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
                  tabPanel(HTML("Step 2:<br/>PD relation<br/>Results"), tableOutput("results_PD")),
                  tabPanel(HTML("Step 2:<br/>PD relation<br/>Fit Summary"), verbatimTextOutput("fit_step2"))
                  )
      )
    )
  )


#--------------------------------------------------------------------------------------------------------------------------------------
#Server function ----------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------
server <- function(input, output,session) {

#____________________event and input functions_____________________________________________________________________________________________
  
  observeEvent(input$file1, {
    df_all <- getData_all_R()
    updateSelectInput(session, "selected_strain", label = "Select a strain", choices = levels(df_all$strain_name))
    updateSelectInput(session, "selected_drug", label = "Select a drug", choices = levels(df_all$drug_name))
    updateSelectInput(session, "selected_media", label = "Select a media", choices = levels(df_all$media_name))
    updateSelectInput(session, "selected_time", label = 'Exclude time points before:', choices = levels(as.factor((df_all$time))))
    
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
  sliderValues_exp <- eventReactive(input$update_2,{
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
  
#_____________________________reactive wrapper functions_(everything input or app specific goes here)___________________________________________________________________________________
  
  getData_all_R <- reactive({
    req(input$file1)
    df_all <- getData_all(input$file1$datapath)
    return(df_all)
  })
  
  getData_subset_R <- reactive({
    df_all <- getData_all_R()
    df_subset <- getData_subset(df_all)
  return(df_subset)
    })

  # step 1 fitting/ not modularized yet
  Fit_Growthrate_R <- eventReactive(input$update_1,{
    df <- getData_subset_R()
    many_fits <- Fit_Growthrate(df)
    return(many_fits)
  })

  #step 2 fitting
  Fit_PD_R <- eventReactive(input$update_2,{
    comb <- combinations_R()
    many_fits <- Fit_Growthrate_R() # IS many_fits necessary here? - rewrite Fit_PD()
    m <- Fit_PD(comb,many_fits)
    return(m)
  })

  combinations_R <- reactive({
    df <- getData_subset_R() # why subset
    comb <- combinations(df)
    return(comb)
  })  
  
  results_PD_R <- reactive({
    comb <- combinations_R()
    many_fits <- Fit_Growthrate_R()
    res <- results(many_fits) # replace with Results()???
    m <- Fit_PD_R()
    d <- results_PD(comb,many_fits,res,m)
    return(d)
  })
  
  #____________________functional functions_____________________________________________________________________________________________
  
  getData_all <- function(f) {
    df_all <- read.csv(f, stringsAsFactors=TRUE)
    return(df_all)
  }
  
  getData_subset <- function(df_all) {
    if (input$all=="no") {
      df_subset <- subset(df_all, (df_all$strain_name==input$selected_strain) & (df_all$drug_name==input$selected_drug) & (df_all$media_name==input$selected_media))
    }
    else if (input$all=="yes") {
      df_subset <- df_all
    }
    df_subset$corrected_measurement <- as.numeric(as.character(df_subset$corrected_measurement))
    df_subset <- subset(df_subset, (df_subset$time>=input$selected_time))
    return(df_subset)
  }
  
  Fit_Growthrate <- function(df) {
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
    return(many_fits)
  }
  
  
  # Fit_PD <- function(comb, many_fits){
  #   m <- list()
  #   for(row in 1:nrow(comb)) { # for every combinations of strain/drug/media
  #     part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
  #     res <- results(part_many_fits)
  #     
  #     if (input$model=="exp"){
  #       part_m <- nls(mumax ~ a*exp(b*drug_concentration)+c,
  #                     data = res, start=c(a=sliderValues_exp()$Est[1], b= sliderValues_exp()$Est[2], c=sliderValues_exp()$Est[3]))
  #     }
  #     
  #     else if (input$model=="sigmoid_Emax") {
  #       if (input$guess=="yes"){
  #         Inits <- guessInitials(res)
  #       } else {Inits <-c(E0=sliderValues_Emax()$Est[1],E_max=sliderValues_Emax()$Est[2], EC50=sliderValues_Emax()$Est[3], k=sliderValues_Emax()$Est[4])}
  #       part_m <- nls(mumax ~ E0 + E_max *(((drug_concentration/EC50)**k)/(1+((drug_concentration/EC50)**k))),
  #                     data = res, start=c(Inits[1], Inits[2], Inits[3], Inits[4]))
  #     }
  #     
  #     else if (input$model=="Emax") {
  #       if (input$guess=="yes"){
  #         Inits <- guessInitials(res)
  #       } else {
  #         Inits <-c(E0=sliderValues_Emax()$Est[1],E_max=sliderValues_Emax()$Est[2], EC50=sliderValues_Emax()$Est[3], k=sliderValues_Emax()$Est[4])}
  #       part_m <- nls(mumax ~ E0 + E_max *(((drug_concentration/EC50))/(1+((drug_concentration/EC50)))),
  #                     data = res, start=c(Inits[1], Inits[2], Inits[3]))
  #     }
  #     
  #     else if (input$model=="capacity_Emax") {
  #       if (input$guess=="yes"){
  #         Inits <- guessInitials(res)
  #       } else {
  #         Inits <-c(E0=sliderValues_Emax()$Est[1],E_max=sliderValues_Emax()$Est[2], EC50=sliderValues_Emax()$Est[3], k=sliderValues_Emax()$Est[4])}
  #       part_m <- nls(mumax ~ E0*(1+ E_max * (((drug_concentration/EC50)**k)/(1 + ((drug_concentration/EC50)**k)))),
  #                     data = res, start=c(Inits[1], Inits[2], Inits[3], Inits[4]))
  #       # E0*(1- E_max * (((drug_concentration/EC50)^k)/(1 + ((drug_concentration/EC50)^k)))
  #     }
  #     
  #     m<-append(m,list(part_m))
  #   }
  #   return(m)
  # }
  
  Fit_PD <- function(comb, many_fits){
    if (input$fitting_type=="individual" | (input$fitting_type=="shared" & input$model!="capacity_Emax")) {
      m <- list()
      for(row in 1:nrow(comb)) { # for every combinations of strain/drug/media
        part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
        res <- results(part_many_fits)
        
        if (input$model=="exp"){
          part_m <- nls(mumax ~ a*exp(b*drug_concentration)+c,
                        data = res, start=c(a=sliderValues_exp()$Est[1], b= sliderValues_exp()$Est[2], c=sliderValues_exp()$Est[3]))
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
      
      
    } else if (input$fitting_type=="shared" & input$model=="capacity_Emax") {
      
      
      res <- results(many_fits)
      #set up stacked variables:
      res_1 <- subset(res,res$strain_name==levels(res$strain_name)[[1]]) #mutant
      res_2 <- subset(res,res$strain_name==levels(res$strain_name)[[2]]) #strain1
      mumax_stacked <- c(res_1$mumax,res_2$mumax) #y
      conc_stacked <- c(res_1$drug_concentration,res_2$drug_concentration) #x
      n1 <- length(res_1$drug_concentration)
      n2 <- length(res_2$drug_concentration)
      lcon1 <- rep(c(1,0), c(n1,n2))
      lcon2 <- rep(c(0,1), c(n1,n2))
      mcon1 <- lcon1
      mcon2 <- lcon2
      
      Inits <- guessInitials(res) # maybe do separately for each strain
      Inits_stacked <- c(E0_1=Inits[1],E0_2=Inits[1],E_max=Inits[2],EC50_1=Inits[3], EC50_2=Inits[3], k=Inits[4])
      names(Inits_stacked) <- c("E0_1", "E0_2", "E_max", "EC50_1", "EC50_2", "k")
      m <- nls(mumax_stacked ~ (E0_1*lcon1 + E0_2*lcon2) *(1+ E_max * (((conc_stacked/(EC50_1*mcon1+EC50_2*mcon2))**k)/(1 + ((conc_stacked/(EC50_1*mcon1+EC50_2*mcon2))**k)))),
               data = res, start=Inits_stacked)
    }
    
    
    return(m)
  }
  
  guessInitials <- function(res) {
    #set E0 guess as first one value
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
    Inits <- c(E0=E0_x,E_max=E_max_x,EC50=EC50_x,k=k_x)
    return(Inits)
  }
  
  combinations <- function(df) {
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
  
  # results_PD <- function(comb,many_fits,res,m) {
  #   
  #   columns <- c( "drug_name", "strain_name", "media_name")
  #   
  #   if (input$model=="exp") {
  #     name = sliderValues_Emax()$Name
  #   } else if (input$model=="sigmoid_Emax") {
  #     name = c("E0","E_max","EC50","k")
  #   } else if (input$model=="Emax") {
  #     name = sliderValues_Emax()$Name[c(1,2,3)]
  #   } else if (input$model=="capacity_Emax") {
  #     name = c("E0","E_max","EC50","k")
  #   }
  #   
  #   columns <- c(columns, name)
  #   ncol <- length(columns)
  #   d <- setNames(data.frame(matrix(ncol = ncol, nrow = 0)), columns)
  #   
  #   for (row in 1:nrow(comb)) {
  #     part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
  #     part_res <- results(part_many_fits)
  #     # columns_part_res <- colnames(part_res)
  #     
  #     vec1 <- list(drug_name = as.character(comb[row,2]),strain_name = as.character(comb[row,1]),media_name = as.character(comb[row,3]))
  #     vec2 <- c(coef(m[[row]]))
  #     vec <-c(vec1,vec2)
  #     
  #     #convert factor columns to character
  #     i <- sapply(d, is.factor)
  #     d[i] <- lapply(d[i], as.character)
  #     d <- rbind(d, vec)
  #   }
  #   return(d)
  # }
  
  results_PD <- function(comb,many_fits,res,m) {
    if (input$model=="capacity_Emax" & input$fitting_type=="shared") {
      columns <- c( "drug_name", "strain_name", "media_name")
      name = c("E0_1","E0_2","E_max","EC50_1", "EC50_2","k")
      columns <- c(columns, name)
      ncol <- length(columns)
      d <- setNames(data.frame(matrix(ncol = ncol, nrow = 0)), columns)
      res <- results(many_fits)
      # columns_part_res <- colnames(part_res)
      
      vec1 <- list(drug_name = as.character(res$drug_name[[1]]),strain_name = paste(as.character(levels(res$strain_name)[[1]]),":",as.character(levels(res$strain_name)[[2]]) ),media_name = as.character(res$media_name[[1]]))
      vec2 <- c(coef(m))
      vec <-c(vec1,vec2)
      
      #convert factor columns to character
      i <- sapply(d, is.factor)
      d[i] <- lapply(d[i], as.character)
      d <- rbind(d, vec)
      
    } else {
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
      
    }
    return(d)
  }
  
  
  # REMOVE THIS ?
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
  
  # show fitting results in table
  output$results <- renderTable(Results(), digits=6)
 
  # plot with fitting
  output$plots <- renderPlot({
    #Plot results
    many_fits <- Fit_Growthrate_R()
    
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
    
    comb <- combinations_R()
    for(row in 1:nrow(comb)) {
      part_plot <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      plot(part_plot,ylab="measurement", xlab="time")
    }
  })
  
  # FIX FOR WHOLE PLATE CAPACITY SHARED
  # plot maximal growth rate estimate over AB
  output$Gconc <- renderPlot({
    m <- Fit_PD_R()
    l_plot <- length(m)
    par(mfrow = c(l_plot, 2))
    par(mar = c(4, 4, 2, 1))
    comb <- combinations_R()
    many_fits <- Fit_Growthrate_R()
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
  

  output$results_PD <- renderTable({
    d <- results_PD_R()
    },digits=6)
    
  # FIX FOR WHOLE PLATE CAPACITY SHARED
  output$fit_step2 <- renderPrint({
    m <- Fit_PD_R()
    comb <- combinations_R()
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
      write.csv(Results(), file)
    }
  )
  
  output$downloadData_step2 <- downloadHandler(
    filename = function() {
      paste("Step2_results-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(results_PD_R(), file)
    }
  )
}

# Create Shiny app ----
shinyApp(ui, server)
