# Maximal growth rate estimator - Functions
# Author: Catharina Meyer, Quantitative Pharmacology
# Contact: c.meyer@lacdr.leidenuniv.nl
# Description: This shiny app fits time kill data in two steps. In step 1, 
#              the growth curves over time are fitted and the maximal growth rates 
#              are estimated. In step 2, the PD relation (maximal growth rates over 
#              drug concentration) is fitted.
# Date of last change: 22.04.2021
#----------------------------------------------------------------------------------------------------------

# ------------------retrieve values from manual user input of inital parameters----------------------------
# creates data frame with input values for step 1 (baranyi and huang) from user input
sliderValues <- function(input) { 
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
}

# creates data frame with input values for step 2 (exp) from user input
sliderValues_exp <- function(input) { 
  data.frame(
    Name = c("a", "b","c"),
    Est_s = c(input$a,
              input$b,
              input$c),
    stringsAsFactors = FALSE)
}

# creates data frame with input values for step 2 (all variants of Emax) from user input
sliderValues_Emax <- function(input) {
  data.frame(
    Name = c("E0", "E_max","EC50","k"),
    Est = c(input$E0,
            input$E_max,
            input$EC50,
            input$k),
    stringsAsFactors = FALSE)
}

#------------wrapped functions---------------------------------------------------------------------------------------------
# These functions are called in the app and wrapped in reactive functions
# reads in data from input file
getData_all <- function(f) {
  df_all <- read.csv(f, stringsAsFactors=TRUE)
  return(df_all)
}

# this selects the subset of the data selected by user
getData_subset <- function(input,df_all) {
  #either the rows containing selected strain/drug/media is chosen
  if (input$all=="no") {
    df_subset <- subset(df_all, (df_all$strain_name==input$selected_strain) & (df_all$drug_name==input$selected_drug) & (df_all$media_name==input$selected_media))
  }
  # or the whole dataset is chosen, in that case the rows with NA entries for strain/drug and media are removed 
  # because these rows are only necessary for the correction which is not done in this app but in the preprocessing app
  else if (input$all=="yes") {
    df_subset <- subset(df_all, (!is.na(df_all$strain_name)) & (!is.na(df_all$drug_name)) & (!is.na(df_all$media_name)))
  }
  df_subset$corrected_measurement <- as.numeric(as.character(df_subset$corrected_measurement)) # makes sure the corrected measurements values are numeric
  # data from timepoints is removed as selected in input
  df_subset <- subset(df_subset, (df_subset$time>=as.numeric(input$selected_time))) 
  
  return(df_subset)
}

# sets inital parameters and ranges for fitting step 1
getInitials_step1 <- function(input) {
  if (input$method == "baranyi") {
    # take inital parameters and ranges from user input
    p   <- c(y0 = sliderValues(input)$InitialEstimate[1], mumax = sliderValues(input)$InitialEstimate[2], K = sliderValues(input)$InitialEstimate[3], h0 = sliderValues(input)$InitialEstimate[4])
    lower   <- c(y0 = sliderValues(input)$LowerEst[1], mumax = sliderValues(input)$LowerEst[2], K = sliderValues(input)$LowerEst[3], h0 = sliderValues(input)$LowerEst[4])
    upper   <- c(y0 = sliderValues(input)$UpperEst[1], mumax = sliderValues(input)$UpperEst[2], K = sliderValues(input)$UpperEst[3], h0 = sliderValues(input)$UpperEst[4])
  }
  else if (input$method == "huang") {
    # take inital parameters and ranges from user input
    p   <- c(y0 = sliderValues(input)$InitialEstimate[1], mumax = sliderValues(input)$InitialEstimate[2], K = sliderValues(input)$InitialEstimate[3], alpha = sliderValues(input)$InitialEstimate[5],lambda = sliderValues(input)$InitialEstimate[6] )
    lower   <- c(y0 = sliderValues(input)$LowerEst[1], mumax = sliderValues(input)$LowerEst[2], K = sliderValues(input)$LowerEst[3], alpha = sliderValues(input)$LowerEst[5],lambda = sliderValues(input)$LowerEst[6])
    upper   <- c(y0 = sliderValues(input)$UpperEst[1], mumax = sliderValues(input)$UpperEst[2], K = sliderValues(input)$UpperEst[3], alpha = sliderValues(input)$UpperEst[5],lambda = sliderValues(input)$UpperEst[6])
  }
  
  return(list(p=p,lower=lower,upper=upper))
}

# fits growthcurves using functions from the growthrates package
Fit_Growthrate <- function(input,df) {
  # select measurement type according to user input
  # changes column name df, sorry for the bad style,
  if (input$selected_measurement=="raw_measurement"){
    colnames(df)[which((colnames(df)=="raw_measurement"))] <- "selected_measurement"
  } else if (input$selected_measurement=="corrected_measurement") {
    colnames(df)[which((colnames(df)=="corrected_measurement"))] <- "selected_measurement"
  }
  # for pooled fitting, this means all the replicas are grouped together 
  if (input$pool=="pooled") {
    if (input$method == "smoothing") {
      many_fits <- all_splines(selected_measurement ~ time | strain_name + drug_concentration + drug_name + media_name ,
                               data = df, spar = 0.5)
    }
    else if (input$method == "baranyi") {
      Inits <- getInitials_step1(input)
      
      many_fits <- all_growthmodels(
        selected_measurement ~ grow_baranyi(time, parms) | strain_name + drug_concentration + drug_name + media_name,
        data = df,
        p = Inits[["p"]], lower = Inits[["lower"]], upper = Inits[["upper"]],
        transform = "log", ncores = 2)
    }
    
    else if (input$method == "huang") {
      Inits <- getInitials_step1(input)
      
      many_fits <- all_growthmodels(
        selected_measurement ~ grow_huang(time, parms) |  strain_name + drug_concentration + drug_name + media_name,
        data = df,
        p = Inits[["p"]], lower = Inits[["lower"]], upper = Inits[["upper"]],
        transform = "log", ncores = 2)
    }
  }
  
  # for non-pooled fitting, this means all the replicas are fitted individually 
  else if (input$pool=="per replicate"){
    if (input$method == "smoothing") {
      # Smoothing for mutliple data sets 
      many_fits <- all_splines(selected_measurement ~ time | strain_name + drug_concentration + replicate_ID + drug_name + media_name,
                               data = df, spar = 0.5)
    }
    
    else if (input$method == "baranyi") {
      Inits <- getInitials_step1(input)
      
      many_fits <- all_growthmodels(
        selected_measurement ~ grow_baranyi(time, parms) | strain_name + drug_concentration + replicate_ID + drug_name + media_name,
        data = df,
        p = Inits[["p"]], lower = Inits[["lower"]], upper = Inits[["upper"]],
        transform = "log", ncores = 2)
    }
    
    else if (input$method == "huang") {
      Inits <- getInitials_step1(input)
      
      many_fits <- all_growthmodels(
        selected_measurement ~ grow_huang(time, parms) | strain_name + drug_concentration + replicate_ID + drug_name + media_name,
        data = df,
        p = Inits[["p"]], lower = Inits[["lower"]], upper = Inits[["upper"]],
        transform = "log", ncores = 2)
    }
  }
  return(many_fits)
}

# determines initial parameters for step 2 either using the automated guessing function or the manual user input
getInitials_step2 <- function(input,res) {
  # automated guessing is used if for all Emax methods if guess is set to yes
  if (input$guess=="yes" & input$model!="exp") {
    Inits <- guessInitials(input,res)
    # user input is used if for all Emax methods if guess is set to no
  } else if (input$guess=="no" & input$model!="exp") {
    Inits <-c(E0=sliderValues_Emax(input)$Est[1],E_max=sliderValues_Emax(input)$Est[2], EC50=sliderValues_Emax(input)$Est[3], k=sliderValues_Emax(input)$Est[4])
    # user input is used if for exp, no matter the guess input
  } else if (input$model=="exp") {
    Inits <- c(a=sliderValues_exp(input)$Est[1], b= sliderValues_exp(input)$Est[2], c=sliderValues_exp(input)$Est[3])
  }
  return(Inits)
}

# this function fits the pd relation maximal growth rates over drug concentration for the maximal growth rates estimates from many_fits for every combination of drug/strain/media in comb
Fit_PD <- function(input, comb, res){
  # this is done for all methods expect when capacity emax is chosen and shared fitting
  if (input$fitting_type=="individual" | (input$fitting_type=="shared" & input$model!="capacity_Emax")) {
    m <- list()
    for(row in 1:nrow(comb)) { # for every combinations of strain/drug/media
      # part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      # res <- results(part_many_fits)
      part_res <- subset(res, as.character(res$strain_name)==as.character(comb[row,1]) & as.character(res$drug_name)==as.character(comb[row,2]) & as.character(res$media_name)==as.character(comb[row,3]))

      Inits <- getInitials_step2(input,part_res)
      
      if (input$model=="exp"){
        part_m <- nls(mumax ~ a*exp(b*drug_concentration)+c,
                      data = part_res, start=c(a=sliderValues_exp(input)$Est[1], b= sliderValues_exp(input)$Est[2], c=sliderValues_exp(input)$Est[3]))
      }
      
      else if (input$model=="sigmoid_Emax") {
        part_m <- nls(mumax ~ E0 + E_max *(((drug_concentration/EC50)**k)/(1+((drug_concentration/EC50)**k))),
                      data = part_res, start=c(Inits[1], Inits[2], Inits[3], Inits[4]))
      }
      
      else if (input$model=="Emax") {
        part_m <- nls(mumax ~ E0 + E_max *(((drug_concentration/EC50))/(1+((drug_concentration/EC50)))),
                      data = part_res, start=c(Inits[1], Inits[2], Inits[3]))
      }
      
      else if (input$model=="capacity_Emax") {
        part_m <- nls(mumax ~ E0*(1- E_max * (((drug_concentration/EC50)**k)/(1 + ((drug_concentration/EC50)**k)))),
                      data = part_res, start=c(Inits[1], Inits[2], Inits[3], Inits[4]))
      }
      
      m<-append(m,list(part_m))
    }
    
    # this is done only if capacity emax and shared fitting is chosen
    # !THIS ONLY WORKS for a dataset with one drug, one media and two strains!
  } else if (input$fitting_type=="shared" & input$model=="capacity_Emax") {
    
    # res <- results(many_fits)
    #set up stacked variables
    res_1 <- subset(res,res$strain_name==levels(res$strain_name)[[1]]) # first strain
    res_2 <- subset(res,res$strain_name==levels(res$strain_name)[[2]]) # second strain
    mumax_stacked <- c(res_1$mumax,res_2$mumax) #y
    conc_stacked <- c(res_1$drug_concentration,res_2$drug_concentration) #x
    n1 <- length(res_1$drug_concentration)
    n2 <- length(res_2$drug_concentration)
    lcon1 <- rep(c(1,0), c(n1,n2))
    lcon2 <- rep(c(0,1), c(n1,n2))
    mcon1 <- lcon1
    mcon2 <- lcon2
    
    Inits <- getInitials_step2(input,res)
    Inits_stacked <- c(E0_1=Inits[1],E0_2=Inits[1],E_max=Inits[2],EC50_1=Inits[3], EC50_2=Inits[3], k=Inits[4])
    names(Inits_stacked) <- c("E0_1", "E0_2", "E_max", "EC50_1", "EC50_2", "k")
    
    m <- nls(mumax_stacked ~ (E0_1*lcon1 + E0_2*lcon2) *(1- E_max * (((conc_stacked/(EC50_1*mcon1+EC50_2*mcon2))**k)/(1 + ((conc_stacked/(EC50_1*mcon1+EC50_2*mcon2))**k)))),
             data = res, start=Inits_stacked)
  }
  
  
  return(m)
}

# makes a guess of the inital parameters for the emax model variants
# this only works if the data roughly follows the shape of an emax function
# make sure the first value res$mumax[1] is not 0! 
guessInitials <- function(input,res) {
  #set E0 guess as first value, assuming this is the lowest drug concentration in the dataset
  E0_x<-res$mumax[1]
  if (input$model!="capacity_Emax") {
    # E_max guess is the difference between the last value minus E0_x
    E_max_x <- res$mumax[length(res$mumax)]-E0_x
  } else if (input$model=="capacity_Emax") {
    E_max_x <- 1-(res$mumax[length(res$mumax)]/E0_x)
  }
  # EC50 is the drug concentration for which half the maximal growth rate is reached
  EC50_yx <- E0_x-((E0_x-res$mumax[length(res$mumax)])/2) # growth rate value halfway between E0 and Emax
  EC50_yx_dif <- abs(res$mumax-EC50_yx) # find closest growth rate datapoint
  which_x <- which(EC50_yx_dif==min(EC50_yx_dif))
  half_mumax <- res$mumax[which_x]  # lookup which drug concentration this relates to
  EC50_x <-  mean(res$drug_concentration[which_x]) # drug concentration at which growth rate is closest to half maximal effect
  # default for k
  k_x<-5
  Inits <- c(E0=E0_x,E_max=E_max_x,EC50=EC50_x,k=k_x)
  return(Inits)
}

# makes a data frame with all the combinations of strain/drug/media present in df
combinations <- function(input,df) {
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

# creates a data frame d from the results in step 2 res in the right output format
results_PD <- function(input,comb,res,m) {
  # only for shared fitting
  if (input$model=="capacity_Emax" & input$fitting_type=="shared") {
    columns <- c( "drug_name", "strain_name", "media_name")
    name = c("E0","E_max","EC50","k")
    columns <- c(columns, name)
    ncol <- length(columns)
    d <- setNames(data.frame(matrix(ncol = ncol, nrow = 0)), columns)
    
    # for the first strain
    vec1 <- list(drug_name = as.character(comb[1,2]),strain_name = as.character(comb[1,1]),media_name = as.character(comb[1,3]))
    coefs_m <- c(coef(m))
    vec2 <- c(E0=coefs_m[[1]], E_max=coefs_m[[3]],EC50=coefs_m[[4]],k=coefs_m[[6]])
    vec <-c(vec1,vec2)
    #convert factor columns to character
    i <- sapply(d, is.factor)
    d[i] <- lapply(d[i], as.character)
    d <- rbind(d, vec)
    
    # for the second strain
    vec1 <- list(drug_name = as.character(comb[2,2]),strain_name = as.character(comb[2,1]),media_name = as.character(comb[2,3]))
    coefs_m <- c(coef(m))
    vec2 <- c(E0=coefs_m[[2]], E_max=coefs_m[[3]],EC50=coefs_m[[5]],k=coefs_m[[6]])
    vec <-c(vec1,vec2)
    #convert factor columns to character
    i <- sapply(d, is.factor)
    d[i] <- lapply(d[i], as.character)
    d <- rbind(d, vec)
    
  } else {
    # all other methods 
    columns <- c( "drug_name", "strain_name", "media_name")
    
    if (input$model=="exp") {
      name = sliderValues_Emax(input)$Name
    } else if (input$model=="sigmoid_Emax") {
      name = c("E0","E_max","EC50","k")
    } else if (input$model=="Emax") {
      name = sliderValues_Emax(input)$Name[c(1,2,3)]
    } else if (input$model=="capacity_Emax") {
      name = c("E0","E_max","EC50","k")
    }
    
    columns <- c(columns, name)
    ncol <- length(columns)
    d <- setNames(data.frame(matrix(ncol = ncol, nrow = 0)), columns)
    
    for (row in 1:nrow(comb)) {
      # part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      # part_res <- results(part_many_fits)
      part_res <- subset(res, as.character(res$strain_name)==as.character(comb[row,1]) & as.character(res$drug_name)==as.character(comb[row,2]) & as.character(res$media_name)==as.character(comb[row,3]))
      

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

# creates growthcurve plots
plot_Growthrate <- function(input, many_fits, comb) {
  
  # determines number of subplots and sets number of rows l_plot
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
  
  #creates the subplots using the plotting function from the growthrates package
  for(row in 1:nrow(comb)) {
    part_plot <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
    plot(part_plot,ylab="measurement", xlab="time")
  }
}

# plots fitting resulst of step 2 pd relation
plot_PD <- function(input,res,m,comb) {
  # for all methods expect shared fitting of capacity emax
  if (!(input$model=="capacity_Emax" & input$fitting_type=="shared")) {
    l_plot <- length(m)  #sets number of rows l_plot
    par(mfrow = c(l_plot, 2))
    par(mar = c(4, 4, 2, 1))
    
    for(row in 1:nrow(comb)) { # for every combinations of strain/drug/media
      # part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      # part_res <- results(part_many_fits)
      part_res <- subset(res, as.character(res$strain_name)==as.character(comb[row,1]) & as.character(res$drug_name)==as.character(comb[row,2]) & as.character(res$media_name)==as.character(comb[row,3]))
      
      ordered<-order(part_res$drug_concentration)
      
      plot(part_res$drug_concentration[ordered],part_res$mumax[ordered],xlab="drug concentration",ylab="Est. maximal growth rate", main = as.character(part_res$strain_name[[1]]))
      lines(part_res$drug_concentration[ordered],predict(m[[row]])[ordered],lty=2,col="red",lwd=2)
      
      plot(part_res$drug_concentration[ordered],part_res$mumax[ordered], log="x",xlab="Log-scale drug concentration",ylab="Est. maximal growth rate",  main = as.character(part_res$strain_name[[1]]))
      lines(part_res$drug_concentration[ordered],predict(m[[row]])[ordered],lty=2,col="red",lwd=2,log="x")
    }
  } else {
    # only for shared fitting
    l_plot <- nrow(comb) #sets number of rows l_plot
    par(mfrow = c(l_plot, 2))
    par(mar = c(4, 4, 2, 1))
    
    for(row in 1:nrow(comb)) { # for every combinations of strain/drug/media
      # part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      # part_res <- results(part_many_fits)
      part_res <- subset(res, as.character(res$strain_name)==as.character(comb[row,1]) & as.character(res$drug_name)==as.character(comb[row,2]) & as.character(res$media_name)==as.character(comb[row,3]))
      ordered<-order(part_res$drug_concentration)
      
      plot(part_res$drug_concentration[ordered],part_res$mumax[ordered],xlab="drug concentration",ylab="Est. maximal growth rate", main = as.character(part_res$strain_name[[1]]))
      lines(part_res$drug_concentration[ordered],predict(m)[ordered],lty=2,col="red",lwd=2)
      
      
      plot(part_res$drug_concentration[ordered],part_res$mumax[ordered], log="x",xlab="Log-scale drug concentration",ylab="Est. maximal growth rate",  main = as.character(part_res$strain_name[[1]]))
      lines(part_res$drug_concentration[ordered],predict(m)[ordered],lty=2,col="red",lwd=2,log="x")
      
    }}
}



# plots the results from step 1 many-fits$mumax over drug concentration as well as the fit of the chosen model with the inital parameters BEFORE FITTING
plot_PD_init <- function(input, res, comb) {
  
  # if (!(input$model=="capacity_Emax" & input$fitting_type=="shared")) {
    l_plot <- nrow(comb)
    par(mfrow = c(l_plot, 2))
    par(mar = c(4, 4, 2, 1))
    
    for(row in 1:nrow(comb)) { # for every combinations of strain/drug/media
      # part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      # res <- results(part_many_fits)
      part_res <- subset(res, as.character(res$strain_name)==as.character(comb[row,1]) & as.character(res$drug_name)==as.character(comb[row,2]) & as.character(res$media_name)==as.character(comb[row,3]))
      
      ordered<-order(part_res$drug_concentration)
      
      Inits <- getInitials_step2(input,part_res)
      
      if (input$model=="Emax") {
        y_init <- Inits["E0"] + Inits["E_max"] *(((part_res$drug_concentration[ordered]/Inits["EC50"]))/(1+((part_res$drug_concentration[ordered]/Inits["EC50"]))))
      } else if (input$model=="sigmoid_Emax") {
        y_init <- Inits["E0"] + Inits["E_max"] *(((part_res$drug_concentration[ordered]/Inits["EC50"])**Inits["k"])/(1+((part_res$drug_concentration[ordered]/Inits["EC50"])**Inits["k"])))
      } else if (input$model=="capacity_Emax") {
        y_init <- Inits["E0"]*(1- Inits["E_max"] * (((part_res$drug_concentration[ordered]/Inits["EC50"])**Inits["k"])/(1 + ((part_res$drug_concentration[ordered]/Inits["EC50"])**Inits["k"]))))
      } else if (input$model=="exp") {
        y_init <- Inits["a"]*exp(Inits["b"]*part_res$drug_concentration[ordered])+Inits["c"]
      }
      
      plot(part_res$drug_concentration[ordered],part_res$mumax[ordered],xlab="drug concentration",ylab="Est. maximal growth rate", main = as.character(part_res$strain_name[[1]]))
      lines(part_res$drug_concentration[ordered],y_init,lty=2,col="red",lwd=2)
      
      plot(part_res$drug_concentration[ordered],part_res$mumax[ordered], log="x",xlab="Log-scale drug concentration",ylab="Est. maximal growth rate",  main = as.character(part_res$strain_name[[1]]))
      lines(part_res$drug_concentration[ordered],y_init,lty=2,col="red",lwd=2,log="x")
    }
  # } 
}

# prints summary of fitting step 2
Print_PD <- function(input, comb, m) {
  
  if (input$model=="capacity_Emax" & input$fitting_type=="shared") {
    print(paste("Parameters specific for ", comb[1,1]," : E0_1 and EC50_1 " ))[[1]]
    print(paste("Parameters specific for ", comb[2,1]," : E0_2 and EC50_2 " ))[[1]]
    print("Shared parameters: E_max and k")[[1]]
    print(summary(m))
    
  } else {
    
    for (i in 1:length(m)) {
      print(paste("Drug= ",comb[i,2], " , Strain= ", comb[i,1], " Media= ", comb[i,3] ))[[1]]
      print(summary(m[[i]]))
    }
  }
}
