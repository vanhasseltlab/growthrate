#____________________functional functions_____________________________________________________________________________________________

# input values step 1 fitting
#Reactive expression to create data frame of step 1 input values, only updated when action button is pressed
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

# input values for step 2 exp fitting
sliderValues_exp <- function(input) { 
  data.frame(
    Name = c("a", "b","c"),
    Est_s = c(input$a,
              input$b,
              input$c),
    stringsAsFactors = FALSE)
}

sliderValues_Emax <- function(input) {
  data.frame(
    Name = c("E0", "E_max","EC50","k"),
    Est = c(input$E0,
            input$E_max,
            input$EC50,
            input$k),
    stringsAsFactors = FALSE)
}

getData_all <- function(f) {
  df_all <- read.csv(f, stringsAsFactors=TRUE)
  return(df_all)
}

getData_subset <- function(input,df_all) {
  if (input$all=="no") {
    df_subset <- subset(df_all, (df_all$strain_name==input$selected_strain) & (df_all$drug_name==input$selected_drug) & (df_all$media_name==input$selected_media))
  }
  else if (input$all=="yes") {
    df_subset <- subset(df_all, (!is.na(df_all$strain_name)) & (!is.na(df_all$drug_name)) & (!is.na(df_all$media_name)))
    # df_subset <- df_all
  }
  df_subset$corrected_measurement <- as.numeric(as.character(df_subset$corrected_measurement))
  df_subset <- subset(df_subset, (df_subset$time>=as.numeric(input$selected_time)))
  
  # if (any(df_subset$corrected_measurement<=0)){
  #   added_value <- round(abs(min(df_subset$corrected_measurement)), digits=4) 
  #   df_subset$corrected_measurement <- df_subset$corrected_measurement + added_value
  # }
  
  # if (any(df_subset$raw_measurement<=0)){
  #   added_value <- round(abs(min(df_subset$raw_measurement)), digits=4) 
  #   df_subset$raw_measurement <- df_subset$raw_measurement + added_value
  # }
  return(df_subset)
}

Fit_Growthrate <- function(input,df) {
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
      p   <- c(y0 = sliderValues(input)$InitialEstimate[1], mumax = sliderValues(input)$InitialEstimate[2], K = sliderValues(input)$InitialEstimate[3], h0 = sliderValues(input)$InitialEstimate[4])
      lower   <- c(y0 = sliderValues(input)$LowerEst[1], mumax = sliderValues(input)$LowerEst[2], K = sliderValues(input)$LowerEst[3], h0 = sliderValues(input)$LowerEst[4])
      upper   <- c(y0 = sliderValues(input)$UpperEst[1], mumax = sliderValues(input)$UpperEst[2], K = sliderValues(input)$UpperEst[3], h0 = sliderValues(input)$UpperEst[4])
      
      many_fits <- all_growthmodels(
        selected_measurement ~ grow_baranyi(time, parms) | strain_name + drug_concentration + drug_name + media_name,
        data = df,
        p = p, lower = lower, upper = upper,
        transform = "log", ncores = 2)
    }
    
    else if (input$method == "huang") {
      p   <- c(y0 = sliderValues(input)$InitialEstimate[1], mumax = sliderValues(input)$InitialEstimate[2], K = sliderValues(input)$InitialEstimate[3], alpha = sliderValues(input)$InitialEstimate[5],lambda = sliderValues(input)$InitialEstimate[6] )
      lower   <- c(y0 = sliderValues(input)$LowerEst[1], mumax = sliderValues(input)$LowerEst[2], K = sliderValues(input)$LowerEst[3], alpha = sliderValues(input)$LowerEst[5],lambda = sliderValues(input)$LowerEst[6])
      upper   <- c(y0 = sliderValues(input)$UpperEst[1], mumax = sliderValues(input)$UpperEst[2], K = sliderValues(input)$UpperEst[3], alpha = sliderValues(input)$UpperEst[5],lambda = sliderValues(input)$UpperEst[6])
      
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
      p   <- c(y0 = sliderValues(input)$InitialEstimate[1], mumax = sliderValues(input)$InitialEstimate[2], K = sliderValues(input)$InitialEstimate[3], h0 = sliderValues(input)$InitialEstimate[4])
      lower   <- c(y0 = sliderValues(input)$LowerEst[1], mumax = sliderValues(input)$LowerEst[2], K = sliderValues(input)$LowerEst[3], h0 = sliderValues(input)$LowerEst[4])
      upper   <- c(y0 = sliderValues(input)$UpperEst[1], mumax = sliderValues(input)$UpperEst[2], K = sliderValues(input)$UpperEst[3], h0 = sliderValues(input)$UpperEst[4])
      
      many_fits <- all_growthmodels(
        selected_measurement ~ grow_baranyi(time, parms) | strain_name + drug_concentration + replicate_ID + drug_name + media_name,
        data = df,
        p = p, lower = lower, upper = upper,
        transform = "log", ncores = 2)
    }
    
    else if (input$method == "huang") {
      p   <- c(y0 = sliderValues(input)$InitialEstimate[1], mumax = sliderValues(input)$InitialEstimate[2], K = sliderValues(input)$InitialEstimate[3], alpha = sliderValues(input)$InitialEstimate[5],lambda = sliderValues(input)$InitialEstimate[6] )
      lower   <- c(y0 = sliderValues(input)$LowerEst[1], mumax = sliderValues(input)$LowerEst[2], K = sliderValues(input)$LowerEst[3], alpha = sliderValues(input)$LowerEst[5],lambda = sliderValues(input)$LowerEst[6])
      upper   <- c(y0 = sliderValues(input)$UpperEst[1], mumax = sliderValues(input)$UpperEst[2], K = sliderValues(input)$UpperEst[3], alpha = sliderValues(input)$UpperEst[5],lambda = sliderValues(input)$UpperEst[6])
      
      many_fits <- all_growthmodels(
        selected_measurement ~ grow_huang(time, parms) | strain_name + drug_concentration + replicate_ID + drug_name + media_name,
        data = df,
        p = p, lower = lower, upper = upper,
        transform = "log", ncores = 2)
    }
  }
  return(many_fits)
}


Fit_PD <- function(input,comb, many_fits){
  if (input$fitting_type=="individual" | (input$fitting_type=="shared" & input$model!="capacity_Emax")) {
    m <- list()
    for(row in 1:nrow(comb)) { # for every combinations of strain/drug/media
      part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      res <- results(part_many_fits)
      
      if (input$model=="exp"){
        part_m <- nls(mumax ~ a*exp(b*drug_concentration)+c,
                      data = res, start=c(a=sliderValues_exp(input)$Est[1], b= sliderValues_exp(input)$Est[2], c=sliderValues_exp(input)$Est[3]))
      }
      
      else if (input$model=="sigmoid_Emax") {
        if (input$guess=="yes"){
          Inits <- guessInitials(input,res)
        } else {Inits <-c(E0=sliderValues_Emax(input)$Est[1],E_max=sliderValues_Emax(input)$Est[2], EC50=sliderValues_Emax(input)$Est[3], k=sliderValues_Emax(input)$Est[4])}
        part_m <- nls(mumax ~ E0 + E_max *(((drug_concentration/EC50)**k)/(1+((drug_concentration/EC50)**k))),
                      data = res, start=c(Inits[1], Inits[2], Inits[3], Inits[4]))
      }
      
      else if (input$model=="Emax") {
        if (input$guess=="yes"){
          Inits <- guessInitials(input,res)
        } else {
          Inits <-c(E0=sliderValues_Emax(input)$Est[1],E_max=sliderValues_Emax(input)$Est[2], EC50=sliderValues_Emax(input)$Est[3], k=sliderValues_Emax(input)$Est[4])}
        part_m <- nls(mumax ~ E0 + E_max *(((drug_concentration/EC50))/(1+((drug_concentration/EC50)))),
                      data = res, start=c(Inits[1], Inits[2], Inits[3]))
      }
      
      else if (input$model=="capacity_Emax") {
        if (input$guess=="yes"){
          Inits <- guessInitials(input,res)
        } else {
          Inits <-c(E0=sliderValues_Emax(input)$Est[1],E_max=sliderValues_Emax(input)$Est[2], EC50=sliderValues_Emax(input)$Est[3], k=sliderValues_Emax(input)$Est[4])}
        part_m <- nls(mumax ~ E0*(1- E_max * (((drug_concentration/EC50)**k)/(1 + ((drug_concentration/EC50)**k)))),
                      data = res, start=c(Inits[1], Inits[2], Inits[3], Inits[4]))
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
    
    Inits <- guessInitials(input,res) # maybe do separately for each strain
    Inits_stacked <- c(E0_1=Inits[1],E0_2=Inits[1],E_max=Inits[2],EC50_1=Inits[3], EC50_2=Inits[3], k=Inits[4])
    names(Inits_stacked) <- c("E0_1", "E0_2", "E_max", "EC50_1", "EC50_2", "k")
    m <- nls(mumax_stacked ~ (E0_1*lcon1 + E0_2*lcon2) *(1- E_max * (((conc_stacked/(EC50_1*mcon1+EC50_2*mcon2))**k)/(1 + ((conc_stacked/(EC50_1*mcon1+EC50_2*mcon2))**k)))),
             data = res, start=Inits_stacked)
  }
  
  
  return(m)
}

guessInitials <- function(input,res) {
  #set E0 guess as first one value
  E0_x<-res$mumax[1]
  if (input$model!="capacity_Emax") {
    # E_max guess is the difference between the last value minus E0_x
    E_max_x <- res$mumax[length(res$mumax)]-E0_x
  } else if (input$model=="capacity_Emax") {
    E_max_x <- 1-(res$mumax[length(res$mumax)]/E0_x)
  }
  # EC50 is the drug concentration for which half the maximal growth rate is reached
  EC50_yx <- E0_x-((E0_x-res$mumax[length(res$mumax)])/2)
  EC50_yx_dif <- abs(res$mumax-EC50_yx)
  which_x <- which(EC50_yx_dif==min(EC50_yx_dif))
  half_mumax <- res$mumax[which_x]
  EC50_x <-  mean(res$drug_concentration[which_x])
  # for k the default slider value is used
  # k_x<-sliderValues_Emax()$Est[4]
  k_x<-5
  Inits <- c(E0=E0_x,E_max=E_max_x,EC50=EC50_x,k=k_x)
  return(Inits)
}

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

results_PD <- function(input,comb,many_fits,res,m) {
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
      part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      part_res <- results(part_many_fits)

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

plot_Growthrate <- function(input, many_fits, comb) {
  
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
  
  for(row in 1:nrow(comb)) {
    part_plot <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
    plot(part_plot,ylab="measurement", xlab="time")
  }
}

plot_PD <- function(input,many_fits,m,comb) {
  if (!(input$model=="capacity_Emax" & input$fitting_type=="shared")) {
    l_plot <- length(m)
    par(mfrow = c(l_plot, 2))
    par(mar = c(4, 4, 2, 1))
    
    for(row in 1:nrow(comb)) { # for every combinations of strain/drug/media
      part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      part_res <- results(part_many_fits)
      ordered<-order(part_res$drug_concentration)
      
      plot(part_res$drug_concentration[ordered],part_res$mumax[ordered],xlab="drug concentration",ylab="Est. maximal growth rate", main = as.character(part_res$strain_name[[1]]))
      lines(part_res$drug_concentration[ordered],predict(m[[row]])[ordered],lty=2,col="red",lwd=2)
      
      plot(part_res$drug_concentration[ordered],part_res$mumax[ordered], log="x",xlab="Log-scale drug concentration",ylab="Est. maximal growth rate",  main = as.character(part_res$strain_name[[1]]))
      lines(part_res$drug_concentration[ordered],predict(m[[row]])[ordered],lty=2,col="red",lwd=2,log="x")
    }
  } else {
    l_plot <- nrow(comb)
    par(mfrow = c(l_plot, 2))
    par(mar = c(4, 4, 2, 1))
    
    for(row in 1:nrow(comb)) { # for every combinations of strain/drug/media
      part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      part_res <- results(part_many_fits)
      ordered<-order(part_res$drug_concentration)
      
      plot(part_res$drug_concentration[ordered],part_res$mumax[ordered],xlab="drug concentration",ylab="Est. maximal growth rate", main = as.character(part_res$strain_name[[1]]))
      lines(part_res$drug_concentration[ordered],predict(m)[ordered],lty=2,col="red",lwd=2)
      
      
      plot(part_res$drug_concentration[ordered],part_res$mumax[ordered], log="x",xlab="Log-scale drug concentration",ylab="Est. maximal growth rate",  main = as.character(part_res$strain_name[[1]]))
      lines(part_res$drug_concentration[ordered],predict(m)[ordered],lty=2,col="red",lwd=2,log="x")
      
    }}
}


plot_PD_init <- function(input,many_fits,comb) {
  
  # if (!(input$model=="capacity_Emax" & input$fitting_type=="shared")) {
    l_plot <- nrow(comb)
    par(mfrow = c(l_plot, 2))
    par(mar = c(4, 4, 2, 1))
    
    for(row in 1:nrow(comb)) { # for every combinations of strain/drug/media
      part_many_fits <- subset(many_fits, as.character(results(many_fits)$strain_name)==as.character(comb[row,1]) & as.character(results(many_fits)$drug_name)==as.character(comb[row,2]) & as.character(results(many_fits)$media_name)==as.character(comb[row,3]))
      res <- results(part_many_fits)
      ordered<-order(res$drug_concentration)

      
      if (input$guess=="yes" & input$model!="exp") {
        Inits <- guessInitials(input,res)
      } else if (input$guess=="no" & input$model!="exp") {
        Inits <-c(E0=sliderValues_Emax(input)$Est[1],E_max=sliderValues_Emax(input)$Est[2], EC50=sliderValues_Emax(input)$Est[3], k=sliderValues_Emax(input)$Est[4])
      } else if (input$guess=="no" & input$model!="exp") {
        Inits <-c(E0=sliderValues_Emax(input)$Est[1],E_max=sliderValues_Emax(input)$Est[2], EC50=sliderValues_Emax(input)$Est[3], k=sliderValues_Emax(input)$Est[4])
      } else if (input$model=="exp") {
        Inits <- c(a=sliderValues_exp(input)$Est[1], b= sliderValues_exp(input)$Est[2], c=sliderValues_exp(input)$Est[3])
      }
      
      
      if (input$model=="Emax") {
        y_init <- Inits["E0"] + Inits["E_max"] *(((res$drug_concentration[ordered]/Inits["EC50"]))/(1+((res$drug_concentration[ordered]/Inits["EC50"]))))
      } else if (input$model=="sigmoid_Emax") {
        y_init <- Inits["E0"] + Inits["E_max"] *(((res$drug_concentration[ordered]/Inits["EC50"])**Inits["k"])/(1+((res$drug_concentration[ordered]/Inits["EC50"])**Inits["k"])))
      } else if (input$model=="capacity_Emax") {
        y_init <- Inits["E0"]*(1- Inits["E_max"] * (((res$drug_concentration[ordered]/Inits["EC50"])**Inits["k"])/(1 + ((res$drug_concentration[ordered]/Inits["EC50"])**Inits["k"]))))
      } else if (input$model=="exp") {
        y_init <- Inits["a"]*exp(Inits["b"]*res$drug_concentration[ordered])+Inits["c"]
      }
      
      plot(res$drug_concentration[ordered],res$mumax[ordered],xlab="drug concentration",ylab="Est. maximal growth rate", main = as.character(res$strain_name[[1]]))
      lines(res$drug_concentration[ordered],y_init,lty=2,col="red",lwd=2)
      
      plot(res$drug_concentration[ordered],res$mumax[ordered], log="x",xlab="Log-scale drug concentration",ylab="Est. maximal growth rate",  main = as.character(res$strain_name[[1]]))
      lines(res$drug_concentration[ordered],y_init,lty=2,col="red",lwd=2,log="x")
    }
  # } 
}

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
