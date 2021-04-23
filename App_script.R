# Maximal growth rate estimator - Script version
# Author: Catharina Meyer, Quantitative Pharmacology
# Contact: c.meyer@lacdr.leidenuniv.nl
# Description: This shiny app fits time kill data in two steps. In step 1, 
#              the growth curves over time are fitted and the maximal growth rates 
#              are estimated. In step 2, the PD relation (maximal growth rates over 
#              drug concentration) is fitted.
# Date of last change: 22.04.2021
#----------------------------------------------------------------------------------------------------------
library("deSolve")
library("lattice")
library("growthrates") # this is the package which does the main fitting in step 1

# set your working directory, for example choose the folder where this script is located
setwd("~/ownCloud2/cathi/R_files")
source("App/growthrate/Functions.R")

# set path to data file
# Example
f <-"data/ProcessingOutput-new/Preprocessed_defaultFolder.csv"

#-----------------------User input-------------------------------------------------------------------------

# set user input, in the app all of this is in the sidepanel
input <-data.frame(                                                     # OPTIONS
  all =  'yes',                                           # Process the whole plate?: either 'yes' or 'no'
  # if you don't want to process the whole plate:                 
  selected_strain = "STR_A",                              # Select a strain
  selected_drug = "COL",                                  # Select a drug
  selected_media = "MED_A",                               # select a media
  selected_measurement = "raw_measurement",         # Select measurement data: "corrected_measurement" or "raw_measurement"
  selected_time = "0",                                    # Exclude time points before: enter any number in string format, e.g. "0"
  pool = "pooled",                                        # Compute growth rate with data that is: "pooled" or "per replicate"
  method = "smoothing",                                     # Choose a fitting method for step 1: "smoothing" or "baranyi" or "huang"
  model = "sigmoid_Emax",                                 # Choose a model for step 2: "exp", "Emax", "sigmoid_Emax", "capacity_Emax"
  guess = "no",                                           # Guess step 2 parameters: "yes" or "no"
  fitting_type = "individual",                            # for capacity_Emax choose "individual" or "shared"
  
  
  # initial parameter estimates step 1- manual user input, required for baranyi and huang
  y0 = 0.1,#
  mumax = 0.1,
  K = 0.2,
  # for baranyi
  h0 = 1,
  # for huang
  alpha = 1,
  lambda = 0.1,
  
  # ranges for initial parameters baranyi and huang
  y0_lower = 0.01,
  mumax_lower = 0.01, 
  K_lower = 0.01,
  h0_lower = 0,
  alpha_lower = 0, 
  lambda_lower = 0,
  y0_upper = 0.6,
  mumax_upper = 5, 
  K_upper = 2,
  h0_upper = 10,
  alpha_upper = 7, 
  lambda_upper = 4,
  
  # initial parameter estimates step 2- manual user input, always required for exp and for all variants of Emax if automated guess is set to no
  
  # exponential model
  a = -0.4,
  b = -5,
  c = 0.6,
  
  # all variants of Emax
  E0 = 0.4,
  E_max = -0.4, 
  EC50 = 0.2,
  k = 5
)

#----------------------START------------------------------------------------------------------------------------

df_all <- getData_all(f) # read in all data
df_subset <- getData_subset(input, df_all) # select subset according to user input, use this always

#----------------------------------------------------------------------------------------------------------

# if you want you can make changes to the dataset here
# Example:
# df_subset$corrected_measurement <- abs(df_subset$corrected_measurement) # take absolute values of corrected_measurement
# This is not actually useful 
#----------------------------------------------------------------------------------------------------------


#--------------Fit Growthcurves and estimate maximal growth rates - step 1--------------------------------------

many_fits <- Fit_Growthrate(input, df_subset) # run fitting step 1
res <- results(many_fits) # results from step 1 fitting in data frame format
# look at this dataframe by clinking on the variable in you Environment Tab on the right top panel

#--------------Plot Growthcurves and fit - step 1--------------------------------------------------------------

comb <- combinations(input, df_subset) # create data frame with all the combinations of drug/strain/media in subset of data
plot_Growthrate(input, many_fits, comb)

#----------------------------------------------------------------------------------------------------------
# If you want you can make changes to res here, 
# for example if you want to set some maximal growth rates to zero
# or use results from different methods for different drug concentrations etc.
# Example:
# res$mumax[[1]] <- 0 # again not actually useful
#----------------------------------------------------------------------------------------------------------


#--------------Plot fit of initial parameters BEFORE FITTING- PD relation step 2-------------------------------

plot_PD_init(input,res,comb)

#--------------Fit PD relation - step 2------------------------------------------------------------------------

m <- Fit_PD(input, comb, res) # fit pd relation
d <- results_PD(input, comb, res, m) # results from step 2 fitting in data frame format
# look at this dataframe by clinking on the variable in you Environment Tab on the right top panel

#--------------Plot PD relation and model fit - step 2---------------------------------------------------------

plot_PD(input, res, m, comb)

#--------------Print summary of fit - step 2-------------------------------------------------------------------

Print_PD(input, comb, m)

