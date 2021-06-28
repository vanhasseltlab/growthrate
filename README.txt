# growthrate

Shiny app- Maximal growth rate estimator 
Link: https://vanhasseltlab.lacdr.leidenuniv.nl/growthrate/

Author: Catharina Meyer
Contact: c.meyer@lacdr.leidenuniv.nl
Quantitative Pharmacology, LACDR, Leiden University

This shiny app fits time kill data in two steps. In step 1, the growth curves over time are fitted and the maximal growth rates are estimated. In step 2, the relation of maximal growth rates over drug concentration is fitted.

USER INSTRUCTIONS:

- Before using this shiny app, use the shiny app: Growth Curve for 96-Well Plate (http://132.229.100.197:2222/ot2/PlateAnalysis/GrowthCurve/) to prepocess the data.

- Upload a preprocessed data file in csv format to this shiny app

- Either choose to process all the data from the input file or select a set of strain/drug/media and measurement (raw/corrected)

- It is possible to exclude time points from further analysis. For this select the fisrt time point that should be included. 

- The action buttons run step 1 and step 2 separately. Restart estimation when changing input parameters. If changing step 2 methods or parameters, it is only
 necessary to rerun step 2.
 
- Pooling of the data is available. If 'None' is choosen, all single data sets (each row in the input file) are treated individually. If 'replicates' is chosen, the data sets indicated to be replicates are fitted together in step 1. For 'replicates & concentration', replicates and data sets with different concentrations are fitted together. Plots from step 1 match how the data is pooled.

- In step 1, the measurement data ~ time is fitted and the maximal growth rate (mumax) is estimated. Three different methods from the R package "growthrates" (https://tpetzoldt.github.io/growthrates/doc/Introduction.html) are available. Check the documentation of the growthrates package for further explanations of the methods.
	1) Smoothing splines: No inital parameter estimates need to be given
	2) Baranyi model: provide inital estimates for the parameters, for details see https://doi.org/10.1016/0168-1605(94)00121-L)
	3) Huang model: provide inital estimates for the parameters, for details see https://doi.org/10.1111/j.1750-3841.2008.00785.x
	Tip: to find suitable inital parameter estimates, try to use the smoothing splines method first
	
- In step 2, the estimated maximal growth rates (mumax) ~ drug concentration is fitted. Three different model equations are available for fitting:
	1) exponential model equation: mumax ~ a * exp(b * drug_concentration) + c
	   Parameters a, b & c
	2) Emax model equation: mumax ~ E0 + E_max * (((drug_concentration/EC50))/(1 + ((drug_concentration/EC50))))
	   This model requires estimates from three parameters E0, E_max and EC50 which can be automatically guessed by the shiny app.
	3) sigmoidal Emax model equation: mumax ~ E0 + E_max * (((drug_concentration/EC50)^k)/(1 + ((drug_concentration/EC50)^k)))
	   This model requires estimates from four parameters E0, E_max, EC50 and the Hill parameter k which can be automatically guessed by the shiny app.
	4) capacity-limited Emax model: mumax = E0*(1- E_max * (((drug_concentration/EC50)^k)/(1 + ((drug_concentration/EC50)^k)))
	   This model requires estimates from four parameters E0, E_max, EC50 and the Hill parameter k which can be automatically guessed by the shiny app.
	
- shared parameter fitting: It is possible to fit a whole plate (with one drug, two strains and one media) at once using the capacity-limited Emax model with partially shared parameters. E0 and EC50 are strain-specific, Emax and k are shared between strains. It returns values of Emax and k for both strains together and values of E0 and EC50 for each of the two strains. The individual fitting fits parameters for each strain&drug and media separately.
	   
- the results can be downloaded in csv file format using the buttons on the bottom of the side panel



What to do if the parametric fitting fails: Most of the implemented fitting methods (baranyi,huang,exp,Emax,sigmoid_Emax and capacity E_max) are parametric methods. This means that they require an initial estimate of the parameters in order to successfully run a fitting procedure. When the initial values are not close enough,  the fitting procedure will not finish and there will be no results or plots in step 2. This will instead show a red error message.
For the Emax models, the app is able to generate initial parameters when using the guess parameters option. However the automated guess might not always be able to provide good estimates, especially when the data does not follow a Hill shape well or if there is a lot of variability between the replicas. In that case, it is necessary to adjust the initial parameter estimates manually. Choose "no" for the guess parameter option and manually enter values. In the Tab step 2 Initial Parameters, you can see the fit of the model (red line) with the initial parameters before the fit is optimized. Try and get a good fit for the inital estimates and the rerun step 2 (clicking the run button) to see if the fitting procedure is able to finish. It might be necessary to try out different models. The range for the inital parameters estimates for the parameters of the Emax models is between -5 and 5.
There is no automated guess available for the baranyi and huang model. You can try to find inital values when first running the smoothing method which does not require inital estimates.


Negative measurement values: The app cannot process negative measurement values.  
	

