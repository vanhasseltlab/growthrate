# growthrate

Shiny app- Maximal growth rate estimator 
Link: http://132.229.100.197:2222/growthrate/

Author: Catharina Meyer
Contact: c.meyer@lacdr.leidenuniv.nl
Quantitative Pharmacology, LACDR, Leiden University

This shiny app fits time kill data in two steps. In step 1, the growth curves over time are fitted and the maximal growth rates are estimated. In step 2, the relation of maximal growth rates over drug concentration is fitted.

USER INSTRUCTIONS:

- Before using this shiny app, use the shiny app: Growth Curve for 96-Well Plate (http://132.229.100.197:2222/ot2/PlateAnalysis/GrowthCurve/) to prepocess the data.

- Upload a preprocessed data file in csv format to this shiny app

- Select a set of strain/drug/media and measurement (raw/corrected)

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
	   Provide inital estimates for the three parameters a, b & c
	2) Emax model equation: mumax ~ E0 + E_max * (((drug_concentration/EC50))/(1 + ((drug_concentration/EC50))))
	   This model requires estimates from three parameters E0, E_max and EC50 which can be automatically guessed by the shiny app.
	3) sigmoidal Emax model equation: mumax ~ E0 + E_max * (((drug_concentration/EC50)^k)/(1 + ((drug_concentration/EC50)^k)))
	   This model requires estimates from four parameters E0, E_max, EC50 and the Hill parameter k which can be automatically guessed by the shiny app.

- the results can be downloaded in csv file format using the buttons on the bottom of the side panel




	

