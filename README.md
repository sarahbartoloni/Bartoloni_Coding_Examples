# Bartoloni_Coding_Examples
Examples of code written by Sarah Bartoloni

This repository is intended for viewing purposes only. It represents a compilation of projects I have worked on throughout my time as a graduate student. The codes here are not intended to run outside of their private repositories.

 **Argo_SO_visualization**

This script is a compilation of my work analyzing BGC-Argo floats in the Southern Ocean. This script was done in collaboration with both Dr. Nancy William's lab and Dr. Don Chamber's lab and is still a work in progress. First, matchups between Argo floats and mesoscale eddies are identified using the eddymatchup_fxn initially created by Dr. Veronica Tamsitt and then modified by Dr. Jennifer Bonin and myself. The code then examines if a single float encounters the same eddy multiple times throughout its path. The entire BGC-Argo fleet in the Southern Ocean is then visualized in relation to fronts (from Orsie et al 1995) and mesoscale eddy locations. The code then separates the data depending on season, zonal regions, and the presence of anticyclonic/cyclonic eddies. Finally, a function is created to interpolate Argo profiles onto the same pressure levels.

 **CO2_Examples**

Here is an example of code used to process laboratory CO2 system data. This is a collection of code used among a variety of wet lab projects. First, pyCO2SYS is used to convert pH data using equilibrium constants from Dickson and Millero 1987 to total scale pH. Next, CO2Sys is used to convert pH, at in-situ temperatures, to a temperature of 25°C and then saved into a 2D array. Finally, a quadratic model is created to model pH_25 - pH_insitu (∆pH). It should be noted, I perform most modeling and statistical analysis in R.
 
 **Contour_plot_fxn**

A function was created to make a contour plot using any three variables in the Argo float dataset. This function is used in Argo_SO_visualization.

 **eddmatchup_fxn**

This script is a function used to determine when BGC-Argo floats intersect a mesoscale eddy. This function was initially created by Dr. Veronica Tamsitt and then modified by Dr. Jennifer Bonin and myself.

