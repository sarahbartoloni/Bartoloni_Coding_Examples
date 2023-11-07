#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 16:04:05 2023

@author: sarahbartoloni
"""

import PyCO2SYS as pyco2
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib.path as mpath
import gsw
import matplotlib.pyplot as plt 
import matplotlib.dates as mdates 
from datetime import datetime

# ************************************************************* #
#  Converting Dickson and Millero 1987 from SW to Total pH      #
# ************************************************************* #

#calc pH at diff temp
DM87_pHT_change = pyco2.sys(par1=7.35, par2=2200, par1_type=3, par2_type=2, opt_pH_scale = 2,salinity=35,temperature = 25, 
                            temperature_out=0, opt_k_carbonic = 5, opt_k_bisulfate =1, opt_total_borate = 2)

print(DM87_pHT_change['pH_out'])

#calc TA with pHT & DIC
DM87_TA_change = pyco2.sys(par1=DM87_pHT_change['pH_out'], par2=2200, par1_type=3, par2_type=2, opt_pH_scale = 2,salinity=35,temperature = 25, 
                            temperature_out=0, opt_k_carbonic = 5, opt_k_bisulfate =1, opt_total_borate = 2)

print(DM87_TA_change['alkalinity'])

#calc pH total
DM87_pH_total = pyco2.sys(par1=DM87_TA_change['alkalinity'], par2=2200, par1_type=1, par2_type=2, opt_pH_scale = 1,salinity=35,temperature = 25, 
                            temperature_out=25, opt_k_carbonic = 5, opt_k_bisulfate =1, opt_total_borate = 2)
print(DM87_pH_total['pH_out'])

# *************************************************************** #
#  Use CO2Sys to convert pH to in-situ temp and save to 2D array  #
# *************************************************************** #

temp_outputs = [0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35]
pH_inputs = [7.2, 7.3,7.4, 7.5, 7.6, 7.7,7.8, 7.9,8,8.1,8.2]

eq_constants = 10 # Lueker
bisulfate_constant = 1 #Dickson
total_borate = 2 #Lee

#Create empty 2D arrays to store variables
insitu_pH = np.empty((len(pH_inputs),len(temp_outputs)))
dpH = np.empty((len(pH_inputs),len(temp_outputs)))
temp_constant = np.empty((len(pH_inputs),len(temp_outputs)))

temp_25 = 25

for n in range(len(pH_inputs)):
    
    #Create variables and set them to clear for each pH 
    pH = pH_inputs[n]
    output_pH = []
    pH_change = []
    pH_temp_change = []
    
    for T in temp_outputs:
        change_temp = pyco2.sys(par1=pH, 
                                par2=2200, #DIC
                                par1_type=3, 
                                par2_type=2, 
                                opt_pH_scale = 1, #total
                                salinity=35,
                                temperature = 25, 
                                temperature_out=T, 
                                opt_k_carbonic = eq_constants)
        output_pH.append(change_temp['pH_out'])
        pH_change.append(change_temp['pH_out'] - pH)
        pH_temp_change.append((change_temp['pH_out'] - pH)/(T-temp_25))
        
    insitu_pH[n,:] = output_pH
    dpH[n,:] = pH_change
    temp_constant[n,:] = pH_temp_change
        
print(temp_constant)


# *************************************************************** #
#  Create Models for ∆pH                                          #
# *************************************************************** #

#Attempting quatradic first (dpH vs pH)
polyline = np.linspace(7.2,8.2, 100)

model_0 = np.poly1d(np.polyfit(pH_inputs, dpH[:,0],2))
model_5 = np.poly1d(np.polyfit(pH_inputs, dpH[:,5],2))
model_10 = np.poly1d(np.polyfit(pH_inputs, dpH[:,10],2))

#figure(1)
plt.scatter(pH_inputs, dpH[:,0])
plt.scatter(pH_inputs,dpH[:,5])
plt.scatter(pH_inputs,dpH[:,10])
plt.plot(polyline, model_0(polyline), polyline, model_5(polyline), polyline, model_10(polyline))
print(model_0)
print(model_5)
print(model_10)


# *************************************************************** #
#  Create Function to return statistics                           #
# *************************************************************** #

import math

def stats(y_pred, y_true):
    """ 
    - y_pred: Predicted values
    - y_true: Actaul values
    """
    #Residuals 
    R = y_true - y_pred
    
    #Residual sum of squares (RSS)
    RSS = np.sum(np.square(y_pred - y_true))
    
    #Residual standard error (RSE)
    RSE = math.sqrt(RSS/(len(y_true)-2))
    
    #Mean square of errors (MSE)
    MSE = np.square(np.subtract(y_true, y_pred)).mean()
    
    #Root mean squared error (RMSE)
    RMSE = math.sqrt(MSE)
    
    return R, RSS, RSE, MSE, RMSE
