#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 15:02:48 2023

@author: sarahbartoloni
"""

import pandas as pd 
import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime 


def ctr_plot (df, df_eddy, float_ID, x_var, y_var, z_var,duplicates):

    """
    Inputs
    ----------
    df (pandas.dataframe): pkl of SOCCOM data
    df_eddy (pandas.dataframe): plk of eddy data from matchups between SOCCOM data and eddy trajectories
    fload_ID (str): Cruise number of floats
    x_var (str): variable on x-axis of contour plot
    y_var (str): variable on y-axis of contour plot
    z_var (str): variable on z-axis of contour plot
    
    Returns: tricontour plot 

    """
    
    # Create SOCCOM subset for given cruise #
    flt = df.loc[df['Cruise'] == float_ID]
    Subset = flt[pd.notna(df[z_var])]
    
    # Get eddy location
    eddy_ID_cyclonic = df_eddy.loc[
        (df_eddy['eddy_type'] != 0) & 
        (df_eddy['Cruise'] == float_ID) & 
        (df_eddy['eddy_type'] == -1)
    ]
    eddy_ID_anticyclonic = df_eddy.loc[
        (df_eddy['eddy_type'] != 0) & 
        (df_eddy['Cruise'] == float_ID) & 
        (df_eddy['eddy_type'] == 1)
    ]
    
    #adjust length to appear on plot
    eddy_length_cyclonic = len(eddy_ID_cyclonic)
    eddy_y_cyclonic = [10] * eddy_length_cyclonic
    
    eddy_length_anticyclonic = len(eddy_ID_anticyclonic)
    eddy_y_anticyclonic = [10] * eddy_length_anticyclonic
    
    # Get location of eddy IDs that occure multiple times throughout float track
    #mult_eddies = SOCCOM_eddy.loc[
    #    (SOCCOM_eddy['Cruise'] == float_ID) & 
    #    ((SOCCOM_eddy['eddy_ID'] == 356580)  
    #    | (SOCCOM_eddy['eddy_ID'] == 166673 )
    #    | (SOCCOM_eddy['eddy_ID'] == 167847 )
    #    | (SOCCOM_eddy['eddy_ID'] == 171547 ))
    #]
    
    #mult_eddies_length = len(mult_eddies)
    #mult_eddies_y = [10] * mult_eddies_length
    
    
    fig, ax = plt.subplots(1,1)
    
    ct = ax.tricontourf(Subset[x_var], Subset[y_var], Subset[z_var])
    plt.gca().invert_yaxis()
    
    cbar = fig.colorbar(ct)
    cbar.set_label(z_var, rotation=270, labelpad=15)
    ax.set_xlabel(x_var)
    ax.set_ylabel(y_var)
    plt.ylim(bottom = 250)
    
    plt.scatter(
        Subset[x_var], 
        Subset['Depth[m]'], 
        s=0.25, 
        c='black', 
        alpha=0.1,
    ) 
    
    plt.scatter(eddy_ID_cyclonic[x_var], eddy_y_cyclonic, c='red', s =20)
    plt.scatter(eddy_ID_anticyclonic[x_var], eddy_y_anticyclonic, c='blue', s =20)
    #plt.scatter(mult_eddies[x_var],mult_eddies_y,c='yellow', edgecolor = 'yellow', s=20)
    
    return plt.show(), eddy_ID_cyclonic, eddy_ID_anticyclonic

