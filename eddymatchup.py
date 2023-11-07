#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 14:52:49 2023

@author: sarahbartoloni
"""

### Set Paths
from sys import exit
import re
import os
import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.path as mpath
import netCDF4 as nc
import eddymatchup_fxn as ed
import xarray as xr
from datetime import datetime
import cartopy
import cartopy.crs as ccrs
import Contour_plot 
from Contour_plot import ctr_plot
import gsw
import matplotlib.dates as mdates 
import Contour_plot_fxn as ctr_plot

HOMEDIR = os.getcwd() 
MAINDIR = '/Users/sarahbartoloni/Desktop/GitHub/Saildrone/' 
OUTDIR = MAINDIR+'OUTPUT/'
DATADIR = MAINDIR+'Data/'
eddymatchup = MAINDIR

# ************************************************************* #
#  Figure out when SD is in an eddy, based on an eddy database  #
# ************************************************************* #

# Need databases of "eddybase_c" and "eddybase_a"
eddybase_c = MAINDIR+'Data/Eddy_trajectory_nrt_3.2exp_anticyclonic_20180101_20230802.nc'
eddybase_a = MAINDIR+'Data/Eddy_trajectory_nrt_3.2exp_cyclonic_20180101_20230802.nc'

# Load pickle of BGC Argo Data
argo = pd.read_pickle(MAINDIR+'u_stat_SOCCOM.pkl') 
argo = argo.to_numpy()
lons = argo[:,5]
lats = argo[:,6]

# Put into correct datetime format
datetimes = []
for i in range(len(argo[:,3])):
    d_t = argo[i,3] + argo[i,4]
    datetimes.append(datetime.strptime(d_t, '%m/%d/%Y%H:%M'))

datetimes = np.asarray(datetimes)

d_t = []
for i in datetimes:
    dt = np.datetime64(i)
    d_t.append(dt)
time = np.asarray(d_t)


do_eddy_matchup = 1
if do_eddy_matchup == 1:
   print('\n================= EDDY MATCHUP REOUTINE =================== ')
   eddy = ed( lons=lons, lats=lats, datetimes=time ,
                      filename_cyclonic=eddybase_c, 
                      filename_anticyclonic=eddybase_a,
                      latmin=-70, latmax=-30, lonmin=0, lonmax=100,
                      hourrange=24, radiusrange=1.5)
   print('\nEddy matchups complete!')
   eddy.to_netcdf(path=eddymatchup)
   print('Saved to ', eddymatchup)
   eddy.close()
   print('')
else:
   eddy = xr.load_dataset(eddymatchup)

#Convert back to numpy
eddy_type = eddy['eddy_type'].to_numpy()
eddy_ID   = eddy['eddy_ID'].to_numpy()
eddy_lat  = eddy['eddy_lat'].to_numpy()
eddy_lon  = eddy['eddy_lon'].to_numpy()
eddy_time = eddy['eddy_time'].to_numpy()
eddy_amp  = eddy['eddy_amplitude'].to_numpy()
eddy_vmax = eddy['eddy_vmax'].to_numpy()
eddy_rad_to_vmax  = eddy['eddy_rad_to_vmax'].to_numpy()
eddy_age  = eddy['eddy_age'].to_numpy()
eddy_dist_to_ctr  = eddy['eddy_dist_to_ctr'].to_numpy()



# ************************************************************* #
#  Load Orsi fronts                                             #
# ************************************************************* #

def load_orsi_fronts(DATADIR):
    # Import the Southern Ocean fronts for mapping
    stf = pd.read_csv(DATADIR + 'fronts/stf.txt', header=None, sep='\s+', na_values='%', names=['lon','lat'])
    saf = pd.read_csv(DATADIR + 'fronts/saf.txt', header=None, sep='\s+', na_values='%', names=['lon','lat'])
    pf = pd.read_csv(DATADIR + 'fronts/pf.txt', header=None, sep='\s+', na_values='%', names=['lon','lat'])
    saccf = pd.read_csv(DATADIR + 'fronts/saccf.txt', header=None, sep='\s+', na_values='%', names=['lon','lat'])
    sbdy = pd.read_csv(DATADIR + 'fronts/sbdy.txt', header=None, sep='\s+', na_values='%', names=['lon','lat'])
    
    return stf,saf,pf,saccf,sbdy

stf, saf, pf, saccf, sbdy = load_orsi_fronts(DATADIR)

# ************************************************************* #
#  Single Float Check                                           #
# ************************************************************* #

#Define SOCCOM variables in simplier terms
DIC = 'DIC_LIAR[µmol/kg]'
TA = 'TALK_LIAR[µmol/kg]'
salinity = 'Salinity[pss]'
temperature = 'Temperature[°C]'
depth = 'Depth[m]'
oxygen = 'Oxygen[µmol/kg]'
oxygen_sat = 'OxygenSat[%]'
nitrate = 'Nitrate[µmol/kg]'
chl_a = 'Chl_a_corr[mg/m^3]'
pH = 'pHinsitu[Total]'
pCO2 = 'pCO2_LIAR[µatm]'
lat = 'Lat [°N]'
lon = 'Lon [°E]'


#Single float check: Does float intersect a single eddie mlt times?
float_ID = 5904686                                                                                                                                                                                                               
x_var = 'Station'
y_var = 'Depth[m]'
z_var = oxygen
duplicates = []

test = ctr_plot(argo, eddy, float_ID, x_var, y_var, z_var, duplicates ==True)

#Cylonic = red, Anticyclonic = blue
print(eddy['Station'].loc[eddy_ID == 182518])


# ************************************************************* #
#  Visualize Entire SOCCOM BGC-Argo Float dataset               #
# ************************************************************* #

eddy_ID_cyclonic = eddy.loc[(eddy_type != 0) 
                                   & (eddy_type == -1)
                                  ]
eddy_ID_anticyclonic = eddy.loc[(eddy_type != 0) 
                                   & (eddy_type == 1)
                                  ]
    
plt.figure(figsize =(12,12))
ax = plt.axes(projection=ccrs.SouthPolarStereo())
ax.set_extent([-180,180,-90,-40], ccrs.PlateCarree())
ax.add_feature(cartopy.feature.LAND)
ax.add_feature(cartopy.feature.OCEAN)
ax.gridlines(ccrs.PlateCarree(), draw_labels=True)

theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

ax.set_boundary(circle, transform = ax.transAxes)
plt.plot(stf['lon'], stf['lat'], color='Red', transform=ccrs.PlateCarree(), linewidth=.75)
plt.plot(saf['lon'], saf['lat'], color='Orange', transform=ccrs.PlateCarree(), linewidth=.75)
plt.plot(pf['lon'], pf['lat'], color='Yellow', transform=ccrs.PlateCarree(), linewidth=.75)
plt.plot(saccf['lon'], saccf['lat'], color='Green', transform=ccrs.PlateCarree(), linewidth=.75)
plt.plot(sbdy['lon'], sbdy['lat'], color='Blue', transform=ccrs.PlateCarree(), linewidth=.75)


plt.scatter(argo[lon], argo[lat], color = 'Black', edgecolor = 'black', transform=ccrs.PlateCarree(), s=1, zorder=1001)
plt.scatter(eddy_ID_cyclonic[lon], eddy_ID_cyclonic[lat],c='red', edgecolor = 'red', transform=ccrs.PlateCarree(), s=1, zorder=1001)
plt.scatter(eddy_ID_anticyclonic[lon], eddy_ID_anticyclonic[lat],c='blue', edgecolor = 'blue', transform=ccrs.PlateCarree(), s=1, zorder=1001)


# ********************************************************************** #
#  Visualize all float & eddy matchups in region of Saildrone Mission    #
# ********************************************************************** #


# Set Figure and Map Characteristics
plt.figure()
plt.figure(figsize =(12,6))
plt.title('All Eddy Matchups')
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([0,140,-70,-30], ccrs.PlateCarree())
ax.coastlines()
ax.add_feature(cartopy.feature.LAND)
ax.add_feature(cartopy.feature.OCEAN)
ax.gridlines()

# Plot Fronts (uses circle axes from ARGO Float Track)
#ax.set_boundary(circle, transform=ax.transAxes)
plt.plot(stf['lon'], 
         stf['lat'], 
         color='Red', 
         transform=ccrs.PlateCarree(),
        label='Subtropical Front')
plt.plot(saf['lon'], 
         saf['lat'], 
         color='Orange', 
         transform=ccrs.PlateCarree(),
        label='Subantarctic Front')
plt.plot(pf['lon'], 
         pf['lat'], 
         color='Yellow', 
         transform=ccrs.PlateCarree(),
        label='Polar Front')
plt.plot(saccf['lon'], 
         saccf['lat'], 
         color='Green', 
         transform=ccrs.PlateCarree(),
        label='Southern ACC Front')
plt.plot(sbdy['lon'], 
         sbdy['lat'], 
         color='Blue', 
         transform=ccrs.PlateCarree(),
        label='Southern Boundary')
#plt.legend()

# Plot Argo Float Track
plt.scatter(argo[lon], 
            argo[lat], 
            color = 'Black', 
            edgecolor = 'black', 
            transform=ccrs.PlateCarree(), 
            s=1, 
            zorder=1001)
plt.scatter(eddy_ID_cyclonic[lon], 
            eddy_ID_cyclonic[lat],
            c='red', 
            edgecolor = 'red', 
            transform=ccrs.PlateCarree(), 
            s=1, 
            zorder=1001)
plt.scatter(eddy_ID_anticyclonic[lon], 
            eddy_ID_anticyclonic[lat],
            c='blue', 
            edgecolor = 'blue', 
            transform=ccrs.PlateCarree(), 
            s=1, 
            zorder=1001)


# ********************************************************************** #
#  Separate all float and eddie matchups into cyclonic and anticyclonic  #
# ********************************************************************** #

# Separate above into cyclonic and anticyclonic

fig = plt.figure(figsize=(12, 6))
plt.title('Cyclonic Eddies')
ax1 = plt.axes(projection=ccrs.PlateCarree())
ax1.set_extent([0,140,-70,-30], ccrs.PlateCarree())
ax1.coastlines()
ax1.add_feature(cartopy.feature.LAND)
ax1.add_feature(cartopy.feature.OCEAN)
ax1.gridlines()

# Plot Fronts (uses circle axes from ARGO Float Track)
#ax.set_boundary(circle, transform=ax.transAxes)
plt.plot(stf['lon'], 
         stf['lat'], 
         color='Red', 
         transform=ccrs.PlateCarree(),
        label='Subtropical Front')
plt.plot(saf['lon'], 
         saf['lat'], 
         color='Orange', 
         transform=ccrs.PlateCarree(),
        label='Subantarctic Front')
plt.plot(pf['lon'], 
         pf['lat'], 
         color='Yellow', 
         transform=ccrs.PlateCarree(),
        label='Polar Front')
plt.plot(saccf['lon'], 
         saccf['lat'], 
         color='Green', 
         transform=ccrs.PlateCarree(),
        label='Southern ACC Front')
plt.plot(sbdy['lon'], 
         sbdy['lat'], 
         color='Blue', 
         transform=ccrs.PlateCarree(),
        label='Southern Boundary')
#plt.legend()

# Plot ARGO Float Track
plt.scatter(eddy_ID_cyclonic[lon], eddy_ID_cyclonic[lat],c='red', edgecolor = 'red', transform=ccrs.PlateCarree(), s=1, zorder=1001)

fig2 = plt.figure(figsize=(12, 6))
plt.title('Antiyclonic Eddies')
ax2 = plt.axes(projection=ccrs.PlateCarree())
ax2.set_extent([0,140,-70,-30], ccrs.PlateCarree())
ax2.coastlines()
ax2.add_feature(cartopy.feature.LAND)
ax2.add_feature(cartopy.feature.OCEAN)
ax2.gridlines()

# Plot Fronts (uses circle axes from ARGO Float Track)
#ax.set_boundary(circle, transform=ax.transAxes)
plt.plot(stf['lon'], 
         stf['lat'], 
         color='Red', 
         transform=ccrs.PlateCarree(),
        label='Subtropical Front')
plt.plot(saf['lon'], 
         saf['lat'], 
         color='Orange', 
         transform=ccrs.PlateCarree(),
        label='Subantarctic Front')
plt.plot(pf['lon'], 
         pf['lat'], 
         color='Yellow', 
         transform=ccrs.PlateCarree(),
        label='Polar Front')
plt.plot(saccf['lon'], 
         saccf['lat'], 
         color='Green', 
         transform=ccrs.PlateCarree(),
        label='Southern ACC Front')
plt.plot(sbdy['lon'], 
         sbdy['lat'], 
         color='Blue', 
         transform=ccrs.PlateCarree(),
        label='Southern Boundary')
#plt.legend()

# Plot Argo Float Track
plt.scatter(eddy_ID_cyclonic[lon], 
            eddy_ID_cyclonic[lat],
            c='red', 
            edgecolor = 'red', 
            transform=ccrs.PlateCarree(), 
            s=1, 
            zorder=1001)
plt.scatter(eddy_ID_anticyclonic[lon], 
            eddy_ID_anticyclonic[lat],
            c='blue', 
            edgecolor = 'blue', 
            transform=ccrs.PlateCarree(), 
            s=1, 
            zorder=1001)
