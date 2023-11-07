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
import eddymatchup_jen as ed
#from eddymatchup_jen import match_jen
import xarray as xr
from datetime import datetime

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
#eddybase_c = nc.Dataset(MAINDIR+'Data/Eddy_trajectory_nrt_3.2exp_anticyclonic_20180101_20230802.nc')
#eddybase_a = nc.Dataset(MAINDIR+'Data/Eddy_trajectory_nrt_3.2exp_cyclonic_20180101_20230802.nc')

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
   eddy = ed.match_jen( lons=lons, lats=lats, datetimes=time ,
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
