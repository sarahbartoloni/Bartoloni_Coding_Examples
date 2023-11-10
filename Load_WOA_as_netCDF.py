#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 16 11:59:30 2023

@author: sarahbartoloni
"""

# Import operations
#import os
import re
import os
import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.path as mpath
import netCDF4 as nc

### Set the paths
HOMEDIR = os.getcwd() 
MAINDIR = '/Users/sarahbartoloni/Desktop/Bartoloni_Chambers_Argo/' 
OUTDIR = MAINDIR+'OUTPUT/'
DATADIR = MAINDIR+'Data/'
WOA_T1 = nc.Dataset(DATADIR+'WOA2018/t_nc/woa18_A5B7_t01_01.nc')

##check directories exist
if not os.path.isdir(OUTDIR):
    os.mkdir(OUTDIR)
    print("Made new output dir ", OUTDIR) 
if not os.path.isdir(DATADIR):
    os.mkdir(DATADIR)
    print("Made new data dir ", DATADIR) 
    
    
# ************************************************************* #  
#  Manually read in Lat and Lon from WOA                        #                
#  Note:  Data is in Numpy arrays                               # 
# ************************************************************* #      

initlat = WOA_T1.variables['lat'][:]
initlon = WOA_T1.variables['lon'][:]
initdepth = WOA_T1.variables['depth'][:]

def lon_to_ilon(initial_lon):
    nlon = []
    for i in initial_lon:
        if i < 0:
            i = i+ 360
            nlon.append(i)
        elif i >= 0:
            nlon.append(i)
    del initial_lon
    #print(nlon)
    
    ilon = np.zeros([0,1], dtype=float)
    for j in nlon:
        if j > 180:
            j = j - 180.5
            ilon = np.append(ilon,j)
        elif j < 180:
            j = j +179.5
            ilon = np.append(ilon,j)
    return ilon

def lat_to_ilat(initial_lat):
    ilat = []
    for i in initial_lat:
        i = i + 89.5
        ilat.append(i)
    del initial_lat
    return ilat

def depth_to_idepth(initial_depth):
    idepth = np.zeros([0,1], dtype=float)
    for i in initial_depth:
        if i <= 100:
            i = i/5
            idepth = np.append(idepth,i)
        if i >=125 and i <=500:
            i = ((i - 100)/25) + 20
            idepth = np.append(idepth,i)
        if i >= 550 and i <= 1500:
            i = ((i - 500)/50) + 36
            idepth = np.append(idepth,i)
    del initial_depth
    return idepth

lon = lon_to_ilon(initlon)
lat = lat_to_ilat(initlat)
depth = depth_to_idepth(initdepth)

# ************************************************************* #  
# Create WOA grids for Temp, Oxygen, and Nitrate                # 
# ************************************************************* #   

#Create max vars (data type = int)
ntime = 12
ndepth = len(depth) #57
nlat = len(lat)     #180
nlon = len(lon)     #360

temp = np.empty([ntime, ndepth, nlat, nlon])
sal = np.empty([ntime, ndepth, nlat, nlon])
ox = np.empty([ntime, ndepth, nlat, nlon]) 
nitrate = np.empty([ntime, ndepth, nlat, nlon]) 

for imon in range (0,11,1):
    i = imon + 1
    if i < 10:
        #Temperature
        t = nc.Dataset(DATADIR+'WOA2018/t_nc/woa18_A5B7_t0' 
                         + str(i) + '_01.nc')
        itemp = t.variables['t_an'][:]
        temp[i,0:57,0:180,0:360]=itemp[0,:,:,:]
        
        #Oxygen
        o = nc.Dataset(DATADIR+'WOA2018/o_nc/woa18_all_o0' 
                         + str(i) + '_01.nc')
        iox = o.variables['o_an'][:]
        ox[i,0:57,0:180,0:360]=iox[0,:,:,:]
        
        #Salinity
        s = nc.Dataset(DATADIR+'WOA2018/s_nc/woa18_A5B7_s0' 
                         + str(i) + '_01.nc')
        isal = s.variables['s_an'][:]
        sal[i,0:57,0:180,0:360]=isal[0,:,:,:]
        
        #Nitrate
        n = nc.Dataset(DATADIR+'WOA2018/n_nc/woa18_all_n0' 
                         + str(i) + '_01.nc')
        initrate = n.variables['n_an'][:]
        nitrate[i,0:43,0:180,0:360]=initrate[0,:,:,:]
        
    if i >=10:
        #Temperature
        t = nc.Dataset(DATADIR+'WOA2018/t_nc/woa18_A5B7_t' 
                       + str(i) + '_01.nc')
        itemp = t.variables['t_an'][:]
        temp[i,0:57,0:180,0:360]=itemp[0,:,:,:]
        
        #Oxygen
        o = nc.Dataset(DATADIR+'WOA2018/o_nc/woa18_all_o' 
                       + str(i) + '_01.nc')
        iox = o.variables['o_an'][:]
        ox[i,0:57,0:180,0:360]=iox[0,:,:,:]
        
        #Salinity
        s = nc.Dataset(DATADIR+'WOA2018/s_nc/woa18_A5B7_s' 
                       + str(i) + '_01.nc')
        isal = s.variables['s_an'][:]
        sal[i,0:57,0:180,0:360]=isal[0,:,:,:]
        
        #Nitrate
        n = nc.Dataset(DATADIR+'WOA2018/n_nc/woa18_all_n' 
                          + str(i) + '_01.nc')
        initrate = n.variables['n_an'][:]
        nitrate[i,0:43,0:180,0:360]=initrate[0,:,:,:]


# ************************************************************* #  
# Save varaible arrays as npy files                             # 
# ************************************************************* #   

np.save(DATADIR+ 'nitrate.npy', nitrate)
#Use np.load to load into new file


