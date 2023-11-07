#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 14:44:31 2023

@author: sarahbartoloni
"""

# -*- coding: utf-8 -*-
# ### function to match up oceanographic data with Chelton eddy database
# Author: Veronica Tamsitt, August 2022
# Altered by Jennifer Bonin, May 2023  
# Altered by Sarah Bartoloni, August 2023
#
# Reads in a pre-downloaded eddy trajectory atlas, and matches input position (latitude, longitude, time) data from User to eddy database data for cyclonic and anticyclonic eddies. 
# Returns an xarray array "matchups" which contains the relevant eddy info at each time along the track (eddy types 1, 0, -1, etc) 
#
# #### META3.2_allsat (default)
# https://www.aviso.altimetry.fr/en/data/products/value-added-products/global-mesoscale-eddy-trajectory-product/meta3-2-dt.html
# Reference: Pegliasco, C., Delepoulle, A., Mason, E., Morrow, R., Faugère, Y., Dibarboure, G., 2022. 
# META3.1exp: a new global mesoscale eddy trajectory atlas derived from altimetry. Earth Syst. Sci. Data 14, 1087–1107. https://doi.org/10.5194/essd-14-1087-2022 
# Citation: The altimetric Mesoscale Eddy Trajectories Atlas (META3.2 DT) was produced by SSALTO/DUACS and distributed by AVISO+ (https://aviso.altimetry.fr ) 
# with support from CNES, in collaboration with IMEDEA (DOI: 10.24400/527896/a01-2022.005.210802 for the META3.2 DT allsat version  and 10.24400/527896/a01-2022.006.210802  
# for the META3.2 DT twosat version).
#

#import dependent modules
import numpy as np
import xarray as xr
import datetime
import glob, os
from datetime import datetime


def match_jen(lons,lats,datetimes, filename_cyclonic,filename_anticyclonic, latmin=-90,latmax=-35,lonmin=0,lonmax=360, hourrange=24,radiusrange=1):
    
    """
    Inputs:
    lons (array of floats): required; longitude of data to match up with eddy database, -180 to 180 or 0 to 360
    lats (array of floats): required; latitude of data to match up with eddy database, -90 to 90
    datetimes (array of python datetime): required; dates and times of data to match up with eddy database
    database directory (string): location of META3.2_allsat NRT eddy tracking database/algorithm to use 
    latmin/max, lonmin/max (float): optional; lat/lon bounds to consider for matchups, default is entire southern ocean 
    hourrange (int): optional; maximum time difference (in hours) to consider for eddy matchup, default is +/- 24 hours
    radiusrange (int or float): optional; maximum distance from eddy centre to include in matchups, 
        such that the max distance considered is radiusrange*(radius from eddy centre to maximum speed). Default value is 1.
    """
    
    #Load the eddy databases as xarrays 
    ed_c = xr.load_dataset(filename_cyclonic)
    ed_ac = xr.load_dataset(filename_anticyclonic)
    print('Databases read in.')
    


    #reduce eddy database to south of latmax, within time frame of datetimes
    datetimes_n = datetimes
    datemin = datetimes_n.min()
    datemax = datetimes_n.max()

    #check longitude go from 0 to 360
    lons_pos = lons.copy() 
    if np.any(lons_pos<0):
      lons_pos[lons_pos<0] = lons_pos[lons_pos<0]+360.
        
    eddy_c = ed_c.where((ed_c.time>datemin) & (ed_c.time<datemax),drop=True) 
    eddy_c = eddy_c.where((eddy_c.latitude>latmin) & (eddy_c.latitude<latmax),drop=True)
    eddy_c = eddy_c.where((eddy_c.longitude>lonmin) & (eddy_c.longitude<lonmax),drop=True)

    eddy_ac = ed_ac.where((ed_ac.time>datemin) & (ed_ac.time<datemax),drop=True) 
    eddy_ac = eddy_ac.where((eddy_ac.latitude>latmin) & (eddy_ac.latitude<latmax),drop=True)
    eddy_ac = eddy_ac.where((eddy_ac.longitude>lonmin) & (eddy_ac.longitude<lonmax),drop=True)

    #loop through all obs data times and compare locations to eddy database
    dt = np.timedelta64(hourrange,'h')       #set time window to search for eddies   
    match_type = np.zeros((len(datetimes)))
    match_track = np.empty((len(datetimes)))
    match_track[:] = np.nan
    match_time = np.copy(match_track)
    match_lat = np.copy(match_track)
    match_lon = np.copy(match_track)
    match_amp = np.copy(match_track)
    match_speed = np.copy(match_track)
    match_radius = np.copy(match_track)
    match_age = np.copy(match_track)
    match_dist = np.copy(match_track)
    dr = 2*np.pi*(6371.136*1000)/360.0  #define constant for calculating distance 
    
    for t in range(len(datetimes)):
      #find all eddies within +/- dt of the obs time
      eddy_c_day = eddy_c.where((eddy_c.time>datetimes_n[t]-dt) & (eddy_c.time<datetimes_n[t]+dt),drop=True)
      eddy_ac_day = eddy_ac.where((eddy_ac.time>datetimes_n[t]-dt) & (eddy_ac.time<datetimes_n[t]+dt),drop=True)
            
      #compute distance from eddy centre lats/lons
      dy = (eddy_c_day.latitude.values-lats[t])*dr
      dx = (eddy_c_day.longitude.values-lons_pos[t])*dr*np.cos(np.deg2rad(lats[t]))
      dists = np.sqrt(dx*dx + dy*dy)
      eddy_c_day['distance'] = xr.DataArray(dists,dims=['obs'])
        
      #repeat for anticyclonic
      dy = (eddy_ac_day.latitude.values-lats[t])*dr
      dx = (eddy_ac_day.longitude.values-lons_pos[t])*dr*np.cos(np.deg2rad(lats[t]))
      dists = np.sqrt(dx*dx + dy*dy)
      eddy_ac_day['distance'] = xr.DataArray(dists,dims=['obs'])
      
      print(eddy_c_day)
      #find if distance less than radiusrange*eddy radius                 
      if np.any(eddy_c_day.distance < radiusrange*eddy_c_day.speed_radius) or \
         np.any(eddy_ac_day.distance < radiusrange*eddy_ac_day.speed_radius): 
                
         eddy_c_match = eddy_c_day.where(eddy_c_day.distance<(radiusrange*eddy_c_day.speed_radius),drop=True)
         eddy_ac_match = eddy_ac_day.where(eddy_ac_day.distance<(radiusrange*eddy_ac_day.speed_radius),drop=True)
                              
         #if 1 match save it 
         if len(eddy_c_match.obs)+len(eddy_ac_match.obs) == 1:
           if len(eddy_c_match.obs):
             match_type[t] = -1 #-1 is cyclonic
             match_track[t] = eddy_c_match.track
             match_time[t] = eddy_c_match.time
             match_lat[t] = eddy_c_match.latitude
             match_lon[t] = eddy_c_match.longitude
             match_amp[t] = eddy_c_match.amplitude
             match_speed[t] = eddy_c_match.speed_average
             match_radius[t] = eddy_c_match.speed_radius
             match_age[t] = eddy_c_match.observation_number
             match_dist[t] = eddy_c_match.distance
           elif len(eddy_ac_match.obs):
             match_type[t] = 1 #-1 is cyclonic
             match_track[t] = eddy_ac_match.track
             match_time[t] = eddy_ac_match.time
             match_lat[t] = eddy_ac_match.latitude
             match_lon[t] = eddy_ac_match.longitude
             match_amp[t] = eddy_ac_match.amplitude
             match_speed[t] = eddy_ac_match.speed_average
             match_radius[t] = eddy_ac_match.speed_radius
             match_age[t] = eddy_ac_match.observation_number
             match_dist[t] = eddy_ac_match.distance
         #else >1 match find closest match and save it
         elif len(eddy_c_match.obs)>=1 and len(eddy_ac_match.obs)>=1:
           closest_ind_c = eddy_c_match.distance.argmin()
           closest_ind_ac = eddy_ac_match.distance.argmin()
           if eddy_c_match.distance[closest_ind_c]<eddy_ac_match.distance[closest_ind_ac]:
             match_type[t] = -1
             match_track[t] = eddy_c_match.track[closest_ind_c]
             match_time[t] = eddy_c_match.time[closest_ind_c]
             match_lat[t] = eddy_c_match.latitude[closest_ind_c]
             match_lon[t] = eddy_c_match.longitude[closest_ind_c]
             match_amp[t] = eddy_c_match.amplitude[closest_ind_c]
             match_speed[t] = eddy_c_match.speed_average[closest_ind_c]
             match_radius[t] = eddy_c_match.speed_radius[closest_ind_c]
             match_age[t] = eddy_c_match.observation_number[closest_ind_c]
             match_dist[t] = eddy_c_match.distance[closest_ind_c]
           else:
             match_type[t] = 1
             match_track[t] = eddy_ac_match.track[closest_ind_ac]
             match_time[t] = eddy_ac_match.time[closest_ind_ac]
             match_lat[t] = eddy_ac_match.latitude[closest_ind_ac]
             match_lon[t] = eddy_ac_match.longitude[closest_ind_ac]
             match_amp[t] = eddy_ac_match.amplitude[closest_ind_ac]
             match_speed[t] = eddy_ac_match.speed_average[closest_ind_ac]
             match_radius[t] = eddy_ac_match.speed_radius[closest_ind_ac]
             match_age[t] = eddy_ac_match.observation_number[closest_ind_ac]
             match_dist[t] = eddy_ac_match.distance[closest_ind_ac]
         elif len(eddy_c_match.obs)>=1:
             closest_ind_c = eddy_c_match.distance.argmin()
             match_type[t] = -1
             match_track[t] = eddy_c_match.track[closest_ind_c]
             match_time[t] = eddy_c_match.time[closest_ind_c]
             match_lat[t] = eddy_c_match.latitude[closest_ind_c]
             match_lon[t] = eddy_c_match.longitude[closest_ind_c]
             match_amp[t] = eddy_c_match.amplitude[closest_ind_c]
             match_speed[t] = eddy_c_match.speed_average[closest_ind_c]
             match_radius[t] = eddy_c_match.speed_radius[closest_ind_c]
             match_age[t] = eddy_c_match.observation_number[closest_ind_c]
             match_dist[t] = eddy_c_match.distance[closest_ind_c]
         elif len(eddy_ac_match.obs)>=1:
             closest_ind_ac = eddy_ac_match.distance.argmin()
             match_type[t] = 1
             match_track[t] = eddy_ac_match.track[closest_ind_ac]
             match_time[t] = eddy_ac_match.time[closest_ind_ac]
             match_lat[t] = eddy_ac_match.latitude[closest_ind_ac]
             match_lon[t] = eddy_ac_match.longitude[closest_ind_ac]
             match_amp[t] = eddy_ac_match.amplitude[closest_ind_ac]
             match_speed[t] = eddy_ac_match.speed_average[closest_ind_ac]
             match_radius[t] = eddy_ac_match.speed_radius[closest_ind_ac]
             match_age[t] = eddy_ac_match.observation_number[closest_ind_ac]
             match_dist[t] = eddy_ac_match.distance[closest_ind_ac]  
                          
    #index of times with obs matchups
    match_ind = np.flatnonzero(match_type)     #gives the array index values for non-zero match types only 
    nmatches = len(match_ind)
    print('Total eddy matchups = ' + str(nmatches))
    unique_eddies,unique_ind = np.unique(match_track[match_ind],return_index=True)
    neddies = len(unique_eddies)
    print('Total unique eddy IDs = ' + str(neddies))
 
    #return matchups as xarray Dataset to be merged with original data (use data time as dimensions?)
    # define data with variable attributes
    data_vars = dict( eddy_type = (["obs"],match_type),
                      eddy_ID = (["obs"],match_track),
                      eddy_lat = (["obs"],match_lat),
                      eddy_lon = (["obs"],match_lon),
                      eddy_time = (["obs"],match_time),
                      eddy_amplitude = (["obs"],match_amp),
                      eddy_vmax = (["obs"],match_speed),
                      eddy_rad_to_vmax = (["obs"],match_radius),
                      eddy_age = (["obs"],match_age),
                      eddy_dist_to_ctr = (["obs"],match_dist) )
                    
    # define coordinates
    coords = dict( time = (["obs"],datetimes),
                   lat = (["obs"],lats),
                   lon = (["obs"],lons) )

    # define global attributes
    #attrs = dict(creation_date = datetime.datetime.now())

    # create dataset
    matchups = xr.Dataset(data_vars=data_vars, coords=coords)   ##, attrs=attrs)
    
    return matchups     

