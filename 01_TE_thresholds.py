#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find the daily 90th percentile threshold (15-day centered window) of equivalent temperature 
based on the non-irrigated dataset (IRR = 1850) for each grid cell.

References:
1) Matthews et al. (2022) - Latent heat must be visible in climate communications - WIREs Climate Change
2) Pielke et al. (2004) - Assessing "global warming" with surface heat content - Eos, Transactions, AGU

Author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import os
import numpy as np
import xarray as xr

#access the relevant data
datadir = "/Volumes/GISS/IRR project/figshare/"
resultsdir = "/Volumes/GISS/IRR project/results/"

#change the working directory
os.chdir(datadir)

forcing = 'Irrig_1850'

quantileval = 90

#open HI
df_TE = xr.open_dataset(datadir + "TE_day_GISS-E2-1-G_" + forcing + "_185101-198612.nc")  

latarray = df_TE['lat'].values
lonarray = df_TE['lon'].values
dayofyear = np.unique(df_TE['time'].dt.dayofyear)

#create empty array 
TE_threshold = np.empty((len(latarray), len(lonarray), len(dayofyear)))
TE_threshold[:] = np.NAN

#remove first 5 years
rolleddata = df_TE['TE'][1825:49640,:,:].rolling(time= 15, center = True).construct('tmp')


#for each lat, lon
for latpoint in range(0, len(latarray)):
    print(latpoint)
    for lonpoint in range(0, len(lonarray)):
        print(lonpoint)
        
        currentdata = rolleddata[:,latpoint, lonpoint,:]

        for day in range(0, len(dayofyear)):
            
            TE_threshold[latpoint, lonpoint, day] = currentdata[currentdata['time'].dt.dayofyear == day +1].quantile(quantileval/100)
        
#save all characteristic arrays
#save as xarray dataset                    
features = xr.Dataset({'TE_threshold': (['lat', 'lon', 'day'], TE_threshold)},
                      coords = {'lat': latarray, 'lon': lonarray, 'day': dayofyear})  

#export
#in output filename, include: 1) what threshold used, 2) what years included, 3) what variables used
features.to_netcdf(resultsdir + "TE_fixed_" + str(quantileval) + "th_dailythresholds_15dayrollingwindow.nc")  
