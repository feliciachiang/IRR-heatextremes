#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find heatwave event frequency based on TE 

Author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import os
import numpy as np
import xarray as xr
import more_itertools as mit

import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy as cart

#access the relevant data
datadir = "/Volumes/GISS/IRR project/figshare/"
resultsdir = "/Volumes/GISS/IRR project/results/"
figdir = "/Volumes/GISS/IRR project/figures/"

#change the working directory
os.chdir(datadir)

#for tsmax-------------------------------------------------------------------------
forcing = 'Irrig_1850'
#forcing = 'Irrig_2000half'
#forcing = 'Irrig_2000'

#open tmax data
df_TE = xr.open_dataset(datadir + "TE_day_GISS-E2-1-G_" + forcing + "_185101-198612.nc") 
df_TE = df_TE['TE'][1825:49640,:,:]

latlen = len(df_TE.lat)
lonlen = len(df_TE.lon)
yearlen = len(np.unique(df_TE.time.dt.year))

JJAyeararray = df_TE[df_TE['time.season'] == "JJA"].time.dt.year.values
JJAindinyear = (df_TE['time.season']== "JJA")[0:365].values

#import daily percentile data
quantileval = 90
df_TE_percentile = xr.open_dataset(resultsdir + "TE_fixed_" + str(quantileval) + "th_dailythresholds_15dayrollingwindow.nc")  

#create empty arrays to find the number of events for each year
JJAeventfreq = np.empty((latlen, lonlen))
JJAeventfreq[:] = np.NAN


#for each lat, lon
for latpoint in range(0, latlen):
    print(latpoint)
    for lonpoint in range(0, lonlen):
        print(lonpoint)
        #find current JJA time series for lat and lon
        currentJJA = df_TE[df_TE['time.season'] == "JJA"][:,latpoint, lonpoint].values
        
        #find threshold array for lat and lon
        currentdailythreshold = df_TE_percentile['TE_threshold'][latpoint, lonpoint, :].values
        currentJJAthreshold = currentdailythreshold[JJAindinyear]
        currentJJAthreshold = np.tile(currentJJAthreshold, yearlen)

        #find frequency, duration, intensity 
        hotoccurrences_JJA = np.arange(0,len(currentJJA))[currentJJA > currentJJAthreshold]
                
        #find total number of events during the specified time period
        JJAhotevents = [list(group) for group in mit.consecutive_groups(hotoccurrences_JJA)]
        
        #filter total number of events with length of 3 days or  longer
        JJAhoteventlenfilter = [len(list(group)) for group in mit.consecutive_groups(hotoccurrences_JJA)]
        
        JJAhoteventlenfilterindex = [n for n,i in enumerate(JJAhoteventlenfilter) if i > 2]
        
        #final events 
        JJAhoteventsfinal = []
        for lenii in range(0, len(JJAhoteventlenfilterindex)):
            JJAhoteventsfinal.append(JJAhotevents[JJAhoteventlenfilterindex[lenii]])
    
        #save event frequency
        JJAeventfreq[latpoint, lonpoint] = len(JJAhoteventsfinal)
        
#save all characteristic arrays
#save as xarray dataset                    
features = xr.Dataset({'JJAevents': (['lat', 'lon'], JJAeventfreq)},
                      coords = {'lat': df_TE.lat.values, 'lon': df_TE.lon.values})  

#export
#in output filename, include: 1) what threshold used, 2) what years included, 3) what variables used
features.to_netcdf(resultsdir + "JJA_TEevents_GISS-E2-1-G_" + forcing + "_eventfreqandbvars_dailythreshold_final.nc")      

heattype = 'TE'

TE1850mask = xr.open_dataset(resultsdir + "TEmask_Irrig_1850_JJA_day_GISS-E2-1-G_185601-198612.nc")  
TE2000mask = xr.open_dataset(resultsdir + "TEmask_Irrig_2000_JJA_day_GISS-E2-1-G_185601-198612.nc")  
TE2000halfmask = xr.open_dataset(resultsdir + "TEmask_Irrig_2000half_JJA_day_GISS-E2-1-G_185601-198612.nc")  

TE18502000mask = xr.open_dataset(resultsdir + "TEmask_18502000_JJA_day_GISS-E2-1-G_185601-198612.nc")  
TE18502000halfmask = xr.open_dataset(resultsdir + "TEmask_18502000half_JJA_day_GISS-E2-1-G_185601-198612.nc") 

#plot TE event characteristics----------------------------------------------
#import features
Irrig1850features = xr.open_dataset(resultsdir + "JJA_TEevents_GISS-E2-1-G_Irrig_1850_eventfreqandbvars_dailythreshold_final.nc")      
Irrig2000features = xr.open_dataset(resultsdir + "JJA_TEevents_GISS-E2-1-G_Irrig_2000_eventfreqandbvars_dailythreshold_final.nc")      
Irrig2000halffeatures = xr.open_dataset(resultsdir + "JJA_TEevents_GISS-E2-1-G_Irrig_2000half_eventfreqandbvars_dailythreshold_final.nc")


#plot event frequency 
forcing = 'Irrig_1850'
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = Irrig1850features['JJAevents'].where(TE1850mask == 1)['mask'].plot(vmin = 0, vmax = 300, cmap = 'YlOrBr', add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title('Irrigation = 1850')

cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Number of JJA heatwaves')
pl.savefig(figdir + 'JJATEevents_freq_' + forcing + '_dailythreshold.png', dpi = 800, bbox_inches = 'tight')

forcing = 'Irrig_2000'
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = Irrig2000features['JJAevents'].where(TE2000mask == 1)['mask'].plot(vmin =0, vmax = 300, cmap = 'YlOrBr', add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title('Irrigation = 2000')

cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Number of JJA heatwaves')
pl.savefig(figdir + 'JJATEevents_freq_' + forcing + '_dailythreshold.png', dpi = 800, bbox_inches = 'tight')

fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (Irrig2000features['JJAevents']-Irrig1850features['JJAevents']).where(TE18502000mask == 1)['mask'].plot(vmin = -120, vmax = 120, cmap = 'RdBu_r', add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title('Difference in IRR = 2000 and IRR = 1850')

cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Number of JJA heatwaves')
pl.savefig(figdir + 'JJATEevents_freq_2000diff_dailythreshold.png', dpi = 800, bbox_inches = 'tight')

forcing = 'Irrig_2000half'
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = Irrig2000halffeatures['JJAevents'].where(TE2000halfmask == 1)['mask'].plot(vmin =0, vmax = 300, cmap = 'YlOrBr', add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title('Irrigation = 2000 (Half)')

cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Number of JJA heatwaves')
pl.savefig(figdir + 'JJATEevents_freq_' + forcing + '_dailythreshold.png', dpi = 800, bbox_inches = 'tight')


fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (Irrig2000halffeatures['JJAevents']-Irrig1850features['JJAevents']).where(TE18502000halfmask == 1)['mask'].plot(vmin = -120, vmax = 120, cmap = 'RdBu_r', add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title('Difference in IRR = 2000 (Half) and IRR = 1850')

cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Number of JJA heatwaves')
pl.savefig(figdir + 'JJATEevents_freq_2000halfdiff_dailythreshold.png', dpi = 800, bbox_inches = 'tight')
