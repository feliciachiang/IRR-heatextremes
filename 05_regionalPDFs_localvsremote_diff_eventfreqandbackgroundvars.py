#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create individual differenced (between irrigated and non-irrigated runs) regional area-weighted PDFs for local and remote grid cells
    1. Find grid cells within each SREX region
    2. Find data values and weights for each region
    3. Create weighted PDFs for each event frequency and background variables

Author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import os
import numpy as np
import xarray as xr
import pandas as pd
import re
import fnmatch
import scipy.stats

import matplotlib.pylab as pl
import matplotlib.pyplot as plt
from calendar import monthrange
import seaborn as sns


datadir = "/Volumes/GISS/IRR project/data/"
resultsdir = "/Users/fchiang/GISS/IRR project/results/"
resourcesdir = "/Volumes/GISS/resources/"
figdir = "/Volumes/GISS/IRR project/figures/"

noIRRforcing = 'Irrig_1850'
yesIRRforcing = 'Irrig_2000'
#yesIRRhalfforcing = 'Irrig_2000half'

monthweights = [monthrange(2001, 1)[1]/365, monthrange(2001, 2)[1]/365, monthrange(2001, 3)[1]/365, monthrange(2001, 4)[1]/365, monthrange(2001, 5)[1]/365, monthrange(2001, 6)[1]/365, monthrange(2001, 7)[1]/365, monthrange(2001, 8)[1]/365, monthrange(2001, 9)[1]/365, monthrange(2001, 10)[1]/365, monthrange(2001, 11)[1]/365, monthrange(2001, 12)[1]/365]


def kde_preprocessing(data_array, latitude_array, mask_array, sftlf_array):
#function to: 
    #1) input array (e.g. tsmax feature array)
    #2) input latitude information
    #3) input mask to find grid cells 
    #4) output data points and weights

    #establish grid cells to extract
    latvals = latitude_array[np.where(mask_array == True)[0]]
    #lonvals = longitude_array[np.where(mask_array == True)[1]]
    
    #establish weights
    weightarray = np.cos(np.deg2rad(latvals))
    
    #find data using data_array
    datavals = data_array[np.where(mask_array == True)[0], np.where(mask_array == True)[1]]
    
    #find sftlf data 
    sftlfvals = sftlf_array[np.where(mask_array == True)[0], np.where(mask_array == True)[1]]

    weightarray_sftlf = weightarray*sftlfvals
    
    finalweightarray = weightarray_sftlf/np.sum(weightarray_sftlf)

    returndata = [datavals, finalweightarray]
    
    return returndata

#new longitude and latitude
newlon = np.arange(-179.75, 180.25, 0.5)
newlat = np.arange(-89.75,90.25, 0.5)

#open land fraction dataset
sftlf = xr.open_dataset(resourcesdir + "sftlf_fx_GISS-E2-1-G_piControl_r1i1p1f1_gn.nc")
#convert longitude 0-360 to -180-180
sftlf = sftlf.assign_coords(lon = (((sftlf.lon + 180)%360)-180))
sftlf = sftlf.sortby('lon')

#open irrigation dataset
irr = xr.open_dataset(datadir + "Irrig144x90_1848to2100_FixedFuture_v3.nc", decode_times = False)
irr = irr.rename({'latitude': 'lat', 'longitude': 'lon'})

#find JJA average irrigation rate for the year 2000
#units are converted from m/s to mm/day
irr = irr['irrigation_per_m2'][1829:1832, :, :]*60*60*24*1000

irr_mean = np.average(irr, axis = 0, weights = monthweights[5:8])
#≥0.05 mm/day
#from m/s to m/day to mm/day
IRRmask = (irr_mean >= 0.05)
IRRhalfmask = (irr_mean/2 >= 0.05)

#create xarray to select irrigated 
xr_irrigated = xr.Dataset({'IRRmask': (['lat', 'lon'], IRRmask),
                           'IRRhalfmask': (['lat', 'lon'], IRRhalfmask)},
                      coords = {'lat': irr.lat.values, 'lon': irr.lon.values})  

#create mask for "local" and "remote" irrigation for full irrigation run
localIRRmask = (xr_irrigated['IRRmask'] == True)
remoteIRRmask = ((sftlf['sftlf'] > 0)*(xr_irrigated['IRRmask'] == False))
remoteIRRmask[remoteIRRmask['lat'] < -60] = False 
remoteIRRmask[remoteIRRmask['lat'] > 60] = False

localIRRmask = localIRRmask.reindex(lat = newlat, lon = newlon, method = 'nearest')
remoteIRRmask = remoteIRRmask.reindex(lat = newlat, lon = newlon, method = 'nearest')

#for sftlf
sftlf = sftlf.reindex(lat = newlat, lon = newlon, method = 'nearest')


#open IPCC regions
regionfilenames = fnmatch.filter(os.listdir(resourcesdir + 'referenceRegions/'), '*.nc')
#find region names
regionsearchstr = 'region(.*)_highres.nc';
regionnames = [re.search(regionsearchstr, name).group(1) for name in regionfilenames]


#MED---------change this for different region
regionnum = 8
regionname = regionnames[regionnum]


regionds = xr.open_dataset(resourcesdir + 'referenceRegions/' + regionfilenames[regionnum])
regionds = regionds.sortby('Latitude', ascending = True)
#regionds['Region'].plot()
regionds = regionds.rename({'Latitude': 'lat', 'Longitude': 'lon'})

regmask = xr.Dataset({'mask': (['lat', 'lon'], regionds['Region'].values)}, 
    coords = {'lat': regionds['lat'].values, 'lon': regionds['lon'].values})

#----------------------------------------------

#find data array and regrid

#variable
variable = 'JJAtsmaxevents'

#import irrigation, half irrigation, no irrigation tsmax event frequency
noIRRfeatures = xr.open_dataset(resultsdir + variable + "_GISS-E2-1-G_" + noIRRforcing + "_eventfreqandbvars_dailythreshold_final.nc")      
yesIRRfeatures = xr.open_dataset(resultsdir + variable + "_GISS-E2-1-G_" + yesIRRforcing + "_eventfreqandbvars_dailythreshold_final.nc")      

#reindex to match regional mask
noIRRfeatures = noIRRfeatures.reindex(lat = newlat, lon = newlon, method = 'nearest')
yesIRRfeatures = yesIRRfeatures.reindex(lat = newlat, lon = newlon, method = 'nearest')

reldiff = 100*(yesIRRfeatures['JJAevents']-noIRRfeatures['JJAevents'])/noIRRfeatures['JJAevents']

#local

local_reldiff_kderaw = kde_preprocessing(reldiff.values, regmask['mask']['lat'].values, regmask['mask'].values, (localIRRmask*sftlf['sftlf']).values)
local_reldiff_values = local_reldiff_kderaw[0]
local_reldiff_weights = local_reldiff_kderaw[1]

#remote

remote_reldiff_kderaw = kde_preprocessing(reldiff.values, regmask['mask']['lat'].values, regmask['mask'].values, (remoteIRRmask*sftlf['sftlf']).values)
remote_reldiff_values = remote_reldiff_kderaw[0]
remote_reldiff_weights = remote_reldiff_kderaw[1]

bound_min = np.min([local_reldiff_values, remote_reldiff_values])
bound_max = np.max([local_reldiff_values, remote_reldiff_values])
x_pts = np.linspace(bound_min - (bound_max-bound_min)*0.15, bound_max + (bound_max-bound_min)*0.15, np.int(bound_max-bound_min))  


#define local and remote KDE
local_reldiff_kde = scipy.stats.gaussian_kde(local_reldiff_values, weights = local_reldiff_weights)
local_reldiff_estpdf = local_reldiff_kde.evaluate(x_pts)

remote_reldiff_kde = scipy.stats.gaussian_kde(remote_reldiff_values, weights = remote_reldiff_weights)
remote_reldiff_estpdf = remote_reldiff_kde.evaluate(x_pts)

#plot results and save
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pts, local_reldiff_estpdf, color = sns.color_palette()[0], label = 'Local')
ax.plot(x_pts, remote_reldiff_estpdf, color = sns.color_palette()[2], label = 'Remote')
ax.set_title(regionname + ' Relative difference in tsmax event frequency')
plt.axvline(x= 0, color = 'black', ls = '--')
plt.xlabel('Percent difference (%)')
plt.ylabel('Density')
ax.legend(loc = 'upper right')

pl.savefig(figdir + variable + '_' + regionname + 'reldiffinregionalweightedPDF_localvsremote.png', dpi = 800, bbox_inches = 'tight')



#variable
variable = 'JJATEevents'

#import irrigation, half irrigation, no irrigation tsmax event frequency
noIRRfeatures = xr.open_dataset(resultsdir + "JJA_TEevents_GISS-E2-1-G_" + noIRRforcing + "_eventfreqandbvars_dailythreshold_final.nc")      
yesIRRfeatures = xr.open_dataset(resultsdir + "JJA_TEevents_GISS-E2-1-G_" + yesIRRforcing + "_eventfreqandbvars_dailythreshold_final.nc")      

#reindex to match regional mask
noIRRfeatures = noIRRfeatures.reindex(lat = newlat, lon = newlon, method = 'nearest')
yesIRRfeatures = yesIRRfeatures.reindex(lat = newlat, lon = newlon, method = 'nearest')


reldiff = 100*(yesIRRfeatures['JJAevents']-noIRRfeatures['JJAevents'])/noIRRfeatures['JJAevents']

#local

local_reldiff_kderaw = kde_preprocessing(reldiff.values, regmask['mask']['lat'].values, regmask['mask'].values, (localIRRmask*sftlf['sftlf']).values)
local_reldiff_values = local_reldiff_kderaw[0]
local_reldiff_weights = local_reldiff_kderaw[1]

#remote

remote_reldiff_kderaw = kde_preprocessing(reldiff.values, regmask['mask']['lat'].values, regmask['mask'].values, (remoteIRRmask*sftlf['sftlf']).values)
remote_reldiff_values = remote_reldiff_kderaw[0]
remote_reldiff_weights = remote_reldiff_kderaw[1]

bound_min = np.min([local_reldiff_values, remote_reldiff_values])
bound_max = np.max([local_reldiff_values, remote_reldiff_values])
x_pts = np.linspace(bound_min - (bound_max-bound_min)*0.15, bound_max + (bound_max-bound_min)*0.15, np.int(bound_max-bound_min))  


#define local and remote KDE
local_reldiff_kde = scipy.stats.gaussian_kde(local_reldiff_values, weights = local_reldiff_weights)
local_reldiff_estpdf = local_reldiff_kde.evaluate(x_pts)

remote_reldiff_kde = scipy.stats.gaussian_kde(remote_reldiff_values, weights = remote_reldiff_weights)
remote_reldiff_estpdf = remote_reldiff_kde.evaluate(x_pts)

#plot results and save
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pts, local_reldiff_estpdf, color = sns.color_palette()[0], label = 'Local')
ax.plot(x_pts, remote_reldiff_estpdf, color = sns.color_palette()[2], label = 'Remote')
ax.set_title(regionname + ' Relative difference in TE event frequency')
plt.axvline(x= 0, color = 'black', ls = '--')
plt.xlabel('Percent difference (%)')
plt.ylabel('Density')
ax.legend(loc = 'upper right')

pl.savefig(figdir + variable + '_' + regionname + 'reldiffinregionalweightedPDF_localvsremote.png', dpi = 800, bbox_inches = 'tight')

  

#import JJA average climate variables -------------------------------------#

#-------------------change here for other variables---------------------------#
variable = 'qsurf'
longvarname = 'Specific Humidity'
#----------------------------------------------------------------------------#

#open data
noIRRvar = xr.open_dataset(datadir + variable + "_mon_GISS-E2-1-G_" + noIRRforcing + "_185101-198612.nc")
yesIRRvar = xr.open_dataset(datadir + variable + "_mon_GISS-E2-1-G_" + yesIRRforcing + "_185101-198612.nc")

noIRRvar['time'] = pd.date_range('1851-01', '1986-12', freq = 'MS')
yesIRRvar['time'] = pd.date_range('1851-01', '1986-12', freq = 'MS')


#remove first 5 years
noIRRvar = noIRRvar[variable][60:1632,:,:]
yesIRRvar = yesIRRvar[variable][60:1632,:,:]

month_length = noIRRvar.time.dt.days_in_month
weights = (month_length.groupby('time.season')/month_length.groupby('time.season').sum())

novar_weightedavg = (noIRRvar*weights).groupby('time.season').sum(dim = 'time')
yesvar_weightedavg = (yesIRRvar*weights).groupby('time.season').sum(dim = 'time')


#find seasonal averages
noIRR_kdeinput = novar_weightedavg[1]
yesIRR_kdeinput = yesvar_weightedavg[1]

#reindex to match regional mask
noIRR_kdeinput = noIRR_kdeinput.reindex(lat = newlat, lon = newlon, method = 'nearest')
yesIRR_kdeinput = yesIRR_kdeinput.reindex(lat = newlat, lon = newlon, method = 'nearest')

reldiff = 100*(yesIRR_kdeinput-noIRR_kdeinput)/noIRR_kdeinput

#local

local_reldiff_kderaw = kde_preprocessing(reldiff.values, regmask['mask']['lat'].values, regmask['mask'].values, (localIRRmask*sftlf['sftlf']).values)
local_reldiff_values = local_reldiff_kderaw[0]
local_reldiff_weights = local_reldiff_kderaw[1]

#remote

remote_reldiff_kderaw = kde_preprocessing(reldiff.values, regmask['mask']['lat'].values, regmask['mask'].values, (remoteIRRmask*sftlf['sftlf']).values)
remote_reldiff_values = remote_reldiff_kderaw[0]
remote_reldiff_weights = remote_reldiff_kderaw[1]

bound_min = np.min([local_reldiff_values, remote_reldiff_values])
bound_max = np.max([local_reldiff_values, remote_reldiff_values])
x_pts = np.linspace(bound_min - (bound_max-bound_min)*0.15, bound_max + (bound_max-bound_min)*0.15, np.int(bound_max-bound_min)*100)  


#define local and remote KDE
local_reldiff_kde = scipy.stats.gaussian_kde(local_reldiff_values, weights = local_reldiff_weights)
local_reldiff_estpdf = local_reldiff_kde.evaluate(x_pts)

remote_reldiff_kde = scipy.stats.gaussian_kde(remote_reldiff_values, weights = remote_reldiff_weights)
remote_reldiff_estpdf = remote_reldiff_kde.evaluate(x_pts)

#plot results and save
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pts, local_reldiff_estpdf, color = sns.color_palette()[0], label = 'Local')
ax.plot(x_pts, remote_reldiff_estpdf, color = sns.color_palette()[2], label = 'Remote')
ax.set_title(regionname + ' Relative difference in ' + longvarname)
plt.axvline(x= 0, color = 'black', ls = '--')
plt.xlabel('Percent difference (%)')
plt.ylabel('Density')
ax.legend(loc = 'upper right')

pl.savefig(figdir + variable + '_' + regionname + 'reldiffinregionalweightedPDF_localvsremote.png', dpi = 800, bbox_inches = 'tight')
