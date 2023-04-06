#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Differenced plots for local vs remote - use monthly inputs

    1) Define mask to classify irrigated vs non-irrigated cells 
        -IRR = 2000 JJA average >= 0.05 is considered irrigated
        -Use latitude and land fraction data from GISS model to weight the values from the grid cells
    2) Create PDFs for climate variables

Author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import numpy as np
import xarray as xr
import pandas as pd
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
yesIRRhalfforcing = 'Irrig_2000half'


monthweights = [monthrange(2001, 1)[1]/365, monthrange(2001, 2)[1]/365, monthrange(2001, 3)[1]/365, monthrange(2001, 4)[1]/365, monthrange(2001, 5)[1]/365, monthrange(2001, 6)[1]/365, monthrange(2001, 7)[1]/365, monthrange(2001, 8)[1]/365, monthrange(2001, 9)[1]/365, monthrange(2001, 10)[1]/365, monthrange(2001, 11)[1]/365, monthrange(2001, 12)[1]/365]

#new longitude and latitude
newlon = np.arange(-179.75, 180.25, 0.5)
newlat = np.arange(-89.75,90.25, 0.5)

t = pd.date_range('1848-01', '2100-12', freq = 'MS').strftime("%b%Y").tolist()
t = [date.upper() for date in t]

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

#calculate area-averaged weights with cos
weights = np.cos(np.deg2rad(localIRRmask.lat))
weights.name = 'weights'
        
    
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


variable = 'qsurf' #------------revise here----------------------#

#forcing
noforcing = 'Irrig_1850'
yesforcing = 'Irrig_2000'
#yesforcing = 'Irrig_2000half'

title = 'Irrigation = 2000'
#title = 'Irrigation = 2000 (Half)'

#open data
novar = xr.open_dataset(datadir + variable + "_mon_GISS-E2-1-G_" + noforcing + "_185101-198612.nc")
novar['time'] = pd.date_range('1851-01', '1986-12', freq = 'MS')

yesvar = xr.open_dataset(datadir + variable + "_mon_GISS-E2-1-G_" + yesforcing + "_185101-198612.nc")
yesvar['time'] = pd.date_range('1851-01', '1986-12', freq = 'MS')


#remove first 5 years
novar = novar[variable][60:1632,:,:]
yesvar = yesvar[variable][60:1632,:,:]

month_length = novar.time.dt.days_in_month
weights = (month_length.groupby('time.season')/month_length.groupby('time.season').sum())

novar_weightedavg = (novar*weights).groupby('time.season').sum(dim = 'time')
yesvar_weightedavg = (yesvar*weights).groupby('time.season').sum(dim = 'time')

#find seasonal averages
noIRR_kdeinput = novar_weightedavg[1]
yesIRR_kdeinput = yesvar_weightedavg[1]


#----------------------------relative diff-------------------------------------------#
currentfeatures = 100*(yesIRR_kdeinput-noIRR_kdeinput)/noIRR_kdeinput


#conduct kde preprocessing for local
local_kderaw = kde_preprocessing(currentfeatures.values, localIRRmask['lat'].values, localIRRmask.values, sftlf['sftlf'].values)
local_values = local_kderaw[0]
local_weights = local_kderaw[1]

#conduct kde preprocessing for remote
remote_kderaw = kde_preprocessing(currentfeatures.values, remoteIRRmask['lat'].values, remoteIRRmask.values, sftlf['sftlf'].values)
remote_values = remote_kderaw[0]
remote_weights = remote_kderaw[1]

bound_min = np.min(np.concatenate((local_values, remote_values)))
bound_max = np.max(np.concatenate((local_values, remote_values)))
x_pts = np.linspace(bound_min - (bound_max-bound_min)*0.15, bound_max + (bound_max-bound_min)*0.15, 1000)    

#define diff KDE for local
local_kde = scipy.stats.gaussian_kde(local_values, weights = local_weights)
local_estpdf = local_kde.evaluate(x_pts)

#for remote
remote_kde = scipy.stats.gaussian_kde(remote_values, weights = remote_weights)
remote_estpdf = remote_kde.evaluate(x_pts)

#plot results and save
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pts, local_estpdf, color = sns.color_palette()[0], label = 'Local')
ax.plot(x_pts, remote_estpdf, color = sns.color_palette()[2], label = 'Remote')
ax.set_title('Relative difference in JJA average specific humidity')
plt.axvline(x= 0, color = 'black', ls = '--')
plt.xlabel('Percent difference (%)')
plt.ylabel('Density')
ax.legend(loc = 'upper right')
#plt.xlim([-10, 30])
pl.savefig(figdir + variable + '_reldiffinlocalvsremote_weightedPDF.png', dpi = 800, bbox_inches = 'tight')

