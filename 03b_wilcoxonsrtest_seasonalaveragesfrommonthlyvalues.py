#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Conduct Wilcoxon Signed-Rank Test for seasonal averages from monthly climate data

Author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import os
import numpy as np
import xarray as xr
import pandas as pd
import scipy.stats

#access the relevant data
datadir = "/Volumes/GISS/IRR project/figshare/"
resultsdir = "/Volumes/GISS/IRR project/results/"
figdir = "/Volumes/GISS/IRR project/figures/"

#change the working directory
os.chdir(datadir)

#incsw_grnd, HWV, sensht, pcldt, prec, qsurf

variable = 'qsurf' #---replace with variable of interest

#forcing
noforcing = 'Irrig_1850'
#yesforcing = 'Irrig_2000'
yesforcing = 'Irrig_2000half'

#open data
novar = xr.open_dataset(datadir + variable + "_mon_GISS-E2-1-G_" + noforcing + "_185101-198612.nc")
novar['time'] = pd.date_range('1851-01', '1986-12', freq = 'MS')

yesvar = xr.open_dataset(datadir + variable + "_mon_GISS-E2-1-G_" + yesforcing + "_185101-198612.nc")
yesvar['time'] = pd.date_range('1851-01', '1986-12', freq = 'MS')


#remove first 5 years
novar = novar[variable][60:1632,:,:]
yesvar = yesvar[variable][60:1632,:,:]

novar_JJA = novar[novar['time.season'] == "JJA"]
yesvar_JJA = yesvar[yesvar['time.season'] == "JJA"]

month_length_JJA = novar_JJA.time.dt.days_in_month
weights = (month_length_JJA.groupby('time.year')/month_length_JJA.groupby('time.year').sum())

latarray = novar['lat'].values
lonarray = novar['lon'].values

latlen = len(latarray)
lonlen = len(lonarray)

#find seasonal averages for no IRR and yes IRR
novar_seasonal = (novar_JJA*weights).groupby('time.year').sum(dim = 'time')
yesvar_seasonal = (yesvar_JJA*weights).groupby('time.year').sum(dim = 'time')

wilcoxpval = np.empty((latlen, lonlen))
wilcoxpval[:] = np.NAN

#for each lat, lon
for latpoint in range(0, latlen):
    print(latpoint)
    for lonpoint in range(0, lonlen):
        #print(lonpoint)
        
        if not np.any(novar_seasonal[:,latpoint, lonpoint]) and not np.any(novar_seasonal[:,latpoint, lonpoint]):
             wilcoxpval[latpoint, lonpoint] = 1
        else:
        
            wilcoxpval[latpoint, lonpoint] = scipy.stats.wilcoxon(x = novar_seasonal[:,latpoint, lonpoint].values, y = yesvar_seasonal[:,latpoint, lonpoint].values, mode = 'approx')[1]

#find and export csv of lat lon hatch marks of non-significant grid cells
latind = latarray[np.where(wilcoxpval > 0.01)[0]]
lonind = lonarray[np.where(wilcoxpval > 0.01)[1]]

hatch = np.vstack((lonind, latind)).T
df = pd.DataFrame(hatch, columns = ('lon', 'lat'))

df.to_csv(resultsdir + variable + '_wilcoxonsigpval_' + yesforcing + '.csv')

