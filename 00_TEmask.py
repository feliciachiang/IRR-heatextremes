#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create netcdfs to mask grid cells where more than 5% of daily minimum relative humidity data points exceed 100%
Output:
    -Masks for each individual irrigation run (IRR = 1850, IRR = 2000 (half), IRR = 2000)
    -Masks for combined irrigation runs (IRR = 1850 & IRR = 2000 (half); IRR = 1850 & IRR = 2000)

Author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import os
import numpy as np
import xarray as xr

figdir = "/Volumes/GISS/IRR project/figures/"

noIRRforcing = 'Irrig_1850'
yesIRRforcing = 'Irrig_2000'
#yesIRRforcing = 'Irrig_2000half'

#access the relevant data
datadir = "/Volumes/GISS/IRR project/figshare/"

#change the working directory
os.chdir(datadir)

#create mask based on IRR rsmin values for TE FOR JJA
#use values from 1856-1986, removing the first 5 years of the dataset

IRRforcing = 'Irrig_1850'
#IRRforcing = 'Irrig_2000'
#IRRforcing = 'Irrig_2000half'

#open data
rsmin = xr.open_dataset(datadir + "rsmin_day_GISS-E2-1-G_" + IRRforcing + "_185101-198612.nc")['rsmin'][365*5:49640,:,:]

JJAindinyear = (rsmin['time.season']== "JJA")[0:365].values
JJAindvalues = np.arange(1,366)*JJAindinyear
JJAindvalues = JJAindvalues.astype(float)
JJAindvalues[JJAindvalues == 0] = np.nan

JJAdaymin = np.nanmin(JJAindvalues)-7
JJAdaymax = np.nanmax(JJAindvalues)+7

#create heat index mask to limit tsmax values to those above 0 and rsmin values to those below 100  
TEmask = np.sum(rsmin[rsmin['time.dayofyear'].isin(np.arange(JJAdaymin, JJAdaymax + 1))] >100, axis = 0) < 0.05*131*106
TEmask = TEmask.rename("mask")

#save as netcdf
resultsdir = "/Volumes/GISS/IRR project/results/"
TEmask.to_netcdf(resultsdir + "TEmask_" + IRRforcing + "_JJA_day_GISS-E2-1-G_185601-198612.nc")


#combined TE masks
#access the relevant data
datadir = "/Volumes/GISS/IRR project/figshare/"


os.chdir(datadir)

noIRRforcing = 'Irrig_1850'
yesIRRforcing = 'Irrig_2000'
#yesIRRforcing = 'Irrig_2000half'

#open data
noIRRrsmin = xr.open_dataset(datadir + "rsmin_day_GISS-E2-1-G_" + noIRRforcing + "_185101-198612.nc")['rsmin'][365*5:49640,:,:]
yesIRRrsmin = xr.open_dataset(datadir + "rsmin_day_GISS-E2-1-G_" + yesIRRforcing + "_185101-198612.nc")['rsmin'][365*5:49640,:,:]

JJAindinyear = (noIRRrsmin['time.season']== "JJA")[0:365].values
JJAindvalues = np.arange(1,366)*JJAindinyear
JJAindvalues = JJAindvalues.astype(float)
JJAindvalues[JJAindvalues == 0] = np.nan

JJAdaymin = np.nanmin(JJAindvalues)-7
JJAdaymax = np.nanmax(JJAindvalues)+7

#create heat index mask to limit tsmax values to those above 0 and rsmin values to those below 100  
TEmask = (np.sum(noIRRrsmin[noIRRrsmin['time.dayofyear'].isin(np.arange(JJAdaymin, JJAdaymax + 1))] >100, axis = 0) < 0.05*131*106)*(np.sum(yesIRRrsmin[yesIRRrsmin['time.dayofyear'].isin(np.arange(JJAdaymin, JJAdaymax + 1))] >100, axis = 0) < 0.05*131*106)
TEmask = TEmask.rename("mask")


#save as netcdf
resultsdir = "/Volumes/GISS/IRR project/results/"
TEmask.to_netcdf(resultsdir + "TEmask_18502000_JJA_day_GISS-E2-1-G_185601-198612.nc")
#TEmask.to_netcdf(resultsdir + "TEmask_18502000half_JJA_day_GISS-E2-1-G_185601-198612.nc")
