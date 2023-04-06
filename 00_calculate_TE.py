#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Find the Equivalent Temperature (for each grid cell) using maximum temperature and minimum relative humidity

Equivalent temperature references:
1) Matthews et al. (2022) - Latent heat must be visible in climate communications - WIREs Climate Change
2) Pielke et al. (2004) - Assessing "global warming" with surface heat content - Eos, Transactions, AGU

Author: Felicia Chiang, felicia.chiang@nasa.gov

"""

#import libraries
import os
import numpy as np
import xarray as xr
from itertools import product
from cftime import DatetimeNoLeap

#define functions
def calculate_es(tsmax):
    # from Murray 1966, using tsmax in Celsius
    es = 6.1078*np.exp((17.2693882*tsmax)/(tsmax +237.3))
    return es

def calculate_ea(rsmin, tsmax):    
    #calculate saturation vapor pressure
    es = calculate_es(tsmax)
    
    #calculate actual vapor pressure
    ea = (rsmin/100)*es
    
    return ea

#function to calculate specific humidity from relative humidity and tsmax and psurf
def calculate_q(rsmin, tsmax, prsurf):   
    
    ea = calculate_ea(rsmin, tsmax)
    
    qmin = (0.622*ea)/(prsurf-0.378*ea)
      
    return qmin

#function to calculate latent heat of vaporization
def calculate_le(tsmax):
    
    tsmax_k = tsmax + 273.15
    
    le = (((647.3 - tsmax_k)/(647.3-373.15))**0.38)*(2.25e6)
    
    return le

#function to calculate equivalent temperature
def calculate_TE(rsmin, tsmax, prsurf): 
    #calculate Le
    Le = calculate_le(tsmax) #J/kg
    
    Cp = 1005.7 #J/kg-K
    
    #calculate specific humidity kg (of water) / kg (of air)
    q = calculate_q(rsmin, tsmax, prsurf)
    
    #calculate equivalent temperature
    TE = tsmax + (Le*q)/Cp
    
    return TE

#define forcings to examine
noIRRforcing = 'Irrig_1850'
#yesIRRforcing = 'Irrig_2000'
yesIRRforcing = 'Irrig_2000half'

#access the relevant data
datadir = "/Volumes/GISS/IRR project/figshare/"

#change the working directory
os.chdir(datadir)

#open data
noIRR_tsmax = xr.open_dataset(datadir + "tsmax_day_GISS-E2-1-G_" + noIRRforcing + "_185101-198612.nc")
noIRR_rsmin = xr.open_dataset(datadir + "rsmin_day_GISS-E2-1-G_" + noIRRforcing + "_185101-198612.nc")

yesIRR_tsmax = xr.open_dataset(datadir + "tsmax_day_GISS-E2-1-G_" + yesIRRforcing + "_185101-198612.nc")
yesIRR_rsmin = xr.open_dataset(datadir + "rsmin_day_GISS-E2-1-G_" + yesIRRforcing + "_185101-198612.nc")

#rsmin at 100%
noIRR_rsmin_capped = xr.where(noIRR_rsmin['rsmin'] > 100, 100, noIRR_rsmin['rsmin'])
yesIRR_rsmin_capped = xr.where(yesIRR_rsmin['rsmin'] > 100, 100, yesIRR_rsmin['rsmin'])


#open monthly prsurf
noIRR_prsurf = xr.open_dataset(datadir + "prsurf_day_GISS-E2-1-G_" + noIRRforcing + "_185101-198612.nc")
yesIRR_prsurf = xr.open_dataset(datadir + "prsurf_day_GISS-E2-1-G_" + yesIRRforcing + "_185101-198612.nc")

#add monthly dates to prsurf
monthlydates = [DatetimeNoLeap(year, month, 1) for year, month in product(range(1851, 1987), range(1, 13))]
noIRR_prsurf['time'] = monthlydates
yesIRR_prsurf['time'] = monthlydates

#resample to daily
noIRR_prsurf_resamp = noIRR_prsurf.reindex(time = noIRR_tsmax['time'].values, method = 'ffill')
yesIRR_prsurf_resamp = yesIRR_prsurf.reindex(time = noIRR_tsmax['time'].values, method = 'ffill')

#find actual vapor pressure (capped)
noIRR_ea_capped = calculate_ea(noIRR_rsmin_capped, noIRR_tsmax)
yesIRR_ea_capped = calculate_ea(yesIRR_rsmin_capped, yesIRR_tsmax)

#save
noIRR_ea_capped = noIRR_ea_capped['tsmax'].rename("ea")
yesIRR_ea_capped = yesIRR_ea_capped['tsmax'].rename("ea")

noIRR_ea_capped.to_netcdf(datadir + "ea_day_GISS-E2-1-G_" + noIRRforcing + "_185101-198612.nc") 
yesIRR_ea_capped.to_netcdf(datadir + "ea_day_GISS-E2-1-G_" + yesIRRforcing + "_185101-198612.nc")

#find specific humidity (capped)
noIRR_q_capped = calculate_q(noIRR_rsmin_capped, noIRR_tsmax, noIRR_prsurf_resamp['prsurf'])
yesIRR_q_capped = calculate_q(yesIRR_rsmin_capped, yesIRR_tsmax, yesIRR_prsurf_resamp['prsurf'])

#save
noIRR_q_capped = noIRR_q_capped['tsmax'].rename("q")
yesIRR_q_capped = yesIRR_q_capped['tsmax'].rename("q")

noIRR_q_capped.to_netcdf(datadir + "q_day_GISS-E2-1-G_" + noIRRforcing + "_185101-198612.nc") 
yesIRR_q_capped.to_netcdf(datadir + "q_day_GISS-E2-1-G_" + yesIRRforcing + "_185101-198612.nc")


#find equivalent temperature (capped)
noIRR_TE_capped = calculate_TE(noIRR_rsmin_capped, noIRR_tsmax, noIRR_prsurf_resamp['prsurf'])
yesIRR_TE_capped = calculate_TE(yesIRR_rsmin_capped, yesIRR_tsmax, yesIRR_prsurf_resamp['prsurf'])

#save as netcdfs
noIRR_TE_capped = noIRR_TE_capped['tsmax'].rename("TE")
yesIRR_TE_capped = yesIRR_TE_capped['tsmax'].rename("TE")

noIRR_TE_capped.to_netcdf(datadir + "TE_day_GISS-E2-1-G_" + noIRRforcing + "_185101-198612.nc")  
yesIRR_TE_capped.to_netcdf(datadir + "TE_day_GISS-E2-1-G_" + yesIRRforcing + "_185101-198612.nc")  

noIRR_TE_capped = xr.open_dataset(datadir + "TE_day_GISS-E2-1-G_" + noIRRforcing + "_185101-198612.nc")
yesIRR_TE_capped = xr.open_dataset(datadir + "TE_day_GISS-E2-1-G_" + yesIRRforcing + "_185101-198612.nc")
