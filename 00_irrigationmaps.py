#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Irrigation maps:
    
1) Create maps of amount of irrigation for IRR = 1850, IRR = 2000 (half), and IRR = 2000 (Figure 1)
2) Create maps with irrigated (local) vs non-irrigated (remote) grid cells based on a threshold of IRR >= 0.05mm/day (Figure S1)

Author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy as cart
from calendar import monthrange

#set directory locations
datadir = "/Volumes/GISS/IRR project/figshare/"
resourcesdir = "/Volumes/GISS/resources/"

resultsdir = "/Volumes/GISS/IRR project/results/"
figdir = "/Volumes/GISS/IRR project/figures/"

#forcing names
noIRRforcing = 'Irrig_1850'
yesIRRforcing = 'Irrig_2000'
yesIRRhalfforcing = 'Irrig_2000half'

#create weights for months
monthweights = [monthrange(2001, 1)[1]/365, monthrange(2001, 2)[1]/365, monthrange(2001, 3)[1]/365, monthrange(2001, 4)[1]/365, monthrange(2001, 5)[1]/365, monthrange(2001, 6)[1]/365, monthrange(2001, 7)[1]/365, monthrange(2001, 8)[1]/365, monthrange(2001, 9)[1]/365, monthrange(2001, 10)[1]/365, monthrange(2001, 11)[1]/365, monthrange(2001, 12)[1]/365]

#specify longitude and latitude
newlon = np.arange(-179.75, 180.25, 0.5)
newlat = np.arange(-89.75,90.25, 0.5)

#create time series
t = pd.date_range('1848-01', '2100-12', freq = 'MS').strftime("%b%Y").tolist()
t = [date.upper() for date in t]

#open irrigation dataset
irr = xr.open_dataset(datadir + "Irrig144x90_1848to2100_FixedFuture_v3.nc", decode_times = False)
irr = irr.rename({'latitude': 'lat', 'longitude': 'lon'})

#Plot JJA average irrigation at year 1850 
JJA1850irr = irr['irrigation_per_m2'][29:32, :, : ]*60*60*24*1000
JJA1850irr_mean = np.average(JJA1850irr, axis = 0, weights = monthweights[5:8])

#create xarray
xr_irrigated = xr.Dataset({'irrigated': (['lat', 'lon'], JJA1850irr_mean)},
                      coords = {'lat': irr.lat.values, 'lon': irr.lon.values})  

#plot different irrigation thresholds
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = xr_irrigated['irrigated'].plot(vmin = 0, vmax = 1, add_colorbar = False)
#ax.set_ylim([-60,90])
ax.set_title('1850 JJA Average Irrigation')

cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Irrigation [mm/day]')
pl.savefig(figdir + '1850_JJAavg_irrigation.png', dpi = 800, bbox_inches = 'tight')


#Plot JJA average irrigation at year 2000
JJA2000irr = irr['irrigation_per_m2'][1829:1832, :, : ]*60*60*24*1000
JJA2000irr_mean = np.average(JJA2000irr, axis = 0, weights = monthweights[5:8])

#create xarray
xr_irrigated = xr.Dataset({'irrigated': (['lat', 'lon'], JJA2000irr_mean)},
                      coords = {'lat': irr.lat.values, 'lon': irr.lon.values})  

#plot 2000 JJA Average Irrigation
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = xr_irrigated['irrigated'].plot(vmin = 0, vmax = 1, add_colorbar = False)
#ax.set_ylim([-60,90])
ax.set_title('2000 JJA Average Irrigation')

cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Irrigation [mm/day]')
pl.savefig(figdir + '2000_JJAavg_irrigation.png', dpi = 800, bbox_inches = 'tight')


#Plot JJA average irrigation at year 2000 (halved)
JJA2000irr = (irr['irrigation_per_m2'][1829:1832, :, : ]*60*60*24*1000)/2
JJA2000irr_mean = np.average(JJA2000irr, axis = 0, weights = monthweights[5:8])

#create xarray
xr_irrigated = xr.Dataset({'irrigated': (['lat', 'lon'], JJA2000irr_mean)},
                      coords = {'lat': irr.lat.values, 'lon': irr.lon.values})  


#plot 2000 (half) JJA Average Irrigation
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = xr_irrigated['irrigated'].plot(vmin = 0, vmax = 1, add_colorbar = False)
#ax.set_ylim([-60,90])
ax.set_title('2000 (half) JJA Average Irrigation')

cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Irrigation [mm/day]')
pl.savefig(figdir + '2000half_JJAavg_irrigation.png', dpi = 800, bbox_inches = 'tight')


#------------------------------------------------------------------------------
#Plot the local and the remote grid cells

#create local irrigation mask and remote irrigation mask
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
IRRmask = (irr_mean >= 0.05)
IRRhalfmask = (irr_mean/2 >= 0.05)

#irrigated = (irr_halved > 0.01/(60*60*24*1000))
#create xarray to select irrigated 
xr_irrigated = xr.Dataset({'IRRmask': (['lat', 'lon'], IRRmask),
                           'IRRhalfmask': (['lat', 'lon'], IRRhalfmask)},
                      coords = {'lat': irr.lat.values, 'lon': irr.lon.values})  

#create mask for "local" and "remote" irrigation for full irrigation run
localIRRmask = (xr_irrigated['IRRmask'] == True)
remoteIRRmask = ((sftlf['sftlf'] > 0)*(xr_irrigated['IRRmask'] == False))
remoteIRRmask[remoteIRRmask['lat'] < -60] = False 
remoteIRRmask[remoteIRRmask['lat'] > 60] = False

#create xarray
xr_irrigated = xr.Dataset({'localIRRmask': (['lat', 'lon'], localIRRmask),
                           'remoteIRRmask': (['lat', 'lon'], remoteIRRmask)},
                      coords = {'lat': irr.lat.values, 'lon': irr.lon.values})  

#plot
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (xr_irrigated['localIRRmask']).plot(vmin = 0, vmax = 1, add_colorbar = False)
ax.set_title('Local Irrigation Grid Cells')
pl.savefig(figdir + 'localirrigation_JJAavg_threshold005.png', dpi = 800, bbox_inches = 'tight')

fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (xr_irrigated['remoteIRRmask']).plot(vmin = 0, vmax = 1, add_colorbar = False)
ax.set_title('Remote Irrigation Grid Cells')
pl.savefig(figdir + 'remoteirrigation_JJAavg_threshold005.png', dpi = 800, bbox_inches = 'tight')
