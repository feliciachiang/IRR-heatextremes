#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot JJA averages for all relevant background variables (lhflx, shflx, specific humidity, precipitation, incident solar radiation, total cloud cover)
    -Find the JJA average for each year, then find the overall JJA average
    -Find the difference between the two JJA averages (between 1850fix and 2000, between 1850fix and 2000half)
    -Save plots

Author: Felicia Chiang, felicia.chiang@nasa.gov
"""

#import libraries
import os
import xarray as xr
import pandas as pd

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

t = pd.date_range('1851-01', '1986-12', freq = 'MS').strftime("%b%Y").tolist()
t = [date.upper() for date in t]

markersize = 0.05
alphavalue = 1

#Incident Solar Radiation ----------------------------------------------------
variable = 'incsw_grnd'

#forcing
noforcing = 'Irrig_1850'
#yesforcing = 'Irrig_2000'
yesforcing = 'Irrig_2000half'

#title = 'Irrigation = 2000'
title = 'Irrigation = 2000 (Half)'

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

#plot no forcing plot
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = novar_weightedavg[1].plot(cmap = 'YlOrRd', vmin = 0, vmax = 300, add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title('Irrigation = 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('JJA Average Incident Solar Radiation (W/m^2)')
pl.savefig(figdir + variable + '_JJAmean_' + noforcing + '.png', dpi = 800, bbox_inches = 'tight')

#plot yes forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = yesvar_weightedavg[1].plot(cmap = 'YlOrRd', vmin = 0, vmax = 300, add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title(title)
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('JJA Average Incident Solar Radiation (W/m^2)')
pl.savefig(figdir + variable + '_JJAmean_' + yesforcing + '.png', dpi = 800, bbox_inches = 'tight')

hatch = pd.read_csv(resultsdir + variable + '_wilcoxonsigpval_' + yesforcing + '.csv')

#plot difference in yes and no IRR forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (yesvar_weightedavg[1]-novar_weightedavg[1]).plot(cmap = 'RdBu_r', vmin = -10, vmax = 10, add_colorbar = False)
ax.scatter(hatch['lon'], hatch['lat'], s = markersize, marker = 'x', color = 'grey', alpha = alphavalue)
ax.set_ylim([-60,90])
ax.set_title('Difference between ' + title + ' and 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('JJA Average Incident Solar Radiation (W/m^2)')
pl.savefig(figdir + variable + '_JJAmean_differencefrom_' + yesforcing + '_hatched.png', dpi = 800, bbox_inches = 'tight')

#plot relative difference in yes and no IRR forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (100*((yesvar_weightedavg[1]-novar_weightedavg[1])/novar_weightedavg[1])).plot(cmap = 'RdBu_r', vmin = -10, vmax = 10, add_colorbar = False)
ax.scatter(hatch['lon'], hatch['lat'], s = markersize, marker = 'x', color = 'grey', alpha = alphavalue)
ax.set_ylim([-60,90])
ax.set_title('Relative diff between ' + title + ' and 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Relative difference (%)')
pl.savefig(figdir + variable + '_JJAmean_reldifferencefrom_' + yesforcing + '_hatched.png', dpi = 800, bbox_inches = 'tight')

#Total Cloud Cover----------------------------------------------------

variable = 'pcldt'

#forcing
noforcing = 'Irrig_1850'
#yesforcing = 'Irrig_2000'
yesforcing = 'Irrig_2000half'

#title = 'Irrigation = 2000'
title = 'Irrigation = 2000 (Half)'

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

#plot no forcing plot
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = novar_weightedavg[1].plot(cmap = 'YlOrRd', vmin = 0, vmax = 100, add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title('Irrigation = 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('JJA Average Total Cloud Cover (%)')
pl.savefig(figdir + variable + '_JJAmean_' + noforcing + '.png', dpi = 800, bbox_inches = 'tight')

#plot yes forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = yesvar_weightedavg[1].plot(cmap = 'YlOrRd', vmin = 0, vmax = 100, add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title(title)
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('JJA Average Total Cloud Cover (%)')
pl.savefig(figdir + variable + '_JJAmean_' + yesforcing + '.png', dpi = 800, bbox_inches = 'tight')

hatch = pd.read_csv(resultsdir + variable + '_wilcoxonsigpval_' + yesforcing + '.csv')

#plot difference in yes and no IRR forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (yesvar_weightedavg[1]-novar_weightedavg[1]).plot(cmap = 'BrBG', vmin = -5, vmax = 5, add_colorbar = False)
ax.scatter(hatch['lon'], hatch['lat'], s = markersize, marker = 'x', color = 'grey', alpha = alphavalue)
ax.set_ylim([-60,90])
ax.set_title('Difference between ' + title + ' and 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('JJA Average Total Cloud Cover (%)')
pl.savefig(figdir + variable + '_JJAmean_differencefrom_' + yesforcing + '_hatched.png', dpi = 800, bbox_inches = 'tight')

#plot relative difference in yes and no IRR forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (100*((yesvar_weightedavg[1]-novar_weightedavg[1])/novar_weightedavg[1])).plot(cmap = 'BrBG', vmin = -50, vmax = 50, add_colorbar = False)
ax.scatter(hatch['lon'], hatch['lat'], s = markersize, marker = 'x', color = 'grey', alpha = alphavalue)
ax.set_ylim([-60,90])
ax.set_title('Relative diff between ' + title + ' and 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Relative difference (%)')
pl.savefig(figdir + variable + '_JJAmean_reldifferencefrom_' + yesforcing + '_hatched.png', dpi = 800, bbox_inches = 'tight')


#Precipitation---------------------------------------------------------------
variable = 'prec'

#forcing
noforcing = 'Irrig_1850'
#yesforcing = 'Irrig_2000'
yesforcing = 'Irrig_2000half'

#title = 'Irrigation = 2000'
title = 'Irrigation = 2000 (Half)'

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


#plot no forcing plot
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = novar_weightedavg[1].plot(cmap = 'BrBG', vmin = 0, vmax = 5, add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title('Irrigation = 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('JJA Average Precipitation (mm/day)')
pl.savefig(figdir + variable + '_JJAmean_' + noforcing + '.png', dpi = 800, bbox_inches = 'tight')

#plot yes forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = yesvar_weightedavg[1].plot(cmap = 'BrBG', vmin = 0, vmax = 5, add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title(title)
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('JJA Average Precipitation (mm/day)')
pl.savefig(figdir + variable + '_JJAmean_' + yesforcing + '.png', dpi = 800, bbox_inches = 'tight')

hatch = pd.read_csv(resultsdir + variable + '_wilcoxonsigpval_' + yesforcing + '.csv')

#plot difference in yes and no IRR forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (yesvar_weightedavg[1]-novar_weightedavg[1]).plot(cmap = 'BrBG', vmin = -0.5, vmax = 0.5, add_colorbar = False)
ax.scatter(hatch['lon'], hatch['lat'], s = markersize, marker = 'x', color = 'grey', alpha = alphavalue)
ax.set_ylim([-60,90])
ax.set_title('Difference between ' + title + ' and 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('JJA Average Precipitation (mm/day)')
pl.savefig(figdir + variable + '_JJAmean_differencefrom_' + yesforcing + '_hatched.png', dpi = 800, bbox_inches = 'tight')

#plot relative difference in yes and no IRR forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (100*(yesvar_weightedavg[1]-novar_weightedavg[1])/novar_weightedavg[1]).plot(cmap = 'BrBG', vmin = -40, vmax = 40, add_colorbar = False)
ax.scatter(hatch['lon'], hatch['lat'], s = markersize, marker = 'x', color = 'grey', alpha = alphavalue)
ax.set_ylim([-60,90])
ax.set_title('Relative diff between ' + title + ' and 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Relative difference (%)')
pl.savefig(figdir + variable + '_JJAmean_reldifferencefrom_' + yesforcing + '_hatched_check.png', dpi = 800, bbox_inches = 'tight')


#Latent heat flux---------------------------------------------------------------

variable = 'HWV'

multiplier = -1

#forcing
noforcing = 'Irrig_1850'
#yesforcing = 'Irrig_2000'
yesforcing = 'Irrig_2000half'

#title = 'Irrigation = 2000'
title = 'Irrigation = 2000 (Half)'

#open data
novar = xr.open_dataset(datadir + variable + "_mon_GISS-E2-1-G_" + noforcing + "_185101-198612.nc")
novar['time'] = pd.date_range('1851-01', '1986-12', freq = 'MS')

yesvar = xr.open_dataset(datadir + variable + "_mon_GISS-E2-1-G_" + yesforcing + "_185101-198612.nc")
yesvar['time'] = pd.date_range('1851-01', '1986-12', freq = 'MS')


#remove first 5 years
novar = multiplier*novar[variable][60:1632,:,:]
yesvar = multiplier*yesvar[variable][60:1632,:,:]

month_length = novar.time.dt.days_in_month
weights = (month_length.groupby('time.season')/month_length.groupby('time.season').sum())

novar_weightedavg = (novar*weights).groupby('time.season').sum(dim = 'time')
yesvar_weightedavg = (yesvar*weights).groupby('time.season').sum(dim = 'time')

#plot no forcing plot
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = novar_weightedavg[1].plot(cmap = 'RdBu_r', vmin = -200, vmax = 200, add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title('Irrigation = 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Latent Heat Flux (W/m^2)')
pl.savefig(figdir + variable + '_JJAmean_' + noforcing + '.png', dpi = 800, bbox_inches = 'tight')

#plot yes forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = yesvar_weightedavg[1].plot(cmap = 'RdBu_r', vmin = -200, vmax = 200, add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title(title)
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Latent Heat Flux (W/m^2)')
pl.savefig(figdir + variable + '_JJAmean_' + yesforcing + '.png', dpi = 800, bbox_inches = 'tight')

hatch = pd.read_csv(resultsdir + variable + '_wilcoxonsigpval_' + yesforcing + '.csv')

#plot difference in yes and no IRR forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (yesvar_weightedavg[1]-novar_weightedavg[1]).plot(cmap = 'RdBu_r', vmin = -40, vmax = 40, add_colorbar = False)
ax.scatter(hatch['lon'], hatch['lat'], s = markersize, marker = 'x', color = 'grey', alpha = alphavalue)
ax.set_ylim([-60,90])
ax.set_title('Difference between ' + title + ' and 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Latent Heat Flux (W/m^2)')
pl.savefig(figdir + variable + '_JJAmean_differencefrom_' + yesforcing + '.png', dpi = 800, bbox_inches = 'tight')

#plot relative difference in yes and no IRR forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (100*(yesvar_weightedavg[1]-novar_weightedavg[1])/novar_weightedavg[1]).plot(cmap = 'RdBu_r', vmin = -100, vmax = 100, add_colorbar = False)
ax.scatter(hatch['lon'], hatch['lat'], s = markersize, marker = 'x', color = 'grey', alpha = alphavalue)
ax.set_ylim([-60,90])
ax.set_title('Relative diff between ' + title + ' and 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Relative difference (%)')
pl.savefig(figdir + variable + '_JJAmean_reldifferencefrom_' + yesforcing + '_check.png', dpi = 800, bbox_inches = 'tight')


#Sensible heat flux -----------------------------------------------------------
variable = 'sensht'

multiplier = -1

#forcing
noforcing = 'Irrig_1850'
#yesforcing = 'Irrig_2000'
yesforcing = 'Irrig_2000half'

#title = 'Irrigation = 2000'
title = 'Irrigation = 2000 (Half)'

#open data
novar = xr.open_dataset(datadir + variable + "_mon_GISS-E2-1-G_" + noforcing + "_185101-198612.nc")
novar['time'] = pd.date_range('1851-01', '1986-12', freq = 'MS')

yesvar = xr.open_dataset(datadir + variable + "_mon_GISS-E2-1-G_" + yesforcing + "_185101-198612.nc")
yesvar['time'] = pd.date_range('1851-01', '1986-12', freq = 'MS')


#remove first 5 years
novar = multiplier*novar[variable][60:1632,:,:]
yesvar = multiplier*yesvar[variable][60:1632,:,:]

month_length = novar.time.dt.days_in_month
weights = (month_length.groupby('time.season')/month_length.groupby('time.season').sum())

novar_weightedavg = (novar*weights).groupby('time.season').sum(dim = 'time')
yesvar_weightedavg = (yesvar*weights).groupby('time.season').sum(dim = 'time')

#plot no forcing plot
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = novar_weightedavg[1].plot(cmap = 'RdBu_r', vmin = -200, vmax = 200, add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title('Irrigation = 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Latent Heat Flux (W/m^2)')
pl.savefig(figdir + variable + '_JJAmean_' + noforcing + '.png', dpi = 800, bbox_inches = 'tight')

#plot yes forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = yesvar_weightedavg[1].plot(cmap = 'RdBu_r', vmin = -200, vmax = 200, add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title(title)
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Latent Heat Flux (W/m^2)')
pl.savefig(figdir + variable + '_JJAmean_' + yesforcing + '.png', dpi = 800, bbox_inches = 'tight')

hatch = pd.read_csv(resultsdir + variable + '_wilcoxonsigpval_' + yesforcing + '.csv')

#plot difference in yes and no IRR forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (yesvar_weightedavg[1]-novar_weightedavg[1]).plot(cmap = 'RdBu_r', vmin = -40, vmax = 40, add_colorbar = False)
ax.scatter(hatch['lon'], hatch['lat'], s = markersize, marker = 'x', color = 'grey', alpha = alphavalue)
ax.set_ylim([-60,90])
ax.set_title('Difference between ' + title + ' and 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Latent Heat Flux (W/m^2)')
pl.savefig(figdir + variable + '_JJAmean_differencefrom_' + yesforcing + '.png', dpi = 800, bbox_inches = 'tight')

#plot relative difference in yes and no IRR forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (100*(yesvar_weightedavg[1]-novar_weightedavg[1])/novar_weightedavg[1]).plot(cmap = 'RdBu_r', vmin = -100, vmax = 100, add_colorbar = False)
ax.scatter(hatch['lon'], hatch['lat'], s = markersize, marker = 'x', color = 'grey', alpha = alphavalue)
ax.set_ylim([-60,90])
ax.set_title('Relative diff between ' + title + ' and 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Relative difference (%)')
pl.savefig(figdir + variable + '_JJAmean_reldifferencefrom_' + yesforcing + '_check.png', dpi = 800, bbox_inches = 'tight')



#Specific humidity------------------------------------------------------------
variable = 'qsurf'

#forcing
noforcing = 'Irrig_1850'
#yesforcing = 'Irrig_2000'
yesforcing = 'Irrig_2000half'

#title = 'Irrigation = 2000'
title = 'Irrigation = 2000 (Half)'

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

#plot no forcing plot
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = novar_weightedavg[1].plot(vmin = 0, vmax = 0.02, cmap = 'BrBG', add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title('Irrigation = 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Specific Humidity')
pl.savefig(figdir + variable + '_JJAmean_' + noforcing + '.png', dpi = 800, bbox_inches = 'tight')

#plot yes forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = yesvar_weightedavg[1].plot(vmin = 0, vmax = 0.02, cmap = 'BrBG', add_colorbar = False)
ax.set_ylim([-60,90])
ax.set_title(title)
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Specific Humidity')
pl.savefig(figdir + variable + '_JJAmean_' + yesforcing + '.png', dpi = 800, bbox_inches = 'tight')

hatch = pd.read_csv(resultsdir + variable + '_wilcoxonsigpval_' + yesforcing + '.csv')


#plot difference in yes and no IRR forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (yesvar_weightedavg[1]-novar_weightedavg[1]).plot(vmin = -0.002, vmax =0.002, cmap = 'BrBG', add_colorbar = False)
ax.scatter(hatch['lon'], hatch['lat'], s = markersize, marker = 'x', color = 'grey', alpha = alphavalue)
ax.set_ylim([-60,90])
ax.set_title('Difference between ' + title + ' and 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Specific Humidity')
pl.savefig(figdir + variable + '_JJAmean_differencefrom_' + yesforcing + '.png', dpi = 800, bbox_inches = 'tight')

#plot relative difference in yes and no IRR forcing
fig = plt.figure()
ax = fig.add_subplot(111, projection = ccrs.PlateCarree());
ax.set_global()
ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth = 0.3, facecolor='w')
ax.coastlines(resolution = '50m', linewidth = 0.3)
im = (100*(yesvar_weightedavg[1]-novar_weightedavg[1])/novar_weightedavg[1]).plot(cmap = 'BrBG', vmin = -25, vmax = 25, add_colorbar = False)
ax.scatter(hatch['lon'], hatch['lat'], s = markersize, marker = 'x', color = 'grey', alpha = alphavalue)
ax.set_ylim([-60,90])
ax.set_title('Relative diff between ' + title + ' and 1850')
cbar = plt.colorbar(im, orientation="vertical", fraction = 0.02, pad=0.05)
cbar.set_label('Relative difference (%)')
pl.savefig(figdir + variable + '_JJAmean_reldifferencefrom_' + yesforcing + '_check.png', dpi = 800, bbox_inches = 'tight')

