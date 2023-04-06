#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create PDFs to compare differences between local vs regional responses to irrigation

    1) Define mask to classify irrigated vs non-irrigated cells 
        -IRR = 2000 JJA average >= 0.05 is considered irrigated
        -Use land fraction data from GISS model to define non-irrigated land grid cells
    2) Import frequency results from tasmax and TE heat waves 
    3) Create PDFs to show the differences in tasmax and TE frequency in local vs remote grid cells for the globe

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

datadir = "/Volumes/GISS/IRR project/figshare/"
resultsdir = "/Volumes/GISS/IRR project/results/"
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


#import irrigation, half irrigation, no irrigation tsmax event frequency
noIRRfeatures = xr.open_dataset(resultsdir + "JJAtsmaxevents_GISS-E2-1-G_" + noIRRforcing + "_eventfreqandbvars_dailythreshold_final.nc")      
yesIRRfeatures = xr.open_dataset(resultsdir + "JJAtsmaxevents_GISS-E2-1-G_" + yesIRRforcing + "_eventfreqandbvars_dailythreshold_final.nc")      
yesIRRhalffeatures = xr.open_dataset(resultsdir + "JJAtsmaxevents_GISS-E2-1-G_" + yesIRRhalfforcing + "_eventfreqandbvars_dailythreshold_final.nc")

#variable
variable = 'JJAtsmaxevents'

#----------------------------local--------------------------------------------#

#conduct kde preprocessing
noIRR_kderaw = kde_preprocessing(noIRRfeatures['JJAevents'].values, localIRRmask['lat'].values, localIRRmask.values, sftlf['sftlf'].values)
noIRR_values = noIRR_kderaw[0]
noIRR_weights = noIRR_kderaw[1]

yesIRR_kderaw = kde_preprocessing(yesIRRfeatures['JJAevents'].values, localIRRmask['lat'].values, localIRRmask.values, sftlf['sftlf'].values)
yesIRR_values = yesIRR_kderaw[0]
yesIRR_weights = yesIRR_kderaw[1]

yesIRRhalf_kderaw = kde_preprocessing(yesIRRhalffeatures['JJAevents'].values, localIRRmask['lat'].values, localIRRmask.values, sftlf['sftlf'].values)
yesIRRhalf_values = yesIRRhalf_kderaw[0]
yesIRRhalf_weights = yesIRRhalf_kderaw[1]

bound_min = np.min([noIRR_values, yesIRRhalf_values, yesIRR_values])
bound_max = np.max([noIRR_values, yesIRRhalf_values, yesIRR_values])
x_pts = np.linspace(np.min([0, np.abs(bound_min - (bound_max-bound_min)*0.05)]), bound_max + (bound_max-bound_min)*0.05, 200)

#define no IRR KDE
noIRR_kde = scipy.stats.gaussian_kde(noIRR_values, weights = noIRR_weights)
noIRR_estpdf = noIRR_kde.evaluate(x_pts)

#define yes IRR KDE
yesIRR_kde = scipy.stats.gaussian_kde(yesIRR_values, weights = yesIRR_weights)
yesIRR_estpdf = yesIRR_kde.evaluate(x_pts)

#define yes IRR (half) KDE
yesIRRhalf_kde = scipy.stats.gaussian_kde(yesIRRhalf_values, weights = yesIRRhalf_weights)
yesIRRhalf_estpdf = yesIRRhalf_kde.evaluate(x_pts)

#plot results and save
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pts, noIRR_estpdf, color = sns.color_palette()[3], label = 'IRR = 1850')
ax.plot(x_pts, yesIRRhalf_estpdf, color = sns.color_palette()[2], label = 'IRR = 2000 (half)')
ax.plot(x_pts, yesIRR_estpdf, color = sns.color_palette()[0], label = 'IRR = 2000')
ax.set_title('Local tsmax event frequency')
plt.xlabel('Number of tsmax events')
plt.ylabel('Density')
#ax.legend(loc = 'upper left')
pl.savefig(figdir + variable + '_localweightedPDF_all_noleg.png', dpi = 800, bbox_inches = 'tight')

#plot no and yes (IRR = 2000 only) and save
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pts, noIRR_estpdf, color = sns.color_palette()[3], label = 'IRR = 1850')
#ax.plot(x_pts, yesIRRhalf_estpdf, color = sns.color_palette()[2], label = 'IRR = 2000 (half)')
ax.plot(x_pts, yesIRR_estpdf, color = sns.color_palette()[0], label = 'IRR = 2000')
ax.set_title('Local tsmax event frequency')
plt.xlabel('Number of tsmax events')
plt.ylabel('Density')
ax.legend(loc = 'upper left')
pl.savefig(figdir + variable + '_localweightedPDF.png', dpi = 800, bbox_inches = 'tight')

#----------------------------remote-------------------------------------------#

#conduct kde preprocessing
noIRR_kderaw = kde_preprocessing(noIRRfeatures['JJAevents'].values, remoteIRRmask['lat'].values, remoteIRRmask.values, sftlf['sftlf'].values)
noIRR_values = noIRR_kderaw[0]
noIRR_weights = noIRR_kderaw[1]

yesIRR_kderaw = kde_preprocessing(yesIRRfeatures['JJAevents'].values, remoteIRRmask['lat'].values, remoteIRRmask.values, sftlf['sftlf'].values)
yesIRR_values = yesIRR_kderaw[0]
yesIRR_weights = yesIRR_kderaw[1]

yesIRRhalf_kderaw = kde_preprocessing(yesIRRhalffeatures['JJAevents'].values, remoteIRRmask['lat'].values, remoteIRRmask.values, sftlf['sftlf'].values)
yesIRRhalf_values = yesIRRhalf_kderaw[0]
yesIRRhalf_weights = yesIRRhalf_kderaw[1]

bound_min = np.min([noIRR_values, yesIRRhalf_values, yesIRR_values])
bound_max = np.max([noIRR_values, yesIRRhalf_values, yesIRR_values])
x_pts = np.linspace(np.min([0, np.abs(bound_min - (bound_max-bound_min)*0.05)]), bound_max + (bound_max-bound_min)*0.05, 200)

#define no IRR KDE
noIRR_kde = scipy.stats.gaussian_kde(noIRR_values, weights = noIRR_weights)
noIRR_estpdf = noIRR_kde.evaluate(x_pts)

#define yes IRR KDE
yesIRR_kde = scipy.stats.gaussian_kde(yesIRR_values, weights = yesIRR_weights)
yesIRR_estpdf = yesIRR_kde.evaluate(x_pts)

#define yes IRR (half) KDE
yesIRRhalf_kde = scipy.stats.gaussian_kde(yesIRRhalf_values, weights = yesIRRhalf_weights)
yesIRRhalf_estpdf = yesIRRhalf_kde.evaluate(x_pts)

#plot results and save
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pts, noIRR_estpdf, color = sns.color_palette()[3], label = 'IRR = 1850')
ax.plot(x_pts, yesIRRhalf_estpdf, color = sns.color_palette()[2], label = 'IRR = 2000 (half)')
ax.plot(x_pts, yesIRR_estpdf, color = sns.color_palette()[0], label = 'IRR = 2000')
ax.set_title('Remote tsmax event frequency')
plt.xlabel('Number of tsmax events')
plt.ylabel('Density')
ax.legend(loc = 'upper right')
#plt.xlim([0, 300])
pl.savefig(figdir + variable + '_remoteweightedPDF_all.png', dpi = 800, bbox_inches = 'tight')

#plot no and yes (IRR = 2000 only) and save
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pts, noIRR_estpdf, color = sns.color_palette()[3], label = 'IRR = 1850')
#ax.plot(x_pts, yesIRRhalf_estpdf, color = sns.color_palette()[2], label = 'IRR = 2000 (half)')
ax.plot(x_pts, yesIRR_estpdf, color = sns.color_palette()[0], label = 'IRR = 2000')
ax.set_title('Remote tsmax event frequency')
plt.xlabel('Number of tsmax events')
plt.ylabel('Density')
ax.legend(loc = 'upper left')
pl.savefig(figdir + variable + '_remoteweightedPDF.png', dpi = 800, bbox_inches = 'tight')


#import TE event frequency-------------------------------------------------------------------------------------------------------


#import irrigation, half irrigation, no irrigation TE event frequency
noIRRfeatures = xr.open_dataset(resultsdir + "JJA_TEevents_GISS-E2-1-G_" + noIRRforcing + "_eventfreqandbvars_dailythreshold_final.nc")      
yesIRRfeatures = xr.open_dataset(resultsdir + "JJA_TEevents_GISS-E2-1-G_" + yesIRRforcing + "_eventfreqandbvars_dailythreshold_final.nc")      
yesIRRhalffeatures = xr.open_dataset(resultsdir + "JJA_TEevents_GISS-E2-1-G_" + yesIRRhalfforcing + "_eventfreqandbvars_dailythreshold_final.nc")


#variable
variable = 'JJATEevents'

#----------------------------local--------------------------------------------#

#conduct kde preprocessing
noIRR_kderaw = kde_preprocessing(noIRRfeatures['JJAevents'].values, localIRRmask['lat'].values, localIRRmask.values, sftlf['sftlf'].values)
noIRR_values = noIRR_kderaw[0]
noIRR_weights = noIRR_kderaw[1]

yesIRR_kderaw = kde_preprocessing(yesIRRfeatures['JJAevents'].values, localIRRmask['lat'].values, localIRRmask.values, sftlf['sftlf'].values)
yesIRR_values = yesIRR_kderaw[0]
yesIRR_weights = yesIRR_kderaw[1]

yesIRRhalf_kderaw = kde_preprocessing(yesIRRhalffeatures['JJAevents'].values, localIRRmask['lat'].values, localIRRmask.values, sftlf['sftlf'].values)
yesIRRhalf_values = yesIRRhalf_kderaw[0]
yesIRRhalf_weights = yesIRRhalf_kderaw[1]

bound_min = np.min([noIRR_values, yesIRRhalf_values, yesIRR_values])
bound_max = np.max([noIRR_values, yesIRRhalf_values, yesIRR_values])
x_pts = np.linspace(np.min([0, np.abs(bound_min - (bound_max-bound_min)*0.05)]), bound_max + (bound_max-bound_min)*0.05, 200)

#define no IRR KDE
noIRR_kde = scipy.stats.gaussian_kde(noIRR_values, weights = noIRR_weights)
noIRR_estpdf = noIRR_kde.evaluate(x_pts)

#define yes IRR KDE
yesIRR_kde = scipy.stats.gaussian_kde(yesIRR_values, weights = yesIRR_weights)
yesIRR_estpdf = yesIRR_kde.evaluate(x_pts)

#define yes IRR (half) KDE
yesIRRhalf_kde = scipy.stats.gaussian_kde(yesIRRhalf_values, weights = yesIRRhalf_weights)
yesIRRhalf_estpdf = yesIRRhalf_kde.evaluate(x_pts)

#plot results and save
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pts, noIRR_estpdf, color = sns.color_palette()[3], label = 'IRR = 1850')
ax.plot(x_pts, yesIRRhalf_estpdf, color = sns.color_palette()[2], label = 'IRR = 2000 (half)')
ax.plot(x_pts, yesIRR_estpdf, color = sns.color_palette()[0], label = 'IRR = 2000')
ax.set_title('Local TE event frequency')
plt.xlabel('Number of TE events')
plt.ylabel('Density')
#ax.legend(loc = 'upper right')
pl.savefig(figdir + variable + '_localweightedPDF_all_noleg.png', dpi = 800, bbox_inches = 'tight')

#plot no and yes (IRR = 2000 only) and save
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pts, noIRR_estpdf, color = sns.color_palette()[3], label = 'IRR = 1850')
#ax.plot(x_pts, yesIRRhalf_estpdf, color = sns.color_palette()[2], label = 'IRR = 2000 (half)')
ax.plot(x_pts, yesIRR_estpdf, color = sns.color_palette()[0], label = 'IRR = 2000')
ax.set_title('Local TE event frequency')
plt.xlabel('Number of TE events')
plt.ylabel('Density')
ax.legend(loc = 'upper right')
pl.savefig(figdir + variable + '_localweightedPDF.png', dpi = 800, bbox_inches = 'tight')

#----------------------------remote-------------------------------------------#

#conduct kde preprocessing
noIRR_kderaw = kde_preprocessing(noIRRfeatures['JJAevents'].values, remoteIRRmask['lat'].values, remoteIRRmask.values, sftlf['sftlf'].values)
noIRR_values = noIRR_kderaw[0]
noIRR_weights = noIRR_kderaw[1]

yesIRR_kderaw = kde_preprocessing(yesIRRfeatures['JJAevents'].values, remoteIRRmask['lat'].values, remoteIRRmask.values, sftlf['sftlf'].values)
yesIRR_values = yesIRR_kderaw[0]
yesIRR_weights = yesIRR_kderaw[1]

yesIRRhalf_kderaw = kde_preprocessing(yesIRRhalffeatures['JJAevents'].values, remoteIRRmask['lat'].values, remoteIRRmask.values, sftlf['sftlf'].values)
yesIRRhalf_values = yesIRRhalf_kderaw[0]
yesIRRhalf_weights = yesIRRhalf_kderaw[1]

bound_min = np.min([noIRR_values, yesIRRhalf_values, yesIRR_values])
bound_max = np.max([noIRR_values, yesIRRhalf_values, yesIRR_values])
x_pts = np.linspace(np.min([0, np.abs(bound_min - (bound_max-bound_min)*0.05)]), bound_max + (bound_max-bound_min)*0.05, 200)

#define no IRR KDE
noIRR_kde = scipy.stats.gaussian_kde(noIRR_values, weights = noIRR_weights)
noIRR_estpdf = noIRR_kde.evaluate(x_pts)

#define yes IRR KDE
yesIRR_kde = scipy.stats.gaussian_kde(yesIRR_values, weights = yesIRR_weights)
yesIRR_estpdf = yesIRR_kde.evaluate(x_pts)

#define yes IRR (half) KDE
yesIRRhalf_kde = scipy.stats.gaussian_kde(yesIRRhalf_values, weights = yesIRRhalf_weights)
yesIRRhalf_estpdf = yesIRRhalf_kde.evaluate(x_pts)

#plot results and save
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pts, noIRR_estpdf, color = sns.color_palette()[3], label = 'IRR = 1850')
ax.plot(x_pts, yesIRRhalf_estpdf, color = sns.color_palette()[2], label = 'IRR = 2000 (half)')
ax.plot(x_pts, yesIRR_estpdf, color = sns.color_palette()[0], label = 'IRR = 2000')
ax.set_title('Remote TE event frequency')
plt.xlabel('Number of TE events')
plt.ylabel('Density')
#ax.legend(loc = 'upper right')
pl.savefig(figdir + variable + '_remoteweightedPDF_all.png', dpi = 800, bbox_inches = 'tight')

#plot no and yes (IRR = 2000 only) and save
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_pts, noIRR_estpdf, color = sns.color_palette()[3], label = 'IRR = 1850')
#ax.plot(x_pts, yesIRRhalf_estpdf, color = sns.color_palette()[2], label = 'IRR = 2000 (half)')
ax.plot(x_pts, yesIRR_estpdf, color = sns.color_palette()[0], label = 'IRR = 2000')
ax.set_title('Remote TE event frequency')
plt.xlabel('Number of TE events')
plt.ylabel('Density')
ax.legend(loc = 'upper right')
pl.savefig(figdir + variable + '_remoteweightedPDF.png', dpi = 800, bbox_inches = 'tight')
