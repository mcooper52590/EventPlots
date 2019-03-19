#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  8 13:08:18 2018

@author: matthew
"""
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import cdfstuff.ripper_V2 as rip
import mysql.connector as conn
import dbstuff.databaseFunctions as dbFun
from scipy.signal import blackmanharris
from scipy.fftpack import fft
import matplotlib.dates as mdates
import matplotlib.colors as colors
import matplotlib.ticker as ticker
#==============================================================================
#==============================================================================
def format_func(value, tick_number):
  return '{0:.4g}'.format(1/value)
  
def fmt(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${}*10^{{{}}}$'.format(a, b)


#==============================================================================
#==============================================================================
dates = [dt.datetime(2013,5,1,12,30,16), dt.datetime(2013,5,1,12,50,16),.5]
#2013 05/01, 1230:16, 16.5 minute duration
strtDate = dates[0] - dt.timedelta(minutes=180)
stpDate = dates[1] + dt.timedelta(minutes=180) 
#==============================================================================
magDict = rip.get_CDF_Dict('Mag_1Sec_A', strtDate, stpDate)
magPress(magDict, 300)
TOFxEHDict = rip.get_CDF_Dict('TOFxEH_A', strtDate, stpDate)
betaDict = dbFun.get_Beta_Total_For_Span('TOFxEH_A', strtDate, stpDate)
pressDict = dbFun.get_Pressure_Total_For_Span('TOFxEH_A', strtDate, stpDate)
kappaDict = dbFun.get_Kappa_For_Span('TOFxEH_A', strtDate, stpDate)
#==============================================================================
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, sharey=False, 
     figsize=(17,12), dpi=166) 
ax1T = ax1.twinx()
ax2T = ax2.twinx()
ax3T = ax3.twinx()
fig.subplots_adjust(hspace=0.02)
ax3.set_xlabel('Time Stamp')
#------------------------------------------------------------------------------
#==============================================================================
holder = TOFxEHDict['FPDU']
holder[np.where(holder<0)] = 0
binnedByEnergy = np.sum(holder, 2)
binnedByAngle = np.sum(holder, 1)

get_Particle_Heatmap_SubPlot_EnergyBinned(fig, ax1, ax1T, binnedByEnergy, 
                                          TOFxEHDict)
get_Particle_Heatmap_SubPlot_AngleBinned(fig, ax2, ax2T, binnedByAngle, 
                                         'A')
get_Power_Spectra_Subplot(fig, ax3, ax3T, magDict, betaDict)

pos1 = list(ax1.get_position().bounds)
pos1T = list(ax1T.get_position().bounds)
pos2T = list(ax2T.get_position().bounds)
pos3T = list(ax3T.get_position().bounds)
pos1T[2] = pos1[2]
pos2T[2] = pos1[2]
pos3T[2] = pos1[2]
ax1T.set_position(pos1T)
ax2T.set_position(pos2T)
ax3T.set_position(pos3T)

ax1.axvline(x=dates[0], color='black')
ax1.axvline(x=dates[1], color='black')
ax2.axvline(x=dates[0], color='black')
ax2.axvline(x=dates[1], color='black')
ax3.axvline(x=dates[0], color='black')
ax3.axvline(x=dates[1], color='black')
plt.plot()

Dir = '/home/matthew/CSTR/PlotFolder/FixedPlots/'
fig.savefig(Dir + 'TEST.jpeg')

#==============================================================================
