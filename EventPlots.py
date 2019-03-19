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

def get_Power_Spectra_Subplot(fig, ax, axT, magDict, betaDict):
    ax.text(.98, .8, 'Magnetic Field Power Spectra', horizontalalignment='right', transform=ax.transAxes, 
      bbox=dict(facecolor='white', alpha=0.7))
    ax.set_facecolor('black')
    lowest = 3
    cutoff = 1200
    power, frequency, datesNum, span, N = get_FFT_Power_Info(magDict['Magnitude'], magDict['Epoch'])
    newN = len(frequency) - cutoff
    newFreq = frequency[lowest:newN]
    y,x = np.meshgrid(newFreq, datesNum)
    newPower = power[:,lowest:newN]
    im = ax.pcolormesh(x, y, newPower, cmap='RdBu_r', norm=colors.LogNorm(vmin=0.01, vmax=newPower.max()))
    ax.set_ylim(np.max(newFreq), np.min(newFreq))
    ax.set_xlim(datesNum[0], datesNum[span-1])
    Hours = mdates.HourLocator()   
    Minutes = mdates.MinuteLocator()
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(format_func))
    ax.xaxis.set_major_locator(Hours)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.set_xticks(datesNum[::int(np.ceil(span/15))])
    fig.colorbar(im, ax=ax, fraction = .05, pad = .07)   
    ax.set_ylabel('Period [sec]', labelpad=5)
    
    axT.set_ylabel('Beta Value', rotation=270, labelpad=11)
    axT.plot(betaDict['Epoch'], betaDict['Total'], color='black')
    axT.axhline(y=1, color='black')
    
def get_Particle_Heatmap_SubPlot_EnergyBinned(fig, ax, axT, enBin, TOFxEHDict):
  ax.text(.98, .8, 'Particle Count By Energy Channel', horizontalalignment='right', transform=ax.transAxes, 
       bbox=dict(facecolor='white', alpha=0.7))
  ax.set_facecolor('black')
  ax.set_ylabel('Energy [keV]')
  TOFxE_Epoch = TOFxEHDict['Epoch']
  enBin[np.where(enBin == 0)] = 0.01
  energy = np.linspace(10, 1000, 14)
  x,y = np.meshgrid(TOFxEHDict['Epoch'], energy)
  im = ax.pcolormesh(x, y, enBin.transpose(), cmap='jet', norm=colors.LogNorm(vmin=0.1, vmax=enBin.max()+3e12))
  span = np.ceil(int(TOFxE_Epoch[len(TOFxE_Epoch)-1].timestamp() - TOFxE_Epoch[0].timestamp())/11)
  Hours = mdates.HourLocator()   
  Minutes = mdates.MinuteLocator()
  ax.xaxis.set_major_locator(Minutes)
  ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
  ax.set_xticks(TOFxE_Epoch[::int(span/15)])
  fig.colorbar(im, ax=ax, fraction = .05, pad = .07)
  ax.set_yscale('linear')
  
  axT.plot(TOFxEHDict['Epoch'], TOFxEHDict['FPDU_Density'], color='maroon')  
  axT.set_ylabel('Particle Density', rotation = 270, labelpad = 11)
  
def fmt(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${}*10^{{{}}}$'.format(a, b)

def get_Particle_Heatmap_SubPlot_AngleBinned(fig, ax, axT, angleBin, SC):
  TOFxEHDict = rip.get_CDF_Dict('TOFxEH_'+SC, strtDate, stpDate)
  kappaDict = dbFun.get_Kappa_For_Span('TOFxEH_'+SC, strtDate, stpDate)
  ax.text(.98, .8, 'Particle Count By Pitch Angle', horizontalalignment='right',
          transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7))
  ax.set_facecolor('black')
  ax.set_ylabel('Angle [degrees]')  
  angle = TOFxEHDict['PA_Midpoint'][0,:]
  TOFxE_Epoch = TOFxEHDict['Epoch']
  x,y = np.meshgrid(TOFxE_Epoch, angle)
  im = ax.pcolormesh(x, y, angleBin.transpose(), cmap='jet', vmin=0, 
                     vmax=angleBin.max())
  span = np.ceil(int(TOFxE_Epoch[len(TOFxE_Epoch)-1].timestamp() - 
                                 TOFxE_Epoch[0].timestamp())/11)
  Hours = mdates.HourLocator()   
  Minutes = mdates.MinuteLocator()
  ax.xaxis.set_major_locator(Hours)
  ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
  ax.set_xticks(TOFxE_Epoch[::int(span/15)])
  fig.colorbar(im, ax=ax, fraction = .05, pad = .07, 
               format=ticker.FuncFormatter(fmt))
  ax.set_yscale('linear')
  ax.set_yticks(angle[::4])
  ax.set_yticklabels(angle[::4])
  
  axT.plot(kappaDict['Epoch'], kappaDict['Kappa'], color='oldlace')  
  axT.set_ylabel('Kappa', rotation = 270, labelpad = 11)
  axT.set_ylim(0,2)
    
def get_Mag_Pressure(magnitude):
  mu = (4 * np.pi) * 1e-7
  magVals = magnitude*1e-9
  magPress = (magVals**2)/(2*mu)
  return magPress

def magPress(magDict, maximum):
    magDict['MagPress'] = np.zeros([len(magDict['Epoch']), 4])
    for i in range(magDict['MagPress'].shape[0]):
        magDict['MagPress'][i,0] = get_Mag_Pressure(magDict['Mag'][i,0])
        magDict['MagPress'][i,1] = get_Mag_Pressure(magDict['Mag'][i,1])
        magDict['MagPress'][i,2] = get_Mag_Pressure(magDict['Mag'][i,2])
        magDict['MagPress'][i,3] = get_Mag_Pressure(magDict['Magnitude'][i])
    magDict['MagPress'][np.where(magDict['MagPress'] >= maximum)] = 0
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
