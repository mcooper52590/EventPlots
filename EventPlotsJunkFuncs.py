#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 09:45:48 2019

@author: matthew
"""
#==========================================================================================
#Functions that were in old code that didn't seem to have a purpose
#Saved them here in kind of a junk drawer in case their use arose.
#I didn't document shit back then, so I actually don't know what some of these do.
def get_field_Aligned_UnitVecs(MAGDICT):
    unitVec_Zeta = []
    unitVec_Radial = []
    unitVec_Azimuthal = [] 
    for i in range(len(MAGDICT['Epoch'])):  
        unitVec_Z = MAGDICT['Mag'][i,:]/np.linalg.norm(MAGDICT['Mag'][i,:])
        unitVec_Earth = MAGDICT['coordinates'][i,:]/np.linalg.norm(MAGDICT['coordinates'][i,:])
        vec_Azi = -np.cross(unitVec_Z, unitVec_Earth)
        unitVec_Azi = vec_Azi/np.linalg.norm(vec_Azi)
        vec_Rad = np.cross(unitVec_Z, unitVec_Azi) 
        unitVec_Zeta.append(unitVec_Z)
        unitVec_Azimuthal.append(unitVec_Azi)
        unitVec_Radial.append(vec_Rad/np.linalg.norm(vec_Rad))    
    unitVec_Zeta = np.array(unitVec_Zeta)
    unitVec_Radial = np.array(unitVec_Radial)
    unitVec_Azimuthal = np.array(unitVec_Azimuthal)
    return unitVec_Zeta, unitVec_Radial, unitVec_Azimuthal

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
  
  
  
  
  
  