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
import dbstuff.databaseFunctions as dbFun
from scipy.signal import blackmanharris
import EventPlotFunctions as eve
'''
========================================================================================================================
Code section which retrieves a dictionary from the mySQL database on Dgar at NJIT
'''
diff = 240
date = dt.datetime(2013,5,1,21,0,0)
strtDate = date - dt.timedelta(minutes=diff)
stpDate = date + dt.timedelta(minutes=diff) 
#=======================================================================================================================
magDict = rip.get_CDF_Dict('Mag_1Sec_A', strtDate, stpDate)
tofDict = rip.get_CDF_Dict('TOFxEH_A', strtDate, stpDate)
tofDict['FPDU'][np.where(tofDict['FPDU']<0)] = 0
betaDict = dbFun.get_Beta_Total_For_Span('TOFxEH_A', strtDate, stpDate)
kappaDict = dbFun.get_Kappa_For_Span('TOFxEH_A', strtDate, stpDate)

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, sharey=False, figsize=(17,12), dpi=166) 
ax2T = ax2.twinx()
pos1 = list(ax1.get_position().bounds)
fig.subplots_adjust(hspace=0.04)
#=======================================================================================================================
T = 1
N = 2048
window = blackmanharris(N)
eve.get_Power_Spectra_Subplot(fig, ax1, magDict, N, T, window)
ax1.text(-.08, 0.5, 'Period [s]', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, 
         rotation='vertical')

ax2.text(-.08, 0.6, 'Energy [keV]', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, 
         rotation='vertical')
eve.get_Particle_Heatmap_SubPlot_EnergyBinned(fig, ax2, tofDict)  

ax2T.set_ylabel('Particle Density', rotation = 270, labelpad = 11)
eve.get_ParticleDensity_LinePlot(ax2T, tofDict)

ax3.text(-.08, 0.6, 'Angle [degrees]', horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, 
         rotation='vertical')
eve.get_Particle_Heatmap_AngleBinned(fig, ax3, tofDict)

ticks = []
locs = []
ticksep = 2*diff/16
for i in range(0, 17):
    ticks.append((strtDate + dt.timedelta(minutes = i*ticksep)).strftime('%H:%M'))
    locs.append(i*(1/16))

ax3.set_xticks(locs)
ax3.set_xticklabels(ticks, rotation=45)
    


