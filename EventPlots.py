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
date = dt.datetime(2013,5,1,12,0,0)
strtDate = date - dt.timedelta(minutes=360)
stpDate = date + dt.timedelta(minutes=360) 
#==============================================================================
magDict = rip.get_CDF_Dict('Mag_1Sec_A', strtDate, stpDate)
tofDict = rip.get_CDF_Dict('TOFxEH_A', strtDate, stpDate)
tofDict['FPDU'][np.where(tofDict['FPDU']<0)] = 0
betaDict = dbFun.get_Beta_Total_For_Span('TOFxEH_A', strtDate, stpDate)
kappaDict = dbFun.get_Kappa_For_Span('TOFxEH_A', strtDate, stpDate)


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, sharey=False, figsize=(17,12), dpi=166) 
ax1T = ax1.twinx()
ax2T = ax2.twinx()
ax3T = ax3.twinx()
pos1 = list(ax1.get_position().bounds)
fig.subplots_adjust(hspace=0.04)
#===============================================================================================================
T = 1
N = 2048
window = blackmanharris(N)
ax1.set_ylabel('Period [s]', labelpad = 11)
eve.get_Power_Spectra_Subplot(fig, ax1, magDict, N, T, window)

#ax1T.set_ylabel('Beta Value', rotation=270, labelpad=11)
    

ax2.set_ylabel('Energy [keV]', labelpad=11)
eve.get_Particle_Heatmap_SubPlot_EnergyBinned(fig, ax2, tofDict)  

ax2T.set_ylabel('Particle Density', rotation = 270, labelpad = 11)
eve.get_ParticleDensity_LinePlot(ax2T, tofDict)

ax3.set_ylabel('Angle [degrees]', labelpad=11)  
eve.get_Particle_Heatmap_AngleBinned(fig, ax3, tofDict)

