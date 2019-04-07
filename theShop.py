#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 10:44:53 2019

@author: mattcooper
"""
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#This is the Shop Floor.  Place codes here to work on/test them.  Place function calls/libraries 
#in the ShopFunctions.py file in this same directory.

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import datetime as dt
import cdfstuff.ripper_V2 as rip
import dbstuff.databaseFunctions as dbFun
import EventPlotFunctions as eve
from scipy.signal import blackmanharris
import numpy as np

diff = 180
date = dt.datetime(2013,8,23,6,0,0)
strtDate = date - dt.timedelta(minutes=diff)
stpDate = date + dt.timedelta(minutes=diff) 

magDict = rip.get_CDF_Dict('Mag_1Sec_A', strtDate, stpDate)
betaDict = dbFun.get_Beta_Total_For_Span('TOFxEH_A', strtDate, stpDate)
kappaDict = dbFun.get_Kappa_For_Span('TOFxEH_A', strtDate, stpDate)



eve.FAUV_ToMagDict(magDict)
eve.get_field_Aligned_Mag(magDict)
eve.takeFFT_EMFISISMag(1, 2048, blackmanharris(2048), magDict)
eve.mag_AverageRemoved(600, magDict)
    


for k in range(0,3):
    magDict['MagFA_NoBack'][:,k] = eve.lpf(magDict['MagFA_NoBack'][:,k]) 
    magDict['MagFA_NoBack'][:,k] = eve.hpf(magDict['MagFA_NoBack'][:,k])
    
    
gs = gridspec.GridSpec(2, 4)

ax1 = plt.subplot(gs[0:1, 0:1])
ax2 = plt.subplot(gs[1:2, 0:1])
ax3 = plt.subplot(gs[0:1, 1:])
ax4 = plt.subplot(gs[1:2, 1:])

eve.get_Position_Plot_For_Span(ax1, ax2, strtDate, stpDate)
epBin1 = np.linspace(0,1, len(magDict['Epoch_Avg']))
ax3.plot(epBin1, magDict['MagFA_NoBack'][:,0], color = 'red', label='Field-Aligned Mag')
ax3.set_ylim(-20,20)
ax3.legend()
#ax1.set_facecolor('black')

epBin2 = np.linspace(0,1, len(betaDict['Epoch']))
ax4.plot(epBin2, betaDict['Total'], color='green', label='Beta')
ax4.plot(epBin2, kappaDict['Kappa'], color='blue', label='Kappa') 
ax4.set_ylim(-2,3)
ax4.legend()

gs.update(hspace=.25, wspace=.15)

ticks = []
locs = []
ticksep = 2*diff/18
for i in range(0, 19):
    ticks.append((strtDate + dt.timedelta(minutes = i*ticksep)).strftime('%H:%M'))
    locs.append(i*(1/18))

ax3.set_xticks(locs)
ax3.set_xticklabels(ticks, rotation=45)
ax4.set_xticks(locs)
ax4.set_xticklabels(ticks, rotation=45)
ax4.axhline(0, color='blue', linestyle='--')
ax4.axhline(1, color='green', linestyle='--')