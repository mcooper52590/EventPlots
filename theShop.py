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

date = dt.datetime(2013,5,1,12,0,0)
strtDate = date - dt.timedelta(minutes=135)
stpDate = date + dt.timedelta(minutes=180) 
magDict = rip.get_CDF_Dict('Mag_1Sec_A', strtDate, stpDate)
betaDict = dbFun.get_Beta_Total_For_Span('TOFxEH_A', strtDate, stpDate)
kappaDict = dbFun.get_Kappa_For_Span('TOFxEH_A', strtDate, stpDate)

FAUV_ToMagDict(magDict)
get_field_Aligned_Mag(magDict)
takeFFT_EMFISISMag(1, 2048, blackmanharris(2048), magDict)
mag_AverageRemoved(600, magDict)
   
gs = gridspec.GridSpec(2, 4)

ax1 = plt.subplot(gs[0:1, 0:1])
ax2 = plt.subplot(gs[1:2, 0:1])
ax3 = plt.subplot(gs[0:1, 1:])
ax4 = plt.subplot(gs[1:2, 1:])

eve.get_Position_Plot_For_Span(ax1, ax2, strtDate, stpDate)

eve.FAUV_ToMagDict(magDict)
eve.takeFFT_EMFISISMag(1, 2048, blackmanharris(2048), magDict)

ax3.plot(magDict['Epoch_Avg'], magDict['MagFA_NoBack'][:,0], color = 'red', label='Field-Aligned Mag')
ax3.set_ylim(-20,20)
ax3.legend()
#ax1.set_facecolor('black')


ax4.plot(betaDict['Epoch'], betaDict['Total'], color='black', label='Beta')
ax4.plot(kappaDict['Epoch'], kappaDict['Kappa'], color='blue', label='Kappa') 
ax4.set_ylim(-2,2)
ax4.legend()

gs.update(hspace=.25, wspace=.15)

