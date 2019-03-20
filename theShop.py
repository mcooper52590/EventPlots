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

magDictA = rip.get_CDF_Dict('Mag_1Sec_A', strtDate, stpDate)
magDictB = rip.get_CDF_Dict('Mag_1Sec_B', strtDate, stpDate)

betaDictA = dbFun.get_Beta_Total_For_Span('TOFxEH_A', strtDate, stpDate)
betaDictB = dbFun.get_Beta_Total_For_Span('TOFxEH_B', strtDate, stpDate)

kappaDictA = dbFun.get_Kappa_For_Span('TOFxEH_A', strtDate, stpDate)
kappaDictB = dbFun.get_Kappa_For_Span('TOFxEH_B', strtDate, stpDate)


eve.FAUV_ToMagDict(magDictA)
eve.get_field_Aligned_Mag(magDictA)
eve.takeFFT_EMFISISMag(1, 2048, blackmanharris(2048), magDictA)
eve.mag_AverageRemoved(600, magDictA)
    
eve.FAUV_ToMagDict(magDictB)
eve.get_field_Aligned_Mag(magDictB)
eve.takeFFT_EMFISISMag(1, 2048, blackmanharris(2048), magDictB)
eve.mag_AverageRemoved(600, magDictB)


for k in range(0,3):
    magDictA['MagFA_NoBack'][:,k] = eve.lpf(magDictA['MagFA_NoBack'][:,k])
    magDictB['MagFA_NoBack'][:,k] = eve.lpf(magDictB['MagFA_NoBack'][:,k])   
    magDictA['MagFA_NoBack'][:,k] = eve.hpf(magDictA['MagFA_NoBack'][:,k])
    magDictB['MagFA_NoBack'][:,k] = eve.hpf(magDictB['MagFA_NoBack'][:,k])      
    
gs = gridspec.GridSpec(2, 4)

ax1 = plt.subplot(gs[0:1, 0:1])
ax2 = plt.subplot(gs[1:2, 0:1])
ax3 = plt.subplot(gs[0:1, 1:])
ax4 = plt.subplot(gs[1:2, 1:])

eve.get_Position_Plot_For_Span(ax1, ax2, strtDate, stpDate)

ax3.plot(magDictA['Epoch_Avg'], magDictA['MagFA_NoBack'][:,0], color = 'blue', label='Field-Aligned Mag A')
ax3.plot(magDictB['Epoch_Avg'], magDictB['MagFA_NoBack'][:,0], color = 'green', label='Field-Aligned Mag B')
ax3.set_ylim(-20,20)
ax3.legend()
#ax1.set_facecolor('black')


ax4.plot(betaDictA['Epoch'], betaDictA['Total'], color='royalblue', label='Beta A')
ax4.plot(betaDictB['Epoch'], betaDictB['Total'], color='turquoise', label='Beta B')
ax4.plot(kappaDictA['Epoch'], kappaDictA['Kappa'], color='maroon', label='Kappa A') 
ax4.plot(kappaDictB['Epoch'], kappaDictB['Kappa'], color='lightcoral', label='Kappa B') 
ax4.set_ylim(-2,2)
ax4.legend()

gs.update(hspace=.25, wspace=.35)

