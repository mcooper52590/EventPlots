#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Tue Mar 19 10:55:07 2019

@author: matthew
'''
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import datetime as dt
import cdfstuff.ripper_V2 as rip
import dbstuff.databaseFunctions as dbFun
import EventPlotFunctions as eve
from scipy.signal import blackmanharris
import numpy as np

date = dt.datetime(2013,3,1,0,0,0)
strtDate = date - dt.timedelta(minutes=180)
stpDate = date + dt.timedelta(minutes=180) 
magDict = rip.get_CDF_Dict('Mag_1Sec_A', strtDate, stpDate)
betaDict = dbFun.get_Beta_Total_For_Span('TOFxEH_A', strtDate, stpDate)
kappaDict = dbFun.get_Kappa_For_Span('TOFxEH_A', strtDate, stpDate)

eve.FAUV_ToMagDict(magDict)
eve.get_field_Aligned_Mag(magDict)
eve.mag_AverageRemoved(600, magDict)
   
gs = gridspec.GridSpec(3, 4)

ax1 = plt.subplot(gs[0:2, 0:2])
ax2 = plt.subplot(gs[0:2, 2:4])
ax3 = plt.subplot(gs[2, :])
ax3T = ax3.twinx()

eve.get_Position_Plot_For_Span(ax1, ax2, strtDate, stpDate)


eve.takeFFT_EMFISISMag(1, 2048, blackmanharris(2048), magDict)

l1 = ax3.plot(magDict['Epoch_Avg'], magDict['MagFA_NoBack'][:,0], color = 'red', label='Field-Aligned Mag')
ax3.set_ylim(-20,20)
#ax1.set_facecolor('black')

ax3T.set_ylim(-2,2)

l2 = ax3T.plot(betaDict['Epoch'], betaDict['Total'], color='black', label='Beta')
l3 = ax3T.plot(kappaDict['Epoch'], kappaDict['Kappa'], color='blue', label='Kappa') 
#gs.update(hspace=.5, wspace=.15)

lns = l1+l2+l3
labs = [l.get_label() for l in lns]
ax3.legend(lns, labs, loc=0, prop={'size': 6})