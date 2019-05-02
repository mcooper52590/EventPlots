#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 10:44:53 2019

@author: mattcooper
"""
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

def get_Interp_Mag_Press(pressDict, magDict):
    magEpochNum = [entry.timestamp() for entry in magDict['Epoch']]
    TOFxEEpochNum = [entry.timestamp() for entry in pressDict['Epoch']]
    magInterpd = np.interp(TOFxEEpochNum, magEpochNum, magDict['Magnitude'])*1e-9
    mu = (4 * np.pi) * 1e-7
    magPress = (magInterpd**2)/(2*mu)
    return magPress

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
diff = 240
date = dt.datetime(2013,5,1,12,00,0)
strtDate = date - dt.timedelta(minutes=diff)
stpDate = date + dt.timedelta(minutes=diff) 

magDict = rip.get_CDF_Dict('Mag_1Sec_A', strtDate, stpDate)
betaDict = dbFun.get_Beta_Total_For_Span('TOFxEH_A', strtDate, stpDate)
kappaDict = dbFun.get_Kappa_For_Span('TOFxEH_A', strtDate, stpDate)
pressDict = dbFun.get_Pressure_Total_For_Span('TOFxEH_A', strtDate, stpDate)

eve.FAUV_ToMagDict(magDict)
eve.get_field_Aligned_Mag(magDict)
magDict['magPressFA'] = eve.get_Mag_Press(magDict['MagFA'])
#eve.mag_AverageRemoved(600, magDict)


pressDict['InterpMagPress'] = get_Interp_Mag_Press(pressDict, magDict)
#beta = pressDict['Total']/interpMag
N = 100
hw = int(N/2)
pressDict['Epoch_Avg'] = pressDict['Epoch'][hw:pressDict['Epoch'].shape[0] - hw] 
magDict['MagPress_Avg'] = np.zeros([pressDict['Epoch_Avg'].shape[0], 4])
magDict['Perp_Avg'] = np.zeros([pressDict['Epoch_Avg'].shape[0], 4])

#for i in range(hw, magDict['Epoch_Avg'].shape[0] - hw):
#    for j in range(0,3):



plt.plot(pressDict['Epoch'], pressDict['InterpMagPress'], color='red')
plt.plot(pressDict['Epoch'], pressDict['Perpendicular'], color='black', linestyle='dashed')
plt.plot(pressDict['Epoch'], pressDict['Parallel'], color='black')
#plt.plot(pressDict['Epoch'], beta, color='green')