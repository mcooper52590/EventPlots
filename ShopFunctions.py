#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 10:44:53 2019

@author: mattcooper
"""
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#This is the Shop Functions file.  Place function calls/library imports into this file and execute
#it before attempting to run theShop file.

import spacepy.pycdf as cdf
import os
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

def get_Mag_Pressure(magnitude):
  mu = (4 * np.pi) * 1e-7
  magVals = magnitude*1e-9
  magPress = (magVals**2)/(2*mu)
  return magPress

#Field Aligned Unit Vector (FAUV)
#Takes a dict from the EMFISIS data and adds the unit vectors for field-aligned coordinates
def FAUV_ToMagDict(magDict):
    magDict['FAUV'] = np.zeros([len(magDict['Epoch']), 3, 3])
    for i in range(0, len(magDict['Epoch'])):
        #Z Direction 
        magDict['FAUV'][i,0,:] = magDict['Mag'][i,:]/np.linalg.norm(magDict['Mag'][i,:])
        unitVec_Earth = magDict['coordinates'][i,:]/np.linalg.norm(magDict['coordinates'][i,:])
        #Azimuthal Direction
        magDict['FAUV'][i,1,:] = (-np.cross(magDict['FAUV'][i,0,:], unitVec_Earth)/
                                   np.linalg.norm(np.cross(magDict['FAUV'][i,0,:], unitVec_Earth)))
        #Radial Direction
        magDict['FAUV'][i,2,:] = (np.cross(magDict['FAUV'][i,0,:], magDict['FAUV'][i,1,:])/
                                  np.linalg.norm(np.cross(magDict['FAUV'][i,0,:], magDict['FAUV'][i,1,:])))

#Takes a dict from the EMFISIS data and adds the unit vectors for field-aligned coordinates, as well as
#adding the field values in the field-aligned coordinate system.
def get_field_Aligned_Mag(magDict):
    magDict['MagFA'] = np.zeros([len(magDict['Epoch']), 4])
    for i in range(len(magDict['Epoch'])):
        magDict['MagFA'][i,0] = np.dot(magDict['Mag'][i,:], magDict['FAUV'][i,0,:])
        magDict['MagFA'][i,1] = np.dot(magDict['Mag'][i,:], magDict['FAUV'][i,1,:])
        magDict['MagFA'][i,2] = np.dot(magDict['Mag'][i,:], magDict['FAUV'][i,2,:]) 
        magDict['MagFA'][i,3] = magDict['Magnitude'][i]
magCdf = cdf.CDF('/home/mattcooper/EventPlots/CDF_ForDict/rbsp-a_magnetometer_1sec-sm_emfisis-L3_20180823_v1.6.1.cdf')
magDict = dict(magCdf.copy())
magCdf.close()
FAUV_ToMagDict(magDict)
get_field_Aligned_Mag(magDict)

#dates = [dt.datetime(2013,5,1,12,30,16), dt.datetime(2013,5,1,12,50,16),.5]
##2013 05/01, 1230:16, 16.5 minute duration
#strtDate = dates[0] - dt.timedelta(minutes=180)
#stpDate = dates[1] + dt.timedelta(minutes=180) 
##==============================================================================
#magDict = rip.get_CDF_Dict('Mag_1Sec_A', strtDate, stpDate)
#TOFxEHDict = rip.get_CDF_Dict('TOFxEH_A', strtDate, stpDate)
#magPress = get_MagPress_BelowThreshold(300, magDict['Magnitude'])
#betaDict = dbFun.get_Beta_Total_For_Span('TOFxEH_A', strtDate, stpDate)
#pressDict = dbFun.get_Pressure_Total_For_Span('TOFxEH_A', strtDate, stpDate)
#kappaDict = dbFun.get_Kappa_For_Span('TOFxEH_A', strtDate, stpDate)


    

   

































#==========================================================================================
#Functions that were in old code that didn't seem to have a purpose
#Saved them here in kind of a junk drawer in case their use arose.
#I didn't document shit back then, so I actually don't know what some of these do.
#def sort_Mag(magEntries):
#  cutoff = 500
#  mag_Epoch = []
#  mag_Mag = []
#  mag_Coord = []
#  for entry in magEntries:      
#      mag_Epoch.append(entry[0])
#      if entry[1] > cutoff:
#          mag_Mag.append(cutoff)
#      else:
#          mag_Mag.append(entry[1])
#      mag_Coord.append(entry[12]/6371.2)
#  mag_Epoch = np.array(mag_Epoch)
#  mag_Mag = np.array(mag_Mag)
#  mag_Coord = np.array(mag_Coord)
#  return mag_Epoch, mag_Mag, mag_Coord