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

def format_func(value, tick_number):
  return '{0:.4g}'.format(1/value)

def fmt(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${}*10^{{{}}}$'.format(a, b)

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

#Takes as arguments the sample period, T, the number of samples, N, the convolving window, 
#as well as a dictionary generated by using the copy method for EMFISIS cdf files.
#This calculated FFTs for the three field aligned coordinate directions, as well as the overall 
#magnitude.  It can generate a rather large dict if given a long enough span, so be careful.
def takeFFT_EMFISISMag(T, N, window, magDict):
    magDict['Frequency'] = np.fft.fftfreq(N)
    magDict['FFTEpoch'] = magDict['Epoch'][int(N/2):len(magDict['Epoch']) - int(N/2)]
    magDict['FFT_Raw'] = np.zeros([4, len(magDict['FFTEpoch']), magDict['Frequency'].shape[0]], dtype='complex64')
    get_field_Aligned_Mag(magDict)
    hWin = int(N/2)
    for i in range(hWin, len(magDict['Epoch']) - hWin):
        for j in range(0,4):
            magDict['FFT_Raw'][j, i - hWin] = np.fft.fft(window*magDict['MagFA'][i - hWin:i + hWin, j])

def get_Power_Spectra_Subplot(fig, ax, magDict, N, T, window, coord='FieldAligned'):
    #Function which takes a matplotlib fig and ax, along with the magDict from either the database at NJIT for EMFISIS
    #on Dgar or from a EMFISIS CDF file, and plots the power spectra for the dict on the axes provided.
        #coord = ('FieldAligned', 'Azimuthal', 'Radial')
        #N = number of samples to include in the FFT
        #T is the sampling period
    FAUV_ToMagDict(magDict)
    get_field_Aligned_Mag(magDict)
    takeFFT_EMFISISMag(T, N, window, magDict)
    if coord == 'FieldAligned':
        coordNumber = 0
    elif coord == 'Azimuthal':
        coordNumber = 1
    elif coord == 'Radial':
        coordNumber = 2
    else:
        raise ValueError('Not a possible coordinate.') 
    ax.text(.98, .8, 'Magnetic Field Power Spectra', horizontalalignment='right', transform=ax.transAxes, 
      bbox=dict(facecolor='white', alpha=0.7))
    ax.set_facecolor('black')
    newFreq = magDict['Frequency'][3:int(len(magDict['Frequency'])/2)]
    y,x = np.meshgrid(newFreq, mdates.date2num(magDict['FFTEpoch']))
    newPower = (1/T)*np.log(np.abs(magDict['FFT_Raw'][coordNumber][:,3:int(len(magDict['Frequency'])/2)])**2)
    im = ax.pcolormesh(x, y, newPower, cmap='RdBu_r')
    ax.set_ylim(np.max(newFreq), np.min(newFreq))
#    Hours = mdates.HourLocator()   
#    Minutes = mdates.MinuteLocator()
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(format_func))
#    ax.xaxis.set_major_locator(Hours)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
#    ax.set_xticks(datesNum[::int(np.ceil(len(magDict['Epoch'])/15))])
    fig.colorbar(im, ax=ax, fraction = .05, pad = .07)   
  
def get_Beta_LinePlot(ax, betaDict):
    ax.plot(betaDict['Epoch'], betaDict['Total'], color='black', label='Beta')

def get_Particle_Heatmap_SubPlot_EnergyBinned(fig, ax, tofDict):
    #Function which takes a matplotlib fig and ax, along with the tofDict from either the database at NJIT for EMFISIS
    enBin = np.sum(tofDict['FPDU'], 2)
    enBin[np.where(enBin == 0)] = .01
    ax.text(.98, .8, 'Particle Count By Energy Channel', horizontalalignment='right', transform=ax.transAxes, 
            bbox=dict(facecolor='white', alpha=0.7))
    ax.set_facecolor('black')
    energy = np.linspace(10, 1000, 14)
    x,y = np.meshgrid(tofDict['Epoch'], energy)
    im = ax.pcolormesh(x, y, np.log(enBin.transpose()), cmap='jet', vmax=4*(np.log(enBin)).std())
    span = np.ceil(int(tofDict['Epoch'][len(tofDict['Epoch'])-1].timestamp() - tofDict['Epoch'][0].timestamp())/11)
    #Hours = mdates.HourLocator()   
    #Minutes = mdates.MinuteLocator()
    #ax.xaxis.set_major_locator(Minutes)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.set_xticks(tofDict['Epoch'][::int(span/15)])
    fig.colorbar(im, ax=ax, fraction = .05, pad = .07)
    ax.set_yscale('linear')

def get_ParticleDensity_LinePlot(ax, tofDict):
    ax.plot(tofDict['Epoch'], tofDict['FPDU_Density'], color='maroon')  

def get_Particle_Heatmap_AngleBinned(fig, ax, tofDict):
    angleBin = np.sum(tofDict['FPDU'], 1)
    ax.text(.98, .8, 'Particle Count By Pitch Angle', horizontalalignment='right',
            transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7))
    ax.set_facecolor('black')
    angle = tofDict['PA_Midpoint'][0,:]
    TOFxE_Epoch = tofDict['Epoch']
    x,y = np.meshgrid(tofDict['Epoch'], angle)
    im = ax.pcolormesh(x, y, angleBin.transpose(), cmap='jet', vmin=0, vmax=4*angleBin.std())
    span = np.ceil(int(tofDict['Epoch'][len(tofDict['Epoch'])-1].timestamp() - tofDict['Epoch'][0].timestamp())/11)
    Hours = mdates.HourLocator()   
    ax.xaxis.set_major_locator(Hours)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.set_xticks(tofDict['Epoch'][::int(span/15)])
    ax.set_yscale('linear')
    fig.colorbar(im, ax=ax, fraction = .05, pad = .07, format=ticker.FuncFormatter(fmt))

def get_Kappa_Lineplot(ax, kappaDict):
    ax.plot(kappaDict['Epoch'], kappaDict['Kappa'], color='blue', label='Kappa')  