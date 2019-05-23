#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  3 14:32:44 2019

@author: matthew
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import cdfstuff.ripper_V2 as rip
import dbstuff.databaseFunctions as dbFun
import scipy.signal as sig
import EventPlotFunctions as eve

def hpf(data, fmin=.01, T=1, init=0):
    '''
    A simple high pass filter function.
        data = the numpy array containing the data to filter
        fmin = lower cutoff frequency
        T = sample period
        init = inital boundary condition
    ''' 
    alpha = 1/(1 + 2*np.pi*fmin*T)
    filtData = np.zeros([data.shape[0]])
    filtData[0] = 0 
    for i in range(1, filtData.shape[0]):
        filtData[i] = alpha * filtData[i-1] + alpha*(data[i] - data[i-1])
    return filtData

def lpf(data, fmax=.005, T=1, initial=0):
    '''
    A simple low pass filter function.
        data = the numpy array containing the data to filter
        fmax = upper cutoff frequency
        T = sample period
        init = inital boundary condition
    ''' 
    alpha = (2*np.pi*fmax*T)/(1 + 2*np.pi*fmax*T)
    filtData = np.zeros([data.shape[0]])
    filtData[0] = 0
    alpha*(initial + data[1] - data[0])  
    for i in range(1, filtData.shape[0]):
        filtData[i] = alpha*(data[i] - filtData[i-1]) + filtData[i-1]
    return filtData

def get_Mag_Press(magnitude):
  mu = (4 * np.pi) * 1e-7
  return ((magnitude*1e-9)**2)/(2*mu)

pi = np.pi
diff = 200
SC = 'B'
date = dt.datetime(2013,5,1,8,0,0)
strtDate = date - dt.timedelta(minutes=diff)
stpDate = date + dt.timedelta(minutes=diff) 
#=======================================================================================================================
magDict = rip.get_CDF_Dict('Mag_1Sec_' + SC, strtDate, stpDate)
tofDict = rip.get_CDF_Dict('TOFxEH_' + SC, strtDate, stpDate)
tofDict['FPDU'][np.where(tofDict['FPDU']<0)] = 0

eve.FAUV_ToMagDict(magDict)
eve.get_field_Aligned_Mag(magDict)

wp = [2*pi*8e-4, 2*pi*8e-3]
ws = [2*pi*4e-4, 2*pi*1.6e-2]
gpass = .1
gstop = 10
N, Wn = sig.buttord(wp, ws, gpass, gstop, False)
b, a = sig.butter(N, Wn, 'band', False)
w, h = sig.freqs(b, a)
#plt.plot(w, 20 * np.log10(abs(h)))
#plt.xscale('log')
#plt.title('Butterworth filter frequency response')
#plt.xlabel('Frequency [radians / second]')
#plt.ylabel('Amplitude [dB]')
#plt.margins(0, 0.1)
#plt.grid(which='both', axis='both')
#plt.axvline(Wn[0], color='green') # cutoff frequency
#plt.axvline(Wn[1], color='green') # cutoff frequency
#plt.show()

compMag = magDict['MagFA'][:,0]
#lpfMag = compMag
##lpfMag = lpf(compMag, fmax=.005, T=1, initial=0)
##lpfMag = lpf(lpfMag, fmax=.005, T=1, initial=0)
N, Wn = sig.buttord(wp, ws, gpass, gstop, False)
b, a = sig.butter(N, Wn, 'band', False)
w, h = sig.freqs(b, a)
filtMagPress = sig.lfilter(b, a,get_Mag_Press(magDict['MagFA'][:,0]))

N, Wn = sig.buttord(wp, ws, gpass, gstop, False, fs = 1/11)
b, a = sig.butter(N, Wn, 'band', False, 11)
w, h = sig.freqs(b, a)
filtPressPara = sig.lfilter(b, a, tofDict['FPDU_ParaPressure']*1e-9)
filtPressPerp = sig.lfilter(b, a, tofDict['FPDU_PerpPressure']*1e-9)


plt.plot(magDict['Epoch'], filtMagPress)
plt.plot(tofDict['Epoch'], filtPressPara)
plt.plot(tofDict['Epoch'], filtPressPerp)
#plt.ylim(0,10)
#
#
#samps = int(4096)
#window = sig.windows.blackmanharris(samps)
#hs = int(samps/2)
#freqs = np.fft.fftfreq(samps)
#dynPS = np.zeros([hs, len(compMag)-samps], dtype='complex64')
#for i in range(hs, len(compMag)-hs):
#    dynPS[:,i-hs] = np.fft.fft(filtMag[i-hs:i+hs])[:hs]
#    
#fig, ax1 = plt.subplots(1, 1, sharex=True, sharey=False, figsize=(17,12), dpi=166) 
#img = ax1.imshow(np.log(np.abs(dynPS))[:250], aspect=25)
#fig.colorbar(img, ax=ax1)


