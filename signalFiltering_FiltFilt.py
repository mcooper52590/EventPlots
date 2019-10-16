#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 19:52:54 2019

@author: matthew
"""
import scipy.signal as sig
import numpy as np
import matplotlib.pyplot as plt

sin = np.sin
cos = np.cos
pi = np.pi
zeros = np.zeros
rand=np.random.rand
fft = np.fft.fft

def srand():
  return rand()-rand()

def cta(freq):
    return 2*pi*freq

'''
'''
sampleFreq = int(5e4)
fNy=sampleFreq/2
pad = 100
hpad = int(pad/2)
samSpc = 1/sampleFreq
spaceLen = 10
t = np.linspace(0,spaceLen,sampleFreq*spaceLen)
tpad = np.linspace(-hpad*samSpc, spaceLen+(hpad*samSpc), sampleFreq*spaceLen+pad)

'''
'''
Channels = 20
fbase = 5e3
FDMSig = zeros([t.shape[0]+pad, Channels + 1])
FDMSigNoise = zeros([t.shape[0]+pad, Channels + 1])
for i in range(1, 1 + Channels):
  FDMSig[hpad:(t.shape[0]+hpad),i]=(FDMSig[hpad:(t.shape[0]+hpad),i] + 
                                    sin(2*pi*(fbase+100*i)*t))
  FDMSig[hpad:(t.shape[0]+hpad),0]=(FDMSig[hpad:(t.shape[0]+hpad),0] + 
                                    FDMSig[hpad:(t.shape[0]+hpad),i])
  FDMSigNoise[hpad:(t.shape[0]+hpad),i]=(FDMSigNoise[hpad:(t.shape[0]+hpad),i] + 
                                    (1+srand())*sin(2*pi*(fbase*+100*i)*t + 2*pi*rand()))
  FDMSigNoise[hpad:(t.shape[0]+hpad),0]=(FDMSigNoise[hpad:(t.shape[0]+hpad),0] + 
                                    FDMSigNoise[hpad:(t.shape[0]+hpad),i])


'''
'''
fcarry = 1e4
carrySig = np.zeros([t.shape[0] + pad])
carrySig[hpad:(t.shape[0]+hpad)] = cos(2*pi*fcarry*t)
sigDown = carrySig*FDMSig[:,0]


'''
'''
#fig, axes = plt.subplots(3, 2, sharex=False, sharey=False, figsize=(17,12), dpi=166) 
#axes[0,0].plot(tpad, FDMSig[:,0], linewidth=.3, color='black') 
#axes[1,0].plot(tpad, carrySig, linewidth=.3, color='black') 
#axes[2,0].plot(tpad, sigDown, linewidth=.3, color='black')
#axes[0,0].set_xlim(0,.02)
#axes[1,0].set_xlim(0,.005)
#axes[2,0].set_xlim(0,.02)
#
#N = 64*4096
#hN = int(N/2)
#freqs = np.fft.fftfreq(N, samSpc)
#SigFFT = fft(FDMSig[:N,0])
#axes[0,1].plot(freqs[:hN], np.abs(SigFFT)[:hN], linewidth=.1, color='black')
#axes[0,1].set_yticklabels([])
#axes[0,1].tick_params(left="off")
#
#carryFFT = fft(carrySig[:N])
#axes[1,1].plot(freqs[:hN], np.abs(carryFFT)[:hN], linewidth=.1, color='black')
#axes[1,1].set_yticklabels([])
#axes[1,1].tick_params(left="off")
#
#downFFT = fft(sigDown[:N])
#axes[2,1].plot(freqs[:hN], np.abs(downFFT)[:hN], linewidth=.1, color='black')
#axes[2,1].set_yticklabels([])
#axes[2,1].tick_params(left="off")
#axes[0,0].set_title('Time Domain')
#axes[0,1].set_title('Frequency Domain')
#fig.subplots_adjust(hspace=.5)

'''
'''
fig, axes = plt.subplots(3, 2, sharex=False, sharey=False, figsize=(17,12), dpi=166)
fig.subplots_adjust(hspace=.45)

wNy = cta(fNy)
ubp = cta(6500)/wNy
ubs = cta(ubp + 2000)/wNy
N, Wn = sig.buttord(ubp, ubs, .01, 100, False, fs = sampleFreq)
b, a = sig.butter(N, Wn, 'lowpass', False)
filtSig = sig.lfilter(b, a, sigDown)
filtfiltSig = sig.filtfilt(b, a, sigDown)

axes[0,0].set_title('Signal')
axes[0,0].plot(tpad, sigDown, linewidth=.3, color='black')
axes[0,0].set_xlim(5.0775,5.0825)
axes[0,0].set_xticklabels([])
axes[0,0].tick_params(bottom="off")
fig.subplots_adjust(hspace=.45)

axes[1,0].set_title('lfilter()')
axes[1,0].plot(tpad, sigDown, linewidth=.3, color='black')
axes[1,0].plot(tpad, filtSig, linewidth=.3, color='black')
axes[1,0].set_xlim(5.0775,5.0825)
axes[1,0].set_xticklabels([])
axes[1,0].tick_params(bottom="off")

axes[2,0].set_title('filtfilt()')
axes[2,0].plot(tpad, sigDown, linewidth=.3, color='black')
axes[2,0].plot(tpad, filtfiltSig, linewidth=.3, color='black')
axes[2,0].set_xlim(5.0775,5.0825)
#axes[2,0].set_xticklabels([])
#axes[2,0].tick_params(bottom="off")


N = 64*4096
hN = int(N/2)
freqs = np.fft.fftfreq(N, samSpc)

SigFFT = fft(sigDown[:N])
axes[0,1].set_title('Frequency Domain')
axes[0,1].plot(freqs[:hN], np.abs(SigFFT[:hN]), linewidth=.3, color='black')
axes[0,1].set_yticklabels([])
axes[0,1].tick_params(bottom="off")

filtSigFFT = fft(filtSig[:N])
axes[1,1].plot(freqs[:hN], np.abs(filtSigFFT)[:hN], linewidth=.1, color='black')
axes[1,1].set_yticklabels([])
axes[1,1].tick_params(left="off")

filtfiltSigFFT = fft(filtfiltSig[:N])
axes[2,1].plot(freqs[:hN], np.abs(filtfiltSigFFT)[:hN], linewidth=.1, color='black')
axes[2,1].set_yticklabels([])
axes[2,1].tick_params(left="off")
'''
'''
#plt.plot(tpad, FDMSig[:,0], color='black', label='Ideal Signal')  
#plt.plot(tpad, FDMSigNoise[:,0], color='black', linestyle='dashed', label='Noisy Signal')  
#plt.xlabel('Time')
#plt.ylabel('Signal Amplitude')
#plt.legend()
#plt.scatter(tpad, FDMSig[:,0], s=1, color='red')


#plt.plot(tpad, FDMSig[:,0])
#plt.plot(tpad, filtSig)
#plt.xlim(5.5,5.56)

'''
'''
#fig, axes = plt.subplots(2, 1, sharex=False, sharey=False, figsize=(17,12), dpi=166) 
#axes[0].plot(tpad, FDMSig[:,0], color='black')
#axes[1].scatter(tpad[::10], FDMSig[::10,0], s=1, color='black')
#axes[0].set_xlim(0,.025)
#axes[1].set_xlim(0,.025)












