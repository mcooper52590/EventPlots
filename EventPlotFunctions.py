#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 10:44:53 2019

@author: mattcooper
"""
'''
Functions which are utilized in the creation of event plots for the Climatology paper submitted in Spring 2019.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import cdfstuff.ripper_V2 as rip
import os
'''
========================================================================================================================
'''
def freq_to_Period(value, tick_number):
    return '{0:.4g}'.format(1/value)

def fmt(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    return r'${}*10^{{{}}}$'.format(a, b)

def FAUV_ToMagDict(magDict):
    '''
    Takes a dict from the EMFISIS data and adds the Field Aligned Unit Vector (FAUV)
    '''
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

#
def get_field_Aligned_Mag(magDict):
    '''
    Takes a dict from the EMFISIS data and uses FAUV_ToMagDict, as well as adding the field values in the field-aligned 
    coordinate system.
    '''
    magDict['MagFA'] = np.zeros([len(magDict['Epoch']), 4])
    for i in range(len(magDict['Epoch'])):
        magDict['MagFA'][i,0] = np.dot(magDict['Mag'][i,:], magDict['FAUV'][i,0,:])
        magDict['MagFA'][i,1] = np.dot(magDict['Mag'][i,:], magDict['FAUV'][i,1,:])
        magDict['MagFA'][i,2] = np.dot(magDict['Mag'][i,:], magDict['FAUV'][i,2,:]) 
        magDict['MagFA'][i,3] = magDict['Magnitude'][i]

def takeFFT_EMFISISMag(T, N, window, magDict):
    '''
        T = Sample Period
        N = Number of samples in window
        window = Windowing function array for FFT.  Must be the same size as N
        magDict = a dict from an EMFISIS CDF file
    This calculates FFTs for the three field aligned coordinate directions, as well as the overall 
    magnitude.  It can generate a rather large dict if given a long enough span, so be careful.
    '''
    magDict['Frequency'] = np.fft.fftfreq(N)
    magDict['FFTEpoch'] = magDict['Epoch'][int(N/2):len(magDict['Epoch']) - int(N/2)]
    magDict['FFT_Raw'] = np.zeros([4, len(magDict['FFTEpoch']), magDict['Frequency'].shape[0]], dtype='complex64')
    get_field_Aligned_Mag(magDict)
    hWin = int(N/2)
    for i in range(hWin, len(magDict['Epoch']) - hWin):
        for j in range(0,4):
            magDict['FFT_Raw'][j, i - hWin] = np.fft.fft(window*magDict['MagFA'][i - hWin:i + hWin, j])

def get_Power_Spectra_Subplot(fig, ax, magDict, N, T, window, coord='FieldAligned'):
    '''
    Function which takes a matplotlib fig and ax, along with the magDict from either the database at NJIT for EMFISIS
    on Dgar or from a EMFISIS CDF file, and plots the power spectra for the dict on the axes provided.
        #coord = ('FieldAligned', 'Azimuthal', 'Radial')
        #N = number of samples to include in the FFT
        #T is the sampling period
    '''
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
    ax.yaxis.set_major_formatter(plt.FuncFormatter(freq_to_Period))
#    ax.xaxis.set_major_locator(Hours)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
#    ax.set_xticks(datesNum[::int(np.ceil(len(magDict['Epoch'])/15))])
    fig.colorbar(im, ax=ax, fraction = .05, pad = .07)   
  
def get_Beta_LinePlot(ax, betaDict):
    #Plots the beta values from a betaDict provided by Dgar
    ax.plot(betaDict['Epoch'], betaDict['Total'], color='black', label='Beta')

def get_Particle_Heatmap_SubPlot_EnergyBinned(fig, ax, tofDict):
    '''
    Function which takes a matplotlib fig and ax, along with the tofDict from either the database at NJIT or TOF cdf
    and generates a heatmap of the angle summed particles for different energies
    '''
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
    #Plots particle densities from a tofDict provided by Dgar or a TOF cdf
    ax.plot(tofDict['Epoch'], tofDict['FPDU_Density'], color='maroon')  

def get_Particle_Heatmap_AngleBinned(fig, ax, tofDict):
    '''
    Function which takes a matplotlib fig and ax, along with the tofDict from either the database at NJIT or TOF cdf
    and generates a heatmap of the energy summed particles for different angles
    '''
    angleBin = np.sum(tofDict['FPDU'], 1)
    ax.text(.98, .8, 'Particle Count By Pitch Angle', horizontalalignment='right',
            transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.7))
    ax.set_facecolor('black')
    angle = tofDict['PA_Midpoint'][0,:]
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
    
def get_Position_Plot_For_Span(ax1, ax2, strtDt, endDt):
    magDict = rip.get_CDF_Dict('Mag_1Sec_A', strtDt, endDt)
    coords = magDict['coordinates']/6371.5
    
    ws = 6
    ax1.set_xlim(ws,-ws)
    ax1.set_ylim(ws,-ws)
    img = plt.imread(os.getcwd() + "/Background-XY-GREY.jpeg")
    ax1.imshow(img,zorder=0, extent=[6, -6, 6, -6])
    ax1.scatter(coords[:,0], coords[:,1], marker='^', c='tomato', s=.1, label='RBSP-A')
    ax1.set_facecolor('.82')
    
    for i in range(len(coords[:,0])):
        if i % 10000 ==0 and i-len(coords[:,0]) < 200:
            diffX = coords[i+1, 0] - coords[i,0]
            diffY = coords[i+1, 1] - coords[i,1]
            ax1.arrow(coords[i,0], coords[i,1], diffX, diffY, shape='full', 
                      lw=0, length_includes_head=True, head_width=.2, 
                      color='darkred')
    
    ax2.set_xlim(ws,-ws)
    ax2.set_ylim(-ws,ws)
    img = plt.imread(os.getcwd() + "/Background-XZ-GREY.jpeg")
    ax2.imshow(img,zorder=0, extent=[6, -6, -6, 6])
    ax2.scatter(coords[:,0], coords[:,2], marker='^', c='tomato', s=.1, label='RBSP-A')
    ax2.set_facecolor('.82')
    
    for i in range(len(coords[:,0])):
        if i % 10000 ==0 and i-len(coords[:,0]) < 200:
            diffX = coords[i+1, 0] - coords[i,0]
            diffY = coords[i+1, 2] - coords[i,2]
            ax2.arrow(coords[i,0], coords[i,2], diffX, diffY, shape='full', 
                      lw=0, length_includes_head=True, head_width=.2, 
                      color='darkred')
    
    magDict = rip.get_CDF_Dict('Mag_1Sec_B', strtDt, endDt)
    coords = magDict['coordinates']/6371.5
    
    ws = 6
    ax1.set_xlim(ws,-ws)
    ax1.set_ylim(ws,-ws)
    ax1.scatter(coords[:,0], coords[:,1], marker='^', c='green', s=.1, label='RBSP-B')
    
    for i in range(len(coords[:,0])):
        if i % 10000 ==0 and i-len(coords[:,0]) < 200:
            diffX = coords[i+1, 0] - coords[i,0]
            diffY = coords[i+1, 1] - coords[i,1]
            ax1.arrow(coords[i,0], coords[i,1], diffX, diffY, shape='full', 
                      lw=0, length_includes_head=True, head_width=.2, 
                      color='darkgreen')
    
    ax2.set_xlim(ws,-ws)
    ax2.set_ylim(-ws,ws)
    ax2.scatter(coords[:,0], coords[:,2], marker='^', c='green', s=.1, label='RBSP-B')
    
    for i in range(len(coords[:,0])):
        if i % 10000 ==0 and i-len(coords[:,0]) < 200:
            diffX = coords[i+1, 0] - coords[i,0]
            diffY = coords[i+1, 2] - coords[i,2]
            ax2.arrow(coords[i,0], coords[i,2], diffX, diffY, shape='full', 
                      lw=0, length_includes_head=True, head_width=.2, 
                      color='darkgreen')
    
    ax1.minorticks_on()
    ax1.legend(markerscale=20)
    ax1.grid(True, which='major', linewidth=.2, color='black')      
    ax1.grid(True, which='minor', linewidth=.1, color='black')  
    ax1.set_xlabel('X-Axis SM Coordinates')
    ax1.set_ylabel('Y-Axis SM Coordinates')
    
    ax2.legend(markerscale=20)
    ax2.minorticks_on()
    ax2.grid(True, which='major', linewidth=.2, color='black')      
    ax2.grid(True, which='minor', linewidth=.1, color='black')  
    ax2.set_xlabel('X-Axis SM Coordinates')
    ax2.set_ylabel('Z-Axis SM Coordinates')
    
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
    A simple high pass filter function.
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

def mag_AverageRemoved(N, magDict):
    hw = int(N/2)
    magDict['Epoch_Avg'] = magDict['Epoch'][hw:magDict['Epoch'].shape[0] - hw] 
    magDict['MagFA_Avg'] = np.zeros([magDict['Epoch_Avg'].shape[0], 4])
    magDict['MagFA_NoBack'] = np.zeros([magDict['Epoch_Avg'].shape[0], 4])
    for i in range(0, magDict['Epoch'].shape[0] - N):
        for j in range(0,3):
            magDict['MagFA_Avg'][i, j] = np.mean(magDict['MagFA'][i:i + N, j]) 
            magDict['MagFA_NoBack'][i, j] = magDict['MagFA'][i,j] - magDict['MagFA_Avg'][i,j]