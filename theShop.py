#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 10:44:53 2019

@author: mattcooper
"""
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#This is the Shop Floor.  Place codes here to work on/test them.  Place function calls/libraries 
#in the ShopFunctions.py file in this same directory.


#Need to reconfigure to take a dict itead
#Also, need to be able to specify field-aligned coordinate system, as well as which component to
#get the FFT information for.
#def get_FFT_Power_Info(mag_Mag, mag_Epoch): 

T = 1
N = 256
window = blackmanharris(N)
fig, (ax1) = plt.subplots(1, 1, sharex=True, sharey=False, figsize=(17,12), dpi=166) 
ax1T = ax1.twinx()
ax2T = ax2.twinx()
ax3T = ax3.twinx()
fig.subplots_adjust(hspace=0.02)
ax3.set_xlabel('Time Stamp')
coordNumber = 1
ax = ax1
def get_Power_Spectra_Subplot(fig, ax, axT, magDict):
FAUV_ToMagDict(magDict)
get_field_Aligned_Mag(magDict)


takeFFT_EMFISISMag(T, N, window, magDict)
ax.text(.98, .8, 'Magnetic Field Power Spectra', horizontalalignment='right', transform=ax.transAxes, 
  bbox=dict(facecolor='white', alpha=0.7))
ax.set_facecolor('black')
lowest = 3
cutoff = 1200
#power, frequency, datesNum, span, N = get_FFT_Power_Info(magDict['Magnitude'], magDict['Epoch'])
#newN = len(frequency) - cutoff
#newFreq = frequency[lowest:newN]
y,x = np.meshgrid(magDict['Frequency'][3:], mdates.date2num(magDict['FFTEpoch']))
newPower = np.abs(np.log(np.abs(magDict['FFT_Raw'][coordNumber][:,3:])**2))
im = ax.pcolormesh(x, y, newPower, cmap='RdBu_r')
#ax.set_ylim(np.max(newFreq), np.min(newFreq))
#ax.set_xlim(datesNum[0], datesNum[span-1])
#Hours = mdates.HourLocator()   
#Minutes = mdates.MinuteLocator()
#ax.set_yscale('log')
#ax.yaxis.set_major_formatter(plt.FuncFormatter(format_func))
#ax.xaxis.set_major_locator(Hours)
#ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
#ax.set_xticks(datesNum[::int(np.ceil(len(magDict['Epoch'])/15))])
#fig.colorbar(im, ax=ax, fraction = .05, pad = .07)   
#ax.set_ylabel('Period [sec]', labelpad=5)



        


 
#Editing field aligned files for insertion on this line    
#if i > int(N/2) and len(mag_Mag) - i > int(N/2):
#magSlice = 
#if i < int(N/2):
#firstSlice = mag_Mag[0:int(N/2)-i]
#magSlice = np.append(firstSlice, mag_Mag[0:i+int(N/2)])
#if len(mag_Mag) - i < int(N/2):
#firstSlice = mag_Mag[i-int(N/2):]
#magSlice = np.append(firstSlice, mag_Mag[len(mag_Mag) + 
#(len(mag_Mag)-(int(N/2)+i)):])
#secondSlice = datesNum[i-int(N/2):]
#timeSlice = np.append(secondSlice, datesNum[len(mag_Mag) + 
#(len(mag_Mag)-(int(N/2)+i)):])
#  
#fitMag = get_Fit_Values(timeSlice, magSlice)


#==========================================================================================
#Functions that were in old code that didn't seem to have a purpose
#Saved them here in kind of a junk drawer in case their use arose.
#I didn't document shit back then, so I actually don't know what some of these do.
def get_field_Aligned_UnitVecs(MAGDICT):
    unitVec_Zeta = []
    unitVec_Radial = []
    unitVec_Azimuthal = [] 
    for i in range(len(MAGDICT['Epoch'])):  
        unitVec_Z = MAGDICT['Mag'][i,:]/np.linalg.norm(MAGDICT['Mag'][i,:])
        unitVec_Earth = MAGDICT['coordinates'][i,:]/np.linalg.norm(MAGDICT['coordinates'][i,:])
        vec_Azi = -np.cross(unitVec_Z, unitVec_Earth)
        unitVec_Azi = vec_Azi/np.linalg.norm(vec_Azi)
        vec_Rad = np.cross(unitVec_Z, unitVec_Azi) 
        unitVec_Zeta.append(unitVec_Z)
        unitVec_Azimuthal.append(unitVec_Azi)
        unitVec_Radial.append(vec_Rad/np.linalg.norm(vec_Rad))    
    unitVec_Zeta = np.array(unitVec_Zeta)
    unitVec_Radial = np.array(unitVec_Radial)
    unitVec_Azimuthal = np.array(unitVec_Azimuthal)
    return unitVec_Zeta, unitVec_Radial, unitVec_Azimuthal

