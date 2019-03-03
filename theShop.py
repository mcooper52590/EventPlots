#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  2 10:44:53 2019

@author: mattcooper
"""
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#This is the Shop Floor.  Place codes here to work on/test them.  Place function calls/libraries 
#in the ShopFunctions.py file in this same directory.


#Need to reconfigure to take a dict instead
#Also, need to be able to specify field-aligned coordinate system, as well as which component to
#get the FFT information for.
#def get_FFT_Power_Info(mag_Mag, mag_Epoch): 

T = 1
N = 2048
window = blackmanharris(N)



        


 
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

