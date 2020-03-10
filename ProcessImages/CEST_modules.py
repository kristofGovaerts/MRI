# -*- coding: utf-8 -*-
"""
Created on Mon May 26 14:09:18 2014

@author: u0091609
"""

from ReadInput import *
import numpy as np
import pylab as pl
import os
import glob
from scipy.interpolate import interp1d

def order_im(img, freqs):
    '''Orders image in the t-domain according to frequencies, in
    ascending fashion. Necessary for scipy's interpolation functions.'''
    r=len(freqs) 
    orders=sorted(zip(freqs,range(r)))
    simg=np.zeros(img.shape)
    for i in range(r):
        simg[:,:,:,i]=img[:,:,:,orders[i][1]]
    return simg
    
def find_offset(x_freq, y_i, bin=0.1, return_xy=False):
    '''Finds the frequency offset from zero when provided a list of frequencies
    and a list of image intensities. In principle simply finds the minimum
    of a curve found via interpolation between points.
    
    x_freq = Frequencies
    y_i = Intensities
    bin = step size for the interpolation function. The smaller, the more 
    accurate but the more RAM you need.
    return_xy: if True, returns interpolated frequencies and intensities.'''
    ip=interp1d(x_freq, y_i, kind="cubic")
    xn=np.arange(x_freq[0], x_freq[-1], bin)
    yn=ip(xn)
    if return_xy==False:
        output=xn[np.argmin(yn)]
    else:
        output=(xn[np.argmin(yn)],xn,yn)
    return output

def create_b0map(img, freqs, bin=1, thresh=1e9):
    '''Creates a b0map based on the frequency offsets calculated per voxel.'''
    cim_f=np.reshape(img, (np.prod(img.shape[:-1]), img.shape[-1]))    
    arr=np.zeros(len(cim_f))
    for vx in range(len(cim_f)):
        if cim_f[vx,-1]>thresh:
            arr[vx]=find_offset(freqs, cim_f[vx,:], bin=bin)
    b0map=np.reshape(arr, img.shape[:-1])
    return b0map
    
def find_nearest(array,value):
    '''Finds index of nearest value in an array.'''
    idx = (np.abs(array-value)).argmin()
    return idx
    
def correct_b0(img, freqs, b0, offsets, bin=10, thresh=1e5):
    '''Corrects for B0 effects and returns desired offset images.
    
    img = 4D array
    freqs = frequencies corresponding in order to the 4th dim of the img
    offset = desired CEST offset (in a list, allows multiple)
    thresh = Masking threshold
    
    Returns:
    b0map - b0 offset map
    plusmap - +offset map corrected for b0
    minmap - -offset map corrected for b0'''
    cim_f=np.reshape(img, (np.prod(img.shape[:-1]), img.shape[-1]))    
    b0map=np.reshape(b0, (np.prod(img.shape[:-1])))    
    plusmap=np.zeros((len(cim_f),len(offsets)))
    minmap=np.zeros((len(cim_f),len(offsets)))
    
    for vx in range(len(cim_f)):
        if cim_f[vx,-1]>thresh:
            ip=interp1d(freqs, cim_f[vx, :], kind="cubic")
            xn=np.arange(freqs[0], freqs[-1], bin)
            yn=ip(xn)
            off = b0map[vx]
            for i,o in enumerate(offsets):
                plus=yn[find_nearest(xn, o + off)]
                minus=yn[find_nearest(xn, -o + off)]
                plusmap[vx,i]=plus
                minmap[vx,i]=minus
           
    nshape=list(img.shape[:-1])        
    nshape.append(len(offsets))
    plusmap=np.reshape(plusmap, nshape)
    minmap=np.reshape(minmap, nshape)
    return plusmap, minmap 
    
    
def correct_b02(img, freqs, offsets, thresh=1e9):
    '''OLD, ONLY USE IF NECESSARY
    
    Corrects for B0 effects and returns desired offset images.
    
    img = 4D array
    freqs = frequencies corresponding in order to the 4th dim of the img
    offset = desired CEST offset (in a list, allows multiple)
    thresh = Masking threshold
    
    Returns:
    b0map - b0 offset map
    plusmap - +offset map corrected for b0
    minmap - -offset map corrected for b0'''
    cim_f=np.reshape(img, (np.prod(img.shape[:-1]), img.shape[-1]))    
    b0map=np.zeros(len(cim_f))
    plusmap=np.zeros((len(cim_f),len(offsets)))
    minmap=np.zeros((len(cim_f),len(offsets)))
    
    for vx in range(len(cim_f)):
        if cim_f[vx,-1]>thresh:
            off, xn, yn=find_offset(freqs, cim_f[vx,:], return_xy=True)
            for i,o in enumerate(offsets):
                plus=yn[find_nearest(xn, o + off)]
                minus=yn[find_nearest(xn, -o + off)]
                plusmap[vx,i]=plus
                minmap[vx,i]=minus
            b0map[vx]=off
            
    nshape=list(img.shape[:-1])        
    b0map=np.reshape(b0map, nshape)
    nshape.append(len(offsets))
    plusmap=np.reshape(plusmap, nshape)
    minmap=np.reshape(minmap, nshape)
    return b0map, plusmap, minmap 
    
def cest_spectrum(img, freqs, n=150, step=10, thresh=1e9):
    cim_f=np.reshape(img, (np.prod(img.shape[:-1]), img.shape[-1])) 
    b0map=np.zeros(len(cim_f))
    cmap=np.zeros((len(cim_f), n+1))
    
    for vx in range(len(cim_f)):
        if cim_f[vx,-1]>thresh:
            off, xn, yn=find_offset(freqs, cim_f[vx,:], bin=step, return_xy=True)
            zero=find_nearest(xn, off)
            minus=yn[:zero+1][-n-1:][::-1] #have to reverse 
            plus=yn[zero+1:][:n+1]
            try:
                cest=(plus-minus)/plus
                cmap[vx,:]=cest
            except ValueError:
                cmap[vx,:]=-np.ones(n+1)
            b0map[vx]=off
    
    b0map=np.reshape(b0map, img.shape[:-1])
    cmap=np.reshape(cmap, list(img.shape[:-1])+[n+1])
    return b0map, cmap
    
def scansToArray(filelist, start=0, end=None, correct=True, sorted=False):
    '''Returns an array concatenating a series of 3D arrays across the
    t-dimension.

    Inputs:
        filelist = list of files, no extension
        start = where to start concatenating, optional (necessary if matrix size not equal for all images in your filelist)
        end=where to stop'''
    if sorted == False:
        sfl=sort_scans(filelist)
    else:
        sfl=filelist
    if end is None:
        scans=sfl[start:]
    else:
        scans=sfl[start:end]
    data1=nib.load(scans[0]).get_data()
    if correct==True:
        data1=correctSlope(data1, scans[0])
    array=np.zeros([data1.shape[0], data1.shape[1],1, len(scans)])
    array[:,:,:,0]=data1[:,:,:,0]
    for i, scan in enumerate(scans[1:]):
        data=nib.load(scan).get_data()
        if correct==True:
            data=correctSlope(data, scan, pr=False)
        array[:,:,:,i+1]=data[:,:,:,0]
    return array