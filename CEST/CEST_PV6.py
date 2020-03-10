# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 14:54:00 2017

CEST processing batch script. Compatible with PV6 data. 

@author: Kristof Govaerts
"""

import os
import numpy as np 
import nibabel as nib #medical image I/O - http://nipy.org/nibabel/
from scipy.interpolate import PchipInterpolator
import re
from scipy import ndimage as ndi
from Tkinter import Tk
from tkFileDialog import askopenfilename
import pandas as pd
import imreg_dft as ird #image registration - for motion correction. https://pypi.org/project/imreg_dft/
from pandas import rolling_mean

#Globals
FUNC_EXT='7.nii'

B0_BIN = 1 #bin size for the B0 correction. 
B0_RANGE = [10,68]#B0 range of the ordered spectrum. For a symmetric spectrum these should then be the middle freqs
B0_NFREQS = 40 #amount of freqs around center freq that can be used for B0 correction
CEST_BIN = 20 #bin size for the actual CEST spectrum
THRESH = 40.0 #signal threshold (at largest frequency offset)
THRESH_MC = 20.0 #for motion correction
CTRL = [20000, -20000] #control frequencies
MC = False
AUC_RANGE = [80, 800] #range for calculating AUC
T_SMOOTH = 3
SMOOTH = 1

#This is a dictionary with known data types. Important for reading the file correctly
dtypes = {
'_16BIT_SGN_INT': np.int16,
'_32BIT_SGN_INT': np.int32
}

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def temporal_smoothing(data, t_smooth):
    out_arr = np.array(data)
    for x in range(data.shape[0]):
        for y in range(data.shape[1]):
            for z in range(data.shape[2]):
                out_arr[x,y,z,:] = rolling_mean(out_arr[x,y,z,:], t_smooth, min_periods=1)   
    return out_arr

def load_bruker(f):
    '''Loads a Bruker 2dseq file and reformats based on associated visu_pars file.
    Outputs an image array.'''
    
    #Step 1: Extract dimensionality from visu_pars file
    vp = os.path.join(f, 'pdata/1/visu_pars') #load visu_pars
    with open(vp) as vpfile:
        l = vpfile.readlines()
    for i, line in enumerate(l):        
        if '##$VisuCoreSize' in line: #get matrix dimensions
            dims = l[i+1][:-1].split(' ')
            dims = [int(x) for x in dims]
        elif '##$VisuCoreFrameCount' in line: #get frame count
            fc = int(l[i][l[i].find('=') + 1 : -1])
        elif '##$VisuCoreWordType=' in line: #image bits
            dtype = l[i][l[i].find('=') + 1 : -1]
        elif '##$VisuCoreDataSlope' in line: #slopes
            temp_string = str(l[i+1:])
            slopes = temp_string[:temp_string.find('##') - 3]
            slopes = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", slopes)   #Regular expression that extracts all numbers from the string...Not sure how it works!   
            slopes = np.array([float(s) for s in slopes])
            
    if len(dims) < 3:
        dims += [1]
    size = dims + [fc] #create image shape list

    #Step 2: Load 2dseq file and reformat
    im_arr = np.fromfile(os.path.join(f, 'pdata/1/2dseq'), dtype=dtypes[dtype])
    im_arr = im_arr.reshape(size, order='F').astype('float') #Bruker files are ordered Fortran-style ('F')
    im_arr *= slopes #slope correction
    
    return im_arr
    
def get_parameters(f):
    method = os.path.join(f, 'method')
    par_dict = {}
    with open(method) as metfile:
        l = metfile.readlines()
    for i, line in enumerate(l):        
        if '##$PVM_SpatResol' in line: #get matrix dimensions
            resol = l[i+1][:-1].split(' ')
            resol = [float(x) for x in resol]
            par_dict['resolution'] = resol
        if '##$PVM_SliceThick=' in line: #get matrix dimensions
            s_thick = float(l[i][l[i].find('=')+1:-1])
            par_dict['slicethickness'] = s_thick
        if '##$CEST_FrequencyOffset' in line:
            temp_string = str(l[i+1:])
            freqs = temp_string[:temp_string.find('$$') - 3]
            freqs = re.findall("[-+]?\d+[\.]?\d*[eE]?[-+]?\d*", freqs)   #Regular expression that extracts all numbers from the string...Not sure how it works!   
            freqs = np.array([float(s) for s in freqs])
            par_dict['CEST_freqs'] = freqs
    return par_dict
    
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
    ip=PchipInterpolator(x_freq, y_i) #Piecewise spline interpolation
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
            ip=PchipInterpolator(freqs, cim_f[vx, :])
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

def get_edges(arr2d):
    eh = ndi.sobel(arr2d,0)
    ev = ndi.sobel(arr2d, 1)
    return np.hypot(eh, ev)

def brain_mask(a4d):
    im = ndi.gaussian_filter(a4d, (2,2,0,0))
    im[im<THRESH_MC/2.0] = 0.0
    return im

def show_slices(a4d, axis=-1):
    nim = a4d.shape[axis]
    pl.figure()
    for i in range(nim):
        pl.subplot(int(np.sqrt(nim)+1), int(np.sqrt(nim)+1), i+1)
        pl.imshow(a4d[:,:,0,i])

def motion_correction(arr, buffers = 7):
    '''Simple, quick 2D translation for motion correction.'''
    a4d = np.copy(arr)
    a4d[a4d<THRESH_MC] = 0.0
    shape = a4d.shape
    outmat = np.zeros(shape)
    
    mean_intensities = [np.mean(a4d[:,:,:,i]) for i in range(shape[3])]
    m = np.argmin(mean_intensities) #volume where the average intensity is the lowest
    refvol = np.sum(a4d[:,:,:,m-buffers:m+buffers], axis=3)
    refvol_rs = ndi.interpolation.zoom(refvol, (3.0,3.0,1))    
    #refvol_rs = ndi.filters.gaussian_filter(refvol_rs, (1,1,0))
    for i in range(shape[2]):
        refvol_rs[:,:,i]=get_edges(refvol_rs[:,:,i])
    refvol_mean = np.mean(refvol_rs)
    pl.figure()
    pl.imshow(refvol_rs[:,:,0])
    
    for sl in range(shape[2]):
        for v in range(shape[3]):
            print "Motion correction in progress. Volume %s of %s." %(v+1, shape[3])
            if v > m-buffers and v < m+buffers:
                print "Reference slice. Skipping."
                outmat[:,:,sl,v]=a4d[:,:,sl,v]
            else:
                csl =  a4d[:,:,sl,v]
                csl_rs = ndi.interpolation.zoom(csl, (3.0,3.0))
                #csl_rs = ndi.gaussian_filter(csl_rs, sigma=1.0)
                csl_rs = get_edges(csl_rs)
                csl_rs = (csl_rs / np.mean(csl_rs)) * refvol_mean #intensity normalization
                pl.figure()
                pl.imshow(csl_rs)
                tl = ird.imreg.similarity(refvol_rs[:,:,sl], csl_rs, constraints={"angle":[0.0,45.0], "tx":[0.0,20.0], "ty":[0.0,20.0], "scale":[1.0, 0.0]})
                tvec = tl["tvec"].round(4)
                angle = tl["angle"].round(4)
                print tvec
                print angle
                print tl["scale"]
                outmat[:,:,sl,v] = ird.transform_img(csl, angle=angle, tvec=tvec/3.0)
    return outmat
        
Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file

parfile = pd.read_csv(filename, delimiter=';')
sourcefiles = parfile.path
outnames = parfile.out_name
        
for i, im in enumerate(sourcefiles):
    print 'File %s of %s.\n%s' %(i+1, len(sourcefiles), im)
    try:
        arr = load_bruker(im)
    except IOError:
        print "visu_pars file not found. Scan probably not acquired."
        continue
    pars = get_parameters(im)
    if arr.shape[-1] != len(pars['CEST_freqs']):
        arr = arr.reshape((arr.shape[0], arr.shape[1], arr.shape[-1]/len(pars['CEST_freqs']), len(pars['CEST_freqs'])), order='F')
        
    affine = np.diag(pars['resolution'] + [pars['slicethickness'], 1.0])
    arr_o = order_im(arr, pars['CEST_freqs'])
    f_o = sorted(pars['CEST_freqs'])    
    
    if MC == True: #motion corr
        arr_o = motion_correction(arr_o)
        print "Saving motion corrected image..."
        mcnif = nib.Nifti1Image(arr_o, affine)
        nib.save(mcnif, os.path.join(outnames[i] + '_MC.nii'))
        
    try:
        ctrl_i = [f_o.index(c) for c in CTRL]
    except ValueError:
        print 'Control frequency not found. Skipping.'
        continue
    ctrl_arr = arr_o[:,:,:,ctrl_i]
    
    arr_o = np.delete(arr_o, ctrl_i, axis=-1)
    f_o = list(np.delete(f_o, ctrl_i))
    
    print "Saving sorted image..."
    snif = nib.Nifti1Image(arr_o, affine)
    nib.save(snif, os.path.join(outnames[i] + '_Sorted.nii'))
  
    print "Calculating B0 map..."
    if B0_RANGE != None:
        B0_RANGE = [len(f_o)/2 - B0_NFREQS/2, len(f_o)/2 + B0_NFREQS/2]
    b0 = create_b0map(arr_o[:,:,:,B0_RANGE[0]:B0_RANGE[1]], f_o[B0_RANGE[0]:B0_RANGE[1]], bin = B0_BIN, thresh=THRESH)    
    
    print "Saving B0 map..."
    b0nif = nib.Nifti1Image(b0.astype('float32'), affine)
    nib.save(b0nif, os.path.join(outnames[i] + '_B0.nii'))
    
    print "Performing B0 correction..."
    freqs = range(0, int(f_o[-1])-90, CEST_BIN) #interpolated freqs
    plusmap, minmap = correct_b0(arr_o, f_o, b0, freqs, thresh=THRESH)
    
    print "Saving B0 corrected Z-spectrum..."
    plusnif = nib.Nifti1Image(plusmap.astype('float32'), affine)
    nib.save(plusnif, os.path.join(outnames[i] + '_Pos_corr.nii'))
    negnif = nib.Nifti1Image(minmap.astype('float32'), affine)
    nib.save(negnif, os.path.join(outnames[i] + '_Neg_corr.nii'))   
    fdf = pd.DataFrame(freqs)
    fdf['neg'] = -np.array(freqs)
    fdf.columns = ['pos', 'neg']
    fdf.to_csv(outnames[i] + '_interpolated_freqs.csv', sep=';')
    
    print "Saving MTRasym map..."
    MTRasym = minmap-plusmap #no longer absolute value
    if SMOOTH != 0: 
        MTRasym = ndi.filters.gaussian_filter(MTRasym, [SMOOTH, SMOOTH, 0,0])
    ctrl = ctrl_arr[:,:,:,0]
    for j in range(MTRasym.shape[-1]):
        MTRasym[:,:,0,j] = 100*(MTRasym[:,:,0,j]/ctrl[:,:,0])
    mtrnif = nib.Nifti1Image(MTRasym.astype('float32'), affine)
    nib.save(mtrnif, os.path.join(outnames[i] + '_MTRasym.nii'))
    
    print "Saving AUC map..."
    idx = [find_nearest(np.array(freqs), AUC_RANGE[0]), find_nearest(np.array(freqs), AUC_RANGE[1])]
    if T_SMOOTH != 0: 
        auc_sect = temporal_smoothing(MTRasym, T_SMOOTH)[:,:,:,idx[0]:idx[1]]
    mauc_sect = -auc_sect    
    auc_sect[auc_sect<=0.0] = 0.0
    mauc_sect[mauc_sect<=0.0] = 0.0
    auc_arr = np.zeros(auc_sect.shape[:-1])
    for x in range(auc_arr.shape[0]):
        for y in range(auc_arr.shape[1]):
            for z in range(auc_arr.shape[2]):
                if np.sum(auc_sect[x,y,z,:]) != 0.0:
                    auc_arr[x,y,z] = np.trapz(auc_sect[x,y,z,:], dx=CEST_BIN) - np.trapz(mauc_sect[x,y,z,:], dx=CEST_BIN)#trapezoid AUC calculation
    aucnif = nib.Nifti1Image(auc_arr.astype('float32'), affine)
    fn = os.path.join(outnames[i] + '_AUC%s_%s.nii')  %(str(AUC_RANGE[0]), str(AUC_RANGE[1]))
    nib.save(aucnif, fn)
    
    print "Saving Ctrl image..."
    ctrlnif = nib.Nifti1Image(ctrl.astype('float32'), affine)
    nib.save(ctrlnif, os.path.join(outnames[i] + '_Ctrl.nii'))

       
