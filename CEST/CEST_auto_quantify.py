# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 18:56:51 2018

@author: Kristof
"""

import os
from glob import glob
import nibabel as nib
import pandas as pd
import numpy as np
import pylab as pl
from skimage import measure
from scipy import ndimage
from sklearn import cluster
import Tkinter, tkFileDialog
root = Tkinter.Tk()

BMASK_EXT = '_Bmask.nii'
MTR_EXT = '_MTRasym.nii'
FREQ_EXT = '_interpolated_freqs.csv'
F1 = 280 #frequency at which crypto contrast is the best
F2 = 1800 #frequency at which brain contrast is best
INVERT_AXIS = (False, 0) #if mask were drawn in ITK-Snap - need to invert axis. Either False (no flip), 0 (x-axis) or 1 (y-axis) 
SLICE = 0
BLUR = 1.0
    
def get_outline(binmask):
    outline = (binmask - ndimage.morphology.binary_erosion(binmask)).astype('float')
    outline[outline==0] = np.nan
    return outline

print "Enter the scan directory."

while True:
    rootdir = tkFileDialog.askdirectory(initialdir="/",title='Please select a directory')
    if os.path.isdir(rootdir) is True: #Checks if entered dir exists
        os.chdir(rootdir)
        root.destroy()
        break
    else:
        print "Pathname invalid. Try again."
        continue

os.chdir(rootdir)

filelist = glob('*' + MTR_EXT)
masklist = [f.replace(MTR_EXT, BMASK_EXT) for f in filelist]
freqlist = [f.replace(MTR_EXT, FREQ_EXT) for f in filelist]

hdr = 'name;brainvol;posvox;negvox;0p1vox;m0p1vox;ratio_pos;ratio_0p1\n'

for i, f in enumerate(filelist):
    im = nib.load(f)
    aff = im.get_affine()
    mtr = im.get_data()[:,:,SLICE,:]
    mask = nib.load(masklist[i]).get_data()[:,:,SLICE]
    if INVERT_AXIS[0] == True:
        print "Flipping mask."
        mask = np.flip(mask, INVERT_AXIS[1])
    m_outl = get_outline(mask)
    freqs = np.array(pd.read_csv(freqlist[i], sep=';')['pos'])#extract positive frequencies
    
    i1 = np.where(freqs==F1)[0][0]
    i2 = np.where(freqs==F2)[0][0]
    
    ratio_im = ndimage.filters.gaussian_filter(mask * mtr[:,:,i1], sigma=1.0) / ndimage.gaussian_filter(mask * mtr[:,:,i2], sigma=1)
    ratio_im[np.isnan(ratio_im)] = 0.0
    logratio_flat = np.log10(ratio_im[mask != 0.0]) #remove background and make log
    
    out_dict = {'brain_vox' : np.count_nonzero(mask),
                'log_neg': (logratio_flat < 0.0).sum(),
                'log_pos': (logratio_flat > 0.0).sum(),
                'log_0p1': (logratio_flat >= 0.1).sum(),
                'log_m0p1': (logratio_flat <= -0.1).sum()}
    
    out_dict['ratio_pos'] = np.round(np.float(out_dict['log_pos']) / out_dict['brain_vox'], 3)
    out_dict['ratio_0p1'] = np.round(np.float(out_dict['log_0p1']) / out_dict['brain_vox'], 3)
    
    out_txt = '''Brain voxels: %s
    Positive voxels: %s
    Negative voxels: %s
    Voxels > 0.1: %s
    Voxels < 0.1: %s
    Ratio Pos/Brain: %s
    Ratio +0.1/Brain: %s''' %(out_dict['brain_vox'], out_dict['log_pos'], 
    out_dict['log_neg'], out_dict['log_0p1'], out_dict['log_m0p1'], out_dict['ratio_pos'], out_dict['ratio_0p1'])
    
    #write to output txt file
    hdr += f + ';' + str(out_dict['brain_vox']) + ';' + str(out_dict['log_pos']) + ';' + str(out_dict['log_neg']) + ';' + str(out_dict['log_0p1']) + ';' + str(out_dict['log_m0p1']) + ';' + str(out_dict['ratio_pos']) + ';' + str(out_dict['ratio_0p1']) + '\n'
    
    fig = pl.figure()
    pl.subplot(2,2,1)
    pl.title('Frequency 1')
    pl.imshow(mtr[:,:,i1], vmin=0, vmax=15)
    pl.colorbar(label='MTRasym')
    pl.imshow(m_outl, vmin=0, vmax=1)
    
    pl.subplot(2,2,2)
    pl.title('Frequency 2')
    pl.imshow(mtr[:,:,i2], vmin=0, vmax=15)
    pl.colorbar(label='MTRasym')
    pl.imshow(m_outl, vmin=0, vmax=1)
    
    pl.subplot(2,2,3)
    pl.title('Log ratio')
    pl.imshow(np.log10(ratio_im), vmin=-1, vmax=1, cmap='jet')
    pl.colorbar(label='LogRatio F1/F2')
    
    pl.subplot(2,2,4)
    pl.title('Histogram')
    pl.hist(logratio_flat, bins=100)
    pl.axvline(x=0.0, color='black', ls='dotted')
    
    fig.savefig(f.replace(MTR_EXT, 'Quant_fig.png'))
   
with open('CEST_quant.csv', 'w') as f:
    f.write(hdr)