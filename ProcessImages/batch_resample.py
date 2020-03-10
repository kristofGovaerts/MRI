# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 08:28:24 2015

@author: Kristof Govaerts


"""

import os
import nibabel as nib
import numpy as np
import glob
import Tkinter, tkFileDialog
import subprocess
import shutil

#Globals
LABEL_IM = 'BNLref_DSc_lblR'
REFIM = 'BNLref_DSc.hdr'#has to be in same folder
REG_NRR_EXT= '_BFC_BMd1_TregS_BNLref_DSc_NR_invT.nii' #file extension for the nonrigid transform to use
REG_AFF_EXT= None #file extension for the affine transform to use. Can leave blank
FLOAT_EXT = '_BFC_BMd1_TregS_BNLref_DSc_NR_LblR_BNLref_lbl_DSct_ctxR.hdr' #images in floating space to resample

NEUROMORPH_LOC = "C:/Program Files/NeuroMorph/Packages/MIC/Applications/Modules/Macros/KULNeuroMorph"
#default pars
THRESH=7000
FTYPE='float'

reg_ala = os.path.join(NEUROMORPH_LOC, 'reg_aladin.exe')
reg_f3d = os.path.join(NEUROMORPH_LOC, 'reg_f3d.exe')
reg_transform = os.path.join(NEUROMORPH_LOC, 'reg_transform.exe')
reg_resample = os.path.join(NEUROMORPH_LOC, 'reg_resample.exe')

#functions
def RegAndResample(flt, atl, mode='A'):
    '''Registers [atl] into [flt] space and resamples [flt] into [atl] space.
    Inputs should be image paths.
    Can do affine (mode='A') or rigid (mode='R') registration.'''
    omat = os.path.join(os.getcwd(), flt.split('.')[0] + mode + '.txt')     
    out = os.path.join(os.getcwd(), flt.split('.')[0] + mode + '.nii') 
    if mode == 'A': #affine    
        args = '-flo %s -ref %s -noSym' %(os.path.join(os.getcwd(), flt), atl)         
    elif mode == 'R': #rigid
        args = '-flo %s -ref %s -noSym -rigOnly' %(os.path.join(os.getcwd(), flt), atl)   
        
    flirt = [reg_ala] + args.split()    
 
    process = subprocess.call(flirt, shell=True)
    os.rename('outputResult.nii', out) #NiftyReg saves as outputResult by default
    os.rename('outputAffine.txt', omat)

    return out
    
def RegAndResample_NRR(flt, atl, ext=''):
    '''Registers [atl] into [flt] space and resamples [flt] into [atl] space.
    Inputs should be image paths.
    Does nonrigid registration.'''
 
    out = os.path.join(os.getcwd(), flt.split('.')[0] + 'NRR' + str(ext) +  '.nii') 
    dfout = os.path.join(os.getcwd(), flt.split('.')[0] + 'CPP' + str(ext) + '.nii') 
    
    args = '-flo %s -ref %s' %(os.path.join(os.getcwd(), flt), atl)     
    fnirt = [reg_f3d] + args.split()  
    
    subprocess.call(fnirt, shell=True)
    os.rename('outputResult.nii', out) #NiftyReg saves as outputResult by default
    os.rename('outputCPP.nii', dfout)
    return out, dfout
    
def otsu_filter(img, background=0):
    flat = (img.flatten())[img.flatten()>background]
    thresh = filter.threshold_otsu(flat)
    im2 = np.copy(img)
    im2[im2<=thresh]=0
    return im2
    
def reg_batch(refim, filelist, mode, zoom=1, ext=''):
    '''Performs affine or rigid registration of 'filelist' to 'refim' and 
    resamples into 'refim' space. 'mode' should be 'rigid' or 'affine'.'''
    outlist = []
    dfoutlist = []
    for i, f in enumerate(filelist):
        print "file %s of %s." %(i+1, len(filelist))

        if mode == 'NRR':
            out, dfout = RegAndResample_NRR(f, refim, ext)
            dfoutlist.append(dfout)
        else:
            out = RegAndResample(f, refim, mode=mode) #one registration step
        outlist.append(out)
    if mode == 'NRR':
        o=(outlist, dfoutlist)
    else:
        o=outlist
    return o

def combine_transforms(t1, t2, ref, filename, ref2 = None):
    '''Inputs:
    t1, t2: Filenames of transforms to combine (Should point to control point imgs)
    ref: Reference image used for one of the transforms
    filename: Filename for the output control point image
    Output:
    Filename of the combined transform control point image'''
    if t1[-4:] == '.txt':
        args = '-comp %s %s %s -ref %s -ref2 %s' %(t1, t2, filename, ref, ref2)  
    else:
        args = '-comp %s %s %s -ref %s' %(t1, t2, filename, ref)     
    tf = [reg_transform] + args.split()  
    subprocess.call(tf, shell=True)
    return filename, tf
    
def apply_transform(t, ref, flo, filename):
    '''Inputs:
    t: Filename of a control point image or affine transformation matrix
    ref: Reference image
    flo: Floating image to which to apply the transform
    filename: Name to give the transformed floating image.
    Output:
    Filename of the transformed floating image.'''
    args = '-ref %s -flo %s -trans %s -res %s -inter 0' %(ref, flo, t, filename)     
    tf = [reg_resample] + args.split()  
    subprocess.call(tf, shell=True)
    return filename    
    
def save_sum_img(filelist, filename, thresh=7000, ftype='int32'):
    cim = nib.load(filelist[0]).get_data()
    sum_img = np.zeros((cim.shape[0], cim.shape[1], cim.shape[2], len(filelist)))
    
    for i, f in enumerate(filelist):
        cim = np.squeeze(nib.load(f).get_data())
        sum_img[:,:,:,i] = cim  
       
    nif1 = nib.Nifti1Image(sum_img, nib.load(filelist[0]).get_affine())
    nib.save(nif1, 'test.nii')       
     
    adata = (np.mean(sum_img, axis=-1)).astype(ftype) 
    adata[adata<thresh] = 0 #for intensity normalized images, the background is usually not 0. 
    nif = nib.Nifti1Image(adata, nib.load(filelist[0]).get_affine())
    nib.save(nif, filename)

#Main loop
print "Enter the scan directory. This directory should contain all scans except for the reference image, and no other scans."
root = Tkinter.Tk()
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

#rootdir = '/media/sf_host/data/NRR/Raw_data/test_FSL'
#os.chdir(rootdir)
#refim = '/media/sf_host/data/NRR/Raw_data/Bl6_Ref_ER.hdr'
#maskim = '/media/sf_host/data/NRR/Raw_data/Bl6_ref_ER_BNL_lbl_BM.hdr'

filelist=glob.glob('*' + REG_NRR_EXT)    

print "STEP 1: Applying affine + nonrigid transform to all volumes of input image."
for i, f in enumerate(filelist):
    im = f.replace(REG_NRR_EXT, FLOAT_EXT)
    
    if REG_AFF_EXT != None:
        aff=f.replace(REG_NRR_EXT, REG_AFF_EXT)
        tc, tf = combine_transforms(f, aff, REFIM, f[:-4] + 'CPPA.nii')
    else:
        tc = f
    print "Image %s of %s." %(i+1, len(filelist))
    print "Transform:", tc
    print "Reference:", REFIM
    print "Float:", im
    output = f[:-4] + "_InvLblR.nii"
    print "Output:", output
    nrr=apply_transform(tc, REFIM, im, output) #NRR file has affine and nrr transform in it
    
