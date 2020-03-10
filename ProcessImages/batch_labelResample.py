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
REFIM = 'C:\Users\\u0091609\Desktop\\forReg\BNLref_DSc.hdr'
REG_NRR_EXT= '_BFC_BMd1_TregS_BNLref_DSc_NR_invT.nii' #file extension for the nonrigid transform to use
REG_AFF_EXT= None #file extension for the affine transform to use. Can leave blank
FLOAT_EXT = '_BFC_BMd1_TregS_BNLref_DSc_NR_LblR_BNLref_lbl_DSct_ctx.hdr' #images in floating space to resample

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

def batch_combine_apply_tf(ref, flo_list, tf1_list, tf2_list, ext='CPPC.nii'):
    '''Combines transforms from lists 1 and 2 and applies them to images in flo_list. 
    All input lists should be equally long. ref is the reference image used to 
    generate the transforms.'''
    ext2=ext.split('.')[0] + 'T.nii'
    tflist = []
    tfimlist = []
    assert len(flo_list) == len(tf1_list) and len(flo_list) == len(tf2_list) #inputs should be the same size!    
    
    for i, f in enumerate(flo_list):
        print "Combining transforms for image %s of %s..." %(i+1, len(flo_list))
        fn=combine_transforms(tf1_list[i], tf2_list[i], ref, f.split('.')[0] + ext)
        tflist.append(fn)
        print "Applying transform for image %s of %s..." %(i+1, len(flo_list))
        fn2=apply_transform(fn, ref, f, f.split('.')[0] + ext2)
        tfimlist.append(fn2)
    return tflist, tfimlist    
    
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
    
def invert_transform(tf, out, template=None, fl=None):
    '''
    Inputs:
    tf: Control point image of the transform you want to invert. Can also be an affine matrix file.
    fl: Floating image used for the registration.
    template: Template image used for the registration.
    out: Output filename.
    
    Output:
    out: Output filename.
    '''
    if tf[-4:] == '.nii':
        print "Inverting nonrigid transform."
        args = '-invNrr %s %s %s -ref %s' %(tf, fl, out, template)
        reg_t = [reg_transform] + args.split()      
        subprocess.call(reg_t, shell=True)
    elif tf[-4:] == '.txt':
        print "Inverting affine transform."
        args = '-invAff %s %s ' %(tf, out)
        reg_t = [reg_transform] + args.split()      
        subprocess.call(reg_t, shell=True)
    else:
        print "Cannot identify transform type."
    return out

def get_jacobian(df, template, out):
    '''
    Inputs:
    df: Deformation field or control point image.
    template: Template image.
    out: Output filename.
    
    Output:
    out: Output filename.
    '''
    args = '-cpp %s -target %s -jac %s' %(df, template, out)
    reg_j = [reg_jacobian] + args.split()  
    subprocess.call(reg_j, shell=True)
    return out

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
for f in filelist:
    im = f.replace(REG_NRR_EXT, FLOAT_EXT)
    
    if REG_AFF_EXT != None:
        aff=f.replace(REG_NRR_EXT, REG_AFF_EXT)
        tc, tf = combine_transforms(f, aff, REFIM, f[:-4] + 'CPPA.nii')
    else:
        tc = f
    inv = invert_transform(tc, f[:-4] + 'InvCPPA.nii', template=im, fl =REFIM)
    
    nrr=apply_transform(inv, im, LABEL_IM, f[:-4] + "_InvLblR.nii") #NRR file has affine and nrr transform in it
    
    #mask junk outside of brain
    n=nib.load(nrr)    
    a = n.get_affine()
    i=nib.load(im).get_data()
    
    d=n.get_data()
    d[i==0]=0  
    
    n2=nib.Nifti1Image(d, a)
    nib.save(n2, nrr.replace('.nii', 'R.nii'))
    
for f in glob.glob('*_InvLblRR.nii'):
    n=nib.load(REFIM)
    os.remove(f.replace('R.nii','.nii'))
    os.rename(f, f.replace('R.nii','.nii'))
