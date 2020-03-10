# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 08:28:24 2015

@author: Kristof Govaerts

This script iteratively registers all files in a directory you choose to a template image.
It creates a template image using rigid registration, then repeats the process with affine regisration,
updates the template and floating images, and repeats the process with nonrigid registration.
Adjust the N_ITERATIONS parameter to choose the number of images. NEUROMORPH_LOC needs to point
to the installation directory of Neuromorph (or any directory with all executables for NiftyReg).
THRESH can be adjusted if the background value for your images is not 0. 
"""

import os
import nibabel as nib
import numpy as np
import glob
import Tkinter, tkFileDialog
import subprocess
import scipy.ndimage as ndi

#Globals
FUNC_EXT='_basic_Full_V5mb.nii'
FUNC2_EXT = '_Full_V0.nii'
ANAT_PATTERN='*Brain_M4_FS.hdr' #file extension for your processed images to transform
RESHAPE_AX = 1 #axis along which to reshape - None if anatomical and functional imgs are oriented the same way

VOLUME_INDEX = 0 #which volume to register
NEUROMORPH_LOC = "C:/Program Files/NeuroMorph/Packages/MIC/Applications/Modules/Macros/KULNeuroMorph"
#default pars
BENDING_ENERGY = 0.05
GRID_SPACING = -5 #final grid spacing - if negative, is in voxels, if positive in mm
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
    omat = flt.split('.')[0] + mode + '.txt'    
    out = flt.split('.')[0] + mode + '.nii'
    if mode == 'A': #affine    
        args = '-flo %s -ref %s' %(flt, atl)         
    elif mode == 'R': #rigid
        args = '-flo %s -ref %s -rigOnly' %(flt, atl)   
        
    flirt = [reg_ala] + args.split()    
 
    process = subprocess.call(flirt, shell=True)
    os.rename('outputResult.nii', out) #NiftyReg saves as outputResult by default
    os.rename('outputAffine.txt', omat)

    return out, omat
    
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

def combine_transforms(t1, t2, ref, filename):
    '''Inputs:
    t1, t2: Filenames of transforms to combine (Should point to control point imgs)
    ref: Reference image used for one of the transforms
    filename: Filename for the output control point image
    Output:
    Filename of the combined transform control point image'''
    args = '-comp %s %s %s -ref %s' %(t1, t2, filename, ref)     
    tf = [reg_transform] + args.split()  
    subprocess.call(tf, shell=True)
    return filename
    
def apply_transform(t, ref, flo, filename):
    '''Inputs:
    t: Filename of a control point image or affine transformation matrix
    ref: Reference image
    flo: Floating image to which to apply the transform
    filename: Name to give the transformed floating image.
    Output:
    Filename of the transformed floating image.'''
    args = '-ref %s -flo %s -trans %s -res %s' %(ref, flo, t, filename)     
    tf = [reg_resample] + args.split()  
    print tf
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
    
def nrr_multiple_iterations(affine_temp, affine_list, n=5, thresh=7000, ftype='int32'):
    print "ITERATION 1: Registering affine transformed images to affine template."
    if len(glob.glob('*NRR1.nii')) == len(affine_list):
        print "Registration already performed. Skipping iteration 1."
        nrr_list1 = glob.glob('*NRR1.nii')
        tf_list1 = glob.glob('*CPP1.nii')
    else:
        nrr_list1, tf_list1 = reg_batch(affine_temp, affine_list, mode='NRR', ext='1')   
        
    print "Updating template: Averaging images..."
    nrr_temp = 'nrr_template1.nii'
    save_sum_img(nrr_list1, nrr_temp, thresh=thresh, ftype=ftype)
    
    for i in range(n-1):
        print "ITERATION %s: Registering nonrigidly transformed images to template." %(i+2)
        e = '*NRR%s.nii' %(i+2)
        e2 = '*CPP%s.nii' %(i+2)
        if len(glob.glob(e)) == len(affine_list):
            print "Registration already performed. Skipping iteration %s." %(i+2)
            nrr_list2 =  glob.glob(e)
            tf_list2 = glob.glob(e2)
        else:
            nrr_list2, tf_list2 = reg_batch(nrr_temp, nrr_list1, mode='NRR', ext=(i+2))
        
        print "Updating template: Combining iterations."
        tflist, tfimlist = batch_combine_apply_tf(nrr_temp, affine_list, tf_list1, tf_list2)
        
        print "Averaging images..."
        nrr_temp='nrr_template%s.nii' %(i+2)
        save_sum_img(tfimlist, nrr_temp, thresh=thresh, ftype=ftype)
        
        tf_list1 = tflist #Keep updating your average transforms
        nrr_list1 = tfimlist
    return nrr_temp, tflist, tfimlist

    
    
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
    
def resave_with_affine(fn, aff, fn2):
    '''Resaves image with a new affine matrix.'''
    im = nib.load(fn)
    nif = nib.Nifti1Image(im.get_data(), aff)
    nib.save(nif, fn2)
    
def bfc(ndim, sigma=3.5, mode='2D', thresh=1e5):   
    '''Rudimentary bias field correction according to homomorphic unsharp masking.
    Vovk et al., IEEE Trans Med Imaging. 2007 Mar;26(3):405-21.
    
    Input: 
    ndim = 3D numpy array.
    sigma = sigma to use for the low-pass filter.
    mode = 2D or 3D. Depends on the MRI acquisition. Very thin, adjoining slices can count as 3D.
    thresh = How to threshold the input array. Made to ignore zeros.
    
    Output:
    bfc = Output bias field corrected array.'''     
    ndim[ndim<thresh]=0
    #bfc: homomorphic unsharp masking
    if mode == '3D':
        s = (sigma, sigma, sigma)
    elif mode == '2D':
        s = (sigma, sigma, 0)
    bfc = (ndim * np.mean(ndim[np.nonzero(ndim)]))/ndi.filters.gaussian_filter(ndim, sigma=s)
    return bfc
    
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

filelist=glob.glob('*' + FUNC_EXT)
fl_a=[]
fl_b=[]
fl_c=[]

print "STEP 1: Rigid reg of functional image to anatomical image."
for n, f in enumerate(filelist):
    print "File %s of %s." %(n+1, len(filelist))
    n_vols=0
    a=f.replace(FUNC_EXT, '')[:-8]
    a=glob.glob(a + ANAT_PATTERN)#get anatomical img
    f2 = f.replace(FUNC_EXT, FUNC2_EXT)
    if len(a) == 0:
        print "No anatomical img found for", f
        continue
    a=a[0] #take the first if you have multiple anatomical imgs
    if RESHAPE_AX != None:
        f_im=nib.load(f)
        f2_im = nib.load(f2)
        f_data=f_im.get_data()
        f_data = bfc(f_data, thresh=0) #correct for bias field
        f2_data=f2_im.get_data()
        b_data=np.swapaxes(f_data, 1, 2)
        b2_data=np.swapaxes(f2_data, 1, 2)
        aff = f_im.get_affine()
        b_aff = np.zeros(aff.shape)
        b_aff[0,:] = aff[0,:]
        b_aff[1,:] = aff[2,:]
        b_aff[2,:] = aff[1,:]
        b_aff[2,2]=b_aff[2,1]
        b_aff[2,1]=0
        b_aff[1,1]=b_aff[1,2]
        b_aff[1,2]=0
        b_aff[3,:] = aff[3,:]
        #IMPORTANT STEP - Janaki's conversion tool puts the slice thickness and
        #not the interslice dist in the Z-axis -This will make rigid registration fail
        #as the brain will not be as long as it should be!
        if np.round(b_aff[1,1], decimals=2)==0.6:
            b_aff[1,1]=0.72
            print "Z-dimension reshaped from 0.6 to 0.72."
        elif np.round(b_aff[1,1], decimals=2)==0.5:
            b_aff[1,1]=0.6
            print "Z-dimension reshaped from 0.5 to 0.6."
        #b_aff *= np.diag([10,10,10,1])#FSL works with larger matrix sizes! Scale by 10

        b_data=b_data[:,::-1,::-1]
        b2_data=b2_data[:,::-1,::-1]
        b_nif=nib.Nifti1Image(b_data, b_aff)
        b2_nif=nib.Nifti1Image(b2_data, b_aff)
        nib.save(b_nif, f.replace('.nii', 'cor.nii'))
        nib.save(b2_nif, f2.replace('.nii', 'cor.nii'))
        f=f.replace('.nii', 'cor.nii')
        f2=f2.replace('.nii', 'cor.nii')
        
    a1=nib.load(a)
    orig_affine=a1.get_affine()
#    a2 = nib.Nifti1Image(a1.get_data(), orig_affine*np.diag([10,10,10,1]))
#    a=a[:-4] + 'x10.nii'
#    nib.save(a2, a)#resave anatomical img - x10 resolution

    #aff, mat = RegAndResample(f,a, mode='R')
    aff, mat = (f.replace('.nii', 'R.nii'), f.replace('.nii', 'R.txt')) #for just doing resampling
    aff2 = apply_transform(mat, a, f2, f2.replace('.nii', 'R.nii'))
    
    #resave all images with proper resolution
    
    #resave_with_affine(aff, orig_affine, aff.replace('.nii', 'R.nii')) 
    resave_with_affine(aff2, orig_affine, aff2.replace('.nii', 'R.nii'))     
    
    fl_a.append(a)
    fl_b.append(aff)
    fl_c.append(aff2)
#    os.remove(a)
#    os.remove(aff)
#    os.remove(aff2)
    
with open('anat_list.txt', 'w') as f:
    f.write('\n'.join(fl_a))
    
with open('func_list.txt', 'w') as f:
    f.write('\n'.join(glob.glob('*V1corA.nii')))

with open('func2_list.txt', 'w') as f:
    f.write('\n'.join(fl_c))