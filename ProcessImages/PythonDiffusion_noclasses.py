#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:     Diffusion processing
#
# Author:      s0203524
#
# Created:     10/09/2013
# Copyright:   (c) s0203524 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

'''
Diffusion datasets are 4-dimensional matrices, with x,y, and z being the spatial
dimensions and t being the dimension holding b-values and diffusion dirs.

bvals is a 1-dimensional array with length t. It holds the b-value for each
point in the t direction. Although this b-value can vary depending on the
direction due to the acquisition, the algorithm requires that each cluster of
diffusion directions has a single homogenous b-value (since you need at least 6
dirs per b-value to calculate a tensor, and otherwise you will only have 1).

bvecs is a 2-dimensional array with width 3 and length t. It holds three
diffusion directions for each point in the t direction.

NOTE: 3D images are formatted as (x,y,z,t), whereas 2D images are formatted as
(x,y,t,z).
'''

def main():
    os.chdir('/media/sf_host/data/kristof/Diffusion')
    file='FVL_23225_1L_a_iI1_4_1'

    bvs=[2500.0]


if __name__ == '__main__':
    main()

import os
import dipy
from dipy.core.gradients import *
from dipy.align import aniso2iso
import numpy as np
import scipy
from ReadInput import *
import nibabel as nib
from datetime import datetime
from dipy.core.gradients import gradient_table
import subprocess

def main():
    os.chdir('/media/sf_host/data/kristof/Diffusion')
    file='FVL_25117_NC_2m_a_lf1_7_1'

def processDiffusion(file, ds=False, ec=False, bvs=None):
    '''
    Process a diffusion-weighted dataset.
    
    file=filename of Analyze or Nifti file (without extension)
    bvs = list of b-values (optional)
    ec = Whether or not eddy current correction should be applied (takes a while,
         and does not work  inside Spyder IDE but only from command line)
    ds = Whether or not the image should be downsampled to an isotropic voxel size
    
    2D images should be formatted as (x,y,t,z) and 3D images as (x,y,z,t)
    which is standard for Bruker files.
    This protocol automatically determines how many diffusion
    dirs there are and how many b-values according to what is saved in the
    accompanying text file. You can provide a list of exact b-values, but if
    you do not the program will calculate mean b-values for each cluster of
    diffusion directions based on the text file.
    '''
    
    dims=list_values(read_line('VisuCoreSize=', file))
    ext=checkFileType(file)
    img=nib.load(file+ext)
    data=img.get_data()
    affine=img.get_affine()
    bvals, avbvals, dwgrad, dwdir, nA0, nbvals, ndirs=getDiffusionPars(file)  
    if len(dims)==2: 
        #2D arrays are arranged differently. scipy.swapaxes is not sufficient for 
        #Paravision's Fortran-style column-major ordering as the t-axis is ordered differently.
        newshape=(data.shape[0],data.shape[1],data.shape[3],data.shape[2])
        print '2D array with shape %r. Reshaping to %r in Fortran-style column major.' %(data.shape,newshape)
        data=np.reshape(data,newshape,order='F')
    rescaleImage(file,data,nbvals,dims)
    img=nib.Nifti1Image(data, affine)
    
    if ds:
        print "Voxel size nonisotropic. Downsampling..."
        data, affine=downsampleImage(img)
        img=nib.Nifti1Image(data, affine)
    else:
        affine=img.get_affine()
        data=img.get_data()
        
    thresh=np.mean(data[:5,:5,:,0])
    mask=data[...,0] > 2.5*thresh
    for i in range(data.shape[3]):
        data[:,:,:,i]*=mask
    if ec:
        starttime=datetime.now()
        print "Applying eddy current correction."
        img=eddyCorrection(img, file+'_Eddy.nii')
        data=img.get_data()
        affine=img.get_affine()
        time=datetime.now()-starttime
        print "Eddy current correction completed in %r seconds." %time.seconds
        
    if bvs==None:
        bvalmat=np.array(avbvals)
        bvalmat[bvalmat<10]=0
    else:
        bvalmat=np.zeros([nA0+(ndirs*len(bvs))]) #entered pars
        for i,b in enumerate(bvs):    #unfortunately the ideal b-vals(not effective b-vals) are not in the text file. Have to enter manually and convert to appropriate matrix
            bvalmat[nA0+ndirs*i:]=b
            
    bvecmat=np.zeros([nA0+ndirs*nbvals, 3])
    for i in range(nbvals):   
        bvecmat[nA0+ndirs*i:nA0+ndirs*(i+1),:]=dwdir #fills b-vector matrix with the different diffusion dirs
    
    if len(bvecmat) != len(bvals):
        print "Error. Cannot process this image."
        
    print bvalmat.shape
    print dwgrad.shape
    gtab=gradient_table(bvalmat, bvecmat) #creates a gradient table with b-vals and diffusion dirs for processing
    
    from dipy.reconst.dti import TensorModel
    
    starttime=datetime.now()
    print "Fitting tensor model."
    ten = TensorModel(gtab)
    tenfit = ten.fit(data, mask)
    time=datetime.now()-starttime
    print "Tensor fit completed in %r seconds." %time.seconds
    
    from dipy.reconst.dti import fractional_anisotropy
    evecs=tenfit.evecs #eigenvectors
    fa = fractional_anisotropy(tenfit.evals)
    fa=np.clip(fa,0,1) #removes voxels where fit failed by thresholding at 0 and 1
    md=tenfit.md
    md[np.isnan(md)] = 0 #removes voxels where fit failed
    print "Calculated eigenvectors, MD and FA."
    
    from dipy.reconst.dti import color_fa
    cfa = color_fa(fa, tenfit.evecs)
    
    return tenfit, cfa, bvalmat, dwgrad, bvecmat


def getDiffusionPars(file):
    '''Reads text files accompanying images and extracts B-value matrices,
    diffusion dirs and number of A0 images.
    
    Input:
    file=filename, no extension.
    
    Output:
    bvals=B-values, ordered as in the .txt file
    avbvals=averaged b-vals so that b-vals are homogenous for each cluster of diffusion dirs
    dwgrad=gradient vectors. Essentially non-normalized diffusion dirs
    dwdir=unit vectors representing diffusion dirs. This is the input dipy accepts
    nA0=number of A0 volumes
    nbvals = number of b-vals
    ndirs = number of diffusion dirs'''
    bvals=np.array(list_values(read_line('DwEffBval=', file)))
    nA0=sum(i<10 for i in bvals)
    dwdir=list_values(read_line('DwDir=', file)) #Gets all diffusion directions. The length of this parameter should be 3*amount of dirs (3 coordinates per dir) Note that A0 images do not have a diffusion direction
    dwdir=vector_array(dwdir, 3)
    dwgrad=list_values(read_line('DwGradVec=', file))
    dwgrad=vector_array(dwgrad, 3)
    ndirs=len(dwdir)
    
    nbvals=int((len(bvals)-nA0)/ndirs)    #there are as many b-values as there are values in the t-dimension.
    avbvals=averageBvals(bvals,nA0,ndirs,nbvals)
    print "Selected image has %r A0 scans. Found %r nonzero b-value(s) and %r diffusion directions." %(nA0, nbvals, ndirs)
    return bvals, avbvals, dwgrad, dwdir, nA0, nbvals, ndirs
    
def averageBvals(bvals, nA0, ndirs, nbvals):
    '''Returns average b-values by removing the nA0 first values, reshaping into
    an ndirsxnbvals matrix and calculating the mean bvals before reshaping back
    into an ndirs*nbvals 1-D array.'''
    shortbvals=bvals[nA0:]
    rbvals=shortbvals.reshape(ndirs, nbvals)
    bvals2=np.zeros(bvals.shape)
    rbvals2=np.zeros(rbvals.shape)
    for i in range(nbvals):
        rbvals2[:,i]=np.mean(rbvals[:,i])
    rbvals2=rbvals2.reshape(ndirs*nbvals)
    bvals2[nA0:]=rbvals2
    return bvals2

def rescaleImage(file, data, nbvals,dims):
    '''Rescales image according to parameters in accompanying text file. File
    is the filename(without extension), data is the data matrix and nbvals is the
    amount of b-values. Does not return anything because this modifies the array.'''
    slopes=np.array(list_values(read_line('VisuCoreDataSlope=', file)))
    
    if len(dims)==3:
        for i in range(data.shape[3]):
            data[:,:,:,i]=data[:,:,:,i]*slopes[i]
    elif len(dims)==2:
        rslopes=np.reshape(slopes,data.shape[:-3:-1])
        rslopes=scipy.swapaxes(rslopes,0,1)
        for j in range(data.shape[2]):
            for i in range(data.shape[3]):
                data[:,:,j,i]=data[:,:,j,i]*rslopes[j,i]
    else:
        print "No b-values found or cannot process slope formatting."

def downsampleImage(img):
    '''Returns an image downsampled to its largest voxel dimension.

    Inputs:
        img = Nibabel Image object (Nifti1, Analyze,...)
    Outputs:
        data2 = Downsampled image data
        affine2 = New affine matrix'''
    d=img.get_data()
    zooms=img.get_header().get_zooms()[:3]  #only need 3 zooms, fourth one is 1 since this is the b-val dir
    affine=img.get_affine()
    nzooms=[np.max(zooms) for i in range(3)]
    
    data2, affine2=aniso2iso.resample(d, affine, zooms, nzooms)
    return data2, affine2

def unitVector(vecmat):
    '''Normalizes 2D input array of vectors into unit vectors (e.g. vectors
    with length 1). This is necessary because this is the only type of
    vector this software accepts.'''
    output=np.zeros(vecmat.shape)
    for i, vec in enumerate(vecmat):
        unitvec=vec/np.linalg.norm(vec)
        output[i,:]=unitvec
    output[np.isnan(output)] = 0
    return output

def tensorTractography(tenfit, img, fa=0.3):
    fa=tenfit.fa
    fa[np.isnan(fa)] = 0
    evecs=tenfit.evecs
    hdr = nib.trackvis.empty_header()
    hdr['voxel_size'] = img.get_header().get_zooms()[:3]
    hdr['voxel_order'] = 'LAS'
    hdr['dim'] = fa.shape
    
    from dipy.data import get_sphere #EuDX needs to map voxel dirs on a sphere
    sphere = get_sphere('symmetric724')
    from dipy.reconst.dti import quantize_evecs
    peak_indices = quantize_evecs(evecs, sphere.vertices)
    
    from dipy.tracking.eudx import EuDX
    eu = EuDX(fa, peak_indices, odf_vertices = sphere.vertices, a_low=fa)
    tensor_streamlines = [streamline for streamline in eu]
    
    tensor_streamlines_trk = ((sl, None, None) for sl in tensor_streamlines)
    ten_sl_fname = 'tensor_streamlines_tests.trk'
    nib.trackvis.write(ten_sl_fname, tensor_streamlines_trk, hdr, points_space='voxel')

def eddyCorrection(img, filename, ref=0):
    '''
    eddyCorrection(img, filename, ref=0)
    
    Uses FSL's eddy_correct program to affinely register your image slices
    to a reference volume.
    
    Unfortunately, this protocol DOES NOT WORK in the Spyder IDE. Run it from
    a Python console in the command line.
    
    Inputs:
        img: Nibabel Image class (Nifti1Image, AnalyzeImage). Should be
             4-D (x,y,z,t)
        filename: Preferred output filename. Should end with .nii
        ref: Reference volume. Image is a composite of t volumes with (x,y,z)
             dimensions. Default reference volume is at t=0.
             
    Outputs:
        nimg = Eddy-corrected image.'''
    nifti=nib.Nifti1Image(img.get_data(), img.get_affine())
    temp='fsl_temp.nii'
    nib.save(nifti, temp)
    infile=os.path.join(os.getcwd(), temp)
    outfile=os.path.join(os.getcwd(), filename)
    pars="%s %s %r" %(infile, outfile, ref)
    command=["eddy_correct "+pars]
    subprocess.call(command, shell=True)
    os.remove('fsl_temp.nii')
    nimg=nib.load(outfile+'.gz')
    return nimg

