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
    import os
    import nibabel as nib
    os.chdir('/media/sf_host/data/kristof/Diffusion/multibv')
    os.chdir('C:/data/kristof/Diffusion/multibv')
    ecf=nib.load('FVL_45368_1L1R_1y_a_ms1_8_1_EC.nii.gz').get_data()
    file='FVL_45368_1L1R_1y_a_ms1_8_1'
    from ProcessImages.imageTypes import *
    from ProcessImages.Saveoutput import *
    dif=DiffusionImg(file)
    dif.pdata=ecf


#    dif.processDiffusion(ec=True)
#
#    bvs=[2500.0]


if __name__ == '__main__':
    main()

import os
import dipy
from dipy.core.gradients import *
from dipy.align import aniso2iso
import numpy as np
import scipy
import scipy.ndimage as ndi
from ReadInput import *
import nibabel as nib
from datetime import datetime
from dipy.core.gradients import gradient_table
import subprocess
import warnings
import pickle
from tensorMaths import construct_bvecs, calculate_b_matrix
from DSI_studio import *
from ImProcessing import brain_mask, largest_component
#import nipype.interfaces.fsl as fsl
#import nipy.algorithms.registration as nar

#def main2():
#    os.chdir('/media/sf_host/data/kristof/Diffusion')
#    file='FVL_25117_NC_2m_a_lf1_7_1'

class DKIFit:
    def __init__(self, adc, akc, si0, adcstd, akcstd, si0std, bvecs):
        self.adc = adc
        self.akc = akc
        self.si0 = si0
        self.adcstd = adcstd
        self.akcstd = akcstd
        self.si0std = si0std
        self.bvecs = bvecs
        
    def calc_evals(self):
        self.adc_evals, self.adc_evecs=calculate_tensor(self.adc, self.bvecs)
        self.akc_evals, self.akc_evecs=calculate_tensor(self.akc, self.bvecs)
        
    def fa(self):
        return fractional(self.adc_evals)
        
    def fk(self):
        return fractional(self.akc_evals)
        
    def md(self):
        return np.mean(self.adc_evals, -1)
        
    def mk(self):
        return np.mean(self.akc_evals, -1)
        
    def ad(self):
        return self.adc_evals[...,0]
        
    def ak(self):
        return self.akc_evals[...,0]
        
    def rd(self):
        return np.mean(self.adc_evals[...,1:], -1)
        
    def rk(self):
        return np.mean(self.akc_evals[...,1:], -1)
        
    def color_fa(self):
        fa=self.fa
        evecs=self.evecs
        if (fa.shape != evecs[..., 0, 0].shape) or ((3, 3) != evecs.shape[-2:]):
            raise ValueError("Wrong number of dimensions for evecs")
    
        fa = np.clip(fa, 0, 1)
        rgb = np.abs(evecs[..., 0]) * fa[..., None]
        return rgb

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
    mask=largest_component(brain_mask(data))
    for i in range(data.shape[3]):
        data[:,:,:,i]*=mask

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
        
#    if self.nrep>1:
#        print "Image has %r repetitions. Correcting appropriately."
#        avbv=np.array(bvalmat)
#        tbvec=np.array(bvecmat)
#        for c in range(self.nrep-1):
#            np.concatenate(bvalmat,avbv)
#            np.concatenate(bvecmat,tbv)

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

def getDiffusionPars(self):
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
    from PythonDiffusion import averageBvals
    self.bvals=np.array(list_values(read_line('DwEffBval=', self.name)))
    self.nA0=list_values(read_line('DwAoImages=', self.name))
    if len(self.nA0)==1:
        self.nA0=int(self.nA0[0])
    else:
        self.nA0=0
    self.dwdir=list_values(read_line('DwDir=', self.name)) #Gets all diffusion directions. The length of this parameter should be 3*amount of dirs (3 coordinates per dir) Note that A0 images do not have a diffusion direction
    if len(self.dwdir)>0:
        self.dwdir=vector_array(self.dwdir, 3)
    else:
        self.dwdir=[0.0]
    self.dwgrad=list_values(read_line('DwGradVec=', self.name))
    if len(self.dwgrad)>0 and len(self.dwgrad)%3.0==0.0: #DTI scans with very large amounts of volumes often have truncated GradVec lines - make sure this is divisible by 3
        self.dwgrad=vector_array(self.dwgrad, 3)
    else:
        self.dwgrad=[0.0]
    self.ndirs=len(self.dwdir)

    self.nbvals=int((len(self.bvals)-self.nA0)/self.ndirs)    #there are as many b-values as there are values in the t-dimension.
    self.avbvals=averageBvals(self.bvals,self.nA0,self.ndirs,self.nbvals)
    if self.nA0 and self.nbvals>0 and self.ndirs>0:
        print "Selected image has %r A0 scans. Found %r nonzero b-value(s) and %r diffusion directions." %(self.nA0, self.nbvals, self.ndirs)
    else:
        warnings.warn("Error: Could not find all diffusion parameters. Further processing will not be possible.")


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

def rescaleImage(self):
    '''Rescales image according to parameters in accompanying text file. File
    is the filename(without extension), data is the data matrix and nbvals is the
    amount of b-values. Does not return anything because this modifies the array.''' 
    
    if self.dim==3:
        for i in range(self.shape[3]):
            self.pdata[:,:,:,i]=self.pdata[:,:,:,i]*self.slopes[i]
    elif self.dim==2:
        rslopes=np.reshape(self.slopes,self.shape[:-3:-1])
        rslopes=scipy.swapaxes(rslopes,0,1)
        for j in range(self.shape[2]):
            for i in range(self.shape[3]):
                self.pdata[:,:,j,i]=self.pdata[:,:,j,i]*rslopes[j,i]
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

def tensorTractography(tenfit, img, falow=0.3):
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
    eu = EuDX(fa, peak_indices, odf_vertices = sphere.vertices, a_low=falow, ang_thr=30.0)
    tensor_streamlines = [streamline for streamline in eu]
    
    tensor_streamlines_trk = ((sl, None, None) for sl in tensor_streamlines)
    ten_sl_fname = 'tensor_streamlines_tests.trk'
    nib.trackvis.write(ten_sl_fname, tensor_streamlines_trk, hdr, points_space='voxel')

def eddyCorrection(self, filename, protocol='eddy_correct', ref=0):
    '''
    eddyCorrection(img, filename, ref=0)

    Uses FSL's eddy_correct program to affinely register your image slices
    to a reference volume.

    Unfortunately, this protocol DOES NOT WORK in the Spyder IDE (permission errors)
    unless you are using the eddy_v protocol. Run it from
    a Python console in the command line.

    Inputs:
        img: Nibabel Image class (Nifti1Image, AnalyzeImage). Should be
             4-D (x,y,z,t)
        filename: Preferred output filename. Should end with .nii
        ref: Reference volume. Image is a composite of t volumes with (x,y,z)
             dimensions. Default reference volume is at t=0.

    Outputs:
        nimg = Eddy-corrected image.'''
    temp='fsl_temp.nii'
    infile=os.path.join(os.getcwd(), temp)
    outfile=os.path.join(os.getcwd(), filename)
    if protocol == 'eddy_correct':
        save_mask(self)
        mask=self.mask
        for bd in range(self.shape[3]):
            self.pdata[:,:,:,bd]*=mask
        nifti=nib.Nifti1Image(self.pdata, self.affine)
        nib.save(nifti, temp)
        pars="%s %s %r" %(infile, outfile, ref)
        command=["eddy_correct "+pars]
        subprocess.call(command, shell=True)
        os.remove('fsl_temp.nii')
        nimg=nib.load(outfile+'.nii.gz')
        self.pdata=nimg.get_data()
    
    elif protocol == 'eddy':
        print "Running FSL's 'eddy' protocol."
        create_acqp(self)
        print "ACPQ parameters are:", self.acpq
        save_bvals_bvecs(self)
        save_mask(self)
        temp='fsl_temp.nii'
        nifti=nib.Nifti1Image(self.pdata, self.affine)
        nib.save(nifti, temp)
        infile=os.path.join(os.getcwd(), temp)
        mask=os.path.join(os.getcwd(), self.name+"_mask.nii")
        acqp=os.path.join(os.getcwd(), self.name+"_acqp.txt")
        index=os.path.join(os.getcwd(), self.name+"_index.txt") 
        bvecs=os.path.join(os.getcwd(), self.name+"_bvecs.txt") 
        bvals=os.path.join(os.getcwd(), self.name+"_bvals.txt")
        outfile=os.path.join(os.getcwd(), filename)
        command="eddy --imain=%s --mask=%s --acqp=%s --index=%s --bvecs=%s --bvals=%s --out=%s" %(infile, mask, acqp, index, bvecs, bvals, outfile)        
        subprocess.call(command, shell=True)
        os.remove('fsl_temp.nii')
        nimg=nib.load(outfile+'.nii.gz')
        self.pdata=nimg.get_data()
    
    elif protocol == 'eddy_v':
        shifts,simila=eddy_correct_vertical(self)
        shst=["Diffusion experiment: 0 Shift (px): 0"]
        with open(self.name +"_shifts.txt", 'w') as txt:
            for i in range(1,self.shape[3]):
                shst.append("Diffusion experiment: "+str(shifts[i, 0])+" Shift (px): " +str(shifts[i, 1]))
            st='\n'.join(shst)
            txt.write(st)   
        corrimg=nib.Nifti1Image(self.pdata, self.affine)
        nib.save(corrimg, filename)
    
def eddy_correct_vertical(self, dim=1):
    '''Runs eddy current correction by sampling 10 vertical 1-pixel shifts
    and settling on the shift with the highest similarity to the reference 
    image (b0 image 1).
    
    Inputs: 
        self = DiffusionImg object
        dim = Vertical dimension (typically 1)
        mode = Whether to do correction in 2D or 3D. 3D is faster.
    
    Outputs:
        simila = shape[4]-length list of similarity values for each of the shifts
        px = shape[4]-length list of shift magnitude'''
    simila=[]
    px=[]
    now=datetime.now()
    if self.shape[2]>10:
        c=int(self.shape[2]/2)
        r=range(c-3, c+3)
    elif self.shape[2]<=10 and self.shape[2] >2:
        c=int(self.shape[2]/2)
        r=range(c-1, c+1)
    elif self.shape[2] <= 2:
        r=[self.shape[2]]
    pxa=np.zeros([self.shape[2], self.shape[3]])
    for sli in r:
        print "Calculating eddy current shift for slice %r." %sli
        ref=self.pdata[:,:,sli,0]
        for bd in range(1, self.shape[3]):
            arr=self.pdata[:,:,sli,bd]
            #slices=[]
            similarities=[]
            pxr=range(-10,11)
            for i in pxr:
                curslice=np.roll(arr, i, dim)
                #slices.append(curslice)
                similarities.append(simil(self, ref, curslice))
            simila.append((sli, similarities))
            indices=[i for i,j in enumerate(similarities) if j==max(similarities)]
            ind=min(np.abs(indices))
            px.append(pxr[ind])
            pxa[sli, bd]=pxr[ind]
            #self.pdata[:,:,sli,bd]=slices[ind]
    shifts=np.zeros([self.shape[3],2])
    for bd in range(1, self.shape[3]):
        shift=int(np.median(pxa[r,bd]))
        self.pdata[:,:,:,bd]=np.roll(self.pdata[:,:,:,bd], shift, dim)
        shifts[bd,0]=bd
        shifts[bd,1]=shift

    time=datetime.now()-now
    print "Eddy current correction completed in %r seconds." %time.seconds
    return shifts, simila                 

def simil(self, ar1, ar2, bins=256.0, similarity='crl1'):
    '''Calculates similarity between two image arrays.'''
    import nipy.algorithms.registration as nar    
    if ar1.shape != ar2.shape:
        print "Input images do not have the same shape."
    if len(ar1.shape)==2:
        ar1=np.reshape(ar1,(ar1.shape[0], ar1.shape[1], 1))
        ar2=np.reshape(ar2,(ar2.shape[0], ar2.shape[1], 1))
    im1=nib.Nifti1Image(ar1, self.affine)
    im2=nib.Nifti1Image(ar2, self.affine)
    
    reg=nar.HistogramRegistration(im1,im2, bins=bins, similarity=similarity)
    simil=reg.eval(nar.Affine())
    return simil
    
       
#    elif protocol == 'flirt':
#        if sli==None:
#            r=range(1, self.shape[3])
#        else:
#            r=[sli]
#        for sl in r:
#            refdata=self.pdata[:,:,sl,ref]
#            refnifti=nib.Nifti1Image(refdata, self.affine)
#            nib.save(refnifti, temp)
#            for bd in range(1, self.shape[3]):
#                pdata=np.array(self.pdata)
#                print "Registering slice", sl+1, ", diffusion experiment", bd+1
#                temp2='fsl_temp_2.nii'
#                temp3='fsl_temp_3.nii.gz'
#                slicedata=self.pdata[:,:,sl,bd]
#                slicenifti=nib.Nifti1Image(slicedata, self.affine)
#                nib.save(slicenifti, temp2)
#                command = "flirt -in %s -ref %s -out %s -nosearch -cost mutualinfo -2D" %(temp2, temp, temp3)
#                subprocess.call(command, shell=True)
#                rslice=nib.load(temp3).get_data()
#                pdata[:,:,sl,bd]=rslice
#                os.remove('fsl_temp_2.nii')
#                os.remove('fsl_temp_3.nii.gz')
#            os.remove('fsl_temp.nii')
#        self.pdata=pdata
#        outnifti=nib.Nifti1Image(pdata, self.affine)
#        nib.save(outnifti, self.name+'_EC_f.nii')
    
def create_acqp(self):
    self.readdir="L_R"
    self.echospacing=0.256
    self.epifact=32
    if self.readdir=='L_R':
        par1="0 1 0 " #specifies what dir the PE direction is. In this case A-P
    #par2=str(self.echospacing*0.001*self.epifact) #time between center of first and last echoes in s
    par2=str(0.01*10.052)    
    st=par1+par2
    with open(self.name+"_acqp.txt", 'w') as txt:
        txt.write(st)
    st2="1 "*self.shape[3]
    with open(self.name+"_index.txt", 'w') as txt:
        txt.write(st2)
    self.acpq=st

def save_bvals_bvecs(self):
    bvstr=" ".join([str(b) for b in self.bvals])
    with open(self.name+"_bvals.txt", 'w') as txt:
        txt.write(bvstr)
    A0dir="0 "*self.nA0
    dirx=(" ".join([str(b) for b in self.dwdir[:,0]])+" ")*self.nbvals
    diry=(" ".join([str(b) for b in self.dwdir[:,1]])+" ")*self.nbvals
    dirz=(" ".join([str(b) for b in self.dwdir[:,2]])+" ")*self.nbvals
    joiner='\n'
    dirstr=A0dir + dirx + joiner + A0dir + diry + joiner + A0dir + dirz
    with open(self.name+"_bvecs.txt", 'w') as txt:
        txt.write(dirstr)

def save_mask(self):
    thresh=np.mean(self.pdata[:5,:5,:,0])
    mask=self.pdata[...,0] > 7*thresh
    self.mask=mask
    nif=nib.Nifti1Image(mask.astype(int),self.affine)
    nib.save(nif, self.name+'_mask.nii')

def tensorFit(self, m='dti', bv=None, removea0=0, mask=True):
    if bv != None:
        self.backups=[self.pdata, self.avbvals, self.nbvals]
        self.pdata, self.avbvals=selectBvals(self, bv)
        self.avbvals=np.array(self.avbvals)
        self.nbvals=len(bv)
        
    if removea0!=0:
        print "Removing %s A0 images." %removea0
        self.pdata=self.pdata[:,:,:,removea0:]
        self.nA0-=removea0
        self.avbvals=self.avbvals[removea0:]
        self.bvals=self.bvals[removea0:]
        
#    thresh=np.mean(self.pdata[:5,:5,:,0])
    if mask:
        print "Generating brain mask: Erosion>Dilation>MedianOtsu>LargestComponent"
        self.mask=ndi.binary_dilation(ndi.binary_erosion(ndi.binary_dilation(brain_mask(self.pdata))))
        self.mask=largest_component(self.mask)
    else:
        print "Not using brain mask."
        self.mask=np.zeros(self.shape[:-1]) + 1
    
    bvalmat=np.array(self.avbvals)
    self.bvecmat=construct_bvecs(self)
          
    if len(self.bvecmat) != len(self.avbvals):
        print "Error. Cannot process this image."
        
    print bvalmat.shape
    self.gtab=gradient_table(bvalmat, self.bvecmat) #creates a gradient table with b-vals and diffusion dirs for processing
    
    if m == 'dti':
        from dipy.reconst.dti import TensorModel
        
        starttime=datetime.now()
        print "Fitting tensor model."
        ten = TensorModel(self.gtab)
        self.tenfit = ten.fit(self.pdata, self.mask)
        time=datetime.now()-starttime
        tenfile=self.name+"_tenfit.pickle"
        print "Tensor fit completed in %r seconds. Saving file %s." %(time.seconds, tenfile)
#        with open(tenfile, 'w') as f: TOO INEFFICIENT
#            pickle.dump([self, self.tenfit], f)
        
        from dipy.reconst.dti import color_fa
        self.cfa = color_fa(self.tenfit.fa, self.tenfit.evecs)
        
    if m == 'RESTORE':
        from dipy.reconst.dti import TensorModel
        mean_std = 2.5 * np.mean(np.std(self.pdata[..., self.gtab.b0s_mask], -1)) #conservative thresh
        starttime=datetime.now()
        print "Fitting tensor model using the RESTORE method. Sigma=", mean_std
        ten = TensorModel(self.gtab, fit_method="RESTORE", sigma=mean_std)
        self.tenfit = ten.fit(self.pdata, self.mask)
        time=datetime.now()-starttime
        tenfile=self.name+"_tenfit.pickle"
        print "Tensor fit completed in %r seconds. Saving file %s." %(time.seconds, tenfile)
#        with open(tenfile, 'w') as f: TOO INEFFICIENT
#            pickle.dump([self, self.tenfit], f)
        
        from dipy.reconst.dti import color_fa
        self.cfa = color_fa(self.tenfit.fa, self.tenfit.evecs)
        
    elif m == 'dki':
        from dipy.reconst.dki import DiffusionKurtosisModel
        starttime=datetime.now()
        print "Fitting kurtosis tensor."
        kurt = DiffusionKurtosisModel(self.gtab)
        self.tenfit = kurt.fit(self.pdata)
        time=datetime.now()-starttime
        print "Kurtosis fit completed in %r seconds." %time.seconds
    
    if bv != None:
        self.pdata=self.backups[0]
        self.avbvals=self.backups[1]
        self.nbvals=self.backups[2]

def selectBvals(self,bv):
    '''Isolates only the b-values you would like to fit. 'bv' needs to be a list.'''
    realbvs=[0]*self.nA0
    for i in range(len(bv)):
        realbvs.append(min(self.avbvals, key=lambda x:abs(x-bv[i])))
    print "Isolating the following b-values:\n",realbvs
    bvb=np.zeros(self.shape, dtype=bool) #Creates array of negative booleans to decide which b-values to include or omit
    bvb[:,:,:,:self.nA0]=True #make sure you include A0 imgs    
    for i in range(len(self.avbvals)):
        if self.avbvals[i] in realbvs:
            bvb[:,:,:,i]=True
    simg = self.pdata[bvb] #masks original image
    sbvs=[0]*self.nA0
    for i in range(self.ndirs):
        sbvs.extend(realbvs[self.nA0:])
    print len(sbvs)
    simg=np.reshape(simg,[self.shape[0], self.shape[1], self.shape[2], len(sbvs)])
    return simg, sbvs
    
def calc_diffusionpars(dmap):
    '''Helper function for class DKIFit. Calculates mean value over last dimension,
    eliminating values under zero.'''
    c=np.array(dmap)
    c[c<0]=nan
    ma=np.ma.masked_array(c, np.isnan(c))
    mean = np.mean(ma,-1)
    return mean 
    
def fractional(evals):
    ev1 = evals[..., 0]
    ev2 = evals[..., 1]
    ev3 = evals[..., 2]
    
    f = np.sqrt(0.5 * ((ev1 - ev2)**2 + (ev2 - ev3)**2 + (ev3 - ev1)**2)
                      / (ev1*ev1 + ev2*ev2 + ev3*ev3))                      
    return f
    
def mode_of_anisotropy(tenfit):
    l1 = tenfit.evals[:,:,:,0]
    l2 = tenfit.evals[:,:,:,1]
    l3 = tenfit.evals[:,:,:,2]
    k1 = l1+l2+l3
    m1 = k1/3
    k2 = np.sqrt((l1-m1)**2 + (l2 - m1)**2 + (l3-m1)**2)
    m2 = (k2**2)/3
    k3 = (l1*l2*l3)/(k2**3)