# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 18:27:14 2014

@author: u0091609

Kristof Govaerts
Batch script - Convert DTI data to DSI studio-compatible .src data
"""
import numpy as np
import scipy.io
import os
from datetime import datetime
import glob
import Tkinter, tkFileDialog
root = Tkinter.Tk()
import warnings
import subprocess
import numpy as np
import gzip
import nibabel as nib

#global variables. change these to change the variables in the program
FILENAME=0 #shift-corrected image filenames automatically generated if 0
REF=0  #reference number of scan to register other images to
PROTOCOL='eddy_correct' #can be eddy_v, or FSL's eddy or eddy_correct. eddy_v for speed, eddy_correct for more complicated registration. 
SRC_FILENAME=0 #change if you want to edit the filename of the .src imgs
DSI_STUDIO_PROCESS=0 #Whether or not to process this data with DSI studio. Can also be done manually, takes very little time.
REMOVE_A0=2

def checkFileType(file):
    '''Checks the current folder to see what extension the input file has.
    Returns this extension so that it can be easily appended to your filename.

    Inputs:
        file: A string without extension.
    Outputs:
        ext: A string containing the file's extension.'''
    files=glob.glob(file+'.*')
    if file + '.img' in files:
        return '.img'
    elif file + '.nii.gz' in files:
        return '.nii.gz'
    elif file + '.nii' in files:
        return '.nii'
		
def getResolution(f):
    '''Returns x,y,z resolution in mm.'''
    res=list_values(read_line('VisuCoreResolution=',f))
    if len(res)==0:
        try:
            res=list_values(read_line('VisuCoreResolutionMS=',f))
        except AttributeError:
            pass
    return res
	
def correctSlope(img, filename, pr=True):
    st="VisuCoreDataSlope="
    cimg=np.array(img)
    slopes=list_values(read_line(st, filename))
    if len(slopes)==img.shape[-1]:
        if pr:
            print "Correcting slopes..."
        for i, slope in enumerate(slopes):
            cimg[:,:,:,i]*=slope
    return cimg
    
def vector_array(list, m):
    '''Reshapes list of n*m values into an n x m array of vectors.'''
    va=np.array(list).reshape(len(list)/m, m)
    return va
	
def read_line(str, file):
    '''Reads file.txt and looks for str.
    If str is in a line, returns rest of line.'''
    try:
        f = open(file + '.txt')
        for l in f:
            if  str in l:
                value = l.replace(str, '')
                return value[:-1] #to remove /n char
        f.close()
    except IOError:
        pass

def list_values(str):
    '''Puts values between [] in list< '''
    li=str.replace('[','').replace(']','').split()
    return [float(x) for x in li]
    
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

class Br2AInfo:
    def __init__(self, name):
        self.name=name
        self.filetype=checkFileType(self.name)
        self.protocol=read_line("VisuAcquisitionProtocol=",self.name)
        self.date=read_line("Subject_date=",self.name)
        try:
            self.nrep=int(list_values(read_line("NRepetitions=",self.name))[0])
        except (IndexError, AttributeError) as e:
            self.nrep=1    
        try:
            self.te=[list_values(read_line("EffectiveEchoTime=",self.name))]
        except AttributeError:
            self.te=[None]
        if len(self.te)==1:
            self.te=self.te[0]
        self.tr=list_values(read_line("RepetitionTime=",self.name))
        if len(self.tr)==1:
            self.tr=self.tr[0]
        self.flipangle=list_values(read_line("FlipAngle=",self.name))
		
class Bruker2AnalyzeImg(Br2AInfo):
    '''Class for handling converted Bruker images. It loads various parameters
    for easy access. Important parameters:
        self.data = original data array
        self.pdata = duplicate array, to use for processing without losing original data
        self.img = Nibabel image class
    Also consider daughter classes for specialized images:
        DiffusionImg(Bruker2AnalyzeImg)'''
    def __init__(self, name, repetition_avg=True):
        #metadata
        Br2AInfo.__init__(self, name)
        self.img=nib.load(self.name+self.filetype)
        #arrays
        try:
            self.data=self.img.get_data()
        except MemoryError:
            self.data=None
            warnings.warn('MemoryError. Image array too large for available memory.')
        self.affine=self.img.get_affine()
        #geometric
        self.resolution=getResolution(self.name)
        self.dim=int(list_values(read_line("VisuCoreDim=",self.name))[0])
#        if self.dim==2:
#            self.data=scipy.swapaxes(self.data,2,3)
        self.shape=self.data.shape
        self.position=list_values(read_line("VisuCorePosition=",self.name))
        self.zposition=self.position[2]
        if self.dim==3:
            self.zthickness=self.resolution[-1]
        elif self.dim==2:
            self.zthickness=list_values(read_line("VisuCoreFrameThickness=",self.name))[0]
        self.slopes=list_values(read_line("VisuCoreDataSlope=",self.name))
        self.slicegap=list_values(read_line("SliceGap=",self.name))[0]
        self.slicedistance=list_values(read_line("SliceDistance=",self.name))[0]
        #processing output
        self.pdata=np.array(self.data)
        if repetition_avg:
            if self.nrep>1:
                print "Multiple repetitions. Correcting slope first."
                rescaleImage(self)
                print "Averaging across %r repetitions..." %self.nrep
                interv=self.pdata.shape[3]/self.nrep
                nz=np.zeros([self.shape[0],self.shape[1],self.shape[2],interv])
                for i in range(self.nrep):
                    nz+=self.pdata[:,:,:,i*interv:(i+1)*interv]
                self.pdata=nz/self.nrep
                self.shape=self.pdata.shape

    def listPars(self):
        print """
        Image %s is a %r-D image with a resolution of %r mm. The protocol used
        is %s. This image was recorded at %r mm from the center of the magnet. The slice
        thickness or resolution in the z-dimension is %r mm.
        """ %(self.name+self.filetype, int(self.dim), self.resolution, self.protocol, self.zposition,self.zthickness)

    def resampleImage(self,refimg):
        nzoom=refimg.resolution[:]
        if refimg.dim==2:
            nzoom.append(refimg.zthickness)
        self.data,self.affine=resample(self.data, self.affine,self.resolution, nzoom)
        self.zthickness=nzoom[-1]
        self.resolution=nzoom

    def correctSlope(self):
        if len(self.slopes)==1:
            pass
        else:
            if self.dim==3:
                for i in range(self.shape[3]):
                    self.pdata[:,:,:,i]=self.pdata[:,:,:,i]*self.slopes[i]

            elif self.dim==2:
                rslopes=np.reshape(self.slopes,self.shape[2:])
                for j in range(self.shape[2]):
                    for i in range(self.shape[3]):
                        self.pdata[:,:,j,i]=self.pdata[:,:,j,i]
                        #self.pdata[:,:,j,i]=self.pdata[:,:,j,i]*rslopes[j,i]
		
class DiffusionImg(Bruker2AnalyzeImg):
    '''Bruker2AnalyzeImg subclass. Incorporates various diffusion parameters
    necessary for diffusion processing which are not included in the parent
    class to reduce load.

    Does not throw an error when you try to load a non-diffusion image,
    but will tell you.'''
    def __init__(self, name):
        Bruker2AnalyzeImg.__init__(self,name)
        print '%s-D array with shape %s' %(str(self.dim), str(self.pdata.shape))
        getDiffusionPars(self)      
        
    def save_dsi_btable(self, filename=None):
        if filename==None:
            filename=self.name+'_btable.txt'
        else:
            filename+='.txt'
        btab=dsi_btable(self)
        with open(filename,'w') as f:
            f.write(btab)
        print "b-table saved to", filename
        
    def save_mask(self, filename=None):
        if filename==None:
            filename=self.name+'_bmask.nii'
        else:
            filename+='.nii'
        mimg=nib.Nifti1Image(self.mask, self.affine, header=self.img.get_header())
        mimg.to_filename(filename)   
        
    def save_mrtrix_grads(self, filename=None):
        if filename==None:
            filename=self.name+'_grads.txt'
        else:
            filename+='.txt'
        btab=np.zeros((self.gtab.bvecs.shape[0],4))
        btab[:,:3]=self.gtab.bvecs
        btab[:,2]=-btab[:,2] #for MRtrix
        btab[:,3]=self.gtab.bvals
        np.savetxt(filename, btab)
        print "Gradients saved to", filename

class dsi_fib:
    '''Container class for .fib files.'''
    def __init__(self, filename):
        self.vdict=scipy.io.loadmat(gzip.open(filename))
        
    def fa(self):
        return np.reshape(self.vdict['fa0'], tuple(self.vdict['dimension'][0]), order='F')
    
    def adc(self):
        return np.reshape(self.vdict['adc'], tuple(self.vdict['dimension'][0]), order='F')
        
    def axial_dif(self):
        return np.reshape(self.vdict['axial_dif'], tuple(self.vdict['dimension'][0]), order='F')

    def radial_dif(self):
        return np.reshape(self.vdict['radial_dif'], tuple(self.vdict['dimension'][0]), order='F')
        
    def masked_tract_pars(self, tract_filename, thresh=10):
        '''Creates a mask for each of the diffusion images and stores these in a dictionary.'''
        tract=nib.load(tract_filename).get_data()[::-1,::-1,:]
        tract[tract<thresh]=0
        tract[tract>=thresh]=1
        masked_diffusion_params={"fa":np.array(self.fa()),
                          "adc":np.array(self.adc()),
                          "axial_dif":np.array(self.axial_dif()),
                          "radial_dif":np.array(self.radial_dif()),
                          "volume_pixels":len(tract[tract>0]),
                          "volume_ml":len(tract[tract>0])*np.product(self.vdict["voxel_size"])}
        for item in masked_diffusion_params.keys():
            masked_diffusion_params[item]*=tract
        return masked_diffusion_params
        
    def tract_params(self, tract_filename, thresh=10):
        pdict={}
        masked=self.masked_tract_pars(tract_filename, thresh)
        for key in masked.keys():
            area=masked[key]
            pdict.update({key:np.mean(area[area>0])})
        return pdict


def save_dsi_src(self, filename=0, remove_a0=0):
    '''DSI src file consists of a b-table, the geometrical image dimensions,
    the voxel size, and a number of image volumes equal to nA0+ndir*nbvals.
    
    Input is a DiffusionImg instance.'''
    if filename == 0:
        filename=self.name+'_dsi_mat.src'
    if remove_a0 > 0:
        if remove_a0 > self.nA0-2:
            remove_a0=self.nA0-2 #keep 2 A0 images at least. Also a check that you don't accidentally start removing diffusion weighted imgs
        print "Removing %s A0 images." %remove_a0
        self.pdata=self.pdata[:,:,:,remove_a0:]
        self.shape=self.pdata.shape
        self.nA0-=remove_a0
        self.avbvals=self.avbvals[remove_a0:]
        self.bvals=self.bvals[remove_a0:]
        
    b_table=dsi_btable(self)
    dimension=np.array(self.shape[:-1]) #only geometrical dims
    dimension=np.reshape(dimension,(1,3)) #needs to be in this order
    voxel_size=self.resolution
    if len(voxel_size)<3:
        voxel_size.append(self.zthickness)
    voxel_size=np.reshape(np.array(voxel_size),(1,3))
    vdict={'b_table':b_table,
           'dimension':dimension,
           'voxel_size':voxel_size} #dict of variables to save
           
    #rescale to int16. Have to check if there is a more elegant way to keep numbers positive
    div=np.max(self.pdata)/32768
    self.pdata/=div
           
    for i in range(self.shape[-1]):
        print "Saving volume %s of %s." %(i,self.shape[-1])
        vol=np.array(self.pdata[:,:,:,i],dtype=np.int16)
        if self.dim==3:
            vol=vol[::-1,::-1,:] #flips image appropriately. Such a stupidly trivial thing can impact results as diffusion dirs will no longer be correct
        elif self.dim==2:
            vol=vol[::-1,::-1,:] #flips image appropriately. Such a stupidly trivial thing can impact results as diffusion dirs will no longer be correct
        vol=np.reshape(vol,np.prod(dimension),order='F')
        vol[vol<0]=0
        
        if i == 0:
            print "Saving A0 image for visualization purposes."
            nib.save(nib.Nifti1Image(np.reshape(vol, self.shape[:-1], order='F'), self.affine), self.name+"_A0.nii")        
        
        name='image%s'%i
        vdict.update({name:vol})
        
    print "Saving DSI .src image as %s." %filename
    scipy.io.savemat(filename, vdict,appendmat=False,format='4')
    del vdict
    
def dsi_btable(self):
    '''Returns a b-table, compatible with DSI studio.
    b-table is formatted as n x 4, with n being the t-dimension of the image.'''    
    btable=np.zeros((self.shape[-1],4))
    btable[:,0]=self.bvals
    for i in range(self.ndirs):
        btable[self.nA0+self.nbvals*i:self.nA0+self.nbvals*(i+1),1:]=self.dwdir[i] #fills b-vector matrix with the different diffusion dirs
    btable=np.swapaxes(btable,0,1)
    return btable
    
def main():
    print "Enter the scan directory where the diffusion-weighted images are stored."
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
    print "Found %r Analyze images in directory." %len(glob.glob('*.img'))
    
    filelist = [x.replace('.img', '') for x in glob.glob('*.img')] #cool list comprehension that gets all files
    
    for ind, filen in enumerate(filelist):
        startTime = datetime.now()
        print "\nFile %r of %r. Filename: %s" %(ind +1, len(filelist), filen)
        try:
            im=Br2AInfo(filen)
        except AttributeError:
            print "Error. No text file found."
            continue
        except IOError:
            warnings.warn("Cannot find .hdr or .img file.")
            continue
    
        print "Scan type: ", im.protocol
        
        if ("dti" in im.protocol.lower() or "dki" in im.protocol.lower() or "diff" in im.protocol.lower()) and im.name[-1]=='1':
            dtimg=DiffusionImg(filen)
            if not (len(dtimg.bvals)>=1 and len(dtimg.dwdir)>1):
                print "Diffusion image found, but insufficient diffusion dirs or b-vals."
                continue          
		try:
			if dtimg.nrep==1:
				print "Applying slope correction..."
				rescaleImage(dtimg)
		except ValueError:
			print "Matrices not aligned."
                    
            #saving output
            save_dsi_src(dtimg, SRC_FILENAME, REMOVE_A0)
            
            #performing command line stuff     
            if DSI_STUDIO_PROCESS:
                command=""    
                subprocess.call(command, shell=True)
            del dtimg
        
        else:
            print "Can't process this image, sorry."
            
if __name__ == '__main__':
    '''main loop. Does not activate upon importing this module.'''
    main()