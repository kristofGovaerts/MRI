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
from imageTypes import *
from datetime import datetime
import glob
import Tkinter, tkFileDialog
root = Tkinter.Tk()
import warnings
import subprocess
import numpy as np
import gzip

#global variables. change these to change the variables in the program
FILENAME=0 #shift-corrected image filenames automatically generated if 0
REF=0  #reference number of scan to register other images to
PROTOCOL='eddy_correct' #can be eddy_v, or FSL's eddy or eddy_correct. eddy_v for speed, eddy_correct for more complicated registration. 
SRC_FILENAME=0 #change if you want to edit the filename of the .src imgs
DSI_STUDIO_PROCESS=0 #Whether or not to process this data with DSI studio. Can also be done manually, takes very little time.
REMOVE_A0=2

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
    print "Scan types found: ", list_scans()
    
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
            #eddy current correction step
            if os.path.isfile(dtimg.name+'_EC.nii.gz') or os.path.isfile(dtimg.name+'_EC.nii'):
                print "Found EC file. Not running new eddy current correction."
                if os.path.isfile(dtimg.name+'_EC.nii.gz'):
                    dtimg.pdata=nib.load(dtimg.name+'_EC.nii.gz').get_data()
                elif os.path.isfile(dtimg.name+'_EC.nii'):
                    dtimg.pdata=nib.load(dtimg.name+'_EC.nii').get_data()
            else:
                try:
                    if dtimg.nrep==1:
                        from PythonDiffusion import rescaleImage
                        print "Applying slope correction..."
                        rescaleImage(dtimg)
                    dtimg.eddyCorrection(FILENAME, REF, PROTOCOL)
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