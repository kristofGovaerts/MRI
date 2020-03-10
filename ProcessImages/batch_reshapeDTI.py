# -*- coding: utf-8 -*-
"""
Created on Wed Oct 09 21:06:31 2013

@author: Gebruiker
"""
from ReadInput import *
import nibabel as nib
import numpy as np
import os
import Tkinter, tkFileDialog

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
                from PythonDiffusion import rescaleImage
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

print "Enter the scan directory."
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

filelist = glob.glob('*.hdr')
for f in filelist:
    i=Bruker2AnalyzeImg(f[:-4])
    rescaleImage(i)
    im = nib.AnalyzeImage(i.pdata, i.affine)
    nib.save(im, f[:-4] + '_RS.hdr')