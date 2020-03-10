# -*- coding: utf-8 -*-
"""
Created on Wed Oct 09 21:06:31 2013

@author: Gebruiker
"""
from ReadInput import *
import nibabel as nib
import numpy as np
import warnings
import scipy.ndimage as ndimage

class Bruker2AnalyzeImg(object):
    '''Class for handling converted Bruker images. It loads various parameters
    for easy access. Important parameters:
        self.data = original data array
        self.pdata = duplicate array, to use for processing without losing original data
        self.img = Nibabel image class
    Also consider daughter classes for specialized images:
        DiffusionImg(Bruker2AnalyzeImg)'''
    def __init__(self, name):
        #metadata
        self.name=name
        self.filetype=checkFileType(self.name)
        self.img=nib.load(self.name+self.filetype)
        self.hdr=self.img.get_header()
        self.protocol=read_line("VisuAcquisitionProtocol=",self.name)
        self.date=read_line("Subject_date=",self.name)
        try:
            self.nrep=int(list_values(read_line("NRepetitions=",self.name))[0])
        except IndexError:
            self.nrep=1
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
            self.zthickness=self.resolution[2]
        elif self.dim==2:
            self.zthickness=list_values(read_line("VisuCoreFrameThickness=",self.name))[0]
        self.slopes=list_values(read_line("VisuCoreDataSlope=",self.name))
        self.slicegap=list_values(read_line("SliceGap=",self.name))[0]
        self.slicedistance=list_values(read_line("SliceDistance=",self.name))[0]
        #scan parameters
        self.te=[list_values(read_line("EffectiveEchoTime=",self.name))]
        if len(self.te)==1:
            self.te=self.te[0]
        self.tr=list_values(read_line("RepetitionTime=",self.name))
        if len(self.tr)==1:
            self.tr=self.tr[0]
        self.flipangle=list_values(read_line("FlipAngle=",self.name))
        #processing output
        self.pdata=np.array(self.data)
        if self.nrep>1:
            from PythonDiffusion import rescaleImageC:\Users\s0203524\Desktop\test
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
                        self.pdata[:,:,j,i]=self.pdata[:,:,j,i]*rslopes[j,i]

class DiffusionImg(Bruker2AnalyzeImg):
    '''Bruker2AnalyzeImg subclass. Incorporates various diffusion parameters
    necessary for diffusion processing which are not included in the parent
    class to reduce load.

    Does not throw an error when you try to load a non-diffusion image,
    but will tell you.'''
    def __init__(self, name):
        from PythonDiffusion import getDiffusionPars
        Bruker2AnalyzeImg.__init__(self,name)
        print '2D array with shape', self.pdata.shape
        getDiffusionPars(self)

    def eddyCorrection(self,filename, ref=0):
        from PythonDiffusion import eddyCorrection
        starttime=datetime.now()
        print "Applying eddy current correction."
        eddyCorrection(self, filename, ref)
        time=datetime.now()-starttime
        print "Eddy current correction completed in %r seconds." %time.seconds

    def tensorFit(self, bv=None):
        from PythonDiffusion import tensorFit
        tensorFit(self, bv=bv)

    def processDiffusion(self,ec=False, bv=None):
        print "Processing diffusion, full pipeline."
        from PythonDiffusion import rescaleImage, eddyCorrection
        from datetime import datetime
        if self.nrep==1:
            print "Applying slope correction..."
            rescaleImage(self)
        if ec:
            print "ec=True. Applying eddy current correction."
            eddyCorrection(self,self.name+'_EC', protocol='eddy_v')
        self.tensorFit(bv=bv)

    def compose_image(self):
        '''Creates a 4D matrix containing the various intermediate outputs from the
        CBF processing. The 4th dimension contains the outputs:
            [:,:,:,0]=FA map
            [:,:,:,1]=MD map
            [:,:,:,2]=AD map
            [:,:,:,3]=RD map'''

        print "Composing image."
        shape=list(self.shape[:3])
        shape.append(4)
        mat=np.zeros(shape)
        li=[self.tenfit.fa, self.tenfit.md, self.tenfit.ad, self.tenfit.rd]
        for i in range(4):
            mat[:,:,:,i]=li[i]
        self.outdata=mat

    def save_output(self, filename):
        hdr=self.img.get_header()
        try:
            nimg=nib.Nifti1Image(self.outdata, self.affine, header=hdr)
        except AttributeError:
            self.compose_image()
        hdr.set_data_shape(self.outdata.shape)
        nimg=nib.Nifti1Image(self.outdata, self.affine, header=hdr)
        nimg.set_data_dtype('float32')
        nimg.to_filename(filename)

class T2Img(Bruker2AnalyzeImg):
    """T2 images that can be processed."""
    def __init__(self, name):
        Bruker2AnalyzeImg.__init__(self,name)
        print "Reshaping self.pdata."
        self.pdata=np.swapaxes(self.pdata, 2, 3)
        self.pdata=np.reshape(self.pdata, self.shape)

    def processImage(self):
        '''Runs the full T2 processing pipeline.
        Inputs:
            self.pdata = image array
            self.te = echo times
        Outputs:
            self.t2map = Calculated T2 values
            self.stdevmap = Standard deviation of the fit
            self.si0map = Calculated SI0 values
            self.soffmap = Calculated noise values
            self.t2fig = figure for saving later'''
        from ImProcessing import T2regression,t2r
        self.t2map, self.stdevmap, self.si0map, self.soffmap, self.t2fig=T2regression(self.pdata, np.array(self.te))

    def processImage_basic(self):
        '''Runs the full T2 processing pipeline.
        Inputs:
            self.pdata = image array
            self.te = echo times
        Outputs:
            self.t2map = Calculated T2 values
            self.stdevmap = Standard deviation of the fit
            self.si0map = Calculated SI0 values
            self.soffmap = Calculated noise values
            self.t2fig = figure for saving later'''
        from ImProcessing import T2regression_basic,t2r_basic
        self.t2map_basic, self.stdevmap_basic, self.si0map_basic, self.t2fig_basic=T2regression_basic(self.pdata, np.array(self.te))

    def compose_image(self):
        '''Creates a 4D matrix containing the various intermediate outputs from the
        CBF processing. The 4th dimension contains the outputs:
            [:,:,:,0]=T2 map
            [:,:,:,1]=SI0 map
            [:,:,:,2]=noise map
            [:,:,:,3]=Standard deviation of selective T1 fit
            [:,:,:,4]=T2 data calculated with Paravision (zero if not present)
            [:,:,:,5]=Original dataset at TE0'''

        get_parafile(self,2,proc=2)
        print "Composing image."
        shape=list(self.shape[:3])
        shape.append(6)
        mat=np.zeros(shape)
        li=[self.t2map, self.si0map, self.soffmap, self.stdevmap, self.paraprocdata, self.pdata[:,:,:,0]]
        for i in range(6):
            mat[:,:,:,i]=li[i]
        self.outdata=mat

    def compose_image_basic(self):
        '''Creates a 4D matrix containing the various intermediate outputs from the
        CBF processing. The 4th dimension contains the outputs:
        [:,:,:,0]=T2 map
        [:,:,:,1]=SI0 map
        [:,:,:,2]=Standard deviation of selective T1 fit
        [:,:,:,3]=T2 data calculated with Paravision (zero if not present)
        [:,:,:,4]=Original dataset at TE0'''
        get_parafile(self,2,proc=2)
        print "Composing image."
        shape=list(self.shape[:3])
        shape.append(6)
        mat=np.zeros(shape)
        li=[self.t2map_basic, self.si0map_basic, self.stdevmap_basic, self.paraprocdata, self.pdata[:,:,:,0]]
        for i in range(5):
            mat[:,:,:,i]=li[i]
        self.outdata_basic=mat

    def save_output(self, filename):
        hdr=self.img.get_header()
        try:
            nimg=nib.Nifti1Image(self.outdata, self.affine, header=hdr)
        except AttributeError:
            self.compose_image()
        hdr.set_data_shape(self.outdata.shape)
        nimg=nib.Nifti1Image(self.outdata, self.affine, header=hdr)
        nimg.set_data_dtype('float32')
        nimg.to_filename(filename)

    def save_output_basic(self, filename):
        hdr=self.img.get_header()
        try:
            nimg=nib.Nifti1Image(self.outdata_basic, self.affine, header=hdr)
        except AttributeError:
            self.compose_image_basic()
        hdr.set_data_shape(self.outdata_basic.shape)
        nimg=nib.Nifti1Image(self.outdata_basic, self.affine, header=hdr)
        nimg.set_data_dtype('float32')
        nimg.to_filename(filename)

class ASLImg(Bruker2AnalyzeImg):
    """ASL images that can be processed."""
    def _init(self,name):
        Bruker2AnalyzeImg.__init__(self,name)

    def splitimage(self):
        '''Splits image into selective and nonselective parts. Necessary for
        further processing.'''
        from ImProcessing import ASLimg_split
        self.sel, self.nsel = ASLimg_split(self.pdata)

    def calcdeltaM(self, gauss=None):
        '''Calculates magnetization difference between selective and nonselective
        images. Gauss can be a 4-D tuple used for smoothing the data.'''
        self.deltaM = np.abs(self.sel-self.nsel)
        if gauss != None:
            self.deltaM=ndimage.gaussian_filter(self.deltaM, sigma=(gauss))

    def t1processing(self):
        from ImProcessing import T1Regression,t1r
        if len(self.sel) > 0 and len(self.nsel) > 0:
            print "\nProcessing selective inversion:"
            self.selt1, self.evsi0, self.selstd, self.selt1fig = T1Regression(self.sel, self.ti, self.sel.shape)
            print '\n\nProcessing nonselective inversion:'
            self.nselt1, self.oddsi0, self.nselstd, self.nselt1fig = T1Regression(self.nsel, self.ti, self.nsel.shape)
            print '\n'
        else:
            print "Please use self.splitimage() to split image first."

    def CBFcalc(self):
        '''Calculates CBF, incorporating the blood/tissue partition coefficient,
        T1 of blood and also corrects for global T1.'''
        self.cbfmap=np.zeros(self.selt1.shape)
        #if len(self.cbfmap.shape) == 2:
        iter=self.cbfmap.shape[0]*self.cbfmap.shape[1]
        print "\nCalculating CBF..."
        for px in range(0,iter-1):
            x1, y1=np.unravel_index(px, self.cbfmap.shape[:2])
            if bool(self.selt1[x1,y1] > 0 and self.nselt1[x1,y1] > 0):
                self.cbfmap[x1,y1] = ((60000*83*self.nselt1[x1,y1])/self.t1bl)*((1/self.selt1[x1,y1])-(1/self.nselt1[x1,y1]))[0] #otherwise it makes array objects instead of floats. 60000 is because CBF is in ml/mg/min and TI/T1 in ms.
            else:#ignore values below thresh
                self.cbfmap[x1,y1] = 0

    def processimage(self):
        self.splitimage()
        self.t1processing()
        self.CBFcalc()

    def compose_image(self):
        '''Creates a 4D matrix containing the various intermediate outputs from the
        CBF processing. The 4th dimension contains the outputs:
            [:,:,:,0]=CBF map
            [:,:,:,1]=Selective T1 map
            [:,:,:,2]=Standard deviation of selective T1 fit
            [:,:,:,3]=Nonselective T1 map
            [:,:,:,4]=Standard deviation of nonselective T1 fit
            [:,:,:,5]=Perfusion data calculated with Paravision (zero if not present)
            [:,:,:,6]=Original dataset at TI0'''
        get_parafile(self,4)
        print "Composing image."
        shape=list(self.shape[:3])
        shape.append(7)
        mat=np.zeros(shape)
        li=[self.cbfmap,self.selt1, self.selstd, self.nselt1, self.nselstd, self.paraprocdata, self.pdata[:,:,:,0]]
        for i in range(7):
            mat[:,:,:,i]=li[i]
        self.outdata=mat

    def save_output(self, filename):
        hdr=self.img.get_header()
        try:
            nimg=nib.Nifti1Image(self.outdata, self.affine, header=hdr)
        except AttributeError:
            self.compose_image()
        hdr.set_data_shape(self.outdata.shape)
        nimg=nib.Nifti1Image(self.outdata, self.affine, header=hdr)
        nimg.set_data_dtype('float32')
        nimg.to_filename(filename)

def get_parafile(self, num, proc=0):
    '''Gets paravision processed file with reco number 'num'.
    proc is location where paravision stores the actual T2 map in the fourth dimension.'''
    try:
        para=Bruker2AnalyzeImg(self.name[:-1]+str(num))
        if self.dim==2:
            shape=para.shape
            para.pdata=np.swapaxes(para.data, 2, 3)
            para.pdata=np.reshape(para.pdata, shape)
        para.correctSlope()
        self.paraprocdata=para.pdata[:,:,:,proc]
    except TypeError:
        print "No Paravision processed file."
        self.paraprocdata=np.zeros(self.shape[:3])
