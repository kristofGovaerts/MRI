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
from datetime import datetime
import os
import subprocess
#from skimage.restoration import nl_means_denoising #may need to install scikit-image

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

class DiffusionImg(Bruker2AnalyzeImg):
    '''Bruker2AnalyzeImg subclass. Incorporates various diffusion parameters
    necessary for diffusion processing which are not included in the parent
    class to reduce load.

    Does not throw an error when you try to load a non-diffusion image,
    but will tell you.'''
    def __init__(self, name):
        from PythonDiffusion import getDiffusionPars
        Bruker2AnalyzeImg.__init__(self,name)
        print '%s-D array with shape %s' %(str(self.dim), str(self.pdata.shape))
        getDiffusionPars(self)

    def eddyCorrection(self,filename=0, ref=0,protocol='eddy_v'):
        from PythonDiffusion import eddyCorrection
        starttime=datetime.now()
        print "Applying eddy current correction."
        if filename==0:
            filename=self.name+'_EC'
        eddyCorrection(self, filename=filename, ref=ref, protocol=protocol)
        time=datetime.now()-starttime
        print "Eddy current correction completed in %r seconds." %time.seconds

    def tensorFit(self, bv=None, removea0=0, m="dti", mask=True):
        from PythonDiffusion import tensorFit
        tensorFit(self, bv=bv, removea0=0, m=m, mask=mask)

    def processDiffusion(self,ec=False, bv=None, mode='dti', removea0=0, mask=True):
        print "Processing diffusion, full pipeline."
        from PythonDiffusion import rescaleImage, eddyCorrection
        from datetime import datetime
        if self.nrep==1:
            print "Applying slope correction..."
            rescaleImage(self)
        if ec:
            print "ec=True. Applying eddy current correction."
            eddyCorrection(self,self.name+'_EC', protocol='eddy_correct')
        if mode == 'dti' or mode == 'RESTORE':
            self.tensorFit(bv=bv, removea0=removea0, m=mode, mask=mask)
        elif mode == 'dki':
            self.processKurtosis()
            
    def trace_1b(self, sc=True):
        '''Calculates diffusion trace. ln(S/S0)=exp(-b*ADC)'''
        from PythonDiffusion import rescaleImage
        if sc:
            print "Applying slope correction..."
            rescaleImage(self)
        shape = list(self.shape)
        shape[-1] = self.ndirs+2 #ndirs+2 volumes: A0, trace_average, and ndirs*traces
        trace = np.zeros(shape)
        a0 = np.average(self.pdata[:,:,:,:self.nA0], axis=-1) #take mean of all A0s
        for i in range(self.ndirs):
            dwi = self.pdata[:,:,:,self.nA0 + i]
            trace[:,:,:,2 +i] = -np.log(dwi/a0)/self.bvals[self.nA0+i]
        trace[:,:,:,0] = a0 #to draw rois
        trace[:,:,:,1] = np.average(trace[:,:,:,2:], axis=-1)
        self.trace = trace
            
    def trace_multib(self, sc=True, mask=True):
        from PythonDiffusion import rescaleImage
        from ImProcessing import T2regression_basic
        if sc:
            print "Applying slope correction..."
            rescaleImage(self)
        shape = list(self.shape)
        shape[-1] = self.ndirs+3 #ndirs+3 volumes: A0, trace_total, trace_average, and ndirs*traces
        trace = np.zeros(shape)
        a0 = np.average(self.pdata[:,:,:,:self.nA0], axis=-1) #take mean of all A0s
        #fit all directions at once
        self.adc, self.stdevmap, self.si0map, self.adcfig=T2regression_basic(self.pdata, np.array(self.bvals), mask=mask)
        trace[:,:,:,1] = 1.0/self.adc
        #and separate dirs
        for i in range(self.ndirs):
            print "dir", i+1
            bv = list(self.bvals[:self.nA0]) + list(self.bvals[self.nA0+i*self.ndirs : self.nA0+(i+1)*self.nbvals])
            si = np.concatenate((self.pdata[:,:,:,:self.nA0], self.pdata[:,:,:,self.nA0+i*self.ndirs : self.nA0+(i+1)*self.nbvals]), axis=-1)
            self.adcn, self.stdevmapn, self.si0mapn, self.adcfign=T2regression_basic(si, np.array(bv), mask=mask)
            trace[:,:,:,3+i] = 1.0/self.adcn
        trace[:,:,:,2] = np.average(trace[:,:,:,3:], axis=-1)
        trace[:,:,:,0] = a0
        self.trace = trace
            
    def tensorFitAllBvs(self, removea0=0):
        for b in range(self.nbvals):
            bval = self.avbvals[self.nA0 + b]
            print "Current b-value:", str(bval)
            self.tensorFit(bv=[bval], removea0=0)
            self.compose_image()
            self.save_output(self.name+'_EC_b' + str(int(bval)))
        
    def processKurtosis(self):
        from ImProcessing import dki_fit, calc_diffusionpars
        print "Processing diffusion kurtosis."
        time=datetime.now()
        adcmap, akcmap, si0map, adcstdevmap, akcstdevmap=dki_fit(self)
        elapsed=datetime.now()-time
        print "Kurtosis processing took %r seconds. Saving..." %elapsed.seconds
        np.save(self.name+'_adcmap', adcmap)
        np.save(self.name+'_akcmap', akcmap)
        np.save(self.name+'_si0map', si0map)
        np.save(self.name+'_adcstdevmap', adcstdevmap)
        np.save(self.name+'_akcstdevmap', akcstdevmap)
        self.madc, self.fmap=calc_diffusionpars(adcmap)
        self.makc, self.fmap=calc_diffusionpars(akcmap)
        self.madcstd, self.fmap=calc_diffusionpars(adcstdevmap)
        self.makcstd, self.fmap=calc_diffusionpars(akcstdevmap)
        self.msi0, self.fmap=calc_diffusionpars(si0map)

    def compose_image(self, mode='dti'):
        '''Creates a 4D matrix containing the various intermediate outputs from the
        CBF processing. The 4th dimension contains the outputs:
            [:,:,:,0]=FA map
            [:,:,:,1]=MD map
            [:,:,:,2]=AD map
            [:,:,:,3]=RD map
            [:,:,:,4]=Mode of anisotropy map (-1 = linear, +1 = flat)
            [:,:,:,5]=B0 image'''

        print "Composing image."
        shape=list(self.shape[:3])
        if mode =='dti':
            shape.append(6)
            mat=np.zeros(shape)
            li=[self.tenfit.fa, self.tenfit.md, self.tenfit.ad, self.tenfit.rd, self.tenfit.mode, self.pdata[:,:,:,0]]
            for i in range(6):
                mat[:,:,:,i]=li[i]
        elif mode == 'dki':
            shape.append(6)
            mat=np.zeros(shape)
            li=[self.madc, self.makc, self.msi0, self.madcstd, self.makcstd, self.fmap]
            for i in range(6):
                mat[:,:,:,i]=li[i]
        self.outdata=mat

    def save_output(self, filename, savebm=True):
        hdr=self.img.get_header()
        try:
            nimg=nib.Nifti1Image(self.outdata, self.affine, header=hdr)
        except AttributeError:
            self.compose_image()
        hdr.set_data_shape(self.outdata.shape)
        nimg=nib.Nifti1Image(self.outdata, self.affine, header=hdr)
        nimg.set_data_dtype('float32')
        nimg.to_filename(filename)            
        
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
        
    def processDiffusionTrace(self):
        from PythonDiffusion import processDiffusionTrace
        processDiffusionTrace(self)

class T2Img(Bruker2AnalyzeImg):
    """T2 images that can be processed."""
    def __init__(self, name, repetition_avg=True, nlm_denoise=False, slopecorr=True, reshape=True):
        Bruker2AnalyzeImg.__init__(self,name, repetition_avg)
        if reshape: #necessary for PV5 data
            print "Reshaping self.pdata."
            self.pdata=np.swapaxes(self.pdata, 2, 3)
            self.pdata=np.reshape(self.pdata, self.shape)
        if slopecorr:
            print "Correcting for slope..."
            self.correctSlope()
        if nlm_denoise:
            print "Denoising image using nonlocal means."
            for i in range(self.shape[-1]):
                for j in range(self.shape[-2]): #for 2D scans
                    print "Volume %s, slice %s." %(i+1,j+1)
                    self.pdata[:,:,j,i] = nl_means_denoising(self.pdata[:,:,j,i], patch_size=7, patch_distance=11, h=25000, multichannel=False)
                

    def processImage(self, mask=True):
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
        from ImProcessing import T2regression            
        self.t2map, self.stdevmap, self.si0map, self.soffmap, self.t2fig=T2regression(self.pdata, np.array(self.te), mask)

    def processImage_basic(self, mask=True):
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
        from ImProcessing import T2regression_basic
        self.t2map_basic, self.stdevmap_basic, self.si0map_basic, self.t2fig_basic=T2regression_basic(self.pdata, np.array(self.te), mask)

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

    def compose_image_basic(self, filter_im=1e3):
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
        shape.append(5)
        mat=np.zeros(shape)
        li=[self.t2map_basic, self.si0map_basic, self.stdevmap_basic, self.paraprocdata, self.pdata[:,:,:,0]]
        for i in range(5):
            mat[:,:,:,i]=li[i]
            mat[:,:,:,i][li[2]>filter_im]=0 #If the std dev of the fit is larger than 1000ms its probably wrong
        self.outdata_basic=mat
        
    def get_te1(self, mask=True):
        '''Gets the image at the first echo time. Also masks if necessary.'''
        im = self.outdata_basic
        d = im[:,:,:,-1]
        if mask:
            m = im[:,:,:,0]
            m[m>0]=1
            m[m<0]=0
            d *= m
        return d

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
        
class T1Img(Bruker2AnalyzeImg):
    """T2 images that can be processed."""
    def __init__(self, name, repetition_avg=True):
        Bruker2AnalyzeImg.__init__(self,name, repetition_avg)

    def processImage(self, mask=True):
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
        from ImProcessing import T2regression,t1r_vtr
        self.t2map, self.stdevmap, self.si0map, self.soffmap, self.t2fig=T2regression(self.pdata, np.array(self.te), mask, mode=t1r_vtr)

    def processImage_basic(self, mask=True):
        '''Runs the full T2 processing pipeline.
        Inputs:
            self.pdata = image array
            self.te = echo times
        Outputs:
            self.t2map = Calculated T1 values
            self.stdevmap = Standard deviation of the fit
            self.si0map = Calculated SI0 values
            self.soffmap = Calculated noise values
            self.t2fig = figure for saving later'''
        from ImProcessing import T2regression_basic,t1r_vtr_basic
        self.t2map_basic, self.stdevmap_basic, self.si0map_basic, self.t2fig_basic=T2regression_basic(self.pdata, np.array(self.te), mask, mode=t1r_vtr_basic)

    def compose_image(self):
        '''Creates a 4D matrix containing the various intermediate outputs from the
        CBF processing. The 4th dimension contains the outputs:
            [:,:,:,0]=T1 map
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

    def compose_image_basic(self, filter_im=1e3):
        '''Creates a 4D matrix containing the various intermediate outputs from the
        CBF processing. The 4th dimension contains the outputs:
        [:,:,:,0]=T1 map
        [:,:,:,1]=SI0 map
        [:,:,:,2]=Standard deviation of selective T1 fit
        [:,:,:,3]=T1 data calculated with Paravision (zero if not present)
        [:,:,:,4]=Original dataset at TE0'''
        get_parafile(self,2,proc=2)
        print "Composing image."
        shape=list(self.shape[:3])
        shape.append(5)
        mat=np.zeros(shape)
        li=[self.t2map_basic, self.si0map_basic, self.stdevmap_basic, self.paraprocdata, self.pdata[:,:,:,0]]
        for i in range(5):
            mat[:,:,:,i]=li[i]
            mat[:,:,:,i][li[2]>filter_im]=0 #If the std dev of the fit is larger than 1000ms its probably wrong
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
        print self.sel.shape, self.nsel.shape

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

    def compose_image(self, filter_im=5e3):
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
            mat[:,:,:,i][li[2]>filter_im]=0 #If the std dev of the fit is larger than 5000ms its probably wrong
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
        para=Bruker2AnalyzeImg(Br2AInfo(self.name[:-1]+str(num)))
        if self.dim==2:
            shape=para.shape
            para.pdata=np.swapaxes(para.data, 2, 3)
            para.pdata=np.reshape(para.pdata, shape)
        para.correctSlope()
        self.paraprocdata=para.pdata[:,:,:,proc]
    except (TypeError, AttributeError) as e:
        print "No Paravision processed file."
        self.paraprocdata=np.zeros(self.shape[:3])
        
def dke_process(self):
    with open(self.name+'_bvecs.dat', 'w') as f:
        f.write(diffusion_dirs_str(self))
    with open(self.name+'_dke_params.dat', 'w') as f:
        f.write(dke_paramfile(self))
    data=self.pdata[:,:,:,self.nA0-1:]
    a0s=np.mean(self.pdata[:,:,:,:self.nA0],axis=3)
    data[:,:,:,0]=a0s
    data=data[::-1,::-1,:]
    
    temp='dke_temp.nii.gz'
    nifti=nib.Nifti1Image(data, self.affine)
    nib.save(nifti, temp)
    command=["dke " + self.name+'_dke_params.dat']
    subprocess.call(command, shell=True)
    #os.remove('dke_temp.nii.gz')


def dke_paramfile(self, outdir=0, name=0, form='nifti', kmin=0, kmax=3, kmin_f=0, kmax_f=3, brainmask=1, robust=1, noise_tol=0.9):
    '''For details of parameters, refer to: http://www.nitrc.org/docman/view.php/652/1263/DKE%20User's%20Guide '''
    if outdir==0:
        outdir=os.getcwd()+'\\'
    date='Scan date:'+self.date    
    if name==0:
        name='dke_temp.nii.gz'
    bv='[0 '
    for i in range(self.nbvals):
        bv+=str(int(round(self.avbvals[self.nA0+i])))+' '
    bv+=']'
    b0_thresh=2.5*np.mean(self.pdata[:5,:5,:,0]) #brain mask
    gvecs=os.getcwd()+'\\'+self.name+'_bvecs.dat'
    idxes=''
    for i in range(self.nbvals):
        idxes+='idx_gradients{%s} = [];\n'%str(i+1)
    dinds=str(range(1,self.ndirs+1))
    dinds=dinds.replace(',','')
    idxes=idxes.replace('[]',dinds)
    params = '''%% %s
studydir = '%s';
subject_list = {''};
preprocess_options.format = '%s';
preprocess_options.fn_nii = '%s';
fn_img_prefix = 'rdki';
bval = %s;
ndir = %s;
idx_1st_img = 1;
Kmin  = %s;
NKmax = %s;
Kmin_final = %s;
Kmax_final = %s;
T = %s;
find_brain_mask_flag = %s;
dki_method.no_tensor = 0;
dki_method.linear_weighting = 1;
dki_method.linear_constrained = 1;
dki_method.nonlinear = 0;
dki_method.linear_violations = 0;
dki_method.robust_option = %s;
dki_method.noise_tolerance = %s;
dti_method.dti_flag = 0;
dti_method.dti_only = 0;
dti_method.no_tensor = 0;
dti_method.linear_weighting = 1;
dti_method.b_value = %s;
dti_method.directions{1} = %s;
dti_method.robust_option = 0;
dti_method.noise_tolerance = 0.09;
fn_noise = '';
fwhm_img   = 0;
fwhm_noise = 0;
median_filter_method = 1;
map_interpolation_method.flag = 0;
map_interpolation_method.order = 1;
map_interpolation_method.resolution =      1;
fn_gradients = '%s';\n%s''' %(date,outdir, form, name, bv, str(self.ndirs), str(kmin), str(kmax),str(kmin_f), str(kmax_f), str(b0_thresh),str(brainmask), str(robust), str(noise_tol), bv, dinds, gvecs, idxes)
    
    params=params.replace('\\','/')    
    return params
    
def dsi_btable(self):
    '''Returns a b-table, compatible with DSI studio.'''    
    btablestr='0 0 0 0\n'*self.nA0
    for i in range(len(self.avbvals)):
        for d in range(self.ndirs):
            b=self.avbvals[self.nA0+d]
            vec=str(self.dwdir[d])[1:-1]
            btablestr+='%s%s\n' %(b,vec)  
    return btablestr
    
def diffusion_dirs_str(self):
    '''Returns a delimited string with all diffusion directions including 1 b0.
    Accepted by many types of software.'''
    dirs=self.dwdir
    st=''
    for ddir in dirs:
        d=[str(i) for i in ddir]
        st+=' '.join(d) + '\n'
        
    return st