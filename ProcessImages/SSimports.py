#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      s0203524
#
# Created:     28/08/2013
# Copyright:   (c) s0203524 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import os
from ReadInput import *
from ImProcessing import *
from Saveoutput import *
from imageTypes import *
import nibabel as nib
import pylab as pl
import numpy as np
from scipy.optimize import curve_fit
from datetime import datetime
import glob
import Tkinter, tkFileDialog
import warnings

def main():
    pass

if __name__ == '__main__':
    main()

import os, numpy as np, pylab as pl, nibabel as nib, math, glob, Image
import wx

def loadDir(self,str):
    dirname = ''
    f1 = wx.DirDialog(self, str, dirname)
    if f1.ShowModal() == wx.ID_OK:
        curdir = f1.GetPath()
        return curdir
    f1.Destroy()
    
def processAll(filelist, t1_cbf, t1bl, ti, t2noise, t2nonoise):
    '''Batch processing manager. Only works inside dedicated environment.
    Depends on a number of global variables.'''
    for ind, filen in enumerate(filelist):
        startTime = datetime.now()
        print "\nFile %r of %r. Filename: %s" %(ind +1, len(filelist), filen)
        try:
            im=Bruker2AnalyzeImg(filen)
        except AttributeError:
            print "Error. No text file found."
            continue
        except IOError:
            warnings.warn("Cannot find .hdr or .img file.")
            continue
    
        print "Scan type: ", im.protocol
    
        print "Echo times found: ", im.te
    
        print "Image shape: ", im.shape
    
        if (t2noise or t2nonoise) and len(im.te)>2 and filen[-2:] == '_1' and os.path.isdir(filen) is False:  #makes sure it is not a processed image, either by Paravision or by this script
            print "Can do T2 relaxometry on this image. Proceeding..."
            t2img=T2Img(filen)
            print "Correcting for slope..."
            t2img.correctSlope()
            if t2noise:
                t2img.processImage()
                if os.path.isdir(filen) is False:
                    os.mkdir(filen)
                np.save(os.path.join(filen, filen + '_t2map'), t2img.t2map)
                np.save(os.path.join(filen, filen + '_si0map'), t2img.si0map)
                np.save(os.path.join(filen, filen + '_soffmap'), t2img.soffmap)
                np.save(os.path.join(filen, filen + '_stdevmap'), t2img.stdevmap)
                t2nif=nib.Nifti1Image(t2img.t2map,t2img.affine)
                nib.save(t2nif, os.path.join(filen,filen+'_T2Map.nii'))
                t2img.t2fig.savefig(os.path.join(filen, filen + "_T2fit"))
                t2img.save_output(os.path.join(filen, filen + '_Full'))
            if t2nonoise:
                if t2noise:
                    print "\n" #for cleanness. Kind of messy
                t2img.processImage_basic()
                if os.path.isdir(filen) is False:
                    os.mkdir(filen)            
                np.save(os.path.join(filen, filen + '_t2map'), t2img.t2map_basic)
                np.save(os.path.join(filen, filen + '_si0map'), t2img.si0map_basic)
                np.save(os.path.join(filen, filen + '_stdevmap'), t2img.stdevmap_basic)
            #saveslices3d(t2img.t2map, [0.0,75.0], os.path.join(filen, filen + "_T2map")) #saves .png of all slices in one figure
            #saveslices3d(t2img.stdevmap, [0.0,10.0], os.path.join(filen, filen + "_Std Devs"))
                t2img.t2fig_basic.savefig(os.path.join(filen, filen + "_T2fit_basic"))
                t2img.save_output_basic(os.path.join(filen, filen + '_basic_Full'))
            pl.close()
        
        elif t1_cbf and 'fair' in im.protocol.lower() and filen[-2:] == '_1' and os.path.isdir(filen) is False: 
            print "Can do FAIR-RARE ASL processing on this image. Proceeding..."
            aslimg=ASLImg(filen)
            print "Correcting for slope..."
            aslimg.correctSlope()
            aslimg.ti=ti
            aslimg.t1bl = t1bl
            aslimg.processimage()
            if os.path.isdir(filen) is False:
                os.mkdir(filen)
            np.save(os.path.join(filen, filen + '_selt1'), aslimg.selt1)
            np.save(os.path.join(filen, filen + '_nselt1'), aslimg.nselt1)
            np.save(os.path.join(filen, filen + '_selstd'), aslimg.selstd)
            np.save(os.path.join(filen, filen + '_nselstd'), aslimg.nselstd)
            np.save(os.path.join(filen, filen + '_CBF'), aslimg.cbfmap)
            #saveslices3d(aslimg.cbfmap, [-100.0, 600.0], file + "_CBFmap")
            aslnif=nib.Nifti1Image(aslimg.cbfmap, aslimg.affine)
            nib.save(aslnif,os.path.join(filen, filen+'_CBFmap.nii'))
            aslimg.selt1fig.savefig(os.path.join(filen, filen + "_selT1fit"))
            aslimg.nselt1fig.savefig(os.path.join(filen, filen + "_nonselT1fit"))
            aslimg.save_output(os.path.join(filen, filen + '_Full'))
            pl.close()
        else:
            print "Can't process this image, sorry."
    
        timecomp = datetime.now()-startTime
    
        print "Process took %rs to complete." %timecomp.seconds
    print "Processed all images."