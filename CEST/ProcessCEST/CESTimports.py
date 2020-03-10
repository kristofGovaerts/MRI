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

def main():
    pass

if __name__ == '__main__':
    main()

import os, numpy as np, pylab as pl, nibabel as nib, math, glob, Image
import wx
from ReadInput import *
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from scipy.optimize import curve_fit

def loadFile(self,str):
    dirname = ''
    f1 = wx.FileDialog(self, str, dirname, "", "*.*", wx.OPEN)
    if f1.ShowModal() == wx.ID_OK:
        l1f = f1.GetFilename()
        l1dir = f1.GetDirectory()
        l1=os.path.join(l1dir,l1f)
        return l1
    f1.Destroy()

def correctslope(img,file):
    '''Extracts the slope value from the .txt file and multiplies array
    by this value. First input is a numpy array, second is the file path.
    This function also thresholds by removing noise voxels (anything below
    7x the noise level).'''
    f,ext=os.path.splitext(file)
    slope=list_values(read_line("VisuCoreDataSlope=", f))
    sl=list(slope)
    for i in range(len(slope)):
        if i == 0:
            pass
        elif slope[i] == slope[i-1]:
            sl.remove(slope[i])  #removes duplicates from list
    if len(sl)>1:
        pass
    else:
        print "\nCorrecting image by slope ", sl[0]
        slimg=img*float(sl[0])
    noise=np.mean(slimg[:10,:10,0,0])
    slimg[slimg<noise*7]=noise*7
    return slimg

def plotGraph(img):
    figure=pl.figure(num=None, figsize=(1.2, 1), dpi=380, facecolor='w', edgecolor='k')#size fig at 420x420px to fit window
    figure.add_axes([0,0,1,1])#removes gray stuff from edges of figure
    pl.imshow(img)
    pl.axis('off')
    cb=pl.colorbar(shrink=.7, pad=0, format='%.0e') #format applies scientific notation. .0 = significant digits, e = exponential
    cb.outline.set_linewidth(0.2)
    cb.ax.tick_params(size=1,labelsize=2, pad=0)
    return figure, cb

def computeDifference(plusimg, negimg):
    '''Simple subtraction of images.'''
    return np.abs(negimg-plusimg)

def computePercentDiff(plusimg,negimg):
    negnozeros=np.array(negimg)
    negnozeros[negnozeros==0]=1
    percentdiff=100*np.abs((negimg-plusimg)/negnozeros)
    return percentdiff
    

def correctB0(arr,fr,offs,stepsize,il=True,B0=None):
    '''Corrects 2D or 3D image stack for B0 inhomogeneity.
    Inputs:
        arr=Array of 3D CEST images, with the t-axis being the frequency offset.
        fr=Array of frequencies matching the order of the CEST images.
        offs=Desired offset to quantify.
        stepsize=Stepsize taken between elements of arr.
        il=True if images were recorded in an interleaved fashion (e.g. +2000/-2000,+1975/-1975,...
        B0=Optional fieldmap. If not provided, algorithm extrapolates from the image array.
        Note that enough frequencies are necessary to construct a fieldmap.
        
    Outputs:
        offsimg: 4-D image containing positive and negative offset images 
        corrected for B0 inhomogeneity. Dimensions are (x,y,z,3), with t=0
        being the positive offset image and t=1 the negative. t=2 is the image
        at zero offset, corrected for B0.
        
        pseudofieldmap=Fieldmap calculated from the signal-to-frequency offset 
        curves based on local minima.
        
        indexerrors=x,y,z indices of pixels where the range of frequencies was insufficient.
        '''
    x,y,z,t=arr.shape
    offsimg=np.zeros([x,y,z,3])
    array=np.array(arr)
    freqs=np.array(fr)
    pseudofieldmap=np.zeros([x,y,z])
    indexerrors=[]
    if il==True:
        array=unInterleave(array)
        freqs=np.arange(np.min(freqs),np.max(freqs)+1,stepsize)
        pl.plot(freqs,array[30,30,0,:]) #check if the un-interleaving worked
    thresh=2*np.mean(array[0:5,0:5,:,0])
    if B0 == None:
        for ind in range(x*y*z):
            x1,y1,z1 = np.unravel_index(ind, array.shape[:-1])
            if array[x1,y1,z1,0] > thresh:
                yvals=array[x1,y1,z1,:]
                zero=np.argmin(yvals)
                pseudofieldmap[x1,y1,z1]=freqs[zero]
                offsimg[x1,y1,z1,2]=yvals[zero]
                try:    
                    offsimg[x1,y1,z1,0]=yvals[zero+int(offs/stepsize)]
                    offsimg[x1,y1,z1,1]=yvals[zero-int(offs/stepsize)]
                except IndexError:
                    if zero+int(offs/stepsize)>=yvals.shape[0]:
                        offsimg[x1,y1,z1,0]=yvals[-1]
                        indexerrors.append([x1,y1,z1])
                    elif zero-int(offs/stepsize)<0:
                        offsimg[x1,y1,z1,0]=yvals[0]
                        indexerrors.append([x1,y1,z1])
                    else:
                        raise IndexError('Index out of range at position (%r,%r,%r).'%(x1,y1,z1))
    else:
        for ind in range(x*y*z):
            x1,y1,z1 = np.unravel_index(ind, array.shape[:-1])
            if np.abs(B0[x1,y1,z1,0])>5:
                yvals=array[x1,y1,z1,:]
                zero=np.argmin(yvals)
                pseudofieldmap[x1,y1,z1]=freqs[zero]
                offsimg[x1,y1,z1,2]=yvals[zero]                
                
                freq=roundbase(offs-B0[x1,y1,z1,0],base=stepsize) #rounds to nearest number divisible by step size
                ind1=list(freqs).index(freq) #finds index of frequency of interest, positive
                ind2=list(freqs).index(-freq)
                offsimg[x1,y1,z1,0]=array[x1,y1,z1,ind1]
                offsimg[x1,y1,z1,1]=array[x1,y1,z1,ind2]
    return offsimg, pseudofieldmap, indexerrors

def unInterleave(array):
    narray=np.zeros(array.shape)
    evens=[i for i in range(array.shape[-1]) if i%2==0] 
    odds=[i for i in range(array.shape[-1]) if i%2!=0]            
    narray[:,:,:,:array.shape[-1]/2]=array[:,:,:,odds]
    positives=array[:,:,:,evens][:,:,:,::-1]
    narray[:,:,:,array.shape[-1]/2:]=positives
    return narray

def invgauss(x, *p):
    A, mu, sigma, offs = p
    return -A*np.exp(-(x-mu)**2/(2.*sigma**2))+offs

def fitPeak(xdata, ydata, il=True):
    if il:
        narr=zip(xdata,ydata)
        narr.sort()
        narr=np.array(narr)
        xdata=narr[:,0]
        ydata=narr[:,1]

    p0 = [np.min(ydata), 30, 10, 20000]
    coeff, var_matrix = curve_fit(invgauss, xdata, ydata, p0=p0)
    return xdata, ydata, coeff, var_matrix
    
def fillEdges(fieldmap,cestimg):
    newarray=np.zeros([cestimg.data.shape[0],cestimg.data.shape[1], fieldmap.data.shape[2],fieldmap.data.shape[3]])
    xdist=(cestimg.position[1]-fieldmap.position[1])/cestimg.resolution[0]
    ydist=(cestimg.position[0]-fieldmap.position[0])/cestimg.resolution[1]
    xstart=newarray.shape[1]+xdist-fieldmap.data.shape[1]+1
    ystart=newarray.shape[0]+ydist-fieldmap.data.shape[0]+1
    newarray[ystart:ydist,xstart:xdist,:,:]=fieldmap.data
    return newarray

def getCorrespondingSlice(fieldmap,cestimg):
    '''Isolates a slice of interest corresponding to the slice in a 2-D image. 
    Inputs are Bruker2AnalyzeImg.'''
    fieldmapendpos=(fieldmap.zposition-cestimg.zposition)/fieldmap.zthickness
    fieldmapstartpos=(fieldmap.zposition-(cestimg.zposition-cestimg.zthickness))/fieldmap.zthickness
    startslice=fieldmap.data.shape[2]-int(fieldmapstartpos)
    endslice=fieldmap.data.shape[2]-int(math.ceil(fieldmapendpos))+1
    return fieldmap.data[:,:,startslice:endslice,:]
    
def roundbase(x, base=25):
    return int(base * round(float(x)/base))