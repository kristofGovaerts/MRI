#Copyright  Copyright 2013 by Kristof Govaerts <kristof88@gmail.com>

'''Modules that process Analyze images.'''

import pylab as pl
import numpy as np
from scipy.optimize import *
import scipy.ndimage as ndi
import sys
import warnings
import matplotlib
import math

def t2r(x, si0, t2, soff):
    """Basic T2 regression function incorporating noise offset."""
    return si0*np.exp(-x/t2) + soff

def t2r_basic(x, si0, t2):
    """Basic T2 regression function (no noise offset)."""
    return si0*np.exp(-x/t2)

def t1r(x, si0, t1):
    """Basic inversion recovery T1 regression function."""
    return si0*abs(1-2*np.exp(-x/t1))
    
def t1r_vtr(x, si0, t1, soff):
    """Basic VTR recovery T1 regression function incorporating noise offset."""
    return si0*(1-np.exp(-x/t1)) + soff
    
def t1r_vtr_basic(x, si0, t1):
    """Basic VTR T1 regression function."""
    return si0*(1-np.exp(-x/t1))

def biexp_aslr(ti, m0, alpha, cbf, lambd, t1a, t1bl):
    """Biexponential ASL signal function."""
    return 2*m0*alpha*(cbf/lambd)*((np.exp(-ti/t1a)-np.exp(-ti/t1bl))/((1/t1bl)-(1/t1a)))

def dkifunc(b, si0, adc, akc):
    """Expanded Stejskal-Tanner equation."""
    return np.log(si0) - b*adc + (b**2 * adc**2 * akc)/6.0

def T2regression(img, te, mask=True, mode=t2r):
    ''''T2 curve fitting. Based on image intensities, it selects which pixels
    should be processed and uses scipy.optimize.curve_fit to fit the data.
    
    'mode' can be 't2r' for T2 regression or 't1r_vtr' for T1 regression with 
    variable TR.'''
    te0 = img[:,:,:,0]
    x,y,z,t=img.shape
    iter = x * y * z
    t2map = np.zeros((x, y, z)) #create arrays beforehand. Fill it with zeros so that you don't need an Else statement for the threshold
    stdevmap = np.zeros((x, y, z))
    si0map = np.zeros((x, y, z))
    soffmap = np.zeros((x, y, z))
    noiselevel=np.mean(img[:10,:10,:,0]) #takes 10x10x0xz area and uses mean signal intensity as noise level for starting values in regression
    
    if mask:
        mask=brain_mask(img,size=4) 
    else:
        mask=np.ones(t2map.shape)
    
    print "Iterating over image..."

    for px in range(0, iter): #iterates over pixels in image
        x1,y1,z1 = np.unravel_index(px, t2map.shape) #index stuff
        steps=iter/100
        if float(px)%steps == 0:
            print '\r', round(float(px*100)/(iter-1)), 'percent complete', #fancy progress display. First and last few percent are faster because less pixels are processed here
            sys.stdout.flush()
        if mask[x1,y1,z1] > 0: 
            si=np.array(img[x1, y1, z1,:])
            if mode == t2r:
                x0 = np.array([1.5*si[0], 33, noiselevel])
            elif mode == t1r_vtr:
                x0 = np.array([si[0], 1000, noiselevel])
            try:
                with warnings.catch_warnings(record=True) as w: #makes things cleaner, otherwise it throws overflow warnings occasionally. Warnings also stop progress bar from working
                    coeffs, vmat = curve_fit(mode,te,si, x0)
                    if type(vmat) == np.ndarray:
                        t2map[x1,y1,z1]=coeffs[1]
                        stdevmap[x1,y1,z1]=np.sqrt(vmat[1,1])
                        si0map[x1,y1,z1]=coeffs[0]
                        soffmap[x1,y1,z1]=coeffs[2]
                    else: #means that covar matrix is Inf
                        t2map[x1,y1,z1]=0.0
                        stdevmap[x1,y1,z1]=0.0
                        si0map[x1,y1,z1]=0.0
                        soffmap[x1,y1,z1]=0.0
            except RuntimeError:  #regression occasionally does not converge. This is fairly rare and the -1.0 value can be visualized in images and is easy to find as well.
                t2map[x1,y1,z1]=-1.0
                si0map[x1,y1,z1]=-1.0
                soffmap[x1,y1,z1]=-1.0
                stdevmap[x1,y1,z1]=-1.0

    sigraph = [float(np.mean(img[(x/2-2):(x/2+2), (y/2-2):(y/2+2), z/2,i])) for i in range(img.shape[3])] #selects 4x4 region in middle of image for test fit. "float" is because it is otherwise a memmap
    err=[float(np.std(img[(x/2-2):(x/2+2), (y/2-2):(y/2+2), z/2,i])) for i in range(img.shape[3])]
    x0g = np.array([1.5*sigraph[0], 35, noiselevel])
    try:
        gcoeffs, gvmat = curve_fit(t2r,te,sigraph, x0g)
    except RuntimeError:
        sigraph = [float(np.mean(img[(x/2-20):(x/2-16), (y/2-2):(y/2+2), z/2,i])) for i in range(img.shape[3])]
        gcoeffs, gvmat = curve_fit(t2r,te,sigraph, x0g)
    t2fig, (ax1, ax2)=pl.subplots(2)
    pl.errorbar(te, sigraph, fmt='co', yerr=err, label='Raw data')#individual data points + std dev
    ax2.plot(range(0,int(np.max(np.array(te)))+1,1), t2r(np.array(range(0,int(np.max(np.array(te)))+1,1)), gcoeffs[0], gcoeffs[1], gcoeffs[2]), label=r'Fit: $SI=SI_0 * e^{\frac{-TE}{T_2}} + S_{offs}$', color='blue')#Fancy formula display with TEX-like typesetting
    pl.xlabel('Echo time (ms)')
    pl.ylabel('Signal intensity (a.u.)')
    pl.xticks(te) #One tick per echo time
    ax2.legend(loc='upper right')
    ax2.grid()

    ax1.imshow(t2map[:,:,z/2],interpolation='nearest')
    ax1.axis('off')
    ax1.add_patch(matplotlib.patches.Rectangle((x/2-2, y/2-2), 4,4, fill=False, color='yellow'))
    ax1t='Slice %r' %(z/2 +1)
    ax1.text(10,20, ax1t, color='white', family=('serif'))
    pl.tight_layout()

    return t2map, stdevmap, si0map, soffmap, t2fig


def T2regression_basic(img, te, mask=True, mode=t2r_basic):
    ''''T2 curve fitting. Based on image intensities, it selects which pixels
    should be processed and uses scipy.optimize.curve_fit to fit the data.'''
    te0 = img[:,:,:,0]
    x,y,z,t=img.shape
    iter = x * y * z
    t2map = np.zeros((x, y, z)) #create arrays beforehand. Fill it with zeros so that you don't need an Else statement for the threshold
    stdevmap = np.zeros((x, y, z))
    si0map = np.zeros((x, y, z))
    noiselevel=np.mean(img[:10,:10,:,0]) #takes 10x10x0xz area and uses mean signal intensity as noise level for starting values in regression

    if mask:
        mask=brain_mask(img,size=4) 
    else:
        mask=np.ones(t2map.shape)

    print "Iterating over image..."

    for px in range(0, iter): #iterates over pixels in image
        x1,y1,z1 = np.unravel_index(px, t2map.shape) #index stuff
        steps=iter/100
        if float(px)%steps == 0:
            print '\r', round(float(px*100)/(iter-1)), 'percent complete', #fancy progress display. First and last few percent are faster because less pixels are processed here
            sys.stdout.flush()
        if mask[x1,y1,z1] > 0: 
            si=np.array(img[x1, y1, z1,:])
            if mode == t2r_basic:
                x0 = np.array([1.5*si[0], 36])
            elif mode == t1r_vtr_basic:
                x0 = np.array([si[0], 1000])
            try:
                with warnings.catch_warnings(record=True) as w: #makes things cleaner, otherwise it throws overflow warnings occasionally. Warnings also stop progress bar from working
                    coeffs, vmat = curve_fit(mode,te,si, x0)
                    if type(vmat) == np.ndarray:
                        t2map[x1,y1,z1]=coeffs[1]
                        stdevmap[x1,y1,z1]=np.sqrt(vmat[1,1])
                        si0map[x1,y1,z1]=coeffs[0]
                    else: #means that covar matrix is Inf
                        t2map[x1,y1,z1]=0.0
                        stdevmap[x1,y1,z1]=0.0
                        si0map[x1,y1,z1]=0.0
            except RuntimeError:  #regression occasionally does not converge. This is fairly rare and the -1.0 value can be visualized in images and is easy to find as well.
                t2map[x1,y1,z1]=-1.0
                si0map[x1,y1,z1]=-1.0
                stdevmap[x1,y1,z1]=-1.0

    sigraph = [float(np.mean(img[(x/2-2):(x/2+2), (y/2-2):(y/2+2), z/2,i])) for i in range(img.shape[3])] #selects 4x4 region in middle of image for test fit. "float" is because it is otherwise a memmap
    err=[float(np.std(img[(x/2-2):(x/2+2), (y/2-2):(y/2+2), z/2,i])) for i in range(img.shape[3])]
    x0g = np.array([1.5*sigraph[0], 35])
    gcoeffs, gvmat = curve_fit(t2r_basic,te,sigraph, x0g)
    t2fig, (ax1, ax2)=pl.subplots(2)
    pl.errorbar(te, sigraph, fmt='co', yerr=err, label='Raw data')#individual data points + std dev
    ax2.plot(range(0,int(np.max(np.array(te)))+1,1), t2r_basic(np.array(range(0,int(np.max(np.array(te)))+1,1)), gcoeffs[0], gcoeffs[1]), label=r'Fit: $SI=SI_0 * e^{\frac{-TE}{T_2}}$', color='blue')#Fancy formula display with TEX-like typesetting
    pl.xlabel('Echo time (ms)')
    pl.ylabel('Signal intensity (a.u.)')
    pl.xticks(te) #One tick per echo time
    ax2.legend(loc='upper right')
    ax2.grid()

    ax1.imshow(t2map[:,:,z/2],interpolation='nearest')
    ax1.axis('off')
    ax1.add_patch(matplotlib.patches.Rectangle((x/2-2, y/2-2), 4,4, fill=False, color='yellow'))
    ax1t='Slice %r' %(z/2 +1)
    ax1.text(10,20, ax1t, color='white', family=('serif'))
    pl.tight_layout()

    return t2map, stdevmap, si0map, t2fig

def ASLimg_split(img):
    '''Splits img in two images, where output 1 is where the input image
    had even numbers in the Z (3rd) dimension, and output 2 had odd numbers.
    The first frame has index 0 (not 1 like in Matlab).'''
    print "Splitting image..."
    evens = [x for x in range(img.shape[3]) if x%2 == 0]
    odds = [x for x in range(img.shape[3]) if x%2 != 0]
    evimg = img[:,:,:,evens]
    oddimg = img[:,:,:,odds]
    return evimg, oddimg


def T1Regression(img, ti, dims):
    '''T1 relaxometry via inversion recovery. It selects pixels to be
    processed based on image intensity.'''
    ti0 = img[:,:,:,0]
    x,y,z,t=dims
    iter = x * y * z
    t1map = np.zeros((x, y, z))
    si0map = np.zeros((x, y, z))
    stdevmap = np.zeros((x,y,z))
    noiselevel=np.mean(img[:10,:10,0,0])
    
    mask=brain_mask(img) 

    print "Iterating over image..."

    for px in range(0, iter -1): #iterates over pixels in image
        x1,y1,z1 = np.unravel_index(px, t1map.shape) #index stuff
        steps=iter/100
        if float(px)%steps == 0:
            print '\r', round(float(px*100)/(iter-1)), 'percent complete', #fancy progress display. First and last few percent are faster because less pixels are processed here
            sys.stdout.flush()
        if mask[x1,y1,z1] > 0: #could be different for other types of T1 map. Seems to work.
            si=np.array(img[x1, y1, z1, range(t)])
            x0 = np.array([si[0], 1500])
            try:
                with warnings.catch_warnings(record=True) as w: #makes things cleaner, otherwise it throws overflow warnings occasionally. Warnings also stop progress bar from working
                    coeffs, vmat = curve_fit(t1r,np.array(ti),si, x0)
                    if type(vmat) == np.ndarray:
                        t1map[x1,y1,z1]=coeffs[1]
                        stdevmap[x1,y1,z1]=np.sqrt(vmat[1,1])
                        si0map[x1,y1,z1]=coeffs[0]
            except RuntimeError:  #regression occasionally does not converge. This is fairly rare and the -1.0 value can be visualized in images and is easy to find as well.
                t1map[x1,y1,z1]=-1.0
                si0map[x1,y1,z1]=-1.0

    sigraph = [float(np.mean(img[(x/2-2):(x/2+2), (y/2-2):(y/2+2), z/2,i])) for i in range(img.shape[3])]
    err=[float(np.std(img[range(x/2-2, x/2+2), range(y/2-2, y/2+2), z/2,i])) for i in range(img.shape[3])]
    x0g = np.array([sigraph[0], 1500])
    gcoeffs, gvmat = curve_fit(t1r,np.array(ti),sigraph, x0g)
    t1fig, (ax1, ax2)=pl.subplots(2)
    pl.errorbar(ti, sigraph, fmt='co', yerr=err, label='Raw data')#individual data points + std dev

    ax2.plot(range(0,int(ti[-1])+1,1), t1r(np.array(range(0,int(ti[-1])+1,1)), gcoeffs[0], gcoeffs[1]), label=r'Fit: $SI=SI_0 * |1-2e^{\frac{-TI}{T_1}}|$', color='blue')#Fancy formula display with TEX-like typesetting
    pl.xlabel('Inversion time (ms)', fontsize=12)
    pl.ylabel('Signal intensity (a.u.)', fontsize=12)
    pl.yticks(fontsize=11)
    pl.xticks(ti, rotation=45) #One tick per echo time
    ax2.legend(loc='lower right')
    ax2.grid()

    ax1.imshow(t1map[:,:,z/2],interpolation='nearest')
    ax1.axis('off')
    ax1.add_patch(matplotlib.patches.Rectangle((x/2-2, y/2-2), 4,4, fill=False, color='yellow'))
    ax1t='Slice %r' %(math.ceil(z/2) + 1)
    ax1.text(10,20, ax1t, color='white', family=('serif'))
    pl.tight_layout()
    return t1map, si0map, stdevmap, t1fig
    
def within_bounds(x, l=0., u=1.):
    if x>u or x<l:
        return False
    else:
        return True
     

def biexp_asl(self):
    mask=self.mask
    ti=self.ti
    t1amap=self.selt1
    deltaM=np.array(self.deltaM)
    lambd=0.83
    t1bl=self.t1bl

    cbfmap=np.zeros(self.shape[:3])
    stdevmap=np.zeros(self.shape[:3])
    m0map=np.zeros(self.shape[:3])
    alphamap=np.zeros(self.shape[:3])
    iter=np.prod(self.shape[:3])
    for i in range(deltaM.shape[3]):
        deltaM[:,:,:,i]*=mask

    print "Performing biexponential ASL regression..."

    for px in range(0, iter -1): #iterates over pixels in image
        x1,y1,z1 = np.unravel_index(px, cbfmap.shape) #index stuff
        steps=iter/100
        if float(px)%steps == 0:
            print '\r', round(float(px*100)/(iter-1)), 'percent complete', #fancy progress display. First and last few percent are faster because less pixels are processed here
            sys.stdout.flush()
        if mask[x1,y1,z1] > 0: #Using mask instead of just threshold. Hard to mask deltaM image due to low SNR
            si=np.array(deltaM[x1, y1, z1, range(deltaM.shape[3])])
            tsi=list(si)
            tti=list(ti)
            indices=[]
            si=np.array(tsi)
            tti=np.array(tti)
            x0 = np.array([0.95, 0.00150])
            try:
                with warnings.catch_warnings(record=True) as w: #makes things cleaner, otherwise it throws overflow warnings occasionally. Warnings also stop progress bar from working
                    #as T1app varies from voxel to voxel, need to define custom function for every voxel. Lambda da bes
                    func=lambda ti, alpha, cbf: 2*self.oddsi0[x1,y1,z1]*alpha*(cbf/lambd)*((np.exp(-ti/t1amap[x1,y1,z1])-np.exp(-ti/t1bl))/((1/t1bl)-(1/t1amap[x1,y1,z1]))) if within_bounds(alpha) else 1e8
                    #func=lambda ti, m0, alpha, cbf: 2*m0*alpha*(cbf/lambd)*((np.exp(-ti/1700.0)-np.exp(-ti/2400.0))/((1/2400.0)-(1/1600.0)))
                    coeffs, vmat = curve_fit(func, tti, si, x0, maxfev=100000)
                    if type(vmat) == np.ndarray:
                        cbfmap[x1,y1,z1]=coeffs[1]
                        stdevmap[x1,y1,z1]=np.sqrt(vmat[1,1]) #standard dev of CBF fit
                        m0map[x1,y1,z1]=self.oddsi0[x1,y1,z1]
                        alphamap[x1,y1,z1]=coeffs[1]
            except RuntimeError:  #regression occasionally does not converge. This is fairly rare and the -1.0 value can be visualized in images and is easy to find as well.
                    cbfmap[x1,y1,z1]=coeffs[1]
                    m0map[x1,y1,z1]=self.oddsi0[x1,y1,z1]
                    alphamap[x1,y1,z1]=coeffs[1]

    return cbfmap, alphamap, stdevmap, deltaM

def dki_fit(self):
    from PythonDiffusion import DKIfit
    b=self.avbvals[:self.nA0+self.nbvals]
    x,y,z,t=self.shape
    iter = x * y * z
    adcmap = np.zeros((x, y, z, self.ndirs)) #create arrays beforehand. Fill it with zeros so that you don't need an Else statement for the threshold
    akcmap = np.zeros((x, y, z, self.ndirs))
    si0map = np.zeros((x, y, z, self.ndirs))
    akcstdevmap = np.zeros((x, y, z, self.ndirs))
    adcstdevmap = np.zeros((x, y, z, self.ndirs))
    si0stdevmap = np.zeros((x, y, z, self.ndirs))
    noiselevel=np.mean(self.pdata[:10,:10,:,0]) #takes 10x10x0xz area and uses mean signal intensity as noise level for starting values in regression
    
    print "Iterating over image..."
    for ddir in range(self.ndirs):
        print "Now processing direction %r of %r." %(ddir+1, self.ndirs)
        for px in range(0, iter): #iterates over pixels in image
            x1,y1,z1 = np.unravel_index(px, adcmap.shape[:-1]) #index stuff
            steps=iter/100
            if float(px)%steps == 0:
                print '\r', round(float(px*100)/(iter-1)), 'percent complete', #fancy progress display. First and last few percent are faster because less pixels are processed here
                sys.stdout.flush()
            if self.pdata[x1,y1,z1,0] > 10*noiselevel: #could be different for other types of DW img.
                a0=np.log(self.pdata[x1, y1, z1, :self.nA0])
                si=np.log(self.pdata[x1, y1, z1,self.nA0+self.nbvals*ddir:self.nA0+self.nbvals*(ddir+1)]) #logarithmic fit
                si=np.append(a0, si)
                x0 = np.array([np.mean(si[0:self.nA0]), 0.001, 2])
                try:
                    with warnings.catch_warnings(record=True) as w: #makes things cleaner, otherwise it throws overflow warnings occasionally. Warnings also stop progress bar from working            
                        coeffs, vmat = curve_fit(dkifunc, b, si, x0)
                        if type(vmat) == np.ndarray:
                            adcmap[x1,y1,z1, ddir]=coeffs[1]
                            akcmap[x1,y1,z1, ddir]=coeffs[2]
                            si0map[x1,y1,z1, ddir]=coeffs[0]
                            adcstdevmap[x1,y1,z1, ddir]=np.sqrt(vmat[1,1])
                            akcstdevmap[x1,y1,z1, ddir]=np.sqrt(vmat[2,2])
                            si0stdevmap[x1,y1,z1, ddir]=np.sqrt(vmat[3,3])
                        else: #means that covar matrix is Inf
                            pass #points are zero anyways
                except RuntimeError:  #regression occasionally does not converge. This is fairly rare and the -1.0 value can be visualized in images and is easy to find as well.
                    adcmap[x1,y1,z1, ddir]=-1.0
                    akcmap[x1,y1,z1, ddir]=-1.0
                    si0map[x1,y1,z1, ddir]=-1.0
                    adcstdevmap[x1,y1,z1, ddir]=-1.0
                    akcstdevmap[x1,y1,z1, ddir]=-1.0
                    si0stdevmap[x1,y1,z1, ddir]=-1.0
    
    dkifit = DKIFit(adcmap, akcmap, si0map, adcstdevmap, akcstdevmap, si0stdevmap, self.bvecs) 
    return dkifit
        
def brain_mask(mat4d, size=2, numpass=2):
    from skimage.filter import threshold_otsu
    from scipy.ndimage import median_filter
    bm=np.array(mat4d)
    if bm.shape[3]>1:
        bm=np.mean(bm, axis=3)
    for i in range(numpass):
        bm=median_filter(bm, size=size)
    t=threshold_otsu(bm)
    bm[bm<t]=0
    bm[bm>=t]=1    
    return bm.astype(int)
    
def largest_component(array):
    labeled_array, numpatches = ndi.label(array) # labeling
    sizes = ndi.sum(array,labeled_array,range(1,numpatches+1)) 
    # To get the indices of all the min/max patches. Is this the correct label id?
    map = np.where(sizes==sizes.max())[0] + 1 
    mip = np.where(sizes==sizes.min())[0] + 1
    
    # inside the largest, respecitively the smallest labeled patches with values
    max_index = np.zeros(numpatches + 1, np.uint8)
    max_index[map] = 1
    max_feature = max_index[labeled_array]
    return max_feature