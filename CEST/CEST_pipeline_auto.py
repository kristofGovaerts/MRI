# -*- coding: utf-8 -*-
"""
Created on Mon May 26 14:09:18 2014

@author: u0091609
"""
import numpy as np
import os
from scipy.interpolate import interp1d, PchipInterpolator
import nibabel as nib
import Tkinter, tkFileDialog
root = Tkinter.Tk()

def order_im(img, freqs):
    '''Orders image in the t-domain according to frequencies, in
    ascending fashion. Necessary for scipy's interpolation functions.'''
    r=len(freqs) 
    orders=sorted(zip(freqs,range(r)))
    simg=np.zeros(img.shape)
    for i in range(r):
        simg[:,:,:,i]=img[:,:,:,orders[i][1]]
    return simg
    
def find_offset(x_freq, y_i, bin=0.1, return_xy=False):
    '''Finds the frequency offset from zero when provided a list of frequencies
    and a list of image intensities. In principle simply finds the minimum
    of a curve found via interpolation between points.
    
    x_freq = Frequencies
    y_i = Intensities
    bin = step size for the interpolation function. The smaller, the more 
    accurate but the more RAM you need.
    return_xy: if True, returns interpolated frequencies and intensities.'''
    ip=PchipInterpolator(x_freq, y_i) #Piecewise spline interpolation
    xn=np.arange(x_freq[0], x_freq[-1], bin)
    yn=ip(xn)
    if return_xy==False:
        output=xn[np.argmin(yn)]
    else:
        output=(xn[np.argmin(yn)],xn,yn)
    return output

def create_b0map(img, freqs, bin=1, thresh=1e9):
    '''Creates a b0map based on the frequency offsets calculated per voxel.'''
    cim_f=np.reshape(img, (np.prod(img.shape[:-1]), img.shape[-1]))    
    arr=np.zeros(len(cim_f))
    for vx in range(len(cim_f)):
        if cim_f[vx,-1]>thresh:
            arr[vx]=find_offset(freqs, cim_f[vx,:], bin=bin)
    b0map=np.reshape(arr, img.shape[:-1])
    return b0map
    
def create_alphamap(si_a, si_2a, in_degrees = True):
    '''Creates a flip angle map. Inputs are signal intensity at flip angle 
    alpha and flip angle 2alpha.'''
    alpha=np.arccos(si_2a/(2*si_a))
    if in_degrees:
        alpha *= (180/np.pi)
    return alpha
    
def find_nearest(array,value):
    '''Finds index of nearest value in an array.'''
    idx = (np.abs(array-value)).argmin()
    return idx
    
def correct_b0(img, freqs, b0, offsets, bin=10, thresh=1e5):
    '''Corrects for B0 effects and returns desired offset images.
    
    img = 4D array
    freqs = frequencies corresponding in order to the 4th dim of the img
    offset = desired CEST offset (in a list, allows multiple)
    thresh = Masking threshold
    
    Returns:
    b0map - b0 offset map
    plusmap - +offset map corrected for b0
    minmap - -offset map corrected for b0'''
    cim_f=np.reshape(img, (np.prod(img.shape[:-1]), img.shape[-1]))    
    b0map=np.reshape(b0, (np.prod(img.shape[:-1])))    
    plusmap=np.zeros((len(cim_f),len(offsets)))
    minmap=np.zeros((len(cim_f),len(offsets)))
    
    for vx in range(len(cim_f)):
        if cim_f[vx,-1]>thresh:
            ip=PchipInterpolator(freqs, cim_f[vx, :])
            xn=np.arange(freqs[0], freqs[-1], bin)
            yn=ip(xn)
            off = b0map[vx]
            for i,o in enumerate(offsets):
                plus=yn[find_nearest(xn, o + off)]
                minus=yn[find_nearest(xn, -o + off)]
                plusmap[vx,i]=plus
                minmap[vx,i]=minus
           
    nshape=list(img.shape[:-1])        
    nshape.append(len(offsets))
    plusmap=np.reshape(plusmap, nshape)
    minmap=np.reshape(minmap, nshape)
    return plusmap, minmap 
    
def cest_spectrum(img, freqs, n=150, step=10, thresh=1e9):
    cim_f=np.reshape(img, (np.prod(img.shape[:-1]), img.shape[-1])) 
    b0map=np.zeros(len(cim_f))
    cmap=np.zeros((len(cim_f), n+1))
    
    for vx in range(len(cim_f)):
        if cim_f[vx,-1]>thresh:
            off, xn, yn=find_offset(freqs, cim_f[vx,:], bin=step, return_xy=True)
            zero=find_nearest(xn, off)
            minus=yn[:zero+1][-n-1:][::-1] #have to reverse 
            plus=yn[zero+1:][:n+1]
            try:
                cest=(plus-minus)/plus
                cmap[vx,:]=cest
            except ValueError:
                cmap[vx,:]=-np.ones(n+1)
            b0map[vx]=off
    
    b0map=np.reshape(b0map, img.shape[:-1])
    cmap=np.reshape(cmap, list(img.shape[:-1])+[n+1])
    return b0map, cmap
    
def scansToArray(filelist, start=0, end=None, sorted=False):
    '''Returns an array concatenating a series of 3D arrays across the
    t-dimension.

    Inputs:
        filelist = list of files, no extension
        start = where to start concatenating, optional (necessary if matrix size not equal for all images in your filelist)
        end=where to stop'''
    if sorted == False:
        sfl=sort_scans(filelist)
    else:
        sfl=filelist
    if end is None:
        scans=sfl[start:]
    else:
        scans=sfl[start:end]
    data1=nib.load(scans[0]).get_data()
    array=np.zeros([data1.shape[0], data1.shape[1],1, len(scans)])
    array[:,:,:,0]=data1[:,:,:,0]
    for i, scan in enumerate(scans[1:]):
        data=nib.load(scan).get_data()
        array[:,:,:,i+1]=data[:,:,:,0]
    return array
    
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

print "Enter the scan directory."

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

with open("CESTparams.txt") as f:
    pfile=f.readlines()
    
wrange = pfile[2].split('\t')[-1][1:-2].split(',') #WASSR range from second line in file
full_wrange = range(int(wrange[0]), int(wrange[1]) + 1)
crange = pfile[3].split('\t')[-1][1:-2].split(',') #CEST range from third line in file
full_crange = range(int(crange[0]), int(crange[1]) + 1)
coffs = pfile[4].split('\t')[-1][1:-2].split(',') #CEST offset from fourth line in file
coffs = [int(c) for c in coffs]

saveB0 = pfile[5].split('\t')[1][:-1]
saveZ = pfile[6].split('\t')[1][:-1]
savecorr = pfile[7].split('\t')[1][:-1]
thresh = float(pfile[8].split('\t')[1][:-1])
cest_bin = float(pfile[9].split('\t')[1][:-1])
wassr_bin = float(pfile[10].split('\t')[1][:-1])
filename = pfile[11].split('\t')[1][:-1]
save0hz = pfile[12].split('\t')[1][:-1]

wassr_files=[]
wassr_freqs=[]
cest_files=[]
cest_freqs=[]

for l in pfile: #somewhat ugly way of reading in the required files
    pars = l.split('\t')
    if ' ' + pars[0] + ', ' in str(full_wrange) or ' ' + pars[0] + ']' in str(full_wrange) or '[' + pars[0] + ', ' in str(full_wrange) or '[' + pars[0] + ']' in str(full_wrange):
        wassr_files.append(pars[1])
        wassr_freqs.append(int(pars[5][:-1]))
    elif ' ' + pars[0] + ', ' in str(full_crange) or ' ' + pars[0] + ']' in str(full_crange) or '[' + pars[0] + ', ' in str(full_crange) or '[' + pars[0] + ']' in str(full_crange):
        cest_files.append(pars[1])
        cest_freqs.append(int(pars[5][:-1]))

wassr = scansToArray(wassr_files, sorted = True)
cest = scansToArray(cest_files, sorted = True)

wassr_o=order_im(wassr, wassr_freqs)
wassr_f_o=sorted(wassr_freqs)
cest_o=order_im(cest, cest_freqs)
cest_f_o=sorted(cest_freqs)
affine=nib.load(wassr_files[0]).get_affine()

if saveZ == "True":
    print "Saving Z-spectrum..."
    zspec=np.zeros((wassr_o.shape[0], wassr_o.shape[1], wassr_o.shape[2], wassr_o.shape[3] + cest_o.shape[3]))
    cestlen=cest_o.shape[3]
    wassrlen = wassr_o.shape[3]
    zspec[:,:,:,:cestlen/2] = cest_o[:,:,:,:cestlen/2]
    zspec[:,:,:,cestlen/2:cestlen/2+wassrlen] = wassr_o
    zspec[:,:,:,cestlen/2+wassrlen:] = cest_o[:,:,:,cestlen/2:]
    nif=nib.Nifti1Image(zspec, affine)
    nib.save(nif, filename + "_Zspec.nii")

print "Calculating B0 map..."
b0 = create_b0map(wassr_o, wassr_f_o, bin = wassr_bin, thresh=thresh)

if saveB0 == 'True':
    print "Saving B0 map..."
    b0nif = nib.Nifti1Image(-b0, affine)
    nib.save(b0nif, filename + '_B0.nii')
    
print "Performing B0 correction..."
p, m = correct_b0(cest_o, cest_f_o, b0, coffs, bin=cest_bin, thresh=thresh)

if savecorr == 'True': 
    print "Saving B0 corrected images at frequencies", coffs
    nifplus = nib.Nifti1Image(p, affine)
    nifmin = nib.Nifti1Image(m, affine)
    nib.save(nifplus, filename + '_pos_corr.nii')
    nib.save(nifmin, filename + '_neg_corr.nii')

pos_cest=(100*(m-p))/m #save 4D arrays containing CEST images with negative offset image as ref
neg_cest=(100*(p-m))/p #same, but with positive image as reference
nifpos=nib.Nifti1Image(pos_cest, affine)
nifneg=nib.Nifti1Image(neg_cest, affine)
nib.save(nifpos, filename + 'pos_cest.nii')
nib.save(nifneg, filename + 'neg_cest.nii')

for i,f in enumerate(coffs):    
    cest_calc = 100 * (m[:,:,:,i]-p[:,:,:,i])/m[:,:,:,i]
    nifcest = nib.Nifti1Image(cest_calc, affine)
    ext = '%sHz.nii' %str(f) 
    nib.save(nifcest, filename + ext)

    
if save0hz == 'True':
    print "Saving B0 corrected image at 0Hz."
    p0, m0 = correct_b0(wassr_o, wassr_f_o, b0, range(-50,51), bin=1, thresh=thresh)
    nif0hz = nib.Nifti1Image(p0, affine)
    ext = '0Hz_B0corr.nii' 
    nib.save(nif0hz, filename + ext)