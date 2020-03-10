# -*- coding: utf-8 -*-
"""
Created on Fri Jun 06 13:18:13 2014

@author: u0091609
"""

import os
from glob import glob
import nibabel as nib
import numpy as np
import pylab as pl
import Tkinter, tkFileDialog
import pandas
root = Tkinter.Tk()

PAR_EXT="_MTRasym.nii"
POS_EXT="_Pos_corr.nii"
NEG_EXT="_Neg_corr.nii"
LAB_EXT="_LR.nii"
F_EXT = "_interpolated_freqs.csv"
LAB_NAMES=["Left", "Right"]
COL_NAMES = [" MTRasym", " Pos", " Neg"]

#RANGE=np.arange(0,2401,40)
IN_FOLDER=False
UPPER_BOUND=1e8
XRANGE_GRAPH = (-2000,2000) #range for the x-axis of the plots
YRANGE_GRAPH = (0,20) #range for the MTRasym graph

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

def quantify_labels(cim, lim):
    '''Quantifies the means and std devs of parameter image pim for every label
    in label image lim. Label image should consist of ROIs colored with integer
    values.'''
    hist = np.histogram(clab, bins=np.max(lim)) #Check how many labels there are
    values = hist[1][1:]
    print "Found %s labels." %len(values)
    lab_means=[]
    lab_stdevs=[]
    lab_datapoints=[]
    cmask = np.zeros(lim.shape)
    pim = np.array(cim)
    cmask[pim[:,:,:,-1]>0] = 1
    cmask[pim[:,:,:,-1]>UPPER_BOUND] = 0    
    pim[np.isnan(pim)] = 0    
    for i in range(pim.shape[-1]):
        pim[:,:,:,i] *= cmask 
    for v in values:
        print "label", v
        rmask=np.array(lim).astype('float')
        rmask[rmask!=v]=0
        rmask/=v
        par_means=[]
        par_stdevs=[]
        par_datapoints=[]
        for i in range(pim.shape[-1]):
            mpar=pim[:,:,:,i] * rmask
            par_means.append(np.mean(mpar[mpar!=0]))
            par_stdevs.append(np.std(mpar[mpar!=0]))   
            par_datapoints.append(mpar[mpar!=0])                
        lab_means.append(par_means)
        lab_stdevs.append(par_stdevs)
        lab_datapoints.append(par_datapoints)
    return lab_means, lab_stdevs, lab_datapoints
    
def give_plots(means, stdevs, pnames):
    fig=pl.figure(figsize=(40, 30), dpi=80)
    l=len(LAB_NAMES)
    x=np.sqrt(l)
    for p in range(l):
        pl.subplot(np.ceil(x), np.ceil(x), p+1)
        pl.title(pnames[p])
        y=np.array(means[p])
        error=np.array(stdevs[p])
        pl.plot(RANGE, y, 'k', color='#1B2ACC')
        pl.fill_between(RANGE, y-error, y+error,
        alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF',
        linewidth=4, linestyle='dashdot', antialiased=True)
    return fig    

def give_full_plots(lm, ls, pm, ps, nm, ns, pnames):
    fig=pl.figure(figsize=(40, 30), dpi=80)
    l=len(pnames)
    xa=np.sqrt(l+1)
    axc = fig.add_subplot(np.ceil(xa), np.ceil(xa), l+1)
    for p in range(l):
        ax1 = fig.add_subplot(np.ceil(xa), np.ceil(xa), p+1)
        ax2 = ax1.twinx()
        pl.title(pnames[p])
        y=(np.array(lm[p]),np.array(pm[p]),np.array(nm[p]))
        x=(RANGE, RANGE, -RANGE)
        e=(np.array(ls[p]),np.array(ps[p]),np.array(ns[p]))
        ax=(ax2, ax1, ax1)
        ylim=((0,50), None, None)
        ylab=('MTRasym (%)', 'SI (a.u.)', 'SI(a.u.)')   
        colors=('r','b','b')
        legend=("MTRasym", "SI", "")
        for i in range(3):
            ax[i].plot(x[i], y[i], color=colors[i])
            ax[i].fill_between(x[i], y[i]-e[i], y[i]+e[i], color=colors[i], alpha=0.2, linewidth=4, linestyle='dashdot', antialiased=True)
            ax[i].set_ylim(ylim[i])
            ax[i].set_ylabel(ylab[i], color=colors[i])
            ax[i].set_xlabel("Saturation offset (Hz)")
            ax[i].legend(legend[i])
        pl.xlim((XRANGE_GRAPH))
        
        axc.plot(x[0], y[0])
        axc.fill_between(x[0], y[0]-e[0], y[0]+e[0], alpha=0.2, linewidth=4, linestyle='dashdot', antialiased=True)
        axc.set_ylim(ylim[0])
        axc.set_ylabel(ylab[0])
        axc.set_xlabel("Saturation offset (Hz)")
        axc.set_xlim(0, XRANGE_GRAPH[-1])
        axc.set_ylim(0, YRANGE_GRAPH[-1])
        axc.legend(pnames)
        
    return fig    

def to_csv_line(string):
    substr=['[', ']', '(' , ')']
    o=str(string)
    for s in substr:
        o=o.replace(s, '')
    o=o.replace(',', ';')
    o=o.replace(".",",")
    return o

if IN_FOLDER:
    dirlist = [x[0] for x in os.walk(os.getcwd())][1:]
    par_imgs = [os.path.join(d, (d.split('\\')[-1] + PAR_EXT))  for d in dirlist]
    lab_imgs = [os.path.join(d, (d.split('\\')[-1] + LAB_EXT))  for d in dirlist]
else:
    par_imgs=glob("*" + PAR_EXT)
    pos_imgs=glob("*" + POS_EXT)
    neg_imgs=glob("*" + NEG_EXT)
    lab_imgs=glob("*" + LAB_EXT)

for lab_img in lab_imgs:
    try:
        clab=nib.load(lab_img).get_data()
        rname = lab_img[:-len(LAB_EXT)]
        cname = rname + PAR_EXT
        pname = rname + POS_EXT
        nname = rname + NEG_EXT
        fname = rname + F_EXT
        print "Current image:", cname
        
        freqs = pandas.read_csv(fname, sep=';')
        RANGE = np.array(freqs['pos'])
        
        cim=nib.load(cname).get_data()
        pim=nib.load(pname).get_data()
        nim=nib.load(nname).get_data()
        if len(cim.shape) == 3:
            sh=list(cim.shape)
            sh.append(1)
            cim=np.reshape(cim, sh)
        lm, ls, ld = quantify_labels(cim, clab)
        pm, ps, pd = quantify_labels(pim, clab)
        nm, ns, nd = quantify_labels(nim, clab)
        
        means = [lm, pm, nm]
        stdevs = [ls, ps, ns]
        #save csv file
        csv_arr=np.zeros((len(RANGE), 6*len(LAB_NAMES)))
        for a in range(len(LAB_NAMES)):
            for c in range(len(COL_NAMES)):
                csv_arr[:, c*2*len(LAB_NAMES) + 2*a]=means[c][a]
                csv_arr[:, c*2*len(LAB_NAMES) + 2*a + 1]=stdevs[c][a]
            
        header = ''
        for col in COL_NAMES:
            header += ';'.join(sum([[n+ col + " mean", n + col + " stdev"] for n in LAB_NAMES], [])) + ';'
            
        np.savetxt(rname + '_quant.csv', csv_arr, delimiter=';', header=header)
        
        f=give_full_plots(lm,ls,pm,ps,nm,ns,LAB_NAMES)
        f.savefig(rname + "_basic_plot.png")
    except IOError:
        print "Label image for %s not found." %lab_img