import numpy as np
import pylab as pl
import math
from PIL import Image
import os
import matplotlib.animation as animation

def saveslices3d(img, clim, title):
    '''Saves a single PNG image of 3D image array img, equalizes colormap with
    list object clim and saves it under filename title.'''
    nslices=img.shape[2]
    axis = np.ceil(math.sqrt(nslices))
    fig = pl.figure(num=None, figsize=(10, 6), dpi=80)
    fig.suptitle(title, fontsize = 14)
    for n in range(nslices):
        if axis != 1:
            pl.subplot(axis-1, axis+1, n+1, frameon=False, xticks=[], yticks=[])
        else:
            pl.subplot(1, 1, n+1, frameon=False, xticks=[], yticks=[])
        pl.imshow(rotateslice(img[:,:,n]), interpolation='nearest')
        stitle="Slice %r" %(n+1)
        pl.clim(clim)
        cb=pl.colorbar(shrink=.75, ticks=range(0, int(round(clim[1])), int(round(clim[1]))/4))
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(8)
        pl.xlabel(stitle, fontsize=8)
    pl.tight_layout()
    pl.savefig(title + '.png', bbox_inches=0)

def savetiffs(img, filename):
    '''Creates a new dir called filename and saves each slice of
    3-D image array img as filename_s[x]. Tiff preserves intensity
    values, unlike png.'''
    if os.path.isdir(filename) is False:
        os.mkdir(filename)
    for x in range(img.shape[2]):
        im = Image.fromarray(img[:,:,x]).rotate(270)
        st = filename + '_s' + str(x+1) + '.tiff'
        st = os.path.join(filename, st)
        im.save(st)

def rotateslice(sl):
    img=Image.fromarray(sl).rotate(270)
    return np.array(img)

def animateSlices(data, data2=None, sl=0, fps=10):
    '''Animates a 4-D image array, or two arrays side by side. Default animation
    is at 10fps.

    Output is a matplotlib figure which can be saved as a movie by calling
    ani.save('filename.mp4').

    Note that ffmpeg must be installed and available on the PATH before movies
    can be saved.'''
    fig=pl.figure()
    if data2 == None:
        im1=pl.imshow(data[:,:,sl,0])
        pl.xticks(())
        pl.yticks(())
    else:
        pl.subplot(1,2,1)
        im1=pl.imshow(data[:,:,sl,0])
        pl.xticks(())
        pl.yticks(())
        pl.subplot(1,2,2)
        pl.xticks(())
        pl.yticks(())
        im2=pl.imshow(data2[:,:,sl,0])
        pl.tight_layout()
    def update_img(i):
        f=data[:,:,sl,i]
        im1.set_data(f)
        return im1
    def update_imgs(i):
        f1=data[:,:,sl,i]
        f2=data2[:,:,sl,i]
        im1.set_data(f1)
        im2.set_data(f2)
        return im1, im2
    if data2 == None:
        ani=animation.FuncAnimation(fig, update_img, data.shape[3], interval=1000/fps)
    else:
        ani=animation.FuncAnimation(fig, update_imgs, data.shape[3], interval=1000/fps)
    return ani
    
def animate_overlay(vol, overlay, fps=20, clim=None):
    fig=pl.figure()
    ovol=vol+1000*overlay
    im1=pl.imshow(ovol[:,:,0], cmap='gray')
    if clim != None:
        pl.clim(clim)
    pl.xticks(())
    pl.yticks(())
    
    def update_img(i):
        f=ovol[:,:,i]
        im1.set_data(f)
        return im1
    
    ani=animation.FuncAnimation(fig, update_img, vol.shape[2], interval=1000/fps)
    return ani