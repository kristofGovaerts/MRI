# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 17:35:51 2014

@author:Kristof Govaerts
        Biomedical MRI Unit
        KU Leuven 
        u0091609 
        
Container class and functions for blood vessel analysis. Initiating an instance
of the VesselSkeleton class allows you to easily run various functions and 
visualize graphs.
"""

import os
import nibabel as nib
import pylab as pl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from ReadInput import *
import math

class VesselSkeleton:
    def __init__(self, name, tEXT='_VesselPars', thresh=0):
        self.name=name
        Img=nib.load(name+'.img')
        self.img=np.squeeze(Img.get_data())
        if thresh!=0:
            self.img[self.img!=96]=0 #imageJ stores longest most efficient path at value 96 (not 255). Go figure.
        else:
            self.avthickness=np.mean(self.img[self.img>0])
        self.affine=Img.get_affine()
        self.start=list_values(read_line('StartPoint=', name+tEXT))
        self.end=list_values(read_line('EndPoint=', name+tEXT))
        self.res=list_values(read_line('Resolution=', name+tEXT))
        self.world_pos=list_values(read_line('WorldPos=', name+tEXT))
        self.total_thickness=list_values(read_line('MeanThickness=', name+tEXT))[0]
        self.total_stdev=list_values(read_line('ThicknessVariance=', name+tEXT))[0]
        
    def tortuosity(self):
        self.pts=im_to_points(self.img)
        self.coords=pts_to_coords(self.pts, self.res, self.world_pos)
        tdist, edist, fraction, points=curve_length(tuple(map(tuple, self.coords)), start=None, end=None, return_pts=True)
        edist=dist(points[0], points[-1])     #comment out if needed?
        tort=tdist/edist
        self.points=points
        return tort, tdist, edist, fraction, points
        
    def save_tracked_pts(self):
        tp=(np.array(self.points)-self.world_pos)/self.res
        t_array=np.zeros(self.img.shape)
        for i in range(len(tp)):
            t_array[tuple(tp[i].astype('int'))]=i+1
        tpname=self.name+'_points.nii'
        tpim=nib.Nifti1Image(t_array, self.affine)
        nib.save(tpim, tpname)
        print "Tracked points saved in", tpname
            
    def plot2D(self, dim=0, plot_traveled=True):
        '''Plots a 2D projection ofthe vessel as points on a scatterplot.'''
        if plot_traveled:
            pts=np.array(self.points)
        else:
            pts=self.coords
        c=[0,1,2]
        c.pop(dim)
        pl.plot(pts[:,c[0]], pts[:,c[1]], '-o')
        wp1,wp2=self.world_pos[c[0]], self.world_pos[c[1]]
        pl.xlim(wp1,self.res[c[0]]*self.img.shape[c[0]]+wp1)
        pl.ylim(wp2,self.res[c[1]]*self.img.shape[c[1]]+wp2)
        l=['x','y','z']
        l.pop(dim)
        pl.xlabel=l[0]
        pl.ylabel=l[1]
        
    def plot3D(self):
        fig=pl.figure()      
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.coords[:,0],self.coords[:,1],self.coords[:,2], '-o', c='b')
        wp1,wp2,wp3=self.world_pos
        ax.set_xlim(wp1,self.res[0]*self.img.shape[0]+wp1)
        ax.set_ylim(wp2,self.res[1]*self.img.shape[1]+wp2)
        ax.set_zlim(wp3,self.res[2]*self.img.shape[2]+wp3)
        pl.show()

def im_to_points(im):
    '''Converts a binary image to a list of points.'''
    a,b,c=np.where(im>0)
    return zip(a,b,c)
    
def pts_to_coords(pts, res, world_pos):
    '''Converts a list of points to a list of world coordinates.'''
    return np.array(pts) * res + world_pos    
    
def dist(a,b):
    '''Returns the Euclidian distance between two points.'''
    if len(a) != len(b):
        raise ValueError('Input points have an inequal number of coordinates.')
    d=[(b[i] - a[i])**2 for i in range(len(a))]
    return np.sqrt(sum(d))
    
def find_closest(point, pointlist):
    '''Finds the closest point in a given point list to the input point.'''
    d=np.inf
    for p in pointlist:
        nd=dist(point, p)
        if nd<d:
            d=nd
            closest=p
    return closest
    
def curve_length(pointlist, start=None, end=None, return_pts=False):
    '''Traverses all points in a list and measures the cumulative length of 
    the lines connecting each point and its closest neighbour. The endpoint is
    given as a boundary so that the algorithm knows when to stop. Tortuosity
    values can get very large if spurious points past the end point are included.
    
    Inputs:
        pointlist: List or array of x,y,z coordinates. Can be matrix indices or world coordinates.
        start: Tuple of x,y,z coordinates of start point. Algorithm starts here.
        end: Tuple of x,y,z coordinates of end point. Algorithm terminates here.
        return_pts: Returns a list of all points traversed (in order) if True.
        
    Outputs:
        cum_dist: Length of the curve tracked through the input points.
        eucl_dist: Euclidian distance between start and end points.
        fraction_traversed: Fraction of points in the list that were used. If 
                            this fraction is low, there were many points past 
                            the end point and there may have been errors.
        points: (optional): Ordered list of all points that were used.'''
    plist=list(pointlist)#copy list
    if start == None:
        minpoint=np.array(pointlist)[:,2].argmin()
        cpoint=plist[minpoint]
    else:
        cpoint=start

    if end == None:
        endpoint=np.array(pointlist)[:,2].argmax()
        end=plist[endpoint]        
        plist.append(end)
        
    s=np.array(cpoint)
    cum_dist=0
    points=[]

    for x in range(len(plist)):
        repeat=False
        point=find_closest(cpoint, plist)
        fraction_traversed=float(len(points))/(float(len(pointlist)))
        d=dist(cpoint, point)
        if fraction_traversed > .90 and  d >= dist(cpoint, s):
            break        
        cum_dist+=d
        
        if point == end:
            break
        if return_pts:
            points.append(cpoint)
        cpoint=plist.pop(plist.index(point))

    eucl_dist=dist(points[0], points[-1])
    outs=(cum_dist, eucl_dist, fraction_traversed)
    if return_pts:
        outs+=tuple([points])
    return outs
    
def curvature(a,b,c):
    '''cosine rule: AC**2=AB**2 + BC**2 -2*AB*BC*cos(alpha)
    acos((AC**2-AB**2-BC**2)/-2*AB*BC)'''
    a=tuple(a)
    b=tuple(b)
    c=tuple(c)
    d=(dist(a,c)**2-dist(a,b)**2-dist(b,c)**2)/(-2*dist(a,b)*dist(b,c))  
    
    if d <= -1.0:
        c=math.pi/2
    else:
        c=math.acos(d)
    return math.pi - c
        
def cumul_curvature(pointlist, start, end):
    global p1,p2,p3
    plist=list(pointlist)#copy list
    plist.append(end)
    p1=plist.pop(plist.index(find_closest(start, plist)))
    p2=plist.pop(plist.index(find_closest(p1, plist)))
    cum_curv=0
    points=[]
    for x in range(len(plist)):
        p3=plist.pop(plist.index(find_closest(p2, plist)))
        if p3 == end:
            break
        p3=np.array(p3)
        print curvature(p1,p2,p3)
        cum_curv+=curvature(p1,p2,p3)
        p1=np.array(p2)
        p2=np.array(p3)

    eucl_dist=dist(start, p3)
    outs=(cum_curv, eucl_dist)
    return outs    

def main():
    os.chdir('C:\Users\u0091609\Desktop\TS_Stroke_Sham_BFC\LMCA\Skeleton')
    v=VesselSkeleton('BFCTS_stroke_sham_32b_ed1_6_1_BM_M4_FS_LMCA')
    a,b,c,d,e=v.tortuosity()
    v.plot3D()
    

if __name__ == '__main__':
    main()