#-------------------------------------------------------------------------------
# Name:        module2
# Purpose:
#
# Author:      Kristof
#
# Created:     21/09/2013
# Copyright:   (c) Kristof 2013
# Licence:     <your licence>
#-------------------------------------------------------------------------------

# -*- coding: utf-8 -*-

import pylab as pl
import numpy as np
from itertools import product, combinations
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

#draw a vector

def main():
    fig = pl.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect("equal")
    
    drawPoints(points_pos, ax)
    drawSphere(ax, 1)
    pl.show()

class Arrow3D(FancyArrowPatch):
    '''3D arrow object. Plot it with the drawArrow function.'''
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def drawArrow(vec, ax):
    '''Draws a single vector (vec) on axes object ax.'''
    a = Arrow3D([0,vec[0]],[0,vec[1]],[0,vec[2]], mutation_scale=20, lw=1, arrowstyle="->", color="r", alpha=0.5)
    ax.add_artist(a)

def drawVecs(vectab, ax):
    '''Draws multiple vectors. vectab is a (N,3) matrix of N vectors.'''
    for vec in vectab:
        drawArrow(vec, ax)

def drawSphere(ax, r):
    '''Draws a 3D wireframe sphere with radius r on axes object with its centre
    at position (0,0,0).'''
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j] #j = imaginary number (synonym to i) and specifies amount of points to create.
    x=r*np.cos(u)*np.sin(v)
    y=r*np.sin(u)*np.sin(v)
    z=r*np.cos(v)
    ax.plot_wireframe(x, y, z, color="b", alpha=0.4, linestyle='dashed')

def drawEllipsoid(coefs, center):
    '''doesn't work yet'''
    fig = pl.figure(figsize=pl.figaspect(1))  # Square figure
    ax = fig.add_subplot(111, projection='3d')
    
    # Radii corresponding to the coefficients:
    rx, ry, rz = [1/np.sqrt(coef) for coef in coefs]
    
    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    
    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = center[0] + rx * np.outer(np.cos(u), np.sin(v))
    y = center[1] + ry * np.outer(np.sin(u), np.sin(v))
    z = center[2] + rz * np.outer(np.ones_like(u), np.cos(v))
    
    # Plot:
    ax.plot_wireframe(x, y, z, color='b')
    
    # Adjustment of the axes, so that they all have the same span:
    max_radius = max(rx, ry, rz)
    for axis in 'xyz':
        getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))
        
    drawPoints(points, ax)

drawPoints = lambda vectab, ax: ax.scatter(vectab[:,0], vectab[:,1], vectab[:,2])

if __name__ == '__main__':
    main()