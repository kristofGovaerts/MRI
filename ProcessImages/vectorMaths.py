# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 19:07:34 2013

@author: Gebruiker
"""
from __future__ import division
import numpy as np
import numpy.linalg as la

def mvee(points, tol = 0.001, mevals=10000):
    """
    Find the minimum volume ellipse.
    Return A, c where the equation for the ellipse given in "center form" is
    (x-c).T * A * (x-c) = 1
    
    This ellipsoid is in principle equivalent to the diffusion tensor.
    Note that the tolerance should be a smaller number when calculating values
    in the ADC range (around 0.00001)!
    """
    points = np.asmatrix(points)
    N, d = points.shape
    Q = np.column_stack((points, np.ones(N))).T #transpose to vertical
    err = tol+1.0
    u = np.ones(N)/N
    nevals=0
    for i in range(mevals):
        # assert u.sum() == 1 # invariant
        nevals += 1
        X = Q * np.diag(u) * Q.T
        M = np.diag(Q.T * la.inv(X) * Q)
        jdx = np.argmax(M)
        step_size = (M[jdx]-d-1.0)/((d+1)*(M[jdx]-1.0))
        new_u = (1-step_size)*u
        new_u[jdx] += step_size
        err = la.norm(new_u-u)
        u = new_u
        if err < tol: 
            break
    c = u*points
    A = la.inv(points.T*np.diag(u)*points - c.T*c)/d    
    return np.asarray(A), np.squeeze(np.asarray(c)), nevals
    
def unitVector(vecmat):
    '''Normalizes 2D input array of vectors into unit vectors (e.g. vectors
    with length 1). This is necessary because this is the only type of
    vector Dipy accepts, and with good reason, as the diffusion tensor is in 
    principle defined as the product of the b-vector matrix and the ADC values.
    '''
    output=np.zeros(vecmat.shape)
    for i, vec in enumerate(vecmat):
        unitvec=vec/np.linalg.norm(vec)
        output[i,:]=unitvec
    output[np.isnan(output)] = 0
    return output
    
def calculate_tensor(inmat, bvecs):
    x,y,z,d=inmat.shape
    it=x*y*z
    evals=np.zeros((x,y,z,3))  
    evecs=np.array(evals)
    
    for i in range(it):
        x1,y1,z1=np.unravel_index(i, (x,y,z))
        points=inmat[x1,y1,z1,:]
        points_pos=np.array([points[i]*bvecs[i,:]  for i in range(d) if points[i] > 0])
        ten, evec=mvee(points_pos, mevals=1000)
        evals[x1,y1,z1,:]=np.sort(np.diag(ten))[::-1]
        evecs[x1,y1,z1,:]=evec
        
    return evals, evecs
        