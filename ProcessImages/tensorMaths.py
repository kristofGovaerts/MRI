# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 19:07:34 2013

@author: Gebruiker
"""
from __future__ import division
import numpy as np
import numpy.linalg as la

    
def calculate_b_matrix(bvec, b):
    """Calculates a b-matrix when given a table of diffusion directions and
    corresponding b-values.
    
    Inputs: bvec -- nx3-dimensional array of diffusion gradients where n=amount of diffusion dirs
            b    -- n-dimensional array of b-values where n=amount of diffusion dirs
            
    Output: bmat -- nx6-dimensional b-matrix in the order of:
                    b*Gxx, 2b*Gxy, 2b*Gxz, b*Gyy, 2b*Gyz, b*Gzz"""
    bmat=np.zeros((len(b),6))
    
    bmat[:,0]=b * bvec[:,0]**2
    bmat[:,1]=b * 2 * bvec[:,0] * bvec[:,1]
    bmat[:,2]=b * 2 * bvec[:,0] * bvec[:,2]
    bmat[:,3]=b * bvec[:,1]**2
    bmat[:,4]=b * 2 * bvec[:,1] * bvec[:,2]
    bmat[:,5]=b * bvec[:,2]**2
    return bmat

def construct_bvecs(self):
    bvecmat=np.zeros([self.nA0+self.ndirs*self.nbvals, 3])
    for i in range(self.ndirs):
        bvecmat[self.nA0+self.nbvals*i:self.nA0+self.nbvals*(i+1),:]=self.dwdir[i] #fills b-vector matrix with the different diffusion dirs
    return bvecmat        

def calculate_diffusion_tensor(img, bmat):
    '''ln(S/S0) * 1/b_mat_element = D_element'''
    li=[]
    for ind, vec in enumerate(bmat[1:,:]):
        lns=np.log(img[:,:,:,ind]/img[:,:,:,0]) 
        d=[lns*1/(bmat[ind,i]) for i in range(6)]
        li.append(d)
    return li
        

def mvee(points, tol = 0.001, mevals=10000):
    """
    Find the minimum volume ellipse.
    Return A, c where the equation for the ellipse given in "center form" is
    (x-c).T * A * (x-c) = 1
    
    This ellipsoid is in principle equivalent to the diffusion tensor (although not really)
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
    
def unit_vector(vecmat):
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
    
def calculate_tensor(inmat, bvecs, upper = None):
    '''Calculates the tensor based on the input array and a matrix of directional
    vectors by fitting an ellipsoid to one tensor per voxel, with the ADCs corresponding
    to the magnitude of diffusion in every direction.
    
    Inputs:
        inmat = Input matrix. 4D matrix with dimensions (x,y,z,d), where d are 
        the diffusion directions.
        bvecs = (dx3) direction matrix.
    Outputs:
        evals = Eigenvalues of the calculated tensor.
        evecs = Centroid of the tensor.'''
    x,y,z,d=inmat.shape
    it=x*y*z
    evals=np.zeros((x,y,z,3))  
    evecs=np.array(evals)
    
    for i in range(it):
        x1,y1,z1=np.unravel_index(i, (x,y,z))
        if z1 == 7:
            points=inmat[x1,y1,z1,:]
            if upper == None:
                upper = np.max(points)
            points_pos=np.array([points[i]*bvecs[i,:]  for i in range(d) if points[i] > 0 and points[i] <= upper])
            if len(points_pos) > 5:
                ten, evec, nev=mvee(points_pos, tol=0.01, mevals=1000)
                evals[x1,y1,z1,:]=np.sort(np.diag(ten))[::-1]
                evecs[x1,y1,z1,:]=np.sort(evec)
        
    return evals, evecs
    
#def lsqnonneg(C, d, x0=None, tol=None, itmax_factor=3):
#    '''Linear least squares with nonnegativity constraints
#
#    (x, resnorm, residual) = lsqnonneg(C,d) returns the vector x that minimizes norm(d-C*x)
#    subject to x >= 0, C and d must be real
#    '''
#
#    eps = 2.22e-16    # from matlab
#    def norm1(x):
#        return abs(x).sum().max()
#
#    def msize(x, dim):
#        s = x.shape
#        if dim >= len(s):UnitVectors;
#Mprime=321;
#g=g(1:Mprime,:)
#for i=0:order
#    for j=0:order-i
#        pop(i+1,j+1,order-i-j+1)=population(i,j,order-i-j,order);
#    end
#end
#for k=1:length(g)
#    c=1;
#    for i=0:order
#		for j=0:order-i
#			C(k,c)=pop(i+1,j+1,order-i-j+1)*(g(k,1)^i)*(g(k,2)^j)*(g(k,3)^(order-i-j));
#			c=c+1;
#        end
#    end
#end
#
#            return 1
#        else:
#            return s[dim]
#
#    if tol is None:
#        tol = 10*eps*norm1(C)*(max(C.shape)+1)
#
#    C = np.asarray(C)
#
#    (m,n) = C.shape
#    P = np.zeros(n)
#    Z = np.arange(1, n+1)
#
#    if x0 is None:
#        x=P
#    else:
#        if any(x0 < 0):
#            x=P
#        else:
#            x=x0
#
#    ZZ=Z
#
#    resid = d - np.dot(C, x)
#    w = np.dot(C.T, resid)
#
#    outeriter=0
#    it=0
#    itmax=itmax_factor*n
#    exitflag=1
#
#    # outer loop to put variables into set to hold positive coefficients
#    while np.any(Z) and np.any(w[ZZ-1] > tol):
#        outeriter += 1
#
#        t = w[ZZ-1].argmax()
#        t = ZZ[t]
#
#        P[t-1]=t
#        Z[t-1]=0
#
#        PP = np.where(P <> 0)[0]+1
#        ZZ = np.where(Z <> 0)[0]+1
#
#        CP = np.zeros(C.shape)
#
#        CP[:, PP-1] = C[:, PP-1]
#        CP[:, ZZ-1] = np.zeros((m, msize(ZZ, 1)))
#
#        z=np.dot(np.linalg.pinv(CP), d)
#
#        z[ZZ-1] = np.zeros((msize(ZZ,1), msize(ZZ,0)))
#
#        # inner loop to remove elements from the positve set which no longer belong
#        while np.any(z[PP-1] <= tol):
#            if it > itmax:
#                max_error = z[PP-1].max()
#                raise Exception('Exiting: Iteration count (=%d) exceeded\n Try raising the tolerance tol. (max_error=%d)' % (it, max_error))
#
#            it += 1
#
#            QQ = np.where((z <= tol) & (P <> 0))[0]
#            alpha = min(x[QQ]/(x[QQ] - z[QQ]))
#            x = x + alpha*(z-x)
#
#            ij = np.where((abs(x) < tol) & (P <> 0))[0]+1
#            Z[ij-1] = ij
#            P[ij-1] = np.zeros(max(ij.shape))
#            PP = np.where(P <> 0)[0]+1
#            ZZ = np.where(Z <> 0)[0]+1
#
#            CP[:, PP-1] = C[:, PP-1]
#            CP[:, ZZ-1] = np.zeros((m, msize(ZZ, 1)))
#
#            z=np.dot(np.linalg.pinv(CP), d)
#            z[ZZ-1] = np.zeros((msize(ZZ,1), msize(ZZ,0)))
#
#        x = z
#        resid = d - np.dot(C, x)
#        w = np.dot(C.T, resid)
#
#    return (x, sum(resid * resid), resid)