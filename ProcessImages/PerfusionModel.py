# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 23:37:32 2013

@author: Gebruiker
"""

import numpy as np
import pylab as pl
import scipy
import wx

#deltaM(t)=2M0 * f * tau * alpha * exp(-t/T1bl) * qp(t)
# M0 = equilibrium magnetization of blood
# f = local perfusion
# tau = transit time
# alpha = tagging efficiency

m0 = 1.0
cbf = 0.8
tau = 2000.0
alpha = 0.95
t1bl = 2400.0
qpt= 0.99
lambd=83
t1t=5000.0
dt=500.0
k=(1/t1bl)-(1/t1t + cbf/lambd)

arange=np.arange(0,4000, dtype = "float")

def aslfunc1(t, m0, cbf, tau, alpha, t1bl, qpt):   
    return 2*m0 * cbf * (t - dt) * alpha * np.exp(-t/t1bl) * qp(t)

def aslfunc2(t, m0, cbf, tau, alpha, t1bl, qpt):
    return 2 * m0 * cbf * tau * alpha * exp(-t/t1bl) * qp(t)
    
def qp(t):
    return (exp(k*t)*(exp(-k*dt)-exp(-k*t)))/k*tau

yvals=[]    
for t in arange:    
    if dt <= t and t <= tau+dt: 
        yvals.append(aslfunc1(t, m0, cbf, tau, alpha, t1bl, qpt))
    elif tau + dt <= t:
        yvals.append(aslfunc2(t, m0, cbf, tau, alpha, t1bl, qpt))
    else:
        yvals.append(0)
    
pl.plot(arange, yvals)

pl.show()