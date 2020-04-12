#!/usr/bin/env python3
import numpy as np
import ctypes
from multiprocessing import Pool
from itertools import product
from astropy.table import Table
def alpha(fs,ss,ms,conv_factor):
    alp = np.zeros((fs.shape[0],ms.shape[0]))
    for i in range(fs.shape[0]): # looping on all the galaxies
        ms = ms*conv_factor[i] # converting to mJy
        #print(ss[i,:])
        u = (fs[i,:]/(ss[i,:]*ss[i,:]))*ms # broadcasting: (1x7)*(nx7) = nx7
        u = u.sum(axis=1) # n
        d = (1/(ss[i,:]*ss[i,:]))*(ms*ms)
        d = d.sum(axis=1) # n
        alp[i] = u/d
    # 21 galaxies where each row is the alpha array as a consequence of n models
    return alp
def chi2(fs,ss,ms,alpha,conv_factor):
    xx = np.zeros((fs.shape[0],ms.shape[0]))
    for i in range(fs.shape[0]):
        u1 = fs[i,:] # (1x7)
        ms = ms*conv_factor[i] # converting to mJy
        u2 = alpha[i,:]*ms.T # alpha[i,:] = (1xn), ms.shape = (nx7)
        u2 = u2.T # u2.shape = (nx7)
        u = u1 - u2 # (nx7)
        u = u*u #(nx7)
        d = ss[i,:] # (1x7)
        d = d*d # (1x7)
        x = u/d # (nx7)
        x = x.sum(axis=1) # n
        xx[i] = x
    return xx

def alpha_2(fs,ss,ms,conv_factor):
    ms = ms*conv_factor # converting to mJy
    #print(ss[i,:])
    u = (fs/(ss*ss))*ms # broadcasting: (1x7)*(nx7) = nx7
    u = u.sum(axis=1) # n
    d = (1/(ss*ss))*(ms*ms)
    d = d.sum(axis=1) # n
    alp = u/d
    # 21 galaxies where each row is the alpha array as a consequence of n models
    return alp
def chi2_2(fs,ss,ms,conv_factor, alpha):
    u1 = fs # (1x7)
    ms = ms*conv_factor # converting to mJy
    u2 = alpha*ms.T # alpha[i,:] = (1xn), ms.shape = (nx7)
    u2 = u2.T # u2.shape = (nx7)
    u = u1 - u2 # (nx7)
    u = u*u #(nx7)
    d = ss # (1x7)
    d = d*d # (1x7)
    x = u/d # (nx7)
    x = x.sum(axis=1) # n
    return x


p_bands = ['MIPS1','PACS_blue','PACS_green','PACS_red','PLW','PMW','PSW']
data = Table.read('edgar.fits')
d = np.array([dist for dist in data['dist']])
conv_factor = 1e-15/(3.086*d)
data['conv_factor'] = conv_factor

## Computing the chi squared
# First the elements for alpha
fs = np.array([data[band] for band in p_bands])
fs = fs.T
ss = np.array([data[band + '_err'] for band in p_bands])
ss = ss.T
# I need to order then per filter as it is ordered in fs and ss
# ['MIPS1', 'PACS_blue', 'PACS_red', 'PACS_green', 'PSW', 'PMW', 'PLW']
ms = mock_ms(n=15)
## fitting in parallel: I'll adapt my alpha and chi2 for each galaxy
# I need a list for galaxy with f,s and conv_fac and the models

params = [(fs[i],ss[i],conv_factor[i],ms) for i in range(fs.shape[0])]
with Pool() as pool:
    alp2 = pool.starmap(alpha_2, params)
# alp2 is a 21 elements list where each element has n alphas
params = [(fs[i],ss[i],conv_factor[i],ms,alp2[i]) for i in range(fs.shape[0])]
with Pool() as pool:
    x2_2 = pool.starmap(chi2_2, params)
print(len(alp2[0]),len(x2_2[5]))
