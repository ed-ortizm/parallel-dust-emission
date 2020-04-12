#!/usr/bin/env python3
import numpy as np
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
MIPS1= np.ones(10)
PACS_blue = 2* MIPS1
PACS_red = 3* MIPS1
PACS_green= 4*MIPS1
PSW = 5 * MIPS1
PMW= 6* MIPS1
PLW = 7* MIPS1
ms = np.array([MIPS1,PACS_blue,PACS_red,PACS_green,PSW,PMW,PLW])
ms = ms.T # (nx7)
# alpha
alpha = alpha(fs,ss,ms,conv_factor)
print(alpha.shape)
chi2 = chi2(fs,ss,ms,alpha,conv_factor)
print(chi2.shape)
