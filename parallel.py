#!/usr/bin/env python3
from models_dust import *
from main_parallel import *
import ctypes
from multiprocessing.sharedctypes import RawArray
## Computing grid of models in parallel

## Lists of physical parameters for a model
data = Data()
dirpath, dirnames, files_array = data.txt_files()
umins, umaxs = umm(files_array)
# models --> it is already on moddels_dust
gammas = np.logspace(-3.,-0.3,10)
p_bands = ['MIPS1','PACS_blue','PACS_green','PACS_red','PLW','PMW','PSW']
bands = ['./filters/MIPS1.dat','./filters/PACS_blue.dat','./filters/PACS_green.dat',\
'./filters/PACS_red.dat','./filters/PLW.dat','./filters/PMW.dat','./filters/PSW.dat']
filters = [Filter_handler(band) for band in bands]
## Permutations of the physical parameters to generate the tuples for starmap
# https://stackoverflow.com/questions/798854/all-combinations-of-a-list-of-lists
params      = [umins,umaxs,models,gammas,[data],filters]
params_comb = list(product(*params))
# print(len(params_comb))
# 716100 combinations: mamma mia

## Using Pool.starmap to manage the sub-processes
# https://stackoverflow.com/questions/5442910/python-multiprocessing-pool-map-for-multiple-arguments

x_i = 1000
x_f = 1002
params = params_comb[x_i:x_f]
with Pool() as pool:
    models = pool.starmap(Model, params)
    lum_densities = [model.L_density() for model in models]
raw_arrays = []

# I'll need this for the alpha and chi squared
MIPS1, PACS_blue, PACS_red, PACS_green, PSW, PMW, PLW= [],[],[],[],[],[],[]
for i in range(len(params)):
    umin   = float(params[i][0])
    umax   = float(params[i][1])
    q_PAH  = model_q_dic[params[i][2]]
    gamma  = params[i][3]
    lum = lum_densities[i][1:] # --> (L_density,band)
    data = [umin, umax, q_PAH, gamma, lum[0]]
    # when there is no min_max file, do not pay attention to the data
    if lum[0] != 0:
        raw_arrays.append(RawArray(ctypes.c_float,data))
        if lum[1]   == 'MIPS1':
            MIPS1.append(lum[0])
        elif lum[1] == 'PACS_blue':
            PACS_blue.append(lum[0])
        elif lum[1] == 'PACS_red':
            PACS_red.append(lum[0])
        elif lum[1] == 'PACS_green':
            PACS_green.append(lum[0])
        elif lum[1] == 'PSW':
            PSW.append(lum[0])
        elif lum[1] =='PMW':
            PMW.append(lum[0])
        else:
            PLW.append(lum[0])

# ms corresponds to the fluxes in the computed models. I need to order then per
# filter as it is ordered in fs and ss

#ms = np.array([MIPS1,PACS_blue,PACS_red,PACS_green,PSW,PMW,PLW])
#ms = ms.T # (nx7)
### My program is very slow, I'll run on mock data.
ms = mock_ms(n=15)

## Converting monochromatic luminosities to mJy
# I'll store the conversion factor per galaxy and when the time comes,
# I'll use it directly in the Fit class
gals = Table.read('edgar.fits')
d = np.array([dist for dist in gals['dist']])
cfs = 1e-15/(3.086*d)
gals['conv_factors'] = cfs

## Computing the chi squared
# First the elements to initialize the Fit class
# the observed flux
fs = np.array([gals[band] for band in p_bands])
fs = fs.T
# the error associated to the observed flux
ss = np.array([gals[band + '_err'] for band in p_bands])
ss = ss.T
## Fitting in parallel: I'll adapt my alpha and chi2 for each galaxy
# I need a list for galaxy with f,s and cf and the models

params = [(fs[i],ss[i],cfs[i],ms) for i in range(fs.shape[0])]
with Pool() as pool:
    fits = pool.starmap(Fit, params)
# I get a list of 21 elements (galaxies) where each element has n (models) alphas and
# n chi sqared
    alps  = [fit.alpha_2() for fit in fits]
    chi2s = [fit.chi2_2() for fit in fits]

print(len(alps),len(chi2s))
print(len(alps[0]),len(chi2s[0]))

# smallest chi sqare
gals_chi = np.zeros(fs.shape[0])
gals_alp = np.zeros(fs.shape[0])
for i in range(len(chi2s)) :
    gals_chi = chi2s[i].min()
    gals_alp = alps[i]
