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
# print(params_comb[1000])
#m = Model(umin= params_comb[1000][0], umax=params_comb[1000][1], model=params_comb[1000][2],gamma = params_comb[1000][3],data = params_comb[1000][4], filter=params_comb[1000][5])
# if __name__ == '__main__':
x_i = 1000
x_f = 1005
params = params_comb[x_i:x_f]
with Pool() as pool:
    models = pool.starmap(Model, params)
    lum_densities = [model.L_density() for model in models]
raw_arrays = []
# I'll need this for the alpha and chi squared
# ['MIPS1', 'PACS_blue', 'PACS_red', 'PACS_green', 'PSW', 'PMW', 'PLW']
MIPS1, PACS_blue, PACS_red, PACS_green, PSW, PMW, PLW= [],[],[],[],[],[],[]
for i in range(len(params)):
    umin   = float(params[i][0])
    umax   = float(params[i][1])
    q_PAH  = model_q_dic[params[i][2]]
    gamma  = params[i][3]
    lum = lum_densities[i][1:] # --> (L_density,band)
    print(lum)
    data = [umin, umax, q_PAH, gamma, lum[0]]
    # when there is no min_max file, do not pay attention to the data
    if lum[0] == 0:
        pass
    else:
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

## Converting monochromatic luminosities to mJy
# I'll store the conversion factor per galaxy and when the time comes,
# I'll use it
gals = Table.read('edgar.fits')
d = np.array([dist for dist in gals['dist']])
conv_factor = 1e-15/(3.086*d)
gals['conv_factor'] = conv_factor

## Computing the chi squared
# First the elements for alpha
fs = np.array([gals[band] for band in p_bands])
fs = fs.T
ss = np.array([gals[band + '_err'] for band in p_bands])
ss = ss.T
# I need to order then per filter as it is ordered in fs and ss
# ['MIPS1', 'PACS_blue', 'PACS_red', 'PACS_green', 'PSW', 'PMW', 'PLW']
ms = np.array([MIPS1,PACS_blue,PACS_red,PACS_green,PSW,PMW,PLW])
ms = ms.T # (nx7)
# alpha
print(type(conv_factor))
alp = alpha(fs,ss,ms,conv_factor)
x2 = chi2(fs,ss,ms,conv_factor,alp)
#best_x2 = x2.min()
# print(alpha.shape,chi2.shape)
# print(alpha)

## fitting in parallel: I'll adapt my alpha and chi2 for each galaxy
# I need a list for galaxy with f,s and conv_fac and the models

params = [(fs[i],ss[i],conv_factor[i],ms) for i in range(fs.shape[0])]
with Pool() as pool:
    alp2 = pool.starmap(alpha_2, params)
# alp2 is a 21 elements list where each element has n alphas
params = [(fs[i],ss[i],conv_factor[i],ms,alp2[i]) for i in range(fs.shape[0])]
with Pool() as pool:
    x2_2 = pool.starmap(chi2_2, params)
print(x2_2)
