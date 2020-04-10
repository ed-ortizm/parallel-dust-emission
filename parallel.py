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
if __name__ == '__main__':
    x_i = 1000
    x_f = 1002
    params = params_comb[x_i:x_f]
    with Pool() as pool:
        models = pool.starmap(Model, params)
        lum_densities = [model.L_density() for model in models]
    raw_arrays = []
    for i in range(len(params)):
        umin   = float(params[i][0])
        umax   = float(params[i][1])
        q_PAH  = model_q_dic[params[i][2]]
        gamma  = params[i][3]
        lum = lum_densities[i][1]
        data = [umin, umax, q_PAH, gamma, lum]
        # when there is no min_max file, do not pay attention to the data
        if lum == 0:
            pass
        else:
            raw_arrays.append(RawArray(ctypes.c_float,data))
