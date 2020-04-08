#!/usr/bin/env python3
from models_dust import *
from main_parallel import *
## Computing grid of models in parallel

# Lists of physical parameters for a model
data = Data()
dirpath, dirnames, files_array = data.txt_files()
umins, umaxs = umm(files_array)
# models --> it is already on moddels_dust
gammas = np.logspace(-3.,-0.3,10)

# Permutations of the physical parameters to generate the tuples for starmap
# https://stackoverflow.com/questions/798854/all-combinations-of-a-list-of-lists
params      = [umins,umaxs,models,gammas]
params_comb = list(product(*params))
# 102300 combinations: mamma mia
# Using Pool.starmap to manage the sub-processes

# https://stackoverflow.com/questions/5442910/python-multiprocessing-pool-map-for-multiple-arguments
#
# def f(x,y):
#     return x*y
# if __name__ == '__main__':
#     with Pool() as pool:
#         results = pool.starmap(f, [(1,2),(3,4),(5,6)])
#     print(results)
