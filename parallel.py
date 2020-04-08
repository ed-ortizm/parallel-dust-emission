#!/usr/bin/env python3
from models_dust import *
from main_parallel import *
## Computing grid of models in parallel

# Lists of physical parameters for a model
data = Data()
dirpath, dirnames, files_array = data.txt_files()
umin, umax = umm(files_array)
# models --> it is already on moddels_dust
gamma = np.logspace(-3.,-0.3,10)

# Using Pool.starmap to manage xxxxxx xxxxxx iterators
# https://stackoverflow.com/questions/5442910/python-multiprocessing-pool-map-for-multiple-arguments
