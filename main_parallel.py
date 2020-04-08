from models_dust import *
from multiprocessing import Pool
from itertools import product
# Function to obtain umin, umax lists
def umm(files_array):
    umin, umax = [], []
    for files in files_array:
        for file in files:
            umin.append(file.split('_')[0][1:])
            umax.append(file.split('_')[1])
    # removing duplicates
    umin = set(umin)
    umin.remove('U1.00')
    umax = set(umax)
    return umin, umax

#Computing the grid of models
# https://stackoverflow.com/questions/5442910/python-multiprocessing-pool-map-for-multiple-arguments
# if __name__ == '__main__':
#     names = ['Brown', 'Wilson', 'Bartlett', 'Rivera', 'Molloy', 'Opie']
#     with multiprocessing.Pool(processes=3) as pool:
#         results = pool.starmap(merge_names, product(names, repeat=2))
#     print(results)
# model = Model(umin= '0.20', umax='1e3', model='MW3.1_60',gamma = 0.3,data = d, filter=filter)
