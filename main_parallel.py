from models_dust import *
from multiprocessing import Pool
from itertools import product
from astropy.table import Table
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
# class Fit:
#     def __init__(self,fs,ss,ms):
#         self.fs = fs
#         self.ss = ss
#         self.ms = ms
#     def alpha(self):
#         alp = 0
#         return alp
#     def chi2(self):
#         x =0
#         return x

def alpha(fs,ss,ms):
    alp = np.zeros((21,ms.shape[1]))
    for i in range(21):
        u = (fs[i,:]/(ss[i:]*ss[i:]))*ms # broadcasting: (1x7)*(nx7) = nx7
        u = u.sum(axis=1) # n
        d = (1/(ss[i,:]*ss[i,:]))*(ms*ms)
        d = d.sum(axis=1) # n
        alp[i] = u/d
    # 21 galaxies where each row is the alpha array as a consequence of n models
    return alp
def chi2(fs,ss,ms,alpha):
    x = np.zeros(())
    return x

#Computing the grid of models
# https://stackoverflow.com/questions/5442910/python-multiprocessing-pool-map-for-multiple-arguments
# if __name__ == '__main__':
#     names = ['Brown', 'Wilson', 'Bartlett', 'Rivera', 'Molloy', 'Opie']
#     with multiprocessing.Pool(processes=3) as pool:
#         results = pool.starmap(merge_names, product(names, repeat=2))
#     print(results)
# model = Model(umin= '0.20', umax='1e3', model='MW3.1_60',gamma = 0.3,data = d, filter=filter)
