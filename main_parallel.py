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
def chi2_2(fs,ss,ms,alpha,conv_factor):
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


#Computing the grid of models
# https://stackoverflow.com/questions/5442910/python-multiprocessing-pool-map-for-multiple-arguments
# if __name__ == '__main__':
#     names = ['Brown', 'Wilson', 'Bartlett', 'Rivera', 'Molloy', 'Opie']
#     with multiprocessing.Pool(processes=3) as pool:
#         results = pool.starmap(merge_names, product(names, repeat=2))
#     print(results)
# model = Model(umin= '0.20', umax='1e3', model='MW3.1_60',gamma = 0.3,data = d, filter=filter)
