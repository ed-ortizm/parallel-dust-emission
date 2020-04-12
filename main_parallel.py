from models_dust import *
from multiprocessing import Pool
from itertools import product
from astropy.table import Table

class Fit:
    def __init__(self,fs,ss,ms,cf):
        self.fs = fs
        self.ss = ss
        # Here we convert the model flux to mJy
        self.ms = ms * cfs
    def alpha_2(self):
        # Eq. 13 from the paper
        u = (self.fs/(self.ss*self.ss))*self.ms # broadcasting: (1x7)*(nx7) = nx7
        u = u.sum(axis=1) # n
        d = (1/(self.ss*self.ss))*(self.ms*self.ms)
        d = d.sum(axis=1) # n
        alp = u/d
        # an n= self.ms.shape[0] dimensional array is returned
        return alp
    def chi2_2(self):
        # Eq. 14 from the paper
        alpha = self.alpha_2()
        u1 = self.fs # (1x7)
        u2 = self.alpha*self.ms.T # alpha = (n), ms.shape = (nx7)
        u2 = u2.T # u2.shape = (nx7)
        u = u1 - u2 # (nx7)
        u = u*u #(nx7)
        d = self.ss # (1x7)
        d = d*d # (1x7)
        x = u/d # (nx7)
        x = x.sum(axis=1) # n
        # an n= self.ms.shape[0] dimensional array is returned
        return x

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

def mock_ms(n=10):
    # mock model's flux data
    MIPS1      = 1e-24*np.ones(n)
    PACS_blue  = 2.* MIPS1
    PACS_red   = 3.* MIPS1
    PACS_green = 4.*MIPS1
    PSW        = 5.*MIPS1
    PMW        = 6.*MIPS1
    PLW        = 7.*MIPS1
    ms = np.array([MIPS1,PACS_blue,PACS_red,PACS_green,PSW,PMW,PLW])
    ms = ms.T # (nx7)
    return ms
## This ones do the computation with out parallelization
# def alpha(fs,ss,ms,conv_factor):
#     alp = np.zeros((fs.shape[0],ms.shape[0]))
#     for i in range(fs.shape[0]): # looping on all the galaxies
#         ms = ms*conv_factor[i] # converting to mJy
#         #print(ss[i,:])
#         u = (fs[i,:]/(ss[i,:]*ss[i,:]))*ms # broadcasting: (1x7)*(nx7) = nx7
#         u = u.sum(axis=1) # n
#         d = (1/(ss[i,:]*ss[i,:]))*(ms*ms)
#         d = d.sum(axis=1) # n
#         alp[i] = u/d
#     # 21 galaxies where each row is the alpha array as a consequence of n models
#     return alp
# def chi2(fs,ss,ms,conv_factor,alpha):
#     xx = np.zeros((fs.shape[0],ms.shape[0]))
#     for i in range(fs.shape[0]):
#         u1 = fs[i,:] # (1x7)
#         ms = ms*conv_factor[i] # converting to mJy
#         u2 = alpha[i,:]*ms.T # alpha[i,:] = (1xn), ms.shape = (nx7)
#         u2 = u2.T # u2.shape = (nx7)
#         u = u1 - u2 # (nx7)
#         u = u*u #(nx7)
#         d = ss[i,:] # (1x7)
#         d = d*d # (1x7)
#         x = u/d # (nx7)
#         x = x.sum(axis=1) # n
#         xx[i] = x
#     return xx
