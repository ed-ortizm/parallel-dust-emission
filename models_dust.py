import numpy as np
from scipy import interpolate
#https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory
from os import walk # --> dirpath, dirnames, filenames
import matplotlib
import matplotlib.pyplot as plt
#1. It needs to load all the models and keep them in memory
    # This means I'll need to have an array with the parameters that identify
    # a model.
# From README.txt
models = [ 'MW3.1_00', 'MW3.1_10', 'MW3.1_20', 'MW3.1_30', 'MW3.1_40', 'MW3.1_50', 'MW3.1_60', 'LMC2_00', 'LMC2_05', 'LMC2_10', 'smc']
q_PAHs  = [ 0.47, 1.12, 1.77, 2.50, 3.19, 3.90, 4.58, 0.75, 1.49, 2.37, 0.10]
model_q_dic = {models[i]:q_PAHs for i in range(len(models))}
# ftp://ftp.astro.princeton.edu/draine/dust/irem4/

# Class to read the directory structure and pass the data to arrays withing a dictionary

class Data:
    def __init__(self, path='dust_models/'):
        self.path = path
    def txt_files(self):
        # This list store the outputs from walk
        dirpath,dirnames,files = [],[],[]
        # dirpath contains the paths to all the folders that have txt filenames
        # dirnames is a list with only the name of the folders containing the txt filenames
        # files is a list with the names of all the files
        # The map is one to one among each element of these lists, ie, element 1 refer to the same Uumin
        for (a,b,c) in walk(self.path):
            dirpath.append(a)
            dirnames.append(b)
            files.append(c)
        # Doing some make up
        dirpath = dirpath[1:]
        dirnames = dirnames[0]
        files.remove([])
        return dirpath, dirnames,files
    def arr_dat(self,file):
        # This method generates a Numpy array for one txt file of data
        f = open(file,'r')
        vals = []
        # Flag variable
        i = 0
        for line in f:
            if ('lambda' in line):
                i = 1
            elif i == 1:
                i = 2
            elif i == 2:
                val = line.split()
                vals.append(val)
        f.close()
        # if the file is empty, return 3 zeros as a np.array
        if i == 0:
            return np.zeros(3)
        data = np.array(vals)
        data = data.astype(float)
        return data
    def dic_files(self):
        # This method returns a dictionary where the key is the directory path
        # and the values are lists with the files in the corresponding directory.
        dirpath, dirnames,files = self.txt_files()
        all_txts = {dirpath[i]:files[i] for i in range(len(dirpath))}
        return all_txts
    def dic_arr(self):
        dic_arr = {}
        # This method returns a dictionary containing all the data, where a key is 
        # the name of a .txt file and the value is the array containing the relevant
        # data from this file, relevant for the purposes of this project.
        all_txts = self.dic_files()
        for key,values in all_txts.items():
            for val in values:
                arr = self.arr_dat(key + '/' + val)
                dic_arr[val] = arr
        return dic_arr

# Filter handled from last graded pratical, slightly modified for this practical
class Filter_handler():
    def __init__(self,filter):
        self.name = filter
        self.filter = np.loadtxt(filter)
    def lamb_f(self):
        #converting filter wavelength to nm
        lambdas = self.filter[:,0] * 0.1
        # Lines added to make sure interpolate fills with zeros out of the range
        aux1 = np.linspace(lambdas.min() -100,lambdas.min() - 1., 10)
        aux2 = np.linspace(lambdas.max() + 1.,lambdas.max() + 100, 10)
        lambdas = np.concatenate((aux1,lambdas,aux2))
        return lambdas
    # phtons to energy and normalizing the integral of the filter.
    def energy(self):
        # Checking if the filter is defined in photons or energy
        flag = False
        f = open(self.name,'r')
        for line in f:
            if '# photon' in line:
                flag = True
                break
        f.close()
        if flag:
            # Converting photons to energy
            photons = self.filter[:,1]
            # Lines added to make sure interpolate fills with zeros out of the range
            aux = np.zeros(10)
            photons = np.concatenate((aux,photons,aux))
            energies = self.lamb_f() * photons
            norm = np.trapz(energies,self.lamb_f())
            n_energies = energies/norm
            return n_energies
        else:
            energies = self.filter[:,1]
            # Lines added to make sure interpolate fills with zeros out of the range
            aux = np.zeros(10)
            energies = np.concatenate((aux,energies,aux))
            norm = np.trapz(energies,self.lamb_f())
            n_energies = energies/norm
            return n_energies
    def interpolate(self,interval):
        f = interpolate.interp1d(self.lamb_f(),self.energy(),fill_value='extrapolate')
        return f(interval)

def lamb_inter(arr_1,arr_2):
    stack = np.concatenate((arr_1,arr_2))
    # np.unique eliminates the duplicates and returns the array sorted :)
    return np.unique(stack)

# Class to create an emission spectrum for a given combination of umin, umax, model & gamma
class Model:
    def __init__(self, umin, umax, model,gamma,data,filter, model_q_dic = model_q_dic):
        self.umin = umin
        self.umax = umax
        self.model = model
        self.q_PAH = model_q_dic[model]
        if 0 < gamma < 1:
            self.gamma = gamma
        else:
            print("gamma must be a number between 0 and 1")
        self.key_min_min = 'U' + umin + '_' + umin + '_' + model + '.txt'
        self.key_min_max = 'U' + umin + '_' + umax + '_' + model + '.txt'
        self.data = data # This is the data loaded in memory by the last Class
        self.filter = filter # Filter data, cured in the last Class
    # Returning the raw data needed to compute the model
    def raw_model(self):
        return self.key_min_min, self.data.dic_arr()[self.key_min_min], self.key_min_max, self.data.dic_arr()[self.key_min_max]
    # Computing the model (first 4 parameters of the constructor)
    def spectrum(self):
        min_min, j_nu_min_min, min_max, j_nu_min_max = self.raw_model()
        # to be careful about empty txt files
        if (len(j_nu_min_min) == 3):
            print("Impossible to compute the spectra")
            print(min_min + " does not have enough data\n")
            return np.zeros(1), np.zeros(1)
        elif (len(j_nu_min_max) == 3):
            print("Impossible to compute the spectra")
            print(min_max + " does not have enough data\n")
            return np.zeros(1), np.zeros(1)
        else:
            j_nu = (1.-self.gamma)*j_nu_min_min[:,2] + self.gamma*j_nu_min_max[:,2]
            lambdas = j_nu_min_max[:,0] * 1000. # it is in microns, then with this I convert them to nm
            spectrum = j_nu * (0.001/1.66)*4*np.pi*3e7
            spectrum = spectrum/(lambdas*lambdas) # conversion factors for units in W/nm/(Kg of H)
            # This lines of code are to order the data for lambdas in increasing order
            # Otherwise, I found the bolometric integral to be negative
            idx              = np.argsort(lambdas)
            lambdas          = lambdas[idx]
            spectrum         = spectrum[idx]
            return lambdas,spectrum
    def interpolate(self,interval):
        lambdas,spectrum = self.spectrum()
        f = interpolate.interp1d(lambdas,spectrum,fill_value='extrapolate')
        return f(interval)
    def plot_spec(self):
        # it plots the computed spectrum from (1) in a semilog plot (for x), 
        # it shows it and saves it in a PDF file as well.
        lambdas,spectrum = self.spectrum()
        if len(lambdas)== 1:
            print("Impossible to plot the spectrum")
            print("There is no spectral data for this model")
            return None
        plt.semilogx(lambdas,spectrum) # this looks like an emission spectrum
        # plt.loglog do not, neither plt.semilogy
        plt.title(self.key_min_max[:-4] +'_' + str(self.gamma))
        plt.xlabel('$nm$')
        plt.ylabel('$L_{\lambda}$' +' [W/nm/(kg of H)]')
        plt.savefig(self.key_min_max[:-4]+'_' +str(self.gamma)+ '.pdf')
        plt.show()
    def bolometric(self):
        # it returns the bolometric luminosity
        lambdas,spectrum = self.spectrum()
        # Checking if no Data
        if len(lambdas)== 1:
            print("Impossible to compute bolometric luminosity")
            print("There is no spectral data for this model")
            return None
        return np.trapz(spectrum,lambdas)
    def L_density(self):
        # it returns the luminosity density in two formats, per wavelength and per frequency.
        model_lambdas, spectrum = self.spectrum()
        filter_lambdas          = self.filter.lamb_f()
        lambdas                 = lamb_inter(filter_lambdas,model_lambdas)

        filter_energy  = self.filter.interpolate(lambdas)
        spectrum_model = self.interpolate(lambdas)
        TxL = filter_energy * spectrum_model
        lambda2 = 1./np.trapz(filter_energy/(lambdas*lambdas))
        density_w = np.trapz(TxL,lambdas)
        density_f = lambda2*density_w/(3e17)
        return density_w, density_f
