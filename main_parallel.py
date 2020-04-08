from models_dust import *

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
# https://stackoverflow.com/questions/5442910/python-multiprocessing-pool-map-for-multiple-arguments
