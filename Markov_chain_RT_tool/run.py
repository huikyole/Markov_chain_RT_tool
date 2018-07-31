import pandas as pd
import numpy as np
import subprocess
import pickle 


def change_OD_peak_levels_and_write_csv(OD_original, h1=41, h2=44):
    const_type1=0.00801195404875952
    const_type2=0.00200298851218988
    OD_new = np.zeros(OD_original.shape)
    OD_new[:] = OD_original[:]
    nlevs = h2-h1 # number of levels with the peak values
    OD_new[0:-1, 2] = const_type1
    OD_new[0:-1, 3] = const_type2
    OD_new[h1:h2, 2] = np.repeat(OD_original[41, 2]*3./nlevs, nlevs)
    OD_new[h1:h2, 3] = np.repeat(OD_original[41, 3]*3./nlevs, nlevs)
   
    np.savetxt('Inputs/OD.csv', OD_new, fmt='%.17f', delimiter=',') 
    return OD_new
    
def run_March_main(): 
    p = subprocess.Popen('./run', shell=True)
    p.wait()
    return

output1 = 'Case_Jul31/OD_array_380.pickle'
output2 = 'Case_Jul31/output_compiled_pickle_380.pickle'

OD_original = pd.read_csv('Inputs/OD_original.csv', sep=',', header=None).values
ncase = 30
OD_array = np.zeros([ncase, 45, 4])
output = np.zeros([ncase, 214, 3])
for icase in np.arange(ncase):
    print icase
    h2 = np.random.choice(np.arange(10)+35, 1)
    h1 = h2 - np.random.choice(np.arange(5)+2, 1)
    OD_array[icase,:] = change_OD_peak_levels_and_write_csv(OD_original, h1=h1, h2=h2)
    run_March_main()
    p = subprocess.Popen('cp Output/output.csv Output/output%03d.csv' %(icase), shell=True)
    p.wait()
    output[icase, :] = pd.read_csv('Output/output%03d.csv' %(icase), sep=',', header=None).values
pickle.dump(OD_array, open(output1, 'wb'))
pickle.dump(output, open(output2, 'wb'))
