import pandas as pd
import numpy as np
import subprocess
import pickle 

import yaml

def edit_yaml(yaml_file, lamda, file):
    f = open(yaml_file, 'w')
    data = dict(lamda0=lamda, output_file=file)
    yaml.dump(data, f, default_flow_style=False)
    f.close()
    return

def change_OD_peak_levels_and_write_csv(OD_original, h1=41, h2=44):
    const_type1=0.00801195404875952
    const_type2=0.00200298851218988
    OD_new = np.zeros(OD_original.shape)
    OD_new[:] = OD_original[:]
    nlevs = h2-h1 # number of levels with the peak values
    OD_new[13:-1, 2] = const_type1
    OD_new[13:-1, 3] = const_type2
    OD_new[h1:h2, 2] = np.repeat(OD_original[41, 2]*3./nlevs, nlevs)
    OD_new[h1:h2, 3] = np.repeat(OD_original[41, 3]*3./nlevs, nlevs)
   
    np.savetxt('Inputs/OD.csv', OD_new, fmt='%.17f', delimiter=',') 
    return OD_new
    
def run_March_main(): 
    p = subprocess.Popen('./run', shell=True)
    p.wait()
    return

wavelengths = [380, 445, 935]
#wavelengths = [355, 470, 555, 660, 865]
#wavelengths = [355, 380, 445, 470, 555, 660, 865, 935]
for wavelength in wavelengths:
    ncase = 100
    output1 = 'Case_Jul31/OD_array_%dnm.pickle' %wavelength
    output2 = 'Case_Jul31/output_compiled_pickle_%dnm.pickle' %wavelength

    OD_original = pd.read_csv('Inputs/OD_original.csv', sep=',', header=None).values
    OD_array = np.zeros([ncase, 45, 4])
    output = np.zeros([ncase, 214, 4])
    for icase in np.arange(ncase):
        print icase
        h2 = np.random.choice(np.arange(10)+35, 1)
        h1 = h2 - np.random.choice(np.arange(5)+2, 1)
        OD_array[icase,:] = change_OD_peak_levels_and_write_csv(OD_original, h1=h1, h2=h2)
        edit_yaml('choice2.yaml', wavelength*1.e-9, 'Output/output_%dnm_case%03d.csv' %(wavelength, icase))
        run_March_main()
        output[icase, :] = pd.read_csv('Output/output_%dnm_case%03d.csv' %(wavelength, icase), sep=',', header=None).values
    pickle.dump(OD_array, open(output1, 'wb'))
    pickle.dump(output, open(output2, 'wb'))
