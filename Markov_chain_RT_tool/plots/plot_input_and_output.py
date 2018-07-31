import pickle
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

OD_original = pd.read_csv('../Inputs/OD_original.csv', header=None).values
output_original = pd.read_csv('../Output/output.csv', header=None).values

OD = pickle.load(open('../Case_Jul31/OD_array_355.pickle'))
output = pickle.load(open('../Case_Jul31/output_compiled_355.pickle'))

figure1 = 'Random_OD_profiles_355'
figure2 = 'Random_IQ_output_355'

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ncase = OD.shape[0]
y = OD_original[:,0]
for icase in np.arange(ncase):
    ax.plot(OD[icase,:,2], y, lw=0.3, color='grey')
plot = ax.plot(OD_original[:,2], OD_original[:,0], lw=3)
ax.set_ylim(max(OD_original[:,0]), 10)
ax.set_xlabel('OD for Type 1')
ax.set_ylabel('Layer')
fig.savefig(figure1, dpi=600, bbox_inches='tight')

fig = plt.figure()
ax = fig.add_subplot(2, 1, 1)
ncase = output.shape[0]
x = output_original[:,0]
for icase in np.arange(ncase):
    ax.plot(x, output[icase,:,1], lw=0.1, color='grey')
    ax.plot(x, output[icase,:,2], lw=0.1, color='grey')
plot = ax.plot(x, output_original[:,1], lw=0.5, color='blue', label='I')
plot = ax.plot(x, output_original[:,2], lw=0.5, color='red', label = '-10Q')
ax.set_xlim([-90, 90])
ax.set_xlabel('Theta (deg)')
ax.set_ylabel(r'$\pi S/\mu_{0}$')
'''
ax = fig.add_subplot(2, 1, 2)
for icase in np.arange(ncase):
    ax.plot(x, output[icase,:,1]-output_original[:,1], lw=0.1, color='blue')
    ax.plot(x, output[icase,:,2]-output_original[:,2], lw=0.1, color='red')
    print (output[icase,:,1]-output_original[:,1]).max(), (output[icase,:,1]-output_original[:,1]).min()
    print (output[icase,:,2]-output_original[:,2]).max(), (output[icase,:,2]-output_original[:,2]).min()
ax.set_xlim([-90, 90])
ax.set_xlabel('Theta (deg)')
ax.set_ylabel(r'$\pi S/\mu_{0}$ difference')
'''
fig.savefig(figure2, dpi=600, bbox_inches='tight')
