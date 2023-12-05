import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import matplotlib.pyplot as plt
import pandas as pd
import time
import scipy
import scipy.optimize
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import gridspec


save_plots = True
fig, ax = plt.subplots(1, 1, figsize=(9, 5))
name = 'MBR_BurkeSong'

models = {
        #   'mevel':'D:\\Research\\Models\\Mevel\\mevel.cti',
          'LMR-R':"test/data/alzuetamechanism_LMRR.yaml",            
          'Ar':"test/data/alzuetamechanism_LMRR_allAR.yaml",
          r'H$_2$O':"test/data/alzuetamechanism_LMRR_allH2O.yaml",
          }

colors = ['xkcd:purple','r','b']

for i, m in enumerate(list(models.keys())):
    p_list = np.linspace(0,20,5)[1:]
    T = 300.0  # unburned gas temperature [K]

    reactants = 'H2:0.1071, O2:0.1785, He:0.7144'  # premixed gas composition
    width = 0.03  # m
    loglevel = 1  # amount of diagnostic output (0 to 8)
    mbr = []
    for p in p_list:
        gas = ct.Solution(list(models.values())[i])
        gas.TPX = T, p*ct.one_atm, reactants
        f = ct.FreeFlame(gas, width=width)
        f.set_refine_criteria(ratio=2, slope=0.03, curve=0.06)
        f.transport_model = 'mixture-averaged'
        f.solve(loglevel=loglevel, auto=True)
        mbr.append(f.velocity[0]*f.density[0] / 10) # g/cm2*s
    ax.semilogy(p_list, mbr, '-', linestyle='solid', color=colors[i], label=m)


path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"
dataset0 = pd.read_csv(path+'5 FS H2O (Burke)\\black.csv')
ax.plot(dataset0.iloc[:, 0],dataset0.iloc[:, 1],color='k', zorder=2, fillstyle='full', linestyle = 'solid', label="Keromnes et al.")
dataset1 = pd.read_csv(path+'5 FS H2O (Burke)\\exp_pts.csv')
ax.plot(dataset1.iloc[:, 0],dataset1.iloc[:, 1],marker='s',color='k', markersize=6, zorder=3, fillstyle='full', linestyle = 'None', label="Burke et al.")


ax.legend(fontsize=15, frameon=False)#, loc='upper right')  
ax.set_ylabel(r'Mass burning rate (g cm$^\mathdefault{-2}$ s$^\mathdefault{-1}$)', fontsize=18)
ax.set_xlabel(r'Pressure (atm)', fontsize=18)
ax.tick_params(axis='both', direction="in", labelsize=15)
ax.tick_params(axis='both', which='minor', direction="in")#, bottom=False)
ax.set_ylim([0.00,0.04999])


if save_plots == True:
    plt.savefig(name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig(name+'.png', dpi=1000, bbox_inches='tight')
plt.show()     