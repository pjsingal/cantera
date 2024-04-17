
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

import sys, os
import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
import matplotlib as mpl
mpl.rc('font',family='Times New Roman')
mpl.rcParams['mathtext.fontset'] = 'stix'
from matplotlib.legend_handler import HandlerTuple

save_plots = True
fig, ax = plt.subplots(1, 1, figsize=(3.8, 2))


models = {
          'LMR-R':"test/data/alzuetamechanism_LMRR.yaml", 
          'Alzueta':"test/data/alzuetamechanism.yaml",                 
          'Ar':"test/data/alzuetamechanism_LMRR_allAR.yaml",
          r'H$_2$O':"test/data/alzuetamechanism_LMRR_allH2O.yaml",
          }


import csv

# Function to save data to CSV
def save_to_csv(filename, data):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)


# p_list = [50,100,250,760,1500]
p_list=[760]
fuel_list = np.linspace(0.14,0.5,15) #fuel mole fractions

colors = ['xkcd:purple',"xkcd:grey",'r','b']
lines = ['-','-','--',':',(0, (5, 10)),(0, (3, 5, 1, 5, 1, 5))]

Tin = 300.0  # unburned gas temperature [K]
width = 0.03  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

alpha_list = [1.0,0.8,0.6,0.4,0.2,0.0]
a_st = [0.75,0.7,0.65,0.6,0.55,0.5]

# alpha_list = [0]
# a_st = [0.5]

for x, alpha in enumerate(alpha_list):
    for k, m in enumerate(models):
        for i, p in enumerate(p_list):
            mbr = []
            phi_list = []
            for j, fuel_frac in enumerate(fuel_list):
                gas = ct.Solution(list(models.values())[k])
                NH3 = alpha*fuel_frac
                H2 = (1-alpha)*fuel_frac
                ox_frac = 1 - fuel_frac # oxidizer fraction
                O2 = ox_frac*0.21
                N2 = ox_frac*0.79
                phi = np.divide(fuel_frac/O2,1/a_st[x]) # THIS STEP MUST BE REVIEWED
                phi_list.append(phi)
                X = {'NH3':NH3,'H2':H2,'O2':O2,'N2':N2}
                gas.TPX = Tin, (p/760)*ct.one_atm, X
                f = ct.FreeFlame(gas, width=width)
                f.set_refine_criteria(ratio=3, slope=0.06/2, curve=0.12/2)
                # f.transport_model = 'mixture-averaged'
                f.transport_model = 'multicomponent'
                f.solve(loglevel=loglevel, auto=True)
                mbr.append(f.velocity[0] * 100) # cm/s
            # ax.plot(phi_list,mbr,label=m, color=colors[k],linestyle=lines[i])

            # Save phi_list and mbr to CSV
            csv_filename = f'RonneyResults_multicomp\\{m}_{i}_data_{alpha}alpha.csv'
            data = zip(phi_list, mbr)
            save_to_csv(csv_filename, data)
        
        
    # def plotPoints(fname, label, shape,color):
    #     dataset = pd.read_csv(fname)
    #     NH3_list = np.divide(dataset.iloc[:,0],100)
    #     ox_frac_list = np.subtract(1,NH3_list)
    #     O2_list = np.multiply(ox_frac_list, 0.21)
    #     phi_list = np.divide(np.divide(NH3_list,O2_list),np.divide(4,3))
    #     ax.plot(phi_list,dataset.iloc[:,1],shape,linestyle='none',color=color,label=label)
        
    # path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"
    # # # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\50torr.csv','50 torr','o','k')
    # # # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\100torr.csv','100 torr','^','k')
    # # # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\250torr.csv','250 torr','v','k')
    # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv','Ronney','o','k')
    # # # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\1500torr.csv','1500 torr','D','k')

    # # dataset=pd.read_csv(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv')
    # # ax.plot(dataset.iloc[:,0],dataset.iloc[:,1]*100,marker='o',markersize=7,linewidth=3,fillstyle='none',linestyle='none',color='k',label='Ronney')

    # ax.legend(fontsize=15, frameon=False)#, loc='upper right')  
    # ax.set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]', fontsize=18)
    # ax.set_xlabel(r'Equivalence Ratio', fontsize=18)
    # ax.tick_params(axis='both', direction="in", labelsize=15)
    # ax.tick_params(axis='both', which='minor', direction="in")

    # name = f'ronney_flamespeed_{alpha}alpha'
    # if save_plots == True:
    #     plt.savefig(name+'.pdf', dpi=1000, bbox_inches='tight')
    #     plt.savefig(name+'.png', dpi=1000, bbox_inches='tight')
    # # plt.show()     