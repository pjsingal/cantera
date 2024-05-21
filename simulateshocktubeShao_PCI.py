from __future__ import division
from __future__ import print_function
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


import matplotlib as mpl
import sys, os
mpl.rc('font',family='Times New Roman')
mpl.rcParams['mathtext.fontset'] = 'stix'
from matplotlib.legend_handler import HandlerTuple
save_plots = True
fig, ax = plt.subplots(1, 1, figsize=(3.8, 2))
name = 'ShockTubeSpeciesProfile_H2O' #os.path.splitext(os.path.basename(__file__))[0]

refSpecies='H2O'
X_H2O2 = 1163e-6
X_H2O = 1330e-6
X_O2 = 665e-6
X_CO2= 0.2*(1-X_H2O2-X_H2O-X_O2)
X_Ar = 1-X_CO2
def plotXvsTime(fname,pltlabel,pltcolour,lstyle='solid'):
    # gas = ct.Solution('test/data/Burke_H2_ArBath.yaml')
    gas = ct.Solution(fname)
    gas.TPX = 1196, 2.127*101325, {'H2O2':X_H2O2, 'H2O':X_H2O, 'O2':X_O2, 'CO2':X_CO2, 'AR':X_Ar}
    r = ct.Reactor(contents=gas,energy="on")
    reactorNetwork = ct.ReactorNet([r]) # this will be the only reactor in the network
    timeHistory = ct.SolutionArray(gas, extra=['t'])
    estIgnitDelay = 1
    t = 0
    counter = 1
    while t < estIgnitDelay:
        t = reactorNetwork.step()
        if counter % 10 == 0:
            timeHistory.append(r.thermo.state, t=t)
        counter += 1
    tConv = 1e6 #time conversion factor (1e6 converts to microseconds)
    timeShift=0 # [seconds]
    shiftedTime = tConv*(timeHistory.t - timeShift)
    moleFrac = timeHistory(refSpecies).X 
    ax.plot(shiftedTime, moleFrac*100, color=pltcolour,label=pltlabel,linestyle=lstyle,linewidth=0.7)

path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"
# path=os.getcwd()
def plotPoints(filename,mkr='none',mkrw='none',mkrsz='none',line='none',fill='none',colour='k',subplot='off',pltLabel="_hidden"): 
    dataset = pd.read_csv(filename)
    ax.plot(dataset.iloc[:,0],dataset.iloc[:,1]*100,mkr,linewidth=0.7,fillstyle=fill,linestyle=line,color=colour,label=pltLabel,markersize=mkrsz,markeredgewidth=mkrw)

# plotXvsTime("test/data/alzuetamechanism.yaml","Alzueta","xkcd:grey")
plotXvsTime("test/data/alzuetamechanism_LMRR_allAR_PCI.yaml","Ar","r")
plotXvsTime("test/data/alzuetamechanism_LMRR_allH2O_PCI.yaml",r'$\rm H_2O$',"b")
plotXvsTime("test/data/alzuetamechanism_LMRR_PCI.yaml","LMR-R","xkcd:purple")
# plotPoints(path+'\\7 SP H2O X vs t (Shock Tube) (Shao)\\expData.csv',pltLabel='Shao et al.',line=':',colour='k')
plotPoints(path+'\\7 SP H2O X vs t (Shock Tube) (Shao)\\expData.csv',mkr='o',mkrsz=3.5,pltLabel='Shao et al.',mkrw=0.5)
# plotPoints(path+'\\7 SP H2O X vs t (Shock Tube) (Shao)\\troe_k0co2.csv',pltLabel='Troe et al.',line='solid',colour='g')
    
ax.legend(fontsize=7, frameon=False, loc='lower right')  
ax.set_ylabel(r'$\rm H_2O$ mole fraction [%]', fontsize=10) #CHECK UNITS OF Y-AXIS
ax.set_xlabel(r'Time [$\mathdefault{\mu s}$]', fontsize=10)
ax.tick_params(axis='both', direction="in", labelsize=7)
ax.set_xlim([0,300])
ax.set_ylim([0.001*100,0.003*100])

import matplotlib.ticker as ticker
ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}"))

if save_plots == True:
    plt.savefig(name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig(name+'.png', dpi=1000, bbox_inches='tight')
plt.show()     