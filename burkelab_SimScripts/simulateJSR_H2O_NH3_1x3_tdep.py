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
import pandas as pd
import numpy as np
import time
# import cantera as ct
import os.path
from os import path
import matplotlib.pyplot as plt
plt.rcParams.update(plt.rcParamsDefault)
import matplotlib as mpl
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--figwidth', type=float, help="figwidth = ")
parser.add_argument('--figheight', type=float, help="figheight = ")
parser.add_argument('--fsz', type=float, help="mpl.rcParams['font.size'] = ", default=8)
parser.add_argument('--fszxtick', type=float, help="mpl.rcParams['xtick.labelsize'] = ", default=7)
parser.add_argument('--fszytick', type=float, help="mpl.rcParams['ytick.labelsize'] = ", default=7)
parser.add_argument('--fszaxlab', type=float, help="mpl.rcParams['axes.labelsize'] = ", default=8)
parser.add_argument('--lw', type=float, help="lw = ", default=0.7)
parser.add_argument('--mw', type=float, help="mw = ", default=0.5)
parser.add_argument('--msz', type=float, help="msz = ", default=2.5)
parser.add_argument('--lgdw', type=float, help="lgdw = ", default=0.6)
parser.add_argument('--lgdfsz', type=float, help="lgdw = ", default=5)
parser.add_argument('--gridsz', type=int, help="gridsz = ", default=10)
parser.add_argument('--dpi', type=int, help="dpi = ", default=1000)

args = parser.parse_args()
lw=args.lw
mw=args.mw
msz=args.msz
dpi=args.dpi
lgdw=args.lgdw
lgdfsz=args.lgdfsz
gridsz=args.gridsz

mpl.rc('font',family='Times New Roman')
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.size'] = args.fsz
mpl.rcParams['xtick.labelsize'] = args.fszxtick
mpl.rcParams['ytick.labelsize'] = args.fszytick
from matplotlib.legend_handler import HandlerTuple
plt.rcParams['axes.labelsize'] = args.fszaxlab
mpl.rcParams['xtick.major.width'] = 0.5  # Width of major ticks on x-axis
mpl.rcParams['ytick.major.width'] = 0.5  # Width of major ticks on y-axis
mpl.rcParams['xtick.minor.width'] = 0.5  # Width of minor ticks on x-axis
mpl.rcParams['ytick.minor.width'] = 0.5  # Width of minor ticks on y-axis
mpl.rcParams['xtick.major.size'] = 2.5  # Length of major ticks on x-axis
mpl.rcParams['ytick.major.size'] = 2.5  # Length of major ticks on y-axis
mpl.rcParams['xtick.minor.size'] = 1.5  # Length of minor ticks on x-axis
mpl.rcParams['ytick.minor.size'] = 1.5  # Length of minor ticks on y-axis

save_plots = True
f, ax = plt.subplots(3, 7, figsize=(args.figwidth, args.figheight)) 
import matplotlib.ticker as ticker
plt.subplots_adjust(wspace=0.26)
plt.subplots_adjust(hspace=0.15)


ax[1,0].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.4e}"))
ax[1,1].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.4e}"))
ax[1,2].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.4e}"))



# lw=0.7
# mw=0.5
# msz=2.5
# dpi=3000
# lgdw=0.6
# lgdfsz=5
# gridsz=50
# lw=0.7
# mw=0.5
# msz=3.5
# dpi=1000
# lgdw=0.6
# lgdfsz=7
# gridsz=50

name = 'JSR_tdep'
models = {    
        #   'Alzueta':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism.yaml",  
        #   r"$\epsilon_{0,NH_3}(300K)$":"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_epsNH3_T=300K.yaml",  
        #   r"$\epsilon_{0,NH_3}(2000K)$":"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_epsNH3_T=2000K.yaml",            
          'Ar':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR.yaml",
          r'H$_2$O':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allH2O.yaml",
          'LMR-R':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR.yaml", 
          }

# colors = ["xkcd:grey", "xkcd:teal", "orange", 'r', 'b', 'xkcd:purple']
colors = ['r', 'b', 'xkcd:purple']
# lines =['-','-','-','-','-']
# lines =['dotted','dashed','solid']
T_list = np.linspace(800,1050,gridsz)
P = 1.2
tau = 0.5
diluent = 0.94
H2Opercent = 0.20

##############################################################################################################################
# reactorTemperature_list = [800,825,850,875,900,1000]  # Kelvin
reactorTemperature_list = [825,850,875,900,950,1000,1050]  # Kelvin
reactorPressure = P*ct.one_atm  # in atm. This equals 1.06 bars
residenceTime = tau  # s
reactorVolume = 0.000113 #30.5*(1e-2)**3  # m3
heatTransferCoefficient = 79.5 # W/m2/K
reactorRadius = np.cbrt(reactorVolume*3/4/np.pi) # [m3]
reactorSurfaceArea = 4*np.pi*np.square(reactorRadius) # m3
pressureValveCoefficient = 0.01
maxPressureRiseAllowed = 0.01
maxSimulationTime = 500000  # seconds
tempDependence = []

##############################################################################################################################
ax[2, 0].plot(4e4, 2.2, '.', color='white',markersize=0.1,label="20% H$_2$O")  # dummy handle to provide label to lgd column
for k,m in enumerate(models):
    for i, reactorTemperature in enumerate(reactorTemperature_list):
        H2O = diluent * H2Opercent
        Ar = diluent * (1-H2Opercent)
        reactants = {'H2': 0.03, 'O2': 0.03, 'AR': Ar, 'H2O':H2O}    
        
        concentrations = reactants
        gas = ct.Solution(list(models.values())[k])
        gas.TPX = reactorTemperature, reactorPressure, concentrations 
        fuelAirMixtureTank = ct.Reservoir(gas)
        exhaust = ct.Reservoir(gas)
        env = ct.Reservoir(gas)
        stirredReactor = ct.IdealGasReactor(gas, energy='on', volume=reactorVolume)
        massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                    downstream=stirredReactor,
                                                    mdot=stirredReactor.mass/residenceTime)
        pressureRegulator = ct.Valve(upstream=stirredReactor,
                                    downstream=exhaust,
                                    K=pressureValveCoefficient)
        w2 = ct.Wall(stirredReactor, env, A=reactorSurfaceArea, U=heatTransferCoefficient)
        reactorNetwork = ct.ReactorNet([stirredReactor])
        columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
        
        columnNames = ['pressure'] + columnNames
        timeHistory = pd.DataFrame(columns=columnNames)
        # print(columnNames)
        tic = time.time()
        # reactorNetwork.rtol = 1.0e-6
        # reactorNetwork.atol = 1.0e-15
        t = 0
        counter = 1
        while t < maxSimulationTime:
            t = reactorNetwork.step()
            if(counter%10 == 0):
            # if(counter%1 == 0):
                state = np.hstack([stirredReactor.thermo.P, stirredReactor.mass, 
                            stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
                timeHistory.loc[t] = state
            counter += 1
        toc = time.time()

        ax[0,i].plot(timeHistory.index*1e3, np.subtract(timeHistory['temperature'],np.ones(len(timeHistory['temperature']))*reactorTemperature), color=colors[k], linestyle="solid", linewidth=lw, label=m, zorder=k)   
        ax[1,i].plot(timeHistory.index*1e3, timeHistory['O2']*100, color=colors[k], linestyle="solid", linewidth=lw, label=m, zorder=k)   
        ax[2,i].plot(timeHistory.index*1e3, timeHistory['H2']*100, color=colors[k], linestyle="solid", linewidth=lw, label=m, zorder=k) 
        # ax[0,i].set_xscale('log') 
        # ax[0,i].set_yscale('log') 
        # ax[1,i].set_xscale('log') 
        # ax[1,i].set_yscale('linear') 
        # ax[2,i].set_xscale('log') 
        # ax[2,i].set_yscale('linear') 

        

        # ax[1,i].set_title(f"T = {reactorTemperature} K",fontsize=8)
        # ax[0,i].set_xlim([1e-1,5e4])
        # ax[1,i].set_xlim([1e-1,5e4])
        # ax[2,i].set_xlim([1e-1,5e4])
        # ax[0,i].set_ylim([1e-7,500])
        # ax[1,i].set_ylim([0,3.2])
        # ax[2,i].set_ylim([0,3.2])
        
        # ax[0,i].tick_params(axis='both',direction='in')
        # ax[1,i].tick_params(axis='both',direction='in')
        # ax[2,i].tick_params(axis='both',direction='in')
        # ax[0,i].set_ylabel(r'$\Delta$ T [K]') 
        # ax[1,i].set_ylabel('X$_{O_2}$ [%]')
        # ax[2,i].set_ylabel('X$_{H_2}$ [%]')




colors = ['r', 'b', 'xkcd:purple']
P = 1.2
tau = 0.5
diluent = 0.94
NH3percent = 0.1

##############################################################################################################################
# reactorTemperature_list = [800,825,850,875,900,1000]  # Kelvin
reactorTemperature_list = [825,850,875,900,950,1000,1050]  # Kelvin
reactorPressure = P*ct.one_atm  # in atm. This equals 1.06 bars
residenceTime = tau  # s
reactorVolume = 0.000113 #30.5*(1e-2)**3  # m3
reactorRadius = np.cbrt(reactorVolume*3/4*np.pi) # m3
reactorSurfaceArea = 4*np.pi*np.square(reactorRadius) # m3
pressureValveCoefficient = 2e-5
maxPressureRiseAllowed = 0.01
maxSimulationTime = 500000  # seconds
heatTransferCoefficient = 7949.6
heatTransferCoefficient = 7.9496*2.2
tempDependence = []

##############################################################################################################################
ax[2, 0].plot(4e4, 2.2, '.', color='white',markersize=0.1,label="10% NH$_3$")  # dummy handle to provide label to lgd column
for k,m in enumerate(models):
    for i, reactorTemperature in enumerate(reactorTemperature_list):
        NH3 = diluent * NH3percent
        Ar = diluent * (1-NH3percent)
        reactants = {'H2': 0.03, 'O2': 0.03, 'AR': Ar, 'NH3':NH3}    
        
        concentrations = reactants
        gas = ct.Solution(list(models.values())[k])
        gas.TPX = reactorTemperature, reactorPressure, concentrations 
        fuelAirMixtureTank = ct.Reservoir(gas)
        exhaust = ct.Reservoir(gas)
        env = ct.Reservoir(gas)
        stirredReactor = ct.IdealGasReactor(gas, energy='on', volume=reactorVolume)
        massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                    downstream=stirredReactor,
                                                    mdot=stirredReactor.mass/residenceTime)
        pressureRegulator = ct.Valve(upstream=stirredReactor,
                                    downstream=exhaust,
                                    K=pressureValveCoefficient)
        w2 = ct.Wall(stirredReactor, env, A=reactorSurfaceArea, U=heatTransferCoefficient)
        reactorNetwork = ct.ReactorNet([stirredReactor])
        columnNames = [stirredReactor.component_name(item) for item in range(stirredReactor.n_vars)]
        
        columnNames = ['pressure'] + columnNames
        timeHistory = pd.DataFrame(columns=columnNames)
        # print(columnNames)
        tic = time.time()
        # reactorNetwork.rtol = 1.0e-6
        # reactorNetwork.atol = 1.0e-15
        t = 0
        counter = 1
        while t < maxSimulationTime:
            t = reactorNetwork.step()
            if(counter%10 == 0):
            # if(counter%1 == 0):
                state = np.hstack([stirredReactor.thermo.P, stirredReactor.mass, 
                            stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
                timeHistory.loc[t] = state
            counter += 1
        toc = time.time()

        ax[0,i].plot(timeHistory.index*1e3, np.subtract(timeHistory['temperature'],np.ones(len(timeHistory['temperature']))*reactorTemperature), color=colors[k], linestyle="dashed", linewidth=lw, label=m, zorder=k)   
        ax[1,i].plot(timeHistory.index*1e3, timeHistory['O2']*100, color=colors[k], linestyle="dashed", linewidth=lw, label=m, zorder=k)   
        ax[2,i].plot(timeHistory.index*1e3, timeHistory['H2']*100, color=colors[k], linestyle="dashed", linewidth=lw, label=m, zorder=k) 
        ax[0,i].set_xscale('log') 
        ax[0,i].set_yscale('log') 
        ax[1,i].set_xscale('log') 
        ax[1,i].set_yscale('linear') 
        ax[2,i].set_xscale('log') 
        ax[2,i].set_yscale('linear')

        ax[0,i].set_title(f"T = {reactorTemperature} K",fontsize=8)
        ax[2,i].set_xlabel('Time [ms]')
        # ax[0,i].set_xlim([1e-1,5e4])
        # ax[1,i].set_xlim([1e-1,5e4])
        # ax[2,i].set_xlim([1e-1,5e4])
        ax[0,i].set_xlim([1e-1,5e8])
        ax[1,i].set_xlim([1e-1,5e8])
        ax[2,i].set_xlim([1e-1,5e8])
        ax[0,i].set_ylim([1e-7,500])
        ax[1,i].set_ylim([0,3.2])
        ax[2,i].set_ylim([0,3.2])

        ax[0,i].xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=5))
        ax[1,i].xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=5))
        ax[2,i].xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=5))
        
        ax[0,i].tick_params(axis='both',direction='in')
        ax[1,i].tick_params(axis='both',direction='in')
        ax[2,i].tick_params(axis='both',direction='in')


ax[0,0].set_ylabel(r'$\Delta$ T [K]') 
ax[1,0].set_ylabel('X$_{O_2}$ [%]')
ax[2,0].set_ylabel('X$_{H_2}$ [%]')

# ax[1,5].set_xlabel('Time [ms]')
# legend = ax[1,0].legend(fontsize=lgdfsz,frameon=False,loc='center', handlelength=lgdw,ncol=2,columnspacing=0.5)
legend = ax[2,0].legend(fontsize=lgdfsz,frameon=False,loc='center', handlelength=lgdw,ncol=1,columnspacing=0.5)

for text in legend.get_texts():
    if text.get_text() == "10% NH$_3$" or text.get_text() == "20% H$_2$O":
        text.set_fontsize(6)  # Set a larger font size
        text.set_fontweight('bold')  # Make the font bold
# plt.text(0.5, -0.05, 'Time [ms]', ha='center', va='center',fontsize=args.fszaxlab)
if save_plots == True:
    plt.savefig('burkelab_SimScripts/figures/'+name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig('burkelab_SimScripts/figures/'+name+'.png', dpi=1000, bbox_inches='tight')
# plt.show()     