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
parser.add_argument('--date', type=str, help="sim date = ")

args = parser.parse_args()
lw=args.lw
mw=args.mw
msz=args.msz
dpi=args.dpi
lgdw=args.lgdw
lgdfsz=args.lgdfsz
gridsz=args.gridsz
from matplotlib.legend_handler import HandlerTuple
mpl.rc('font',family='Times New Roman')
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.size'] = args.fsz
mpl.rcParams['xtick.labelsize'] = args.fszxtick
mpl.rcParams['ytick.labelsize'] = args.fszytick
plt.rcParams['axes.labelsize'] = args.fszaxlab
mpl.rcParams['xtick.major.width'] = 0.5  # Width of major ticks on x-axis
mpl.rcParams['ytick.major.width'] = 0.5  # Width of major ticks on y-axis
mpl.rcParams['xtick.minor.width'] = 0.5  # Width of minor ticks on x-axis
mpl.rcParams['ytick.minor.width'] = 0.5  # Width of minor ticks on y-axis
mpl.rcParams['xtick.major.size'] = 2.5  # Length of major ticks on x-axis
mpl.rcParams['ytick.major.size'] = 2.5  # Length of major ticks on y-axis
mpl.rcParams['xtick.minor.size'] = 1.5  # Length of minor ticks on x-axis
mpl.rcParams['ytick.minor.size'] = 1.5  # Length of minor ticks on y-axis

plt.subplots_adjust(wspace=0.18)
# colors = ["xkcd:grey"]*3 + ["xkcd:teal"]*3 + ["orange"]*3 + ['r']*3 + ['b']*3 + ['xkcd:purple']*3
lstyles = ["solid","dashed","dotted"]*6
colors = ["xkcd:purple","xkcd:teal","k"]*3
# lstyles = ["solid"]*3 + ["dashed"]*3 + ["dotted"]*3

models = {
    'Alzueta-2023': {
        'a priori': r'test\\data\\alzuetamechanism.yaml',
        'LMRR': r'C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\alzuetamechanism_LMRR.yaml',
        'LMRR-allP': r'C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\alzuetamechanism_LMRR_allP.yaml',
                },
    'Mei-2019': {
        'a priori': r'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Mei-2019\\mei-2019.yaml',
        'LMRR': r'C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\mei-2019_LMRR.yaml',
        'LMRR-allP': r'C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\mei-2019_LMRR_allP.yaml',
                },
    'Zhang-2017': {
        'a priori': r"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Zhang-2017\\zhang-2017.yaml",
        'LMRR': r"C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\zhang-2017_LMRR.yaml",
        'LMRR-allP': r"C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\zhang-2017_LMRR_allP.yaml",
                },
    'Otomo-2018': {
        'a priori': r"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Otomo-2018\\otomo-2018.yaml",
        'LMRR': r"C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\otomo-2018_LMRR.yaml",
        'LMRR-allP': r"C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\otomo-2018_LMRR_allP.yaml",
                },
    'Stagni-2020': {
        'a priori': r"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Stagni-2020\\stagni-2020.yaml",
        'LMRR': r"C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\stagni-2020_LMRR.yaml",
        'LMRR-allP': r"C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\stagni-2020_LMRR_allP.yaml",
                },
}

name = 'JSR_NH3_multimech'
save_plots = True
# figsize=(3.5,7)
f, ax = plt.subplots(3, len(models.keys()), figsize=(args.figwidth,args.figheight)) 
plt.suptitle('JSR_NH3', fontsize=12)
for z, n in enumerate(models):
    mech = n

    import matplotlib.ticker as ticker
    plt.subplots_adjust(hspace=0.3)
    # plt.subplots_adjust(wspace=0.3)
    ax[0,z].yaxis.set_major_locator(ticker.MultipleLocator(5))
    ax[1,z].yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    ax[2,z].yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax[0,z].xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax[1,z].xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax[2,z].xaxis.set_major_locator(ticker.MultipleLocator(100))
    ax[0,z].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax[0,z].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax[1,z].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax[1,z].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
    ax[2,z].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax[2,z].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))

    T_list = np.concatenate((np.linspace(800,829.9,10),np.linspace(830,844,100),np.linspace(848.8,870.9,100),np.linspace(871,1050,50)))
    # T_list = np.linspace(800,1050,50)
    P = 1.2
    tau = 0.5

    # NH3_con_list = [0.02, 0.05, 0.075, 0.10]
    # H2O_con_list = [0, 0.05, 0.10, 0.15, 0.20]
    diluent = 0.94
    # NH3percent_list = [0, 0.05, 0.10, 0.15, 0.20]
    NH3percent_list = [0.10]
    lines =['-','--','-','-','-']


    ##############################################################################################################################

    reactorTemperature = 1000  # Kelvin
    reactorPressure = P*ct.one_atm  # in atm. This equals 1.06 bars
    residenceTime = tau  # s
    reactorVolume = 0.000113 #30.5*(1e-2)**3  # m3
    heatTransferCoefficient = 79.5 # W/m2/K
    reactorRadius = np.cbrt(reactorVolume*3/4/np.pi) # [m3]
    reactorSurfaceArea = 4*np.pi*np.square(reactorRadius) # m3
    pressureValveCoefficient = 2e-5
    maxPressureRiseAllowed = 0.01
    maxSimulationTime = 50  # seconds
    tempDependence = []

    ##############################################################################################################################


    path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\2 JSR NH3\\"
    
    T_10_data = pd.read_csv(path+'JSR_T_NH3_10_data.csv') 
    O2_10_data = pd.read_csv(path+'JSR_O2_NH3_10_data.csv') 
    H2_10_data = pd.read_csv(path+'JSR_H2_NH3_10_data.csv') 
    ax[0,z].plot(T_10_data.iloc[:, 0],T_10_data.iloc[:, 1],marker='o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw, label="Sabia et al.")
    ax[1,z].plot(O2_10_data.iloc[:, 0],O2_10_data.iloc[:, 1],marker='o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw,label="Sabia et al.")
    ax[2,z].plot(H2_10_data.iloc[:, 0],H2_10_data.iloc[:, 1],marker='o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw, label="Sabia et al.")


    for k,m in enumerate(models[n]):
        for i, NH3percent in enumerate(NH3percent_list):
            NH3 = diluent * NH3percent
            Ar = diluent * (1-NH3percent)
            # reactants = {'H2': 0.03, "H":1e-6, 'O2': 0.03, 'AR': Ar, 'NH3':NH3}    
            reactants = {'H2': 0.03, 'O2': 0.03, 'AR': Ar, 'NH3':NH3}    
            
            # for j, T in enumerate(T_list):
            
            concentrations = reactants
            gas = ct.Solution(list(models[n].values())[k])
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
            tic = time.time()
            # reactorNetwork.rtol = 1.0e-6
            # reactorNetwork.atol = 1.0e-15
            t = 0
            counter = 1
            while t < maxSimulationTime:
                t = reactorNetwork.step()
                if(counter%10 == 0):
                    state = np.hstack([stirredReactor.thermo.P, stirredReactor.mass, 
                                stirredReactor.volume, stirredReactor.T, stirredReactor.thermo.X])
                    timeHistory.loc[t] = state
                counter += 1
            toc = time.time()
            
            pressureDifferential = timeHistory['pressure'].max()-timeHistory['pressure'].min()
            if(abs(pressureDifferential/reactorPressure) > maxPressureRiseAllowed):
                print("WARNING: Non-trivial pressure rise in the reactor. Adjust K value in valve")
                
            tempDependence.append(pd.DataFrame(columns=timeHistory.columns))
            tempDependence[i].index.name = 'Temperature'
            
            inletConcentrations = concentrations
            
            for j,T in enumerate(T_list): #temperature in T:
                reactorTemperature = T #temperature  # Kelvin
                gas.TPX = reactorTemperature, reactorPressure, inletConcentrations
                timeHistory = pd.DataFrame(columns=columnNames)
                fuelAirMixtureTank = ct.Reservoir(gas)
                exhaust = ct.Reservoir(gas)
                env = ct.Reservoir(gas)
                # gas.TPX = reactorTemperature, reactorPressure, concentrations
                stirredReactor = ct.IdealGasReactor(gas, energy='on', volume=reactorVolume)
                
                massFlowController = ct.MassFlowController(upstream=fuelAirMixtureTank,
                                                        downstream=stirredReactor,
                                                        mdot=stirredReactor.mass/residenceTime)
                pressureRegulator = ct.Valve(upstream=stirredReactor, 
                                            downstream=exhaust, 
                                            K=pressureValveCoefficient)
                w2 = ct.Wall(stirredReactor, env, A=reactorSurfaceArea, U=heatTransferCoefficient)
                reactorNetwork = ct.ReactorNet([stirredReactor])
                tic = time.time()
                t = 0
                while t < maxSimulationTime:
                    t = reactorNetwork.step()
                state = np.hstack([stirredReactor.thermo.P, 
                                stirredReactor.mass, 
                                stirredReactor.volume, 
                                stirredReactor.T, 
                                stirredReactor.thermo.X])
                toc = time.time()
                concentrations = stirredReactor.thermo.X
                tempDependence[i].loc[T] = state
            ax[0,z].plot(tempDependence[i].index, np.subtract(tempDependence[i]['temperature'],tempDependence[i].index), color=colors[k], linestyle=lstyles[k], linewidth=lw, label=m)   
            ax[1,z].plot(tempDependence[i].index, tempDependence[i]['O2']*100, color=colors[k], linestyle=lstyles[k], linewidth=lw, label=m)   
            ax[2,z].plot(tempDependence[i].index, tempDependence[i]['H2']*100, color=colors[k], linestyle=lstyles[k], linewidth=lw, label=m) 
    ax[0,z].set_title(f"{mech}")
    # ax[0,z].set_xlabel('Temperature [K]')
    ax[0,z].tick_params(axis='both',direction='in')
    # ax[0,z].legend(frameon=False)#,loc='lower right')

    # ax[1,z].set_xlabel('Temperature [K]')

    ax[1,z].tick_params(axis='both',direction='in')
    # ax[1,z].legend(frameon=False)#,loc='upper right')

    ax[2,z].tick_params(axis='both',direction='in')

    ax[2,z].legend(fontsize=lgdfsz,frameon=False,loc='lower left', handlelength=lgdw)

    ax[0,z].set_xlim([780,1070])
    ax[0,z].set_ylim([-1,24.99])
    ax[1,z].set_xlim([780,1070])
    ax[1,z].set_ylim([1.45,3.2])
    ax[2,z].set_xlim([780,1070])
    ax[2,z].set_ylim([0.0001,3.3])

ax[0,0].set_ylabel(r'$\Delta$ T [K]')
ax[1,0].set_ylabel('O$_2$ mole fraction [%]')
ax[2,0].set_ylabel('H$_2$ mole fraction [%]')
ax[2,2].set_xlabel('Temperature [K]')

path=f'burkelab_SimScripts/USSCI_simulations/figures/'+args.date
os.makedirs(path,exist_ok=True)

if save_plots == True:
    plt.savefig(path+f'/{name}.png', dpi=500, bbox_inches='tight')