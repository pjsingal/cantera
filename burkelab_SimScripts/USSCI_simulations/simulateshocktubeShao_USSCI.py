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
import argparse

import matplotlib as mpl
import sys, os
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


import matplotlib.ticker as ticker

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

lstyles = ["solid","dashed","dotted"]*6
colors = ["xkcd:purple","xkcd:teal","k"]*3
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

name = 'ShockTube_H2O_multimech'
save_plots = True
fig, ax = plt.subplots(1, len(models.keys()),figsize=(args.figwidth, args.figheight))

for z, n in enumerate(models):
    mech = n

    import matplotlib.ticker as ticker
    ax[z].xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax[z].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax[z].yaxis.set_major_locator(ticker.MultipleLocator(0.03))
    ax[z].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}"))

    path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"
    shao_data = pd.read_csv(path+'\\7 SP H2O X vs t (Shock Tube) (Shao)\\expData.csv')
    ax[z].plot(shao_data.iloc[:,0],shao_data.iloc[:,1]*100,marker='o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw,label='Shao et al.')

    for k,m in enumerate(models[n]):
        X_H2O2 = 1163e-6
        X_H2O = 1330e-6
        X_O2 = 665e-6
        X_CO2= 0.2*(1-X_H2O2-X_H2O-X_O2)
        X_Ar = 1-X_CO2
        gas = ct.Solution(list(models[n].values())[k])
        X = {'H2O2':X_H2O2, 'H2O':X_H2O, 'O2':X_O2, 'CO2':X_CO2, 'AR':X_Ar}
        gas.TPX = 1196, 2.127*ct.one_atm, X
        r = ct.Reactor(contents=gas,energy="on")
        reactorNetwork = ct.ReactorNet([r])
        timeHistory = ct.SolutionArray(gas, extra=['t'])
        estIgnitDelay = 1
        t = 0
        counter = 1
        while t < estIgnitDelay:
            t = reactorNetwork.step()
            if counter % 10 == 0:
                timeHistory.append(r.thermo.state, t=t)
            counter += 1
        ax[z].plot(timeHistory.t*1e6, timeHistory('H2O').X*100, color=colors[k],linestyle=lstyles[k],linewidth=lw,label=m)


    ax[z].legend(fontsize=lgdfsz,handlelength=lgdw, frameon=False, loc='lower right') 
    ax[z].tick_params(axis='both', direction="in")#, labelsize=7)
    ax[z].set_xlim([0.0001,299.999])
    ax[z].set_ylim([0.120001,0.269999])
    ax[z].set_title(f"{mech}")

ax[0].set_ylabel(r'$\rm H_2O$ mole fraction [%]')
ax[2].set_xlabel(r'Time [$\mathdefault{\mu s}$]')

path=f'burkelab_SimScripts/USSCI_simulations/figures/'+args.date
os.makedirs(path,exist_ok=True)
if save_plots == True:
    plt.savefig(path+f'/{name}.png', dpi=500, bbox_inches='tight')