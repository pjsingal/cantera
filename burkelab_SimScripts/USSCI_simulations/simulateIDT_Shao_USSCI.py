#%%
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import matplotlib.pyplot as plt
import pandas as pd 
import time
import numpy as np
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
parser.add_argument('--dpi', type=int, help="dpi = ", default=500)
parser.add_argument('--LMRtest', type=int, help="LMRtest = ", default=0)

args = parser.parse_args()
lw=args.lw
mw=args.mw
msz=args.msz
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
        'base': r'test\\data\\alzuetamechanism.yaml',
        'LMRR': r'C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\alzuetamechanism_LMRR.yaml',
        'LMRR-allP': r'C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\alzuetamechanism_LMRR_allP.yaml',
                },
    'Mei-2019': {
        'base': r'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Mei-2019\\mei-2019.yaml',
        'LMRR': r'C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\mei-2019_LMRR.yaml',
        'LMRR-allP': r'C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\mei-2019_LMRR_allP.yaml',
                },
    'Zhang-2017': {
        'base': r"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Zhang-2017\\zhang-2017.yaml",
        'LMRR': r"C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\zhang-2017_LMRR.yaml",
        'LMRR-allP': r"C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\zhang-2017_LMRR_allP.yaml",
                },
    'Otomo-2018': {
        'base': r"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Otomo-2018\\otomo-2018.yaml",
        'LMRR': r"C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\otomo-2018_LMRR.yaml",
        'LMRR-allP': r"C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\otomo-2018_LMRR_allP.yaml",
                },
    'Stagni-2020': {
        'base': r"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Stagni-2020\\stagni-2020.yaml",
        'LMRR': r"C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\stagni-2020_LMRR.yaml",
        'LMRR-allP': r"C:\\Users\\pjsin\\Documents\\LMRRfactory\\test\outputs\\Oct22\\stagni-2020_LMRR_allP.yaml",
                },
}


name = 'IDT_shao_multimech'
save_plots = True
f, ax = plt.subplots(4, len(models.keys()), figsize=(args.figwidth, args.figheight)) 
plt.subplots_adjust(wspace=0.18)
plt.suptitle('IDT_shao', fontsize=12)

for z, n in enumerate(models):
    mech = n

    import matplotlib.ticker as ticker
    plt.subplots_adjust(wspace=0.18)
    ax[0,z].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
    ax[1,z].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
    ax[2,z].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
    ax[3,z].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))

    ax[0,z].xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax[1,z].xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax[2,z].xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax[3,z].xaxis.set_major_locator(ticker.MultipleLocator(50))

    ax[0,z].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax[0,z].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0e}"))
    ax[1,z].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax[1,z].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1e}"))
    ax[2,z].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax[2,z].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1e}"))
    ax[3,z].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax[3,z].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1e}"))


    # f.text(0.5, -0.1, r'Temperature [K]', ha='center', va='center')

    def ignitionDelay(states, species):
        # i_ign = states(species).Y.argmax()
        i_ign = np.gradient(states(species).Y.T[0]).argmax()
        return states.t[i_ign]

    ################################################################################################

    path='G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\'
    df = pd.read_csv(path+'Shao_IDT\\1.csv')
    p_df = df['P']
    T_df = df['T']
    IDT_df = df['IDT']

    T_list = np.linspace(1100,1300,gridsz)#[::-1]
    for k, m in enumerate(models[n]):
        estimatedIgnitionDelayTimes = np.ones(len(T_list))
        estimatedIgnitionDelayTimes[:] = 0.05
        ignitionDelays_RG = np.zeros(len(T_list))
        for j, T in enumerate(T_list):
            gas = ct.Solution(list(models[n].values())[k])
            gas.TPX = T, 17*ct.one_atm, {'H2':0.03, 'O2':0.015, 'Ar':1-0.03-0.015}
            r = ct.Reactor(contents=gas)
            reactorNetwork = ct.ReactorNet([r])
            timeHistory = ct.SolutionArray(gas, extra=['t'])
            t0 = time.time()
            t = 0
            counter = 1
            while t < estimatedIgnitionDelayTimes[j]:
                t = reactorNetwork.step()
                if counter % 1 == 0:
                    timeHistory.append(r.thermo.state, t=t)
                counter += 1
            tau = ignitionDelay(timeHistory, 'oh')
            t1 = time.time()
            ignitionDelays_RG[j] = tau
        if colors[k] == 'xkcd:purple':
            zorder_value = 10  # Higher z-order for purple line
        else:
            zorder_value = k  # Default z-order for other lines
        ax[0,z].semilogy(T_list, 1e6*ignitionDelays_RG, '-', linestyle=lstyles[k], linewidth=lw, color=colors[k], label=m, zorder=zorder_value)
        
    ax[0,z].semilogy(T_df,IDT_df,'o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw,label='Shao et al.', zorder=12)
    if z==0:
        ax[0,z].set_ylabel(r'Ignition delay [$\mathdefault{\mu s}$]')
    ax[0,z].tick_params(axis='both', direction="in")
    ax[0,z].tick_params(axis='both', which='minor', direction="in")#, bottom=False)

    ################################################################################################

    df = pd.read_csv(path+'\\Shao_IDT\\2.csv')
    p_df = df['P']
    T_df = df['T']
    IDT_df = df['IDT']
    T_list = np.linspace(1100,1300,gridsz)#[::-1]
    for k, m in enumerate(models[n]):
        estimatedIgnitionDelayTimes = np.ones(len(T_list))
        estimatedIgnitionDelayTimes[:] = 0.05
        ignitionDelays_RG = np.zeros(len(T_list))
        for j, T in enumerate(T_list):
            gas = ct.Solution(list(models[n].values())[k])
            gas.TPX = T, 13*ct.one_atm, {'H2':0.03, 'O2':0.015, 'N2':1-0.03-0.015}
            r = ct.Reactor(contents=gas)
            reactorNetwork = ct.ReactorNet([r])
            timeHistory = ct.SolutionArray(gas, extra=['t'])
            t0 = time.time()
            t = 0
            counter = 1
            while t < estimatedIgnitionDelayTimes[j]:
                t = reactorNetwork.step()
                if counter % 1 == 0:
                    timeHistory.append(r.thermo.state, t=t)
                counter += 1
            tau = ignitionDelay(timeHistory, 'oh')
            t1 = time.time()
            ignitionDelays_RG[j] = tau
        if colors[k] == 'xkcd:purple':
            zorder_value = 10  # Higher z-order for purple line
        else:
            zorder_value = k  # Default z-order for other lines
        ax[1,z].semilogy(T_list, 1e6*ignitionDelays_RG, '-', linestyle=lstyles[k],linewidth=lw, color=colors[k], label=m, zorder=zorder_value)
        
    ax[1,z].semilogy(T_df,IDT_df,'o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw,label='Shao et al.', zorder=12)
    if z==0:
        ax[1,z].set_ylabel(r'Ignition delay [$\mathdefault{\mu s}$]')
    ax[1,z].tick_params(axis='both', direction="in")
    ax[1,z].tick_params(axis='both', which='minor', direction="in")#, bottom=False)

    ################################################################################################

    df = pd.read_csv(path+'\\Shao_IDT\\3.csv')
    p_df = df['P']
    T_df = df['T']
    IDT_df = df['IDT']
    H2O_df = df['H2O']
    T_list = np.linspace(1200,1400,gridsz)#[::-1]
    for k, m in enumerate(models[n]):
        estimatedIgnitionDelayTimes = np.ones(len(T_list))
        estimatedIgnitionDelayTimes[:] = 0.05
        ignitionDelays_RG = np.zeros(len(T_list))
        for j, T in enumerate(T_list):
            gas = ct.Solution(list(models[n].values())[k])
            gas.TPX = T, 15*ct.one_atm, {'H2':0.03, 'O2':0.015, 'H2O':0.09, 'Ar':1-0.03-0.015-0.09}
            r = ct.Reactor(contents=gas)
            reactorNetwork = ct.ReactorNet([r])
            timeHistory = ct.SolutionArray(gas, extra=['t'])
            t0 = time.time()
            t = 0
            counter = 1
            while t < estimatedIgnitionDelayTimes[j]:
                t = reactorNetwork.step()
                if counter % 1 == 0:
                    timeHistory.append(r.thermo.state, t=t)
                counter += 1
            tau = ignitionDelay(timeHistory, 'oh')
            t1 = time.time()
            ignitionDelays_RG[j] = tau
        if colors[k] == 'xkcd:purple':
            zorder_value = 10  # Higher z-order for purple line
        else:
            zorder_value = k  # Default z-order for other lines
        ax[2,z].semilogy(T_list, 1e6*ignitionDelays_RG, '-', linestyle=lstyles[k],linewidth=lw, color=colors[k], label=m, zorder=zorder_value)
    ax[2,z].semilogy(T_df,IDT_df,'o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw,label='Shao et al.', zorder=12)
    if z==0:
        ax[2,z].set_ylabel(r'Ignition delay [$\mathdefault{\mu s}$]')
    ax[2,z].tick_params(axis='both', direction="in")
    ax[2,z].tick_params(axis='both', which='minor', direction="in")#, bottom=False)

    ################################################################################################

    df = pd.read_csv(path+'\\Shao_IDT\\4.csv')
    p_df = df['P']
    T_df = df['T']
    IDT_df = df['IDT']
    T_list = np.linspace(1100,1300,gridsz)#[::-1]
    for k, m in enumerate(models[n]):
        estimatedIgnitionDelayTimes = np.ones(len(T_list))
        estimatedIgnitionDelayTimes[:] = 0.05
        ignitionDelays_RG = np.zeros(len(T_list))
        for j, T in enumerate(T_list):
            gas = ct.Solution(list(models[n].values())[k])
            gas.TPX = T, 12*ct.one_atm, {'H2':0.03, 'O2':0.015, 'CO2':0.20, 'Ar':1-0.2-0.03-0.015}
            r = ct.Reactor(contents=gas)
            reactorNetwork = ct.ReactorNet([r])
            timeHistory = ct.SolutionArray(gas, extra=['t'])
            t0 = time.time()
            t = 0
            counter = 1
            while t < estimatedIgnitionDelayTimes[j]:
                t = reactorNetwork.step()
                if counter % 1 == 0:
                    timeHistory.append(r.thermo.state, t=t)
                counter += 1
            tau = ignitionDelay(timeHistory, 'oh')
            t1 = time.time()
            ignitionDelays_RG[j] = tau
        if colors[k] == 'xkcd:purple':
            zorder_value = 10  # Higher z-order for purple line
        else:
            zorder_value = k  # Default z-order for other lines
        ax[3,z].semilogy(T_list, 1e6*ignitionDelays_RG, '-', linestyle=lstyles[k],linewidth=lw, color=colors[k], label=m, zorder=zorder_value)
    ax[3,z].semilogy(T_df,IDT_df,'o',fillstyle='none',linestyle='none',color='k',markersize=msz,markeredgewidth=mw,label='Shao et al.', zorder=12)
    ax[0,z].legend(fontsize=lgdfsz, frameon=False, loc='lower left',handlelength=lgdw)
    ax[3,z].tick_params(axis='both', direction="in")
    ax[3,z].tick_params(axis='both', which='minor', direction="in")#, bottom=False)
    ax[0,z].set_title(f"{mech}")
    # plt.subplots_adjust(wspace=0.4, hspace=0.4)

    # ax[0,z].set_xlim([1000.1,1499.99])
    # ax[1,z].set_xlim([1000.1,1499.99])
    # ax[2,z].set_xlim([1000.1,1499.99])
    # ax[3,z].set_xlim([1000.1,1499.99])

ax[3,0].set_ylabel(r'Ignition delay [$\mathdefault{\mu s}$]')
ax[3,2].set_xlabel(r'Temperature [K]')

path=f'burkelab_SimScripts/USSCI_simulations/figures/'+args.date
os.makedirs(path,exist_ok=True)

if save_plots == True:
    plt.savefig(path+f'/{name}.png', dpi=500, bbox_inches='tight')