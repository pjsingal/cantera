
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
parser.add_argument('--date', type=str, help="sim date = ",default='May28')
parser.add_argument('--slopeVal', type=float, help="slope value = ",default=-1)
parser.add_argument('--curveVal', type=float, help="curve value = ",default=-1)

args = parser.parse_args()
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

lw=args.lw
mw=args.mw
msz=args.msz
dpi=args.dpi
lgdw=args.lgdw
lgdfsz=args.lgdfsz
date=args.date
fslope=args.slopeVal
fcurve=args.curveVal

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

name = 'Ronney_flamespeed_multimech'
save_plots = True
fig, ax = plt.subplots(2, len(models.keys()),figsize=(args.figwidth, args.figheight))

for z, n in enumerate(models):
    mech = n

    import matplotlib.ticker as ticker
    ax[0,z].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax[0,z].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
    ax[0,z].yaxis.set_major_locator(ticker.MultipleLocator(2))
    ax[0,z].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax[1,z].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax[1,z].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
    ax[1,z].yaxis.set_major_locator(ticker.MultipleLocator(10))
    ax[1,z].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))

    # Plot experimental data
    path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"
    dataset = pd.read_csv(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv')
    NH3_list = np.divide(dataset.iloc[:,0],100)
    ox_frac_list = np.subtract(1,NH3_list)
    O2_list = np.multiply(ox_frac_list, 0.21)
    phi_list = np.divide(np.divide(NH3_list,O2_list),np.divide(4,3))
    ax[0,z].plot(phi_list,dataset.iloc[:,1],marker='o',fillstyle='none',markersize=msz,markeredgewidth=mw,linestyle='none',color='k',label='Ronney',zorder=100)
    dataset = pd.read_csv(path+f'\\Han\\han_0pt6_NH3.csv')
    ax[1,z].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='s',fillstyle='none',markersize=msz,markeredgewidth=mw,linestyle='none',color='k',label='Han',zorder=100)
    dataset = pd.read_csv(path+f'\\Wang\\wang_0pt6_NH3.csv')
    ax[1,z].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='x',fillstyle='none',markersize=msz,markeredgewidth=mw,linestyle='none',color='k',label='Wang',zorder=99)

    # Plot simulation data
    path=f'burkelab_SimScripts\\USSCI_simulations\\data\\Ronney\\'+args.date+'\\'
    alpha_list = [1.0,0.6]
    p_list=[760]
    for x, alpha in enumerate(alpha_list):
      for k, m in enumerate(models[n]):
          for i, p in enumerate(p_list):
            fname = f'{n}_{m}_{p}torr_{alpha}alpha.csv'
            dataset=pd.read_csv(path+fname)
            ax[x,z].plot(dataset.iloc[:,0],dataset.iloc[:,1],color=colors[k],linestyle=lstyles[k],linewidth=lw,label=m)

    ax[0,z].legend(fontsize=lgdfsz, frameon=False, loc='upper right',handlelength=lgdw)
    ax[1,z].legend(fontsize=lgdfsz, frameon=False, loc='upper right',handlelength=lgdw)

    ax[0,z].set_title(f'{mech}')

    ax[0,z].tick_params(axis='both', direction="in")
    ax[0,z].tick_params(axis='both', which='minor', direction="in")
    ax[1,z].tick_params(axis='both', direction="in")
    ax[1,z].tick_params(axis='both', which='minor', direction="in")

    if z==0:
        ax[0,z].set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]')
        ax[1,z].set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]')
    if z==2:
        ax[1,z].set_xlabel(r'Equivalence Ratio')
    ax[0,z].set_xlim([0.6001, 1.7999])
    ax[0,z].set_ylim([0.001, 11.9999])
    ax[1,z].set_xlim([0.6001, 1.7999])
    ax[1,z].set_ylim([0.001, 43])

path=f'burkelab_SimScripts\\USSCI_simulations\\figures\\'+args.date
os.makedirs(path,exist_ok=True)
if save_plots == True:
    plt.savefig(path+f'\\{name}.png', dpi=500, bbox_inches='tight')