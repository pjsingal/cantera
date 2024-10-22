
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import pandas as pd 
import matplotlib as mpl
import numpy as np
import os

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
parser.add_argument('--title', type=str, help="title = ",default='null')


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



save_plots = True
fig, ax = plt.subplots(1,1,figsize=(args.figwidth, args.figheight))

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

name = 'Burke_flamespeed_multimech'
save_plots = True
fig, ax = plt.subplots(1, len(models.keys()),figsize=(args.figwidth, args.figheight))

for z, n in enumerate(models):
    mech = n

    import matplotlib.ticker as ticker
    ax[z].xaxis.set_major_locator(ticker.MultipleLocator(2))
    ax[z].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax[z].yaxis.set_major_locator(ticker.MultipleLocator(0.03))
    ax[z].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}"))

    # Plot experimental data
    path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"
    dataset = pd.read_csv(path+'\\5 FS H2O (Burke)\\exp_pts.csv',header=None)
    ax[z].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='o',fillstyle='none',markersize=msz,markeredgewidth=mw,linestyle='none',color='k',label='Burke et al.',zorder=100)

    # Plot simulation data
    for k, m in enumerate(models[n]):
        path=f'burkelab_SimScripts\\USSCI_simulations\\data\\Burke\\'+args.date+'\\'
        dataset=pd.read_csv(path+f'{n}_{m}.csv')
        ax[z].plot(dataset.iloc[:,0],dataset.iloc[:,1],color=colors[k],linestyle=lstyles[k],linewidth=lw,label=m)
    ax[z].tick_params(axis='both', direction="in")
    ax[z].tick_params(axis='both', which='minor', direction="in")
    ax[z].set_xlim([0.001, 15.999])
    ax[z].set_ylim([-0.005, 0.1299])

ax[0].legend(fontsize=lgdfsz, frameon=False, loc='right', handlelength=lgdw)
ax[0].set_ylabel(r'Mass burning rate [g $\rm cm^{-2}$ $\rm s^{-1}$]',fontsize=args.fszaxlab)
ax[2].set_xlabel(r'Pressure [atm]',fontsize=args.fszaxlab)

path=f'burkelab_SimScripts\\USSCI_simulations\\figures\\'+args.date
os.makedirs(path,exist_ok=True)
if save_plots == True:
    plt.savefig(path+f'\\{name}.png', dpi=500, bbox_inches='tight')