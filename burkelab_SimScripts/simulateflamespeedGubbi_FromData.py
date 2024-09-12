
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

name = 'gubbi_flamespeed'

fig, ax = plt.subplots(1,2,figsize=(args.figwidth, args.figheight))
save_plots = True
lw=args.lw
mw=args.mw
msz=args.msz
dpi=args.dpi
lgdw=args.lgdw
lgdfsz=args.lgdfsz
date=args.date
fslope=args.slopeVal
fcurve=args.curveVal
import matplotlib.ticker as ticker
# for col in range(num_cols_fig-1):
#     ax[0,col].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#     ax[0,col].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
#     ax[0,col].yaxis.set_major_locator(ticker.MultipleLocator(2))
#     ax[0,col].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
#     ax[1,col].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
#     ax[1,col].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
#     ax[1,col].yaxis.set_major_locator(ticker.MultipleLocator(10))
#     ax[1,col].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))




lines = ["solid","dashed",":"]

# ax[0].set_xlabel(r'Equivalence Ratio')

################################ FS-VS-PHI #######################################
# Flame speed across a range of phi, with lines for 1, 10, and 20 bar
numcols=3
col=0
colspacing=0.5
bbval=(0.38,0.01)
lgd_loc='lower center'
P_ls = [1,10,20]
# P_ls = [1,10]
alpha = 1.0
path=f'C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\GubbiResults_vsPhi_'+date+f' (slope={fslope} curve={fcurve})\\'
for i, P in enumerate(P_ls):
    ax[col].plot(0, 0, '.', color='white',markersize=0.1,label=f'{P} bar')  # dummy handle to provide label to lgd column
    label=f'Alzueta'
    dataset=pd.read_csv(path+f'Alzueta_{P}bar_data_{alpha}alpha.csv')
    ax[col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[i],color="xkcd:grey",zorder=30,label=label)
    label=f'LMR-R'
    dataset=pd.read_csv(path+f'LMR-R_{P}bar_data_{alpha}alpha.csv')
    ax[col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[i],color='xkcd:purple',zorder=90,label=label)
    label=f"Mei"
    dataset=pd.read_csv(path+f'Mei_{P}bar_data_{alpha}alpha.csv')
    ax[col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[i],color="xkcd:teal",zorder=60,label=label)

ax[col].set_title(r"$\phi$-dependence for NH$_3$/air",fontsize=8)
ax[col].set_xlim([0.6001, 1.7999])
ax[col].set_xlabel(r'Equivalence ratio',fontsize=7)
ax[col].set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]',fontsize=7)
# ax[col].set_ylim([0.001, 13.9])
# ax[col].annotate('100% NH$_3$\n(760 torr)', xy=(0.97, 0.97), xycoords='axes fraction',ha='right', va='top',fontsize=7)
ax[0].legend(fontsize=lgdfsz, frameon=False, loc=lgd_loc,handlelength=lgdw,ncols=numcols,columnspacing=colspacing,bbox_to_anchor=bbval)



################################ FS-VS-P #######################################
# Flame speed across a range of P, with lines for 1.22 phi
numcols=1
col = 1
colspacing=0.5
lgd_loc='upper right'
P_ls = [1,10,20]
alpha = 1.0
path=f'C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\GubbiResults_vsP_'+date+f' (slope={fslope} curve={fcurve})\\'
# ax[col].plot(0, 0, '.', color='white',markersize=0.1,label=f'{P} bar')  # dummy handle to provide label to lgd column
label=f'Alzueta'
dataset=pd.read_csv(path+f'Alzueta_1.22phi_data.csv')
ax[col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle="solid",color="xkcd:grey",zorder=30,label=label)
label=f'LMR-R'
dataset=pd.read_csv(path+f'LMR-R_1.22phi_data.csv')
ax[col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle="solid",color='xkcd:purple',zorder=30,label=label)
label=f"Mei"
dataset=pd.read_csv(path+f'Mei_1.22phi_data.csv')
ax[col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle="solid",color="xkcd:teal",zorder=30,label=label)

ax[col].set_title(r"P-dependence for NH$_3$/air ($\phi$=1.22)",fontsize=8)
ax[col].set_xlim([0.0001, 50])
ax[col].legend(fontsize=lgdfsz, frameon=False, loc=lgd_loc,handlelength=lgdw,ncols=numcols,columnspacing=colspacing)

# fig.text(.08, 0.5, r'Burning velocity [cm $\rm s^{-1}$]', ha='center', va='center',rotation=90)
ax[col].set_xlabel(r'Pressure [bar]',fontsize=7)

ax[0].tick_params(axis='both',direction='in')
ax[1].tick_params(axis='both',direction='in')

if save_plots == True:
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\"+name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\"+name+'.png', dpi=1000, bbox_inches='tight')

# plt.show()     