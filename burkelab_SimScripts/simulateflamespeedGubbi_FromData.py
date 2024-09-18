
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

# fig, ax = plt.subplots(1,2,figsize=(args.figwidth, args.figheight))
fig, ax = plt.subplots(1,1,figsize=(args.figwidth, args.figheight))
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

ax.xaxis.set_major_locator(ticker.MultipleLocator(4))
ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))

models = {    
          'Alzueta':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism.yaml",
          'LMR-R':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_extraColliders.yaml",
          'Mei':'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Mei-2019\\mei-2019.yaml',
          'Glarborg':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Glarborg-2018\\glarborg-2018.yaml",
          'Zhang':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Zhang-2017\\zhang-2017.yaml",
          'Otomo':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Otomo-2018\\otomo-2018.yaml",
          'Stagni':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Stagni-2020\\stagni-2020.yaml",
          'Shrestha':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Shrestha-2021\\shrestha-2021.yaml",
          'Han':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Han-2021\\han-2021.yaml",
        #   'Cornell':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Cornell-2024\\cornell-2024.yaml",
          }
zorders = [90,100,80,70,60,50,40,30,20,10]
colours = ["xkcd:grey","xkcd:purple", "xkcd:teal", "orange", "r", "b", "xkcd:lime green", "xkcd:magenta", "xkcd:navy blue"]
# colours=[(np.random.rand(), np.random.rand(), np.random.rand()) for _ in range(len(list(models.keys())))]

lines = ["solid","dashed",":"]

# ax[0].set_xlabel(r'Equivalence Ratio')

# ################################ FS-VS-PHI #######################################
# # Flame speed across a range of phi, with lines for 1, 10, and 20 bar
# numcols=3
# col=0
# colspacing=0.5
# bbval=(0.38,0.01)
# lgd_loc='lower center'
# P_ls = [1,10,20]
# # P_ls = [1,10]
# alpha = 1.0
# path=f'C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\GubbiResults_vsPhi_'+date+f' (slope={fslope} curve={fcurve})\\'
# for i, P in enumerate(P_ls):
#     ax[col].plot(0, 0, '.', color='white',markersize=0.1,label=f'{P} bar')  # dummy handle to provide label to lgd column
#     for j, m in enumerate(models):
#         label=f'{m}'
#         dataset=pd.read_csv(path+f'{m}_{P}bar_data_{alpha}alpha.csv')
#         ax[col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[i],color=colours[j],label=label)

# ax[col].set_title(r"$\phi$-dependence for NH$_3$/air",fontsize=8)
# ax[col].set_xlim([0.6001, 1.7999])
# ax[col].set_xlabel(r'Equivalence ratio',fontsize=7)
# ax[col].set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]',fontsize=7)
# # ax[col].set_ylim([0.001, 13.9])
# # ax[col].annotate('100% NH$_3$\n(760 torr)', xy=(0.97, 0.97), xycoords='axes fraction',ha='right', va='top',fontsize=7)
# ax[0].legend(fontsize=lgdfsz, frameon=False, loc=lgd_loc,handlelength=lgdw,ncols=numcols,columnspacing=colspacing,bbox_to_anchor=bbval)



################################ FS-VS-P #######################################
# Flame speed across a range of P, with lines for 1.22 phi
numcols=3
col = 1
colspacing=0.5
lgd_loc='upper right'
alpha = 1.0
path=f'C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\GubbiResults_vsP_'+date+f' (slope={fslope} curve={fcurve})\\'

for j, m in enumerate(models):
    label=f'{m}'
    dataset=pd.read_csv(path+f'{m}_1.22phi_data.csv')
    if m == "Alzueta" or m == "LMR-R":
        ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw*1.3,linestyle="--",color=colours[j],label=label,zorder=zorders[j])
    else:
        ax.plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle="solid",color=colours[j],label=label,zorder=zorders[j])

# ax.set_title(r"P-dependence for NH$_3$/air ($\phi$=1.22)",fontsize=8)
ax.annotate(r'NH$_3$/air'+'\n'+r'$\phi$=1.22'+'\n'+r'T$_{NH_3}$=300 K'+'\n'+r'T$_{air}$=650 K', xy=(0.05, 0.05), xycoords='axes fraction',ha='left', va='bottom',fontsize=6)
ax.set_xlim([0.0001, 19.999])
ax.set_ylim([6, 34])
ax.legend(fontsize=lgdfsz, frameon=False, loc=lgd_loc,handlelength=lgdw,ncols=numcols,columnspacing=colspacing)

# fig.text(.08, 0.5, r'Burning velocity [cm $\rm s^{-1}$]', ha='center', va='center',rotation=90)
ax.set_xlabel(r'Pressure [bar]',fontsize=7)
ax.set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]',fontsize=7)
# ax[0].tick_params(axis='both',direction='in')
ax.tick_params(axis='both',direction='in')

if save_plots == True:
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\"+name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\"+name+'.svg', dpi=500, bbox_inches='tight')

# plt.show()     