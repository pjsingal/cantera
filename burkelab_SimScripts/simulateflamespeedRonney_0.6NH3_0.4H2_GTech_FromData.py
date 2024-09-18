
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
parser.add_argument('--paper', type=str, help="paper = ",default='PCI')

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
fig, ax = plt.subplots(1,2,figsize=(args.figwidth, args.figheight))


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
ax[0].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax[0].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
ax[0].yaxis.set_major_locator(ticker.MultipleLocator(2))
ax[0].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
ax[1].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
ax[1].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
ax[1].yaxis.set_major_locator(ticker.MultipleLocator(10))
ax[1].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))


models = {    
          'Alzueta':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism.yaml",
          'LMR-R':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_extraColliders.yaml",
          'Mei':'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Mei-2019\\mei-2019.yaml',
        #   'Glarborg':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Glarborg-2018\\glarborg-2018.yaml",
          'Zhang':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Zhang-2017\\zhang-2017.yaml",
          'Otomo':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Otomo-2018\\otomo-2018.yaml",
          'Stagni':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Stagni-2020\\stagni-2020.yaml",
        #   'Shrestha':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Shrestha-2021\\shrestha-2021.yaml",
          'Han':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Han-2021\\han-2021.yaml",
        # #   'Cornell':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Cornell-2024\\cornell-2024.yaml",
          }
zorders = [90,100,80,70,60,50,40,30,20,10]
colours = ["xkcd:grey","xkcd:purple", "xkcd:teal", "orange", "r", "b", "xkcd:lime green", "xkcd:magenta", "xkcd:navy blue"]
path=f'C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\GTechResults_'+date+f' (slope={fslope} curve={fcurve})\\'

for j, m in enumerate(models):
    label=f'{m}'
    dataset=pd.read_csv(path+f'{m}_0_data_1.0alpha.csv')
    if m == "Alzueta" or m == "LMR-R":
        ax[0].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw*1.3,linestyle="--",color=colours[j],label=label,zorder=zorders[j])
    else:
        ax[0].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle="solid",color=colours[j],label=label,zorder=zorders[j])

# ax[0].set_xlabel(r'Equivalence Ratio')

path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"

dataset = pd.read_csv(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv')
NH3_list = np.divide(dataset.iloc[:,0],100)
ox_frac_list = np.subtract(1,NH3_list)
O2_list = np.multiply(ox_frac_list, 0.21)
phi_list = np.divide(np.divide(NH3_list,O2_list),np.divide(4,3))
ax[0].plot(phi_list,dataset.iloc[:,1],marker='o',fillstyle='none',markersize=msz,markeredgewidth=mw,linestyle='none',color='k',label='Ronney',zorder=100)
ax[0].legend(fontsize=lgdfsz, frameon=False, loc='upper right',handlelength=lgdw,ncols=2,columnspacing=0.5) 

#### SECOND PLOT
path=f'C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\GTechResults_'+date+f' (slope={fslope} curve={fcurve})\\'
for j, m in enumerate(models):
    label=f'{m}'
    dataset=pd.read_csv(path+f'{m}_0_data_0.6alpha.csv')
    if m == "Alzueta" or m == "LMR-R":
        ax[1].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw*1.3,linestyle="--",color=colours[j],label=label,zorder=zorders[j])
    else:
        ax[1].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle="solid",color=colours[j],label=label,zorder=zorders[j])

# ax.set_title(f'{round(alpha*100)}% NH3/{round((1-alpha)*100)}% H2')

#### ADD DATA POINTS

path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"
ax[0].tick_params(axis='both', direction="in")
ax[0].tick_params(axis='both', which='minor', direction="in")

dataset = pd.read_csv(path+f'\\Han\\han_0pt6_NH3.csv')
ax[1].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='s',fillstyle='none',markersize=msz,markeredgewidth=mw,linestyle='none',color='k',label='Han',zorder=100)
dataset = pd.read_csv(path+f'\\Wang\\wang_0pt6_NH3.csv')
ax[1].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='x',fillstyle='none',markersize=msz,markeredgewidth=mw,linestyle='none',color='k',label='Wang',zorder=99)
ax[1].legend(fontsize=lgdfsz, frameon=False, loc='upper right',handlelength=lgdw,ncols=2,columnspacing=0.5) 

ax[1].tick_params(axis='both', direction="in")
ax[1].tick_params(axis='both', which='minor', direction="in")

ax[0].annotate(r'NH$_3$/air'+'\n1 atm, 296 K', xy=(0.03, 0.95), xycoords='axes fraction',ha='left', va='top',fontsize=7)
ax[1].annotate(r'60% NH$_3$/40% H$_2$/air'+'\n1 atm, 298 K', xy=(0.03, 0.95), xycoords='axes fraction',ha='left', va='top',fontsize=7)

ax[0].set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]')
fig.text(0.5, 0, r'Equivalence Ratio', ha='center', va='center',fontsize=args.fszaxlab)
ax[0].set_xlim([0.6001, 2.1])
ax[0].set_ylim([0.001, 11.9999])
ax[1].set_xlim([0.6001, 2.1])
ax[1].set_ylim([0.001, 43])


name = f'ronney_flamespeed_GTech_'+date

if save_plots == True:
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\"+name+'.pdf', dpi=1000, bbox_inches='tight')
    # plt.savefig('burkelab_SimScripts/figures/'+name+'_ESSCI.png', dpi=1000, bbox_inches='tight')
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\"+name+'.svg', dpi=500, bbox_inches='tight')

# plt.show()     