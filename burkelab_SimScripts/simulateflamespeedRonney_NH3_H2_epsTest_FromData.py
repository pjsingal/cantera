
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

num_cols_fig = 3

fig, ax = plt.subplots(2,num_cols_fig,figsize=(args.figwidth, args.figheight))
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
for col in range(num_cols_fig-1):
    ax[0,col].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax[0,col].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
    ax[0,col].yaxis.set_major_locator(ticker.MultipleLocator(2))
    ax[0,col].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    ax[1,col].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax[1,col].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
    ax[1,col].yaxis.set_major_locator(ticker.MultipleLocator(10))
    ax[1,col].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))


if fslope != -1:
    path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\RonneyResults_"+date+f' (slope={fslope} curve={fcurve})\\'
else:
    path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\RonneyResults_"+date+"\\"

lines = ["solid","dotted","dashed"]

# ax[0].set_xlabel(r'Equivalence Ratio')

################################ EPS-DEPENDENCE #######################################
# if args.plot=="eps-dependence":
numcols=1
colspacing=0.5
bbval=(0.44,-0.03)
lgd_loc='lower center'
P_ls = [760]
alpha_ls = [1.0,0.6]
col=0
for i, alpha in enumerate(alpha_ls):
    for j, P in enumerate(P_ls):
        label=f'Alzueta'
        dataset=pd.read_csv(path+f'Alzueta_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color="xkcd:grey",zorder=30,label=label)

        label=f'LMR-R'
        dataset=pd.read_csv(path+f'LMR-R_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color='xkcd:purple',zorder=30,label=label)

        label=f'LMR-R (extra)'
        dataset=pd.read_csv(path+f'LMR-R-extra_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='o',markersize='2.5',fillstyle='none',markeredgewidth=0.5,linestyle='none',color='xkcd:purple',zorder=30,label=label)

        label=r"$\epsilon_{0,ALL}(300K)$"
        dataset=pd.read_csv(path+f'epsALL-300K_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color="b",zorder=30,label=label)

        label=r"$\epsilon_{0,ALL}(2000K)$"
        dataset=pd.read_csv(path+f'epsALL-2000K_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color="r",zorder=30,label=label)

        label=r"$\epsilon_{0,NH3}(300K)$"
        dataset=pd.read_csv(path+f'epsNH3-300K_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color="xkcd:teal",zorder=30,label=label)

        label=r"$\epsilon_{0,NH3}(2000K)$"
        dataset=pd.read_csv(path+f'epsNH3-2000K_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color="orange",zorder=30,label=label)
# ax[0,col].set_title("Effect of including all "+r"$\epsilon_{0,i}(300K)$"+" and " +r"$\epsilon_{0,i}(2000K)$",fontsize=7)
ax[0,col].set_title("Effect of including all Jasper efficiencies",fontsize=7)
ax[0,col].set_xlim([0.6001, 1.7999])
ax[0,col].set_ylim([0.001, 13.9])
ax[1,col].set_xlim([0.6001, 1.7999])
ax[1,col].set_ylim([0.001, 38])
ax[0,col].annotate('100% NH$_3$\n(760 torr)', xy=(0.97, 0.97), xycoords='axes fraction',ha='right', va='top',fontsize=7)
ax[1,col].annotate('60% NH$_3$/40% H$_2$\n(760 torr)', xy=(0.97, 0.97), xycoords='axes fraction',ha='right', va='top',fontsize=7)
try:
    ax[1,col].legend(fontsize=lgdfsz, frameon=False, loc=lgd_loc,handlelength=lgdw,ncols=numcols,columnspacing=colspacing,bbox_to_anchor=bbval) 
except:
    ax[1,col].legend(fontsize=lgdfsz, frameon=False, loc=lgd_loc,handlelength=lgdw,ncols=numcols,columnspacing=colspacing)

################################ EPS-T-DEPENDENCE #######################################
# if args.plot=="eps-T-dependence":
numcols=1
colspacing=0.5
lgd_loc='lower center'
P_ls = [760]
alpha_ls = [1.0,0.6]
col=1
for i, alpha in enumerate(alpha_ls):
    for j, P in enumerate(P_ls):
        label=f'Alzueta'
        dataset=pd.read_csv(path+f'Alzueta_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color="xkcd:grey",zorder=30,label=label)

        label=f'LMR-R'
        dataset=pd.read_csv(path+f'LMR-R_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color='xkcd:purple',zorder=30,label=label)

        label=r"$\epsilon_{0,NH3}(300K)$"
        dataset=pd.read_csv(path+f'epsNH3-300K_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color="xkcd:teal",zorder=30,label=label)

        label=r"$\epsilon_{0,NH3}(1000K)$"
        dataset=pd.read_csv(path+f'epsNH3-1000K_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color="g",zorder=30,label=label)

        label=r"$\epsilon_{0,NH3}(2000K)$"
        dataset=pd.read_csv(path+f'epsNH3-2000K_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color="orange",zorder=30,label=label)
ax[0,col].set_title("Effect of changing "+r"$\epsilon_{0,NH3}(T)$",fontsize=7)
ax[0,col].set_xlim([0.6001, 1.7999])
ax[0,col].set_ylim([0.001, 13.9])
ax[1,col].set_xlim([0.6001, 1.7999])
ax[1,col].set_ylim([8, 38])
ax[0,col].annotate('100% NH$_3$\n(760 torr)', xy=(0.97, 0.97), xycoords='axes fraction',ha='right', va='top',fontsize=7)
ax[1,col].annotate('60% NH$_3$/40% H$_2$\n(760 torr)', xy=(0.97, 0.97), xycoords='axes fraction',ha='right', va='top',fontsize=7)
try:
    ax[0,col].legend(fontsize=lgdfsz, frameon=False, loc=lgd_loc,handlelength=lgdw,ncols=numcols,columnspacing=colspacing,bbox_to_anchor=bbval) 
except:
    ax[0,col].legend(fontsize=lgdfsz, frameon=False, loc=lgd_loc,handlelength=lgdw,ncols=numcols,columnspacing=colspacing)

# ################################ Tin-DEPENDENCE #######################################
# numcols=1
# colspacing=0.5
# bbval=(1.85,1)
# lgd_loc='upper right'
# P_ls = [760]
# alpha_ls = [1.0,0.6]
# # T_ls = [200,300,400,500,600,700]
# T_ls = [300,400,500]
# col=2
# for i, alpha in enumerate(alpha_ls):
#     for j, P in enumerate(P_ls):
#         for k, T in enumerate(T_ls):
#             label=f'Alzueta '+r'($T_{in}=$'+f'{T}K)'
#             dataset=pd.read_csv(path+f'Alzueta_{P}_{T}K_data_{alpha}alpha.csv')
#             ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[k],color="xkcd:grey",zorder=30,label=label)

#             label=f'LMR-R '+r'($T_{in}=$'+f'{T}K)'
#             dataset=pd.read_csv(path+f'LMR-R_{P}_{T}K_data_{alpha}alpha.csv')
#             ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[k],color='xkcd:purple',zorder=30,label=label)

#             label=r"$\epsilon_{0,NH3}(300K)$ "+r'($T_{in}=$'+f'{T}K)'
#             dataset=pd.read_csv(path+f'epsNH3-300K_{P}_{T}K_data_{alpha}alpha.csv')
#             ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[k],color="xkcd:teal",zorder=30,label=label)

#             label=r"$\epsilon_{0,NH3}(2000K)$ "+r'($T_{in}=$'+f'{T}K)'
#             dataset=pd.read_csv(path+f'epsNH3-2000K_{P}_{T}K_data_{alpha}alpha.csv')
#             ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[k],color="orange",zorder=30,label=label)
# ax[0,col].set_title("Effect of T$_{in}$",fontsize=7)
# ax[0,col].set_xlim([0.6001, 1.7999])
# # ax[0,col].set_ylim([1.6, 13.9])
# ax[1,col].set_xlim([0.6001, 1.7999])
# # ax[1,col].set_ylim([4, 49])
# ax[0,col].annotate('100% NH$_3$\n(760 torr)', xy=(0.97, 0.97), xycoords='axes fraction',ha='right', va='top',fontsize=7)
# ax[1,col].annotate('60% NH$_3$/40% H$_2$\n(760 torr)', xy=(0.97, 0.97), xycoords='axes fraction',ha='right', va='top',fontsize=7)
# name = f'ronney_flamespeed_'+date+f'_0.6NH3_0.4H2'+f'_epsTest'
# try:
#     ax[0,col].legend(fontsize=lgdfsz, frameon=False, loc=lgd_loc,handlelength=lgdw,ncols=numcols,columnspacing=colspacing,bbox_to_anchor=bbval) 
# except:
#     ax[0,col].legend(fontsize=lgdfsz, frameon=False, loc=lgd_loc,handlelength=lgdw,ncols=numcols,columnspacing=colspacing)

################################ P-DEPENDENCE #######################################
# if args.plot=="P-dependence":
lines = ["dotted","solid","dashed"]
numcols=1
colspacing=0.5
bbval=(1.45,1)
# bbval=(1.85,1)
lgd_loc='upper right'
P_ls = [250,760,1500]
alpha_ls = [1.0,0.6]
col=2
for i, alpha in enumerate(alpha_ls):
    for j, P in enumerate(P_ls):
        if i==0:
            ax[0,col].plot(0, 0, '.', color='white',markersize=0.1,label=f'{P} torr')  # dummy handle to provide label to lgd column

        label=f'Alzueta'
        dataset=pd.read_csv(path+f'Alzueta_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color="xkcd:grey",zorder=80,label=label)

        label=f'LMR-R'
        dataset=pd.read_csv(path+f'LMR-R_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color='xkcd:purple',zorder=90,label=label)

        label=r"$\epsilon_{0,NH3}(300K)$"
        dataset=pd.read_csv(path+f'epsNH3-300K_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color="xkcd:teal",zorder=30,label=label)

        label=r"$\epsilon_{0,NH3}(2000K)$"
        dataset=pd.read_csv(path+f'epsNH3-2000K_{P}_data_{alpha}alpha.csv')
        ax[i,col].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,linestyle=lines[j],color="orange",zorder=40,label=label)
ax[0,col].set_title("Effect of pressure",fontsize=7)
ax[0,col].set_xlim([0.6001, 1.7999])
ax[0,col].set_ylim([0.001, 13.9])
ax[1,col].set_xlim([0.6001, 1.7999])
ax[1,col].set_ylim([4, 49])
ax[0,col].annotate('100% NH$_3$', xy=(0.97, 0.97), xycoords='axes fraction',ha='right', va='top',fontsize=7)
ax[1,col].annotate('60% NH$_3$/40% H$_2$', xy=(0.97, 0.97), xycoords='axes fraction',ha='right', va='top',fontsize=7)
name = f'ronney_flamespeed_'+date+f'_0.6NH3_0.4H2'+f'_epsTest'
try:
    legend=ax[0,col].legend(fontsize=lgdfsz, frameon=False, loc=lgd_loc,handlelength=lgdw,ncols=numcols,columnspacing=colspacing,bbox_to_anchor=bbval) 
except:
    legend=ax[0,col].legend(fontsize=lgdfsz, frameon=False, loc=lgd_loc,handlelength=lgdw,ncols=numcols,columnspacing=colspacing)
# pos = ax[0,col].get_position()
# ax[0,col].set_position([pos.x0*1.2, pos.y0, pos.width, pos.height])
# pos = ax[1,col].get_position()
# ax[1,col].set_position([pos.x0*1.2, pos.y0, pos.width, pos.height])

for text in legend.get_texts():
    labels=[]
    for j, P in enumerate(P_ls):
        labels.append(f'{P} torr')
    if text.get_text() in labels:
        text.set_fontsize(6)  # Set a larger font size
        text.set_fontweight('bold')  # Make the font bold

fig.text(.08, 0.5, r'Burning velocity [cm $\rm s^{-1}$]', ha='center', va='center',rotation=90)
ax[1,1].set_xlabel(r'Equivalence Ratio')


if save_plots == True:
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\"+name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig("C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\"+name+'.svg', dpi=1000, bbox_inches='tight')

# plt.show()     