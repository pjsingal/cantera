
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
mpl.rc('font',family='Times New Roman')
mpl.rcParams['mathtext.fontset'] = 'stix'

mpl.rcParams['font.size'] = 8
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 7
from matplotlib.legend_handler import HandlerTuple
plt.rcParams['axes.labelsize'] = 8
mpl.rcParams['xtick.major.width'] = 0.5  # Width of major ticks on x-axis
mpl.rcParams['ytick.major.width'] = 0.5  # Width of major ticks on y-axis
mpl.rcParams['xtick.minor.width'] = 0.5  # Width of minor ticks on x-axis
mpl.rcParams['ytick.minor.width'] = 0.5  # Width of minor ticks on y-axis
mpl.rcParams['xtick.major.size'] = 2.5  # Length of major ticks on x-axis
mpl.rcParams['ytick.major.size'] = 2.5  # Length of major ticks on y-axis
mpl.rcParams['xtick.minor.size'] = 1.5  # Length of minor ticks on x-axis
mpl.rcParams['ytick.minor.size'] = 1.5  # Length of minor ticks on y-axis

save_plots = True
fig, ax = plt.subplots(2,3,figsize=(6.5, 6))

import matplotlib.ticker as ticker


lw=0.7
mw=0.5
msz=2.5
dpi=1000
lgdw=0.6
lgdfsz=5

# alpha_list=[1.0,0.8,0.6,0.4,0.2,0.0]
# idxs=[(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)]
alpha_list=[1]
idxs=[(0,0)]
# pct=["1pt0","0pt8","0pt6","0pt4","0pt2"]
plt.subplots_adjust(wspace=0.3,hspace=0.3)
for x, alpha in enumerate(alpha_list):
    
    # ax[idxs[x]].yaxis.set_major_locator(ticker.MultipleLocator(5))
    # ax[idxs[x]].xaxis.set_major_locator(ticker.MultipleLocator(50))
    # ax[idxs[x]].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    # ax[idxs[x]].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))

    models = {
            'LMR-R':"test/data/alzuetamechanism_LMRR.yaml", 
            'Alzueta':"test/data/alzuetamechanism.yaml",                 
            'Ar':"test/data/alzuetamechanism_LMRR_allAR.yaml",
            r'H$_2$O':"test/data/alzuetamechanism_LMRR_allH2O.yaml",
            }
            
    # def plotPoints(fname, label, shape,color,x):
    #     dataset = pd.read_csv(fname)
    #     NH3_list = np.divide(dataset.iloc[:,0],100)
    #     ox_frac_list = np.subtract(1,NH3_list)
    #     O2_list = np.multiply(ox_frac_list, 0.21)
    #     phi_list = np.divide(np.divide(NH3_list,O2_list),np.divide(4,3))
    #     ax[x].plot(phi_list,dataset.iloc[:,1],marker=shape,fillstyle='none',markersize=3.5,markeredgewidth=0.5,linestyle='none',color=color,label=label)
        
    # path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"
    # # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\50torr.csv','50 torr','o','k')
    # # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\100torr.csv','100 torr','^','k')
    # # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\250torr.csv','250 torr','v','k')
    # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv','Ronney','o','k',x)
    # # plotPoints(path+'\\6 FS NH3 (Stagni-Ronney)\\1500torr.csv','1500 torr','D','k')

    # dataset=pd.read_csv(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv')
    # ax.plot(dataset.iloc[:,0],dataset.iloc[:,1]*100,marker='o',markersize=7,linewidth=3,fillstyle='none',linestyle='none',color='k',label='Ronney')





    path="C:\\Users\\pjsin\\Documents\\cantera\\BurkeSongResults_May28\\"
    dataset=pd.read_csv(path+f'Alzueta_0_data_{alpha}alpha.csv')
    ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color="xkcd:grey",label='Alzueta')

    dataset=pd.read_csv(path+f'Ar_0_data_{alpha}alpha.csv')
    ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color='r',label='Ar')

    dataset=pd.read_csv(path+f'H2O_0_data_{alpha}alpha.csv')
    ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color='b',label=r'H$_2$O')

    dataset=pd.read_csv(path+f'LMR-R_0_data_{alpha}alpha.csv')
    ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],linewidth=lw,color='xkcd:purple',label='LMR-R')

    ax[idxs[x]].set_title(f'alpha,NH3 = {alpha}')

    path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\Graph Reading\\"

    if x==0:
        dataset = pd.read_csv(path+'\\6 FS NH3 (Stagni-Ronney)\\760torr.csv')
        NH3_list = np.divide(dataset.iloc[:,0],100)
        ox_frac_list = np.subtract(1,NH3_list)
        O2_list = np.multiply(ox_frac_list, 0.21)
        phi_list = np.divide(np.divide(NH3_list,O2_list),np.divide(4,3))
        ax[idxs[x]].plot(phi_list,dataset.iloc[:,1],marker='o',fillstyle='none',markersize=msz,markeredgewidth=0.5,linestyle='none',color='k',label='Ronney')

        dataset = pd.read_csv(path+f'\\AlzuetaFig15\\1pt0_NH3.csv')
        ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='s',markersize=msz,markeredgewidth=0.5,linestyle='none',color='b',label='Alz Sim (Graph Read)')
        ax[idxs[x]].legend(fontsize=lgdfsz, frameon=False, loc='lower right') 
    
    if x==2:
        dataset = pd.read_csv(path+f'\\Han\\han_0pt6_NH3.csv')
        ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='s',fillstyle='none',markersize=msz,markeredgewidth=0.5,linestyle='none',color='k',label='Han')
        dataset = pd.read_csv(path+f'\\Wang\\wang_0pt6_NH3.csv')
        ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='x',fillstyle='none',markersize=msz,markeredgewidth=0.5,linestyle='none',color='k',label='Wang')
        ax[idxs[x]].legend(fontsize=lgdfsz, frameon=False, loc='lower right') 

    if x==3:
        dataset = pd.read_csv(path+f'\\Wang\\wang_0pt4_NH3.csv')
        ax[idxs[x]].plot(dataset.iloc[:,0],dataset.iloc[:,1],marker='x',fillstyle='none',markersize=msz,markeredgewidth=0.5,linestyle='none',color='k',label='Wang')
        ax[idxs[x]].legend(fontsize=lgdfsz, frameon=False, loc='lower right') 

    ax[idxs[x]].set_ylabel(r'Burning velocity [cm $\rm s^{-1}$]')
    ax[idxs[x]].set_xlabel(r'Equivalence Ratio')
    ax[idxs[x]].tick_params(axis='both', direction="in")
    ax[idxs[x]].tick_params(axis='both', which='minor', direction="in")
    ax[idxs[x]].set_xlim([0.6, 1.8])

name = f'ronney_flamespeed_allAlpha_Mar04'
if save_plots == True:
    plt.savefig(name+'.pdf', dpi=1000, bbox_inches='tight')
    plt.savefig(name+'.png', dpi=dpi, bbox_inches='tight')

# plt.show()     