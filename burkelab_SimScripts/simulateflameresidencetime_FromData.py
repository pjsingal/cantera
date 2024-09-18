import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import argparse
import csv
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
parser.add_argument('--xscale', type=str)
parser.add_argument('--yscale', type=str)
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

# f, ax = plt.subplots(1, 3, figsize=(args.figwidth, args.figheight))
f, ax = plt.subplots(1, 1, figsize=(args.figwidth, args.figheight))
import matplotlib.ticker as ticker
plt.subplots_adjust(wspace=0.2)

# if args.xscale=='log':
#     ax[0].xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
#     ax[1].xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
#     ax[2].xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
#     ax[0].set_xlim([0.1,1350])
#     ax[1].set_xlim([0.1,1350])
#     ax[2].set_xlim([0.1,1350])
# else:
#     ax[0].xaxis.set_major_locator(ticker.MultipleLocator(300))
#     ax[1].xaxis.set_major_locator(ticker.MultipleLocator(300))
#     ax[2].xaxis.set_major_locator(ticker.MultipleLocator(300))
#     ax[0].set_xlim([0.1,1350])
#     ax[1].set_xlim([0.1,1350])
#     ax[2].set_xlim([0.1,1350])

# if args.yscale=='log':
#     ax[0].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
#     ax[1].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
#     ax[2].yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
#     ax[0].set_ylim([3,1350])
#     ax[1].set_ylim([3,1350])
#     ax[2].set_ylim([3,1350])
# else:
#     ax[0].yaxis.set_major_locator(ticker.MultipleLocator(100))
#     ax[1].yaxis.set_major_locator(ticker.MultipleLocator(100))
#     ax[2].yaxis.set_major_locator(ticker.MultipleLocator(100))
#     ax[0].set_ylim([-80,699])
#     ax[1].set_ylim([-80,699])
#     ax[2].set_ylim([-80,699])

# ax[0].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[0].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[1].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[1].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[2].xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
# ax[2].yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))

if args.xscale=='log':
    ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
    ax.set_xlim([0.1,1350])
else:
    ax.xaxis.set_major_locator(ticker.MultipleLocator(300))
    ax.set_xlim([0.1,1350])

if args.yscale=='log':
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=3))
    ax.set_ylim([3,1350])
else:
    ax.yaxis.set_major_locator(ticker.MultipleLocator(100))
    ax.set_ylim([-80,699])

ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))


models = [
    {'name': 'Alzueta', 'path': 'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism.yaml'},
    {'name': 'LMR-R', 'path': 'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_extraColliders.yaml'},
          {'name': 'Mei', 'path': 'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Mei-2019\\mei-2019.yaml'},
          {'name': 'Glarborg', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Glarborg-2018\\glarborg-2018.yaml"},
          {'name': 'Zhang', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Zhang-2017\\zhang-2017.yaml"},
          {'name': 'Otomo', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Otomo-2018\\otomo-2018.yaml"},
          {'name': 'Stagni', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Stagni-2020\\stagni-2020.yaml"},
          {'name': 'Shrestha', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Shrestha-2021\\shrestha-2021.yaml"},
          {'name': 'Han', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Han-2021\\han-2021.yaml"},
        #   {'name': 'Cornell', 'path': "G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Cornell-2024\\cornell-2024.yaml"},
        ]
colours = ["xkcd:grey","xkcd:purple", "xkcd:teal", "orange", "r", "b", "xkcd:lime green", "xkcd:magenta", "xkcd:navy blue"]
zorders = [90,100,80,70,60,50,40,30,20,10]
# P_list = [1,10,20] # bar
P_list = [20] # bar

expData=['1bar','10bar','20bar']
expData_eq=['1bar_eq','10bar_eq','20bar_eq']
cutoff=int(1)
T_fuel = 300
T_air = 650
phi = 1.22

lw=0.6
mkrw=0.5
mkrsz=3

def getXNOdry(X_NO,X_O2):
    X_O2dry = 0.15 # 15% O2 dry
    return np.multiply(np.multiply(X_NO,np.divide(0.21-X_O2dry,np.subtract(0.21,X_O2))),1e6) # [ppm, 15% O2 dry]

for i,P in enumerate(P_list):
    # path="G:\\Mon disque\\Columbia\\Burke Lab\\09 NOx Mini-Project\\Graph Reading\\"
    # dat = pd.read_csv(path+expData[i]+'.csv',header=None)
    # ax[i].loglog(dat.iloc[:, 0],dat.iloc[:, 1],marker='o',fillstyle='none',linestyle='none',color='k',markersize=mkrsz,markeredgewidth=mkrw, label=f"Gubbi", zorder=110)
    # dat = pd.read_csv(path+expData_eq[i]+'.csv',header=None)
    # x_eq = np.logspace(np.log10(dat.iloc[0, 0]),np.log10(dat.iloc[1, 0]),num=10)
    # y_eq = dat.iloc[0, 1]*np.ones(len(x_eq))
    # ax[i].loglog(x_eq,y_eq,marker='x',fillstyle='none',linestyle='none',color='k',markersize=mkrsz,markeredgewidth=mkrw, label=f"Gubbi_eq", zorder=110)
    path=f'C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\GubbiResidenceTime_'+args.date
    state_eq=pd.read_csv(path+f'\\Mei_{P}bar_data_equilibrium.csv')
    XNOdry_eq = getXNOdry(list(state_eq['X_NO'])[-1], list(state_eq['X_O2'])[-1])
    # print(state_eq['tau [ms]'])
    tau_eq = list(state_eq['tau [ms]'])[:cutoff*(-1)]
    # ax.loglog(np.multiply(tau_eq,P), XNOdry_eq*np.ones(len(tau_eq)), linestyle=":", linewidth=lw,color='k')
    if args.xscale=='log' and args.yscale=='linear':
        ax.semilogx(np.multiply(tau_eq,P), XNOdry_eq*np.ones(len(tau_eq)), linestyle=":", linewidth=lw,color='k')
    elif args.xscale=='linear' and args.yscale=='log':
        ax.semilogy(np.multiply(tau_eq,P), XNOdry_eq*np.ones(len(tau_eq)), linestyle=":", linewidth=lw,color='k')
    elif args.xscale=='log' and args.yscale=='log':
        ax.loglog(np.multiply(tau_eq,P), XNOdry_eq*np.ones(len(tau_eq)), linestyle=":", linewidth=lw,color='k')
    elif args.xscale=='log' and args.yscale=='log':
        ax.plot(np.multiply(tau_eq,P), XNOdry_eq*np.ones(len(tau_eq)), linestyle=":", linewidth=lw,color='k')

    for j, model in enumerate(models):
        # print(f"Adding {model['name']} to plot...")
        state=pd.read_csv(path+f'\\{model['name']}_{P}bar_data.csv')
        XNOdry = getXNOdry(list(state['X_NO']),list(state['X_O2']))
        XNOdry = XNOdry[:cutoff*(-1)]
        tau =list(state['tau [ms]'])[:cutoff*(-1)]

        if model["name"] == "Alzueta" or model["name"] == "LMR-R":
            if args.xscale=='log' and args.yscale=='linear':
                ax.semilogx(np.multiply(tau,P), XNOdry, linestyle='--', linewidth=lw*1.3,color=colours[j],label=f"{model["name"]}",zorder=zorders[j])
            elif args.xscale=='linear' and args.yscale=='log':
                ax.semilogy(np.multiply(tau,P), XNOdry, linestyle='--', linewidth=lw*1.3,color=colours[j],label=f"{model["name"]}",zorder=zorders[j])
            elif args.xscale=='log' and args.yscale=='log':
                ax.loglog(np.multiply(tau,P), XNOdry, linestyle='--', linewidth=lw*1.3,color=colours[j],label=f"{model["name"]}",zorder=zorders[j])
            elif args.xscale=='log' and args.yscale=='log':
                ax.plot(np.multiply(tau,P), XNOdry, linestyle='--', linewidth=lw*1.3,color=colours[j],label=f"{model["name"]}",zorder=zorders[j])
        else:
            if args.xscale=='log' and args.yscale=='linear':
                ax.semilogx(np.multiply(tau,P), XNOdry, linestyle='solid', linewidth=lw,color=colours[j],label=f"{model["name"]}",zorder=zorders[j])
            elif args.xscale=='linear' and args.yscale=='log':
                ax.semilogy(np.multiply(tau,P), XNOdry, linestyle='solid', linewidth=lw,color=colours[j],label=f"{model["name"]}",zorder=zorders[j])
            elif args.xscale=='log' and args.yscale=='log':
                ax.loglog(np.multiply(tau,P), XNOdry, linestyle='solid', linewidth=lw,color=colours[j],label=f"{model["name"]}",zorder=zorders[j])
            elif args.xscale=='log' and args.yscale=='log':
                ax.plot(np.multiply(tau,P), XNOdry, linestyle='solid', linewidth=lw,color=colours[j],label=f"{model["name"]}",zorder=zorders[j])

        

ax.set_xlabel(r'P*$\tau$ [bar*ms]',fontsize=args.fszxtick)
# ax.set_title("1 bar")
if args.yscale=='linear':
    ax.annotate(r'NH$_3$/air'+'\n'+r'$\phi$=1.22'+'\n'+r'T$_{NH_3}$=300 K'+'\n'+r'T$_{air}$=650 K'+'\n'+r'P=20 bar', xy=(0.97, 0.45), xycoords='axes fraction',ha='right', va='top',fontsize=6)
else:
    ax.annotate(r'NH$_3$/air'+'\n'+r'$\phi$=1.22'+'\n'+r'T$_{NH_3}$=300 K'+'\n'+r'T$_{air}$=650 K'+'\n'+r'P=20 bar', xy=(0.05, 0.35), xycoords='axes fraction',ha='left', va='top',fontsize=6)



# ax[1].set_title("10 bar")
# ax.set_title("20 bar")
ax.set_ylabel('NO [ppm, 15% O$_2$ dry]',fontsize=args.fszxtick)
ax.legend(fontsize=6,frameon=False,loc='upper right',handlelength=lgdw,ncols=3,columnspacing=0.5)
ax.tick_params(axis='both',direction='in')
# ax[1].tick_params(axis='both',direction='in')
# ax[2].tick_params(axis='both',direction='in')
ax.tick_params(axis='both', which='minor', direction="in")
# ax[1].tick_params(axis='both', which='minor', direction="in")
# ax[2].tick_params(axis='both', which='minor', direction="in")
path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\"
plt.savefig(path+f'residencetime_NO_{args.date}_{args.xscale}_{args.yscale}.pdf',dpi=1000, bbox_inches='tight')
plt.savefig(path+f'residencetime_NO_{args.date}_{args.xscale}_{args.yscale}.svg',dpi=500, bbox_inches='tight')
# print("Simulation complete!")