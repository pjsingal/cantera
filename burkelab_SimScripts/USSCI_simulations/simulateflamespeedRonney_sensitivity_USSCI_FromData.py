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
parser.add_argument('--transport', type=str, help="flame transport = ")
args = parser.parse_args()

fig, ax = plt.subplots(1,2,figsize=(args.figwidth, args.figheight))
plt.rcParams.update({"axes.labelsize": 10})
save_plots = True
lw=args.lw
mw=args.mw
msz=args.msz
gridsz = args.gridsz
date=args.date
fratio=3
fslope=args.slopeVal
fcurve=args.curveVal
ftransport=args.transport

patterns = ['','x']*9
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

name = 'Ronney_flamespeed_sensitivity_multimech'
save_plots = True
P=760
alpha_list = [1.0,0.6]
a_st = [0.75,0.65]
# fuel_list = np.linspace(0.14,0.4,gridsz)
fuel_frac = [0.21873,0.244168] # fuel fracs needed to get phi=1 for each mixture
width = 0.03  # m
loglevel = 0  # amount of diagnostic output (0 to 8)
Tin = [296,298]  # unburned gas temperature [K]


for x, alpha in enumerate(alpha_list):
    reactionList = []
    for z, n in enumerate(models):
        for k, m in enumerate(models[n]):
            path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\Ronney_Sensitivity\\"+f"{date}_data\\"
            sensitivities_subset=pd.read_csv(path+f"{m}_{alpha}alpha_{Tin[x]}K_1atm_1phi"+'.csv',header=None)
            reactionList += sensitivities_subset.iloc[:,0].to_list()
        reactionList = list(set(reactionList))
        masterDict={}
        for i, reaction in enumerate(reactionList):
            masterDict[reaction] = []

        for k, m in enumerate(models):
            for i, reaction in enumerate(reactionList):
                path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\Ronney_Sensitivity\\"+f"{date}_data\\"
                sensitivities_subset=pd.read_csv(path+f"{m}_{alpha}alpha_{Tin[x]}K_1atm_1phi"+'.csv',header=None)
                sensitivities_subset_sorted = sensitivities_subset.sort_values(by=1, ascending=False)
                rxns = sensitivities_subset_sorted.iloc[:,0].to_list()
                vals = sensitivities_subset_sorted.iloc[:,1].to_list()
                flag = 0
                for z, rxn, in enumerate(rxns):
                    if rxn==reaction:
                        masterDict[reaction].append(vals[z])
                        flag=1
                if flag == 0: # the master list reaction does not exist in the model list
                    masterDict[reaction].append(0)

        sumDict = {}
        for key in masterDict.keys():
            sumDict[key]=abs(sum(masterDict[key]))
        sumDict = dict(sorted(sumDict.items(), key=lambda item: item[1],reverse=True))
        weighted_reactionList = list(sumDict.keys())
        modelData=[]
        for k, m in enumerate(models):
            modelDict={}
            for i, reaction in enumerate(weighted_reactionList):
                path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\Ronney_Sensitivity\\"+f"{date}_data\\"
                sensitivities_subset=pd.read_csv(path+f"{m}_{alpha}alpha_{Tin[x]}K_1atm_1phi"+'.csv',header=None)
                sensitivities_subset_sorted = sensitivities_subset.sort_values(by=1, ascending=False)
                rxns = sensitivities_subset_sorted.iloc[:,0].to_list()
                vals = sensitivities_subset_sorted.iloc[:,1].to_list()
                flag = 0
                for z, rxn, in enumerate(rxns):
                    if rxn==reaction:
                        modelDict[reaction]=vals[z]
                        flag=1
                if flag == 0: # the master list reaction does not exist in the model list
                    modelDict[reaction]=0
            modelData.append(modelDict)
        offsets = [-0.25,-0.15,-0.05,0,0.05,0.15,0.25]
        width = 0.6
        num_models = len(list(models.keys()))
        offsets = np.linspace(-width * (num_models - 1) / 2, width * (num_models - 1) / 2, num_models)

        for k, m in enumerate(models):
            mdl_dict = modelData[k]
            rxns = list(mdl_dict.keys())[:11]
            vals = list(mdl_dict.values())[:11]
            y = np.arange(len(rxns)) # the label locations
            new_y = [5*i for i in y]
            ax[x].barh(new_y+offsets[k], vals, width, label=f"{m}",color=colours[k],hatch=patterns[k])
            ax[x].set_xlabel(r"Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$",fontsize=6)
            ax[x].set_title(f'{alpha}% NH3/{1-alpha}% H2, {Tin[x]}K, 1 atm, '+r'$\phi$=1',fontsize=8)
            ax[x].set_yticks(new_y, rxns, fontsize=6)
            ax[x].invert_yaxis()
                
            print("Subplot successfully generated!")

ax[1].legend(fontsize=6, frameon=True, loc="lower right",handlelength=2.5) 
pos = ax[1].get_position()
ax[1].set_position([pos.x0*1.35, pos.y0, pos.width, pos.height])

            
path="C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\Ronney_Sensitivity\\"
name="sensitivity"+f"{date}"
plt.savefig(path+name+'.pdf', dpi=1000, bbox_inches='tight')
plt.savefig(path+name+'.svg', dpi=1000, bbox_inches='tight')