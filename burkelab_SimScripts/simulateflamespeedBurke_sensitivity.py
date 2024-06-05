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

# Function to save data to CSV
def save_to_csv(filename, data):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)

# # Define the path for the log file
# log_path = 'C://Users//pjsin//Documents//cantera//burkelab_SimScripts//burkesong_log.txt'
# sys.stdout = open(log_path, 'w')
# sys.stderr = open(log_path, 'a')

name = 'MBR_BurkeSong'

# models = {
#           'Alzueta':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism.yaml",            
#           'Ar':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR.yaml",
#           r'H$_2$O':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allH2O.yaml",
#           'LMR-R':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR.yaml", 
#           }

# models = {    
#           'Alzueta':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism.yaml",  
#           'Alzueta-300K':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_epsNH3_T=300K.yaml",  
#           'Alzueta-2000K':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_epsNH3_T=2000K.yaml",            
#           'Ar':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR.yaml",
#           r'H$_2$O':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allH2O.yaml",
#           'LMR-R':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR.yaml", 
#           }

models = {    
          'A1':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR.yaml", 
          'A3':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR_Ax3.yaml", 
          'A5':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR_Ax5.yaml", 
          'A10':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR_Ax10.yaml", 
          'A20':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR_Ax20.yaml", 
          }


parser = argparse.ArgumentParser()
parser.add_argument('--gridsz', type=int, help="gridsz = ", default=10)
parser.add_argument('--date', type=str, help="sim date = ",default='May28')
parser.add_argument('--slopeVal', type=float, help="slope value = ")
parser.add_argument('--curveVal', type=float, help="curve value = ")
parser.add_argument('--transport', type=str, help="flame transport = ")
args = parser.parse_args()

# gridsz = 10
gridsz = args.gridsz
date=args.date
fratio=3
fslope=args.slopeVal
fcurve=args.curveVal
ftransport=args.transport

colors = ['xkcd:purple','r','b']

for i, m in enumerate(list(models.keys())):
    plt.figure(figsize=(6, 4))
    p = 10  # pressure [atm]
    T = 300.0  # unburned gas temperature [K]
    reactants = 'H2:0.1071, O2:0.1785, He:0.7144'  # premixed gas composition
    width = 0.03  # m
    loglevel = 1  # amount of diagnostic output (0 to 8)

    gas = ct.Solution(list(models.values())[i])
    gas.TPX = T, p*ct.one_atm, reactants
    f = ct.FreeFlame(gas, width=width)
    f.set_refine_criteria(ratio=fratio, slope=fslope, curve=fcurve)
    f.transport_model = ftransport
    f.soret_enabled = True
    f.solve(loglevel=loglevel, auto=True)

    sensitivities = pd.DataFrame(index=gas.reaction_equations(), columns=["base_case"])
    sens = f.get_flame_speed_reaction_sensitivities()
    for r in range(gas.n_reactions):
        sensitivities.iloc[r,0]=sens[r]
    
    # Reaction mechanisms can contains thousands of elementary steps. Choose a threshold
    # to see only the top few
    threshold = 0.03

    # For plotting, collect only those steps that are above the threshold
    # Otherwise, the y-axis gets crowded and illegible
    sensitivities_subset = sensitivities[sensitivities["base_case"].abs() > threshold]
    reactions_above_threshold = (
        sensitivities_subset.abs().sort_values(by="base_case", ascending=False).index
    )
    sensitivities_subset.loc[reactions_above_threshold].plot.barh(
        title=f"allAR_{m} at 10 atm", legend=None
    )
    plt.gca().invert_yaxis()

    plt.rcParams.update({"axes.labelsize": 10})
    plt.xlabel(r"Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$");

    plt.tight_layout()  # Adjust the layout to prevent clipping
    # Uncomment the following to save the plot. A higher than usual resolution (dpi) helps
    path=f'C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\BurkeSongResults_'+date+f' (sensitivity)\\'
    os.makedirs(path,exist_ok=True)
    plt.savefig(path+f'allAR_{m}_10atm.png', dpi=300)

# sys.stdout = sys.__stdout__
# sys.stderr = sys.__stderr__