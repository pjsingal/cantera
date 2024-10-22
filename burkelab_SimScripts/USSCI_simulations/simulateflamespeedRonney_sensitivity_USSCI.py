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

def save_to_csv(filename, data):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)

parser = argparse.ArgumentParser()
parser.add_argument('--date', type=str, help="sim date = ",default='May28')
parser.add_argument('--slopeVal', type=float, help="slope value = ",default=-1)
parser.add_argument('--curveVal', type=float, help="curve value = ",default=-1)
parser.add_argument('--transport', type=str, help="flame transport = ")
args = parser.parse_args()

date=args.date
fratio=3
fslope=args.slopeVal
fcurve=args.curveVal
ftransport=args.transport

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

P=760
alpha_list = [1.0,0.6]
a_st = [0.75,0.65]
# fuel_list = np.linspace(0.14,0.4,gridsz)
fuel_frac = [0.21873,0.244168] # fuel fracs needed to get phi=1 for each mixture
width = 0.03  # m
loglevel = 0  # amount of diagnostic output (0 to 8)
Tin = [296,298]  # unburned gas temperature [K]

for z, n in enumerate(models):
    for x, alpha in enumerate(alpha_list):
        for k, m in enumerate(models[n]):
            print(f"{m} ({alpha}% NH3/{1-alpha}% H2, {Tin[x]}K, 1 atm, $\phi$=1)")
            # plt.figure(figsize=(6, 4))
            gas = ct.Solution(list(models[n].values())[k])
            NH3 = alpha*fuel_frac[x]
            H2 = (1-alpha)*fuel_frac[x]
            ox_frac = 1 - fuel_frac[x] # oxidizer fraction
            O2 = ox_frac*0.21
            N2 = ox_frac*0.79
            X = {'NH3':NH3,'H2':H2,'O2':O2,'N2':N2}
            gas.TPX = Tin[x], (P/760)*ct.one_atm, X
            f = ct.FreeFlame(gas, width=width)
            f.set_refine_criteria(ratio=fratio, slope=fslope, curve=fcurve)
            f.transport_model = ftransport
            f.soret_enabled = True
            f.solve(loglevel=loglevel, auto=True)
            # mbr = f.velocity[0] * 100 # cm/s
            Su0 = f.velocity[0] # m/s
            print(f"Su0 = {Su0} m/s")
            sensitivities = pd.DataFrame(index=gas.reaction_equations(), columns=["base_case"])
            dk = 0.1
            for r in range(gas.n_reactions):
                gas.set_multiplier(1.0)
                gas.set_multiplier(1+dk,r)
                f.solve(loglevel=0, refine_grid=False, auto=False)
                Su = f.velocity[0]
                sensitivities.iloc[r,0]=(Su-Su0)/(Su0*dk)
            gas.set_multiplier(1.0)
            threshold = 0.03
            sensitivities_subset = sensitivities[sensitivities["base_case"].abs() > threshold]
            print("Sensitivity analysis complete. Now saving to CSV...")
            path=f'burkelab_SimScripts\\USSCI_simulations\\data\\Ronney\\sensitivity\\'+args.date
            os.makedirs(path,exist_ok=True)
            csv_filename =path+f'\\{n}_{m}_760torr_{Tin[x]}K_{alpha}alpha_1phi.csv'
            data = zip(sensitivities_subset.index,sensitivities_subset["base_case"])
            save_to_csv(csv_filename, data)
            print("Success!")