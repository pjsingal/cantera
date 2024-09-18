# To run file: 
# python "C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\simulateflamespeedGubbi.py"

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
import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--gridsz', type=int, help="gridsz = ", default=10)
parser.add_argument('--date', type=str, help="sim date = ",default='May28')
parser.add_argument('--slopeVal', type=float, help="slope value = ")
parser.add_argument('--curveVal', type=float, help="curve value = ")
parser.add_argument('--transport', type=str, help="flame transport = ")
args = parser.parse_args()

############# CHANGE THESE ####################################################################################
gridsz = args.gridsz
date=args.date



# fuel_list = np.linspace(0.14,0.5,gridsz) #fuel mole fractions
# fuel_list = np.linspace(0.14,0.4,gridsz)
# alpha_list = [1.0,0.8,0.6,0.4,0.2,0.0]
# a_st = [0.75,0.7,0.65,0.6,0.55,0.5]
p_list = np.linspace(1,20,gridsz)

def widthFit(p):
   return round(float(2.0717*np.exp(-0.2586*p)),4)
widths = []
for p in p_list:
   widths.append(widthFit(p))

alpha_list = [1.0]
a_st = [0.75]
# alpha_list = [0.6]
# a_st = [0.65]
Tin = [296]  # unburned gas temperature [K]
# Tin = [298]  # unburned gas temperature [K]

fratio=3
fslope=args.slopeVal #should be low enough that the results don't depend on the value. 0.02 for both is a good place to start. Try 0.01 and 0.05 and see if there are any differences
fcurve=args.curveVal
ftransport=args.transport # 'multicomponent' or 'mixture-averaged'
# models = {    
#           'Alzueta':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism.yaml",
#           'Mei':'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Ammonia\\Mei-2019\\mei-2019.yaml',
#           'LMR-R':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_extraColliders.yaml",
#           }
models = {    
          'Alzueta':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism.yaml",
          'Mei':'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Mei-2019\\mei-2019.yaml',
          'LMR-R':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_extraColliders.yaml",
          'Glarborg':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Glarborg-2018\\glarborg-2018.yaml",
          'Zhang':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Zhang-2017\\zhang-2017.yaml",
          'Otomo':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Otomo-2018\\otomo-2018.yaml",
          'Stagni':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Stagni-2020\\stagni-2020.yaml",
          'Shrestha':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Shrestha-2021\\shrestha-2021.yaml",
          'Han':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Han-2021\\han-2021.yaml",
          'Cornell':"G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Cornell-2024\\cornell-2024.yaml",
          }
###############################################################################################################


phi = 1.22
# width = 3  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

T_fuel = 300
T_air = 650

def save_to_csv(filename, data):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(data)

def cp(T,P,X):
  gas_stream = ct.Solution(list(models.values())[k])
  gas_stream.TPX = T, P*1e5, {X:1}
  return gas_stream.cp_mole # [J/kmol/K]

for x, alpha in enumerate(alpha_list):
    for k, m in enumerate(models):
      mbr = []
      for i, p in enumerate(p_list):
        gas = ct.Solution(list(models.values())[k])
        cp_fuel = cp(T_fuel,p,'NH3') # [J/kmol/K]
        cp_o2 = cp(T_air,p,'O2') # [J/kmol/K]
        cp_n2 = cp(T_air,p,'N2') # [J/kmol/K]
        
        x_fuel = (phi*(1/0.75)*0.21)/(1+phi*(1/0.75)*0.21)
        x_o2 = 0.21*(1-x_fuel)
        x_n2 = 0.79*(1-x_fuel)
        T_mix = (x_fuel*cp_fuel*T_fuel+(x_o2*cp_o2+x_n2*cp_n2)*T_air)/(x_fuel*cp_fuel+ x_o2*cp_o2 + x_n2*cp_n2)
        
        gas.TPX= T_mix, p*1e5, {'NH3':x_fuel,'O2':x_o2,'N2':x_n2}
        f = ct.FreeFlame(gas, width=widths[i])
        f.set_refine_criteria(ratio=fratio, slope=fslope, curve=fcurve)
        f.transport_model = ftransport
        f.soret_enabled = True
        f.solve(loglevel=loglevel, auto=True)
        mbr.append(f.velocity[0] * 100) # cm/s

      # Save phi_list and mbr to CSV
      path=f'C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\GubbiResults_vsP_'+date+f' (slope={fslope} curve={fcurve})'
      os.makedirs(path,exist_ok=True)
      csv_filename =path+f'\\{m}_{phi}phi_data.csv'
      data = zip(p_list, mbr)
      save_to_csv(csv_filename, data)