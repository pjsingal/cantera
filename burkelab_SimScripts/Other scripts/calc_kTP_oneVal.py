import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

T=1000
P=ct.one_atm

# reaction='H2O2 (+M) <=> 2 OH (+M)'
reaction='H + O2 (+M) <=> HO2 (+M)'
# colliders = ["Ar","H2O","CO2", "N2", "H2O2"]
colliders=["Ar","HE","N2", "H2", "CO2", "NH3", "H2O"]
fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\testLmrr.yaml'
for j, X in enumerate(colliders):
    gas = ct.Solution(fname)
    gas.TPX = T, P, {X:1.0}
    kTP = gas.forward_rate_constants[gas.reaction_equations().index(reaction)]
    print("%.5e (T=%.1f K, P=%.1f atm, X=%s)"%(kTP, T, P/ct.one_atm, X))

# python "C:\Users\pjsin\Documents\cantera\burkelab_SimScripts\calc_kTP_oneVal.py"