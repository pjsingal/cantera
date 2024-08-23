
#%%
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

hfont = {'fontname':'sans-serif','fontweight':550,'fontsize':10,'fontstretch':500}


# %%
#STEP 1: GET RATE CONSTANTS FROM TROE PARAMETERS
def getKfromTroe(P_list,T_list,rxn,ref_collider):
    # file = 'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\sandbox.yaml'
    file="C:\\Users\\pjsin\\OneDrive\\Desktop\\sandbox_cheby.yaml"
    reaction = rxn
    gas = ct.Solution(file)
    def getK(Temp,Pres,X) :
        gas.TPX = Temp,Pres,X
        k = gas.forward_rate_constants[gas.reaction_equations().index(reaction)]
        return k
    # plt.figure()
    k_P=[]
    for P in P_list:
        k_T = []
        for T in T_list:
            k_T.append(getK(T,P,{ref_collider:1}))
        k_P.append(k_T)
        # plt.plot(T_list,np.log(k_T),label=str(P/101325)+" atm")
    # plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    # plt.title("Troe Rates for "+rxn)
    # plt.xlabel("Temperature")
    # plt.ylabel("logK [units still undetermined]")
    # plt.show()
    return k_P

Pvals = np.logspace(-6,12,num=140)
# print(Pvals)
# P_list=[0.0001*101325,0.001*101325,0.01*101325,0.1*101325,1*101325,10*101325,100*101325,1000*101325,10000*101325]
P_list=[]
for i in range(len(Pvals)):
    P_list.append(Pvals[i]*101325)

T_list=np.linspace(200,2000,50)
k_P = getKfromTroe(P_list,T_list,'H2O2 (+M) <=> 2 OH (+M)','H2O2')

# AR: 1.0, H2O: 7.5, CO2: 1.6, N2: 1.5, H2O2: 7.7}

# %% STEP 2: PERFORM A PLOG FIT FOR AR COLLIDER
def plogFit(P_list,T_list,k_P,rxn):
    def arrhenius(T, A, n, Ea):
        return np.log(A) + n*np.log(T)+ (-Ea/(1.987*T))
    plt.figure(figsize=(8,6))
    # dataset = pd.read_csv(fname)
    for i,p in enumerate(P_list):
        k_data = k_P[i] #rate constants across a temperature range, for a given P
        popt, pcov = curve_fit(arrhenius, T_list, np.log(k_data),maxfev = 2000)
        print(("- {P: %.3e atm, A: %.5e, b: %.5e, Ea: %.5e}")%(p/101325, popt[0],popt[1],popt[2]))
        lnk_fit = arrhenius(T_list,popt[0],popt[1],popt[2])
        plt.plot(T_list,lnk_fit,label=str(p/101325) + ' fit',linestyle='solid')
        plt.scatter(T_list,lnk_fit,marker='x',label=str(p/101325) + ' data')
    plt.title('PLOG Fit for '+rxn)
    plt.xlabel('Temperature [K]')
    plt.ylabel('ln ( k/[M] [cm^6/molecule^2*s] )')
    plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.show()
plogFit(P_list,T_list,k_P,'H + O2 (+M) <=> HO2 (+M)')
# %%
