# python 'C://Users//pjsin//Documents//cantera//burkelab_SimScripts//simulateK_vs_T.py'

#Validation script
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

models = {    
          'A1':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR.yaml", 
          'A3':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR_Ax3.yaml", 
          'A5':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR_Ax5.yaml", 
          'A10':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR_Ax10.yaml", 
          'A20':"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\alzuetamechanism_LMRR_allAR_Ax20.yaml", 
          }

reactions = ["H + O2 (+M) <=> HO2 (+M)", "H + O2 <=> O + OH"]
pltcolours = ["xkcd:grey", "xkcd:teal", 'xkcd:purple','olive', 'brown','goldenrod']
plt.figure(figsize=(6, 4))
for i, m in enumerate(list(models.keys())):
    gas = ct.Solution(list(models.values())[i])
    Temp = np.linspace(300,3000,70)
    # Temp=[1000]
    # Pres = np.logspace(-2,2,5)
    Pres = [10] # units: Pa
    # Pres=[1e-6,1e-5,1e-4,0.001,0.01,0.1, 1, 10, 100, 1000, 10000, 100000, 1e6, 1e7, 1e8]
    for q, R in enumerate(reactions):
        k_list=[]
        for j, P in enumerate(Pres):
            temp_list = []
            for k,T in enumerate(Temp):
                # gas.TPX = T,P,{'H2O':0.5,'Ar':0.5}
                gas.TPX = T,P*ct.one_atm,'H2:0.1071, O2:0.1785, He:0.7144'
                # print(gas.TPX)
                rc = gas.forward_rate_constants[gas.reaction_equations().index(R)]
                k_list.append(rc)
                # temp_list.append(rc)
            # k_list.append(temp_list)  
            # print(k_list[j]) 
        if q==0:
            plt.semilogy(Temp,k_list,label=m+" "+R,linewidth=0.7,linestyle="-",color=pltcolours[i])
        if q==1:
            plt.semilogy(Temp,k_list,label=m+" "+R,linewidth=0.7,linestyle="--",color=pltcolours[i])
plt.title("10 atm")
plt.legend(fontsize=7,frameon=False,loc='lower right', handlelength=1.0)
plt.xlabel('Temperature [K]')
plt.ylabel('Rate constant (k)')

path=f'C:\\Users\\pjsin\\Documents\\cantera\\burkelab_SimScripts\\figures\\'
plt.savefig(path+f'k_vs_T_allAR.png', dpi=300)