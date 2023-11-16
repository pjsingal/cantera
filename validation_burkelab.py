#Validation script
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as bklabct
import numpy as np
import matplotlib.pyplot as plt
bklabct.print_stack_trace_on_segfault()

file = 'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\kineticsfromscratch_LMRtest.yaml'
reactions = ['H + O2 <=> HO2']
gas = bklabct.Solution(file)
Temp = np.linspace(250,2500,50)
# Pres = np.logspace(-2,2,5)
Pres=[1]
for i, R in enumerate(reactions):
    k_list=[]
    for j, P in enumerate(Pres):
        temp_list = []
        for k,T in enumerate(Temp):
            gas.TPX = T,P,{'H2O':1}
            rc = gas.forward_rate_constants[gas.reaction_equations().index(R)]
            temp_list.append(rc)
        k_list.append(temp_list)  
        print(k_list[j]) 
    plt.figure()
    plt.title(R)
    for j,P in enumerate(Pres):
        plt.semilogy(Temp,k_list[j],label=str(P)+' atm')    
    plt.legend()
    plt.xlabel('Temperature [K]')
    plt.ylabel('k')
    
    plt.savefig("reaction1",bbox_inches="tight")
    plt.show()
    




