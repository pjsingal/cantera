#Validation script
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import matplotlib.pyplot as plt
import cantera as bklabct
import numpy as np

bklabct.print_stack_trace_on_segfault()

def printRateConstant(Temp,Pres,X) :
    # file = 'test/data/kineticsfromscratch_LMRtest.yaml'
    # file = 'test/data/alzuetamechanism_LMRR.yaml'
    file = 'test/data/sandbox.yaml'
    # reactions = ['H + O2 (+M) <=> HO2 (+M)']
    reactions = ['NH2 + NH2 (+M) <=> N2H4 (+M)']
    gas = bklabct.Solution(file)
    for i, R in enumerate(reactions):
        k_list=[]
        for j, P in enumerate(Pres):
            temp_list = []
            for k,T in enumerate(Temp):
                gas.TPX = T,P,X
                rc = gas.forward_rate_constants[gas.reaction_equations().index(R)]
                temp_list.append(rc)
            k_list.append(temp_list)
            print(("%.5e      %s")%(k_list[j][0],str(X)))


# X_list = [{'AR':1.0,'H2O':0.0},{'AR':0.7,'H2O':0.3},{'AR':0.0,'H2O':1.0},
#           {'AR':1.0,'N2':0.0},{'AR':0.7,'N2':0.3},{'AR':0.0,'N2':1.0},
#           {'AR':1.0,'H2':0.0},{'AR':0.7,'H2':0.3},{'AR':0.0,'H2':1.0},
#           {'AR':1.0,'CO2':0.0},{'AR':0.7,'CO2':0.3},{'AR':0.0,'CO2':1.0},
#           {'AR':1.0,'NH3':0.0},{'AR':0.7,'NH3':0.3},{'AR':0.0,'NH3':1.0},
#           {'AR':1.0,'H2O2':0.0},{'AR':0.7,'H2O2':0.3},{'AR':0.0,'H2O2':1.0},
#           {'HE':1.0,'AR':0.0}]
# for X in X_list:
#     printRateConstant([2000],[2*101325],X)

P_list=np.multiply([]
printRateConstant()
