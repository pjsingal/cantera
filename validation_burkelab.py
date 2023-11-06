#Validation script
import sys, os
sys.path.append("C:\\Users\\pjsin\\Documents\\cantera\\build\\python")
import cantera as bklabct

gas = bklabct.Solution('C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\kineticsfromscratch_LMRtest.yaml')
print(gas.species_names)
