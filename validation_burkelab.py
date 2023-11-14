#Validation script
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as bklabct

gas = bklabct.Solution('C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\kineticsfromscratch_LMRtest.yaml')
gas.TPX= 1000, bklabct.one_atm, {"H2O":1}
# reaction = "H + O2 <=> HO2"
# Rc1=gas.forward_rate_constants[gas.reaction_equations().index(reaction)]
# print(Rc1)

# Create a reactor
reactor = bklabct.IdealGasReactor(gas)

# Create a reactor network
reactor_network = bklabct.ReactorNet([reactor])

ratelist=[]
# Advance the simulation to steady state
time = 0.0
end_time = 1  # Set an end time for the simulation
while time < end_time:
    time = reactor_network.step()

    # Output the rate constant
    rate_constant = reactor.kinetics.forward_rate_constants[0]
    ratelist.append(rate_constant)
print(ratelist)




