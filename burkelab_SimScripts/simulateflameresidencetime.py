import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rc('font',family='Times New Roman')

save_fig = True
nominal_model = 'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Ammonia\\Mei-2019\\mei-2019.yaml'
initial_temp = 300

P = 1e5 # Pa

def gas(T,P,X):
    gas_stream = ct.Solution(nominal_model)
    gas_stream.TPX = T, P, X
    return gas_stream

fuel = gas(300,P,{'NH3':1})
air = gas(650,P,{'N2':0.79,'O2':0.21})

phi = 1.22

# Calculate the moles of ammonia and air needed for the desired equivalence ratio
# For phi = (n_fuel/n_air) / (n_fuel/n_air)stoic
n_fuel = 1.0  # Assume 1 mole of ammonia
n_air_stoic = 4.76  # Stoichiometric air requirement for ammonia combustion (3/2 * O2 + 3.76N2)
n_air = n_air_stoic / phi

mixture = ct.Mixture([(fuel, n_fuel), (air, n_air)])
mixture_equivalence_ratio = mixture.equivalence_ratio(gas_ammonia, gas_air)
print(f"Mixture equivalence ratio: {mixture_equivalence_ratio}")


T_f = 400  
pressure = 1
width = 0.3
X = 'H2:0.5, O2:0.5, N2:1.549'
obs = 'NO'
########################################################################################################################################

########################################################################################################################################
ini_cond = (initial_temp, pressure * ct.one_atm, X)
gas = ct.Solution(nominal_model)
gas.TPX = ini_cond
flame = ct.FreeFlame(gas)
flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
flame.solve(loglevel=0, auto=True)         
T_list = flame.T
P = flame.P
velocity_list = flame.velocity
mole = flame.X
distance_list = flame.grid
S_b = list(velocity_list)[-1]
full_time = [(distance)/S_b for distance in distance_list]
distance_index = min(range(len(T_list)), key=lambda i: abs(T_list[i]-T_f))
tau_f = np.array(distance_list[distance_index]) / S_b
tau_list = [ft - tau_f for ft in full_time] 
distance_f = np.array(distance_list[distance_index])
distance_list_new = [dl - distance_f for dl in distance_list]                 
########################################################################################################################################

########################################################################################################################################
plt.figure()
plt.plot(np.multiply(tau_list,1000), np.multiply(mole[gas.species_names.index(obs)],1e6), linestyle='solid', linewidth=3)
plt.xlabel('Time [ms]')
plt.ylabel(obs + ' Mole Fraction [ppm]')
plt.tick_params(axis = 'x', direction='in')
plt.tick_params(axis='x',labelsize=15)
plt.tick_params(axis = 'y', direction='in')
plt.tick_params(axis='y',labelsize=15)          
if save_fig == True:        
    plt.savefig('flame_' + obs + '.pdf',dpi=1000, bbox_inches='tight')
    plt.savefig('flame_' + obs + '.png',dpi=1000, bbox_inches='tight')        
plt.show()
########################################################################################################################################