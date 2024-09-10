import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
mpl.rc('font',family='Times New Roman')

save_fig = True
nominal_model = 'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Ammonia\\Mei-2019\\mei-2019.yaml'
# initial_temp = 300

P = 1e5 # Pa
obs = 'NO'
width = 2

def cp(T,P,X):
    gas_stream = ct.Solution(nominal_model)
    gas_stream.TPX = T, P, {X:1}
    return gas_stream.cp_mole # [J/kmol/K]


T_fuel = 300
T_air = 650
cp_fuel = cp(T_fuel,P,'NH3') # [J/kmol/K]
cp_o2 = cp(T_air,P,'O2') # [J/kmol/K]
cp_n2 = cp(T_air,P,'N2') # [J/kmol/K]

phi = 1.22
x_fuel = (phi*(1/0.75)*0.21)/(1+phi*(1/0.75)*0.21)
x_o2 = 0.21*(1-x_fuel)
x_n2 = 0.79*(1-x_fuel)
x_air=1-x_fuel
T_mix = (x_fuel*cp_fuel*T_fuel+(x_o2*cp_o2+x_n2*cp_n2)*T_air)/(x_fuel*cp_fuel+ x_o2*cp_o2 + x_n2*cp_n2)
# print(T_mix)

# def h(T,P,X):
#     gas_stream = ct.Solution(nominal_model)
#     gas_stream.TPX = T, P, {X:1}
#     return gas_stream.enthalpy_mole 


# T_fuel = 330
# T_air = 650
# h_fuel = h(T_fuel,P,'NH3')
# h_o2 = h(T_air,P,'O2')
# h_n2 = h(T_air,P,'N2')

# phi = 1.22
# x_fuel = (phi*(1/0.75)*0.21)/(1+phi*(1/0.75)*0.21)
# # print(x_fuel)
# x_o2 = 0.21*(1-x_fuel)
# x_n2 = 0.79*(1-x_fuel)
# x_air=1-x_fuel
# print(h_fuel)
# print(h_o2)
# print(h_n2) 
# h_mix = x_fuel*h_fuel + x_o2*h_o2 + x_n2*h_n2
# print(h_mix)
# gas = ct.Solution(nominal_model)
# # gas.HPX = 300,P,{'NH3':x_fuel,'O2':x_o2,'N2':x_n2}
# gas.HP = h_mix,P
# print(gas.T)



# T_mix = (x_fuel*cp_fuel*T_fuel+(x_o2*cp_o2+x_n2*cp_n2)*T_air)/(x_fuel*cp_fuel+ x_o2*cp_o2 + x_n2*cp_n2)
# # print(T_mix)

# # T_explosion = 1000 # combustor exit temperature

mix1 = ct.Solution(nominal_model)
mix1.TPX = T_mix, P,{'NH3':x_fuel,'O2':x_o2,'N2':x_n2}

# mix1.equilibrate("HP")


flame = ct.FreeFlame(mix1,width=width)
flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
flame.solve(loglevel=0, auto=True)
state = {
    'T': flame.T,
    'P': flame.P,
    'vel': flame.velocity,
    'X': flame.X,
    'd': flame.grid #list of distances [m]
}
S_b = list(state['vel'])[-1] # [m/s]
Rjoule = 8.31446261815324 # [J/K/mol]
full_time = [(distance)/S_b for distance in state['d']] # [s]
conc = np.multiply(state['X'][mix1.species_names.index("NH2")], np.divide(P,np.multiply(Rjoule,state['T'])))
distance_index=np.argmax(conc)
tau_f = np.array(state['d'][distance_index])/S_b
tau_list = [(ft - tau_f)*1000 for ft in full_time] # [ms]
X_NO = state['X'][mix1.species_names.index("NO")]
X_O2 = state['X'][mix1.species_names.index("O2")]
X_O2dry = 0.15 # 15% O2 dry
X_NO_15O2dry = np.multiply(np.multiply(X_NO,np.divide(0.21-X_O2dry,np.subtract(0.21,X_O2))),1e6) # [ppm, 15% O2 dry]
# distance_f = np.array(state['d'][distance_index])
# distance_list_new = [dl - distance_f for dl in state['d']]

########################################################################################################################################
plt.figure()
plt.loglog(tau_list, X_NO_15O2dry, linestyle='solid', linewidth=3)
# sims = pd.DataFrame(np.multiply(tau_list,1000), np.multiply(state['X'][mix1.species_names.index(obs)],1e6))
# print(sims)
path="G:\\Mon disque\\Columbia\\Burke Lab\\09 NOx Mini-Project\\Graph Reading\\"
bar1_data = pd.read_csv(path+'1bar.csv',header=None)

plt.loglog(bar1_data.iloc[:, 0],bar1_data.iloc[:, 1],marker='o',fillstyle='none',linestyle='none',color='k',markersize=3,markeredgewidth=1, label="Gubbi et al.", zorder=110)
# plt.loglog(np.multiply(full_time,1000), np.multiply(state['X'][mix1.species_names.index(obs)],1e6), linestyle='solid', linewidth=3)
# plt.loglog(np.multiply(full_time,1000), state['T'], linestyle='solid', linewidth=3)
# print(bar1_data.iloc[:, :2])
plt.xlabel('Time [ms]')
plt.ylim([0.1,5000])
# plt.xlim([0.01,1350])
plt.ylabel(obs + ' Mole Fraction [ppm]')
plt.tick_params(axis = 'x', direction='in')
plt.tick_params(axis='x',labelsize=15)
plt.tick_params(axis = 'y', direction='in')
plt.tick_params(axis='y',labelsize=15)          
plt.savefig('flame_' + obs + '.pdf',dpi=1000, bbox_inches='tight')
plt.savefig('flame_' + obs + '.png',dpi=1000, bbox_inches='tight')
print("Successful simulation!")