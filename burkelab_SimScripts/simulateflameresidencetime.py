import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
mpl.rc('font',family='Times New Roman')
plt.figure()
save_fig = True
nominal_model = 'G:\\Mon disque\\Columbia\\Burke Lab\\07 Mechanisms\\Ammonia\\Mei-2019\\mei-2019.yaml'
T_fuel = 300
T_air = 650
phi = 1.22
# P_list = [1,10,20] # bar
# trimStart = [76,99,112] # number of indices to cut off from the start of the flame simulation
# trimEnd = [1,8,22] # number of indices to cut off from the end of the flame simulation
P_list = [1] # bar
trimStart = [76] # number of indices to cut off from the start of the flame simulation
trimEnd = [1] # number of indices to cut off from the end of the flame simulation

width = 2
colours=["k","xkcd:grey","xkcd:teal"]
for i,P in enumerate(P_list):
    print(f"Simulating {P} bar...")
    def cp(T,P,X):
        gas_stream = ct.Solution(nominal_model)
        gas_stream.TPX = T, P*1e5, {X:1}
        return gas_stream.cp_mole # [J/kmol/K]
    cp_fuel = cp(T_fuel,P,'NH3') # [J/kmol/K]
    cp_o2 = cp(T_air,P,'O2') # [J/kmol/K]
    cp_n2 = cp(T_air,P,'N2') # [J/kmol/K]
    x_fuel = (phi*(1/0.75)*0.21)/(1+phi*(1/0.75)*0.21)
    x_o2 = 0.21*(1-x_fuel)
    x_n2 = 0.79*(1-x_fuel)
    x_air=1-x_fuel
    T_mix = (x_fuel*cp_fuel*T_fuel+(x_o2*cp_o2+x_n2*cp_n2)*T_air)/(x_fuel*cp_fuel+ x_o2*cp_o2 + x_n2*cp_n2)
    mix1 = ct.Solution(nominal_model)
    mix1.TPX = T_mix, P*1e5,{'NH3':x_fuel,'O2':x_o2,'N2':x_n2}
    mix1_eq = ct.Solution(nominal_model)
    mix1_eq.TPX = T_mix, P,{'NH3':x_fuel,'O2':x_o2,'N2':x_n2}
    mix1_eq.equilibrate("HP")
    mixtures = [mix1]
    for j,mix in enumerate(mixtures):
        flame = ct.FreeFlame(mix,width=width)
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
        conc = np.multiply(state['X'][mix.species_names.index("NH2")], np.divide(P*1e5,np.multiply(Rjoule,state['T'])))
        # conc = state['X'][mix.species_names.index("NH2")]
        distance_index=np.argmax(conc)
        tau_f = np.array(state['d'][distance_index])/S_b
        tau_list = [(ft - tau_f)*1000 for ft in full_time] # [ms]
        X_NO = state['X'][mix.species_names.index("NO")]
        X_O2 = state['X'][mix.species_names.index("O2")]
        X_O2dry = 0.15 # 15% O2 dry
        X_NO_15O2dry = np.multiply(np.multiply(X_NO,np.divide(0.21-X_O2dry,np.subtract(0.21,X_O2))),1e6) # [ppm, 15% O2 dry]
        # distance_f = np.array(state['d'][distance_index])
        # distance_list_new = [dl - distance_f for dl in state['d']]
        plt.loglog(np.multiply(tau_list[trimStart[i]:trimEnd[i]*(-1)],P), X_NO_15O2dry[trimStart[i]:trimEnd[i]*(-1)], linestyle='solid', linewidth=1.5,color=colours[i],label=f"{P}bar")
        print("Added to plot!")
print("Adding graph-read data...")
path="G:\\Mon disque\\Columbia\\Burke Lab\\09 NOx Mini-Project\\Graph Reading\\"
expData = [{'f': '1bar','mkr':'o','colour':'k','line':'none'},
        #    {'f': '1bar_eq','mkr':'none','colour':'k','line':':'},
           {'f': '10bar','mkr':'x','colour':"xkcd:grey",'line':'none'},
        #    {'f': '10bar_eq','mkr':'none','colour':"xkcd:grey",'line':':'},
           {'f': '20bar','mkr':'s','colour':"xkcd:teal",'line':'none'},
        #    {'f': '20bar_eq','mkr':'none','colour':"xkcd:teal",'line':':'}
        ]
for exp in expData:
    dat = pd.read_csv(path+exp['f']+'.csv',header=None)
    plt.loglog(dat.iloc[:, 0],dat.iloc[:, 1],marker=exp['mkr'],fillstyle='none',linestyle=exp['line'],color=exp['colour'],markersize=3,markeredgewidth=1, label=f"{exp['f']}", zorder=110)
plt.xlabel(r'P*$\tau$ [bar*ms]')
plt.ylim([5,1350])
# plt.xlim([0.01,1350])
plt.ylabel('NO [ppm, 15% O$_2$ dry]')
plt.legend()
plt.tick_params(axis = 'x', direction='in')
plt.tick_params(axis='x',labelsize=15)
plt.tick_params(axis = 'y', direction='in')
plt.tick_params(axis='y',labelsize=15)          
plt.savefig('flame_NO.pdf',dpi=1000, bbox_inches='tight')
plt.savefig('flame_NO.png',dpi=1000, bbox_inches='tight')
print("Simulation complete!")