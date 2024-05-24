#REACTION 16: H+O2(+M)=HO2(=M)

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

#%%
#STEP 1: PERFORM A PLOG FIT FOR AR COLLIDER, USING GRAPH-READ DATA FROM FIGURE 5 OF SJK PAPER
def arrhenius(T, A, n, Ea):
    return np.log(A) + n*np.log(T)+ (-Ea/(1.987*T))
pltcolours = ['k','g','b','r','olive','orange','brown','indigo']
path="G:\\Mon disque\\Columbia\\Burke Lab\\01 Mixture Rules Project\\"
fname=path+"rateConstantCalculations_SKJFig5.csv"
plt.figure(figsize=(8,6))
p_list = [10,100,300,760,2280,7600,22800,76000]
# p_list = [76000]
dataset = pd.read_csv(fname)
for i,p in enumerate(p_list):
    dataset_p = dataset[dataset['P(torr)']==p]
    p_data = dataset_p['P(atm)'] # atm
    T_data=dataset_p['T'] # K
    N_A = 6.022e23 # molec/mol
    lnk_M_raw = dataset_p['lnK/M (cm6/molec^2/sec)']
    k_M_data=np.multiply(np.exp(lnk_M_raw),np.square(N_A)) # ln(cm6/molec^2/s) converted to cm6/mol^2/s
    R=82.1 #cm3atm/K/mol
    conc = np.divide(np.divide(p,R),T_data) #mol/cm3
    # e^(K/M)*Na^2*conc
    k_data = np.multiply(k_M_data,conc) #cm3/mol/sec
    popt, pcov = curve_fit(arrhenius, T_data, np.log(k_data),maxfev = 2000)
    print(("- {P: %.3e, A: %.5e, b: %.5e, Ea: %.5e}")%(p/760, popt[0],popt[1],popt[2]))
    lnk_fit = arrhenius(T_data,popt[0],popt[1],popt[2])
    lnk_M_fit = np.log(np.divide(np.divide(np.exp(lnk_fit),conc),np.square(N_A))) # ln(cm6/molec^2/s)
    plt.plot(T_data,lnk_M_fit,label=str(p) + ' fit',linestyle='solid',color=pltcolours[i])
    plt.scatter(T_data,lnk_M_raw,marker='x',color=pltcolours[i],label=str(p) + ' data')
plt.title('PLOG Fit for H + O2 (+Ar) <-> HO2 (+Ar) [Data from Graph-Reading]')
plt.xlabel('Temperature [K]')
plt.ylabel('lnk/[M] [cm^6/molecule^2*s] ')
plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.show()
# %%
#STEP 2: MANUALLY PASTE THE PARAMETERS INTO "sandbox.yaml"

#STEP 3: CREATE A PLOT COMPARING THE ARRENHIUS FUNCTION OF THE FITTED PARAMETERS TO THE ORIGINAL GRAPH-READ DATA
fig = plt.figure(figsize=(8, 8))
spec = gridspec.GridSpec(ncols=1, nrows=2,hspace=0, height_ratios=[3, 1])
fig.suptitle(t='PLOG Fit Tester', x=0.5, y=0.95)
fig.supxlabel(t='Temperature [K]', x=0.5, y=0.05)
ax0 = fig.add_subplot(spec[0])
ax1 = fig.add_subplot(spec[1])
ax0.set_ylabel('ln ( k/[M] [cm^6/molecule^2*s] )')
ax1.set_ylabel('Percent Change [%]')
p_list = [10,100,300,760,2280,7600,22800,76000]
T_list = np.linspace(200,2000,50)
reaction_plog = 'H + O2 <=> HO2'
for i, p in enumerate(p_list):
    ratelist_plog = []
    for j, T in enumerate(T_list):
        gas = ct.Solution("C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\sandbox.yaml")
        # print(p)
        gas.TPX = T, p/760, 'H2:0.2,O2:0.2,Ar:0.8'
        rc_plog = gas.forward_rate_constants[gas.reaction_equations().index(reaction_plog)] # m3/kmol/s
        rc_plog = rc_plog*1000 # cm3/mol*s
        rc_plog_M_molec = np.divide(rc_plog, np.square(N_A)) # cm6/molec2*s
        ratelist_plog.append(np.log(rc_plog_M_molec))
    ax0.plot(T_list,ratelist_plog,label=str(p) + ' plog',linestyle='solid',color=pltcolours[i])
    dataset = pd.read_csv(fname) # data from klippy
    dataset_p = dataset[dataset['P(torr)']==p]
    T_data=dataset_p['T']
    k_data_og=dataset_p['lnK/M (cm6/molec^2/sec)']  # cm3/mol/s
    ax0.scatter(T_data,k_data_og,marker='x',color=pltcolours[i],label=str(p) + ' data')
    ax1.plot(T_data, np.multiply(np.divide(np.subtract(np.interp(T_data,T_list,ratelist_plog),k_data_og),k_data_og),100),label=str(p) + ' plog',linestyle='solid',color=pltcolours[i])
ax0.tick_params(direction='in')
ax1.tick_params(direction='in')
ax0.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
plt.show()

# %%
#COMPUTE ARRHENIUS PARAMETERS FOR COLLIDERS USING COLLISION EFFICIENCIES FROM JASPER
from scipy.optimize import least_squares
def plot_ratefit(temperatures, rate_constants, pltcolour, labell): #3-parameter fit (A, b, Ea)
    def arrhenius_rate(T, A, beta, Ea):
        R = 8.314  # Gas constant in J/(mol K)
        return A * T**beta * np.exp(-Ea / (R * T))
    def fit_function(params, T, ln_rate_constants):
        A, beta, Ea = params
        return np.log(arrhenius_rate(T, A, beta, Ea)) - ln_rate_constants
    initial_guess = [3, 0.5, 50.0]  
    result = least_squares(fit_function, initial_guess, args=(temperatures, np.log(rate_constants)))
    A_fit, beta_fit, Ea_fit = result.x
    T_range = np.linspace(200, 2000, 100)
    fit_curve = np.exp(np.log(arrhenius_rate(T_range, A_fit, beta_fit, Ea_fit)))
    print(('- collider: "%s"\n  low-P-rate-constant: {A: %.5e, b: %.5e, Ea: %.5e}')%(labell.upper(),A_fit, beta_fit, Ea_fit))
    plt.plot(T_range, np.log(fit_curve), label=labell, color=pltcolour)
    plt.scatter(temperatures, np.log(rate_constants), label=None, color=pltcolour)

pltcolours = ['k', 'g', 'b', 'r', 'olive', 'orange', 'brown', 'indigo']
a,b,c=1,1,1 #weights of each data point
plt.figure()
plot_ratefit(np.array([300]*a + [1000]*b + [2000]*c), np.array([0.90]*a + [1.17]*b + [1.34]*c), pltcolours[0], "He")
plot_ratefit(np.array([300]*a + [1000]*b + [2000]*c), np.array([1.71]*a+ [1.58]*b + [1.20]*c), pltcolours[1], "N2")
plot_ratefit(np.array([300]*a + [1000]*b + [2000]*c), np.array([3.69]*a+ [3.07]*b + [1.71]*c), pltcolours[2], "H2")
plot_ratefit(np.array([300]*a + [1000]*b + [2000]*c), np.array([13.7]*a+ [8.94]*b + [3.03]*c), pltcolours[3], "CO2")
plot_ratefit(np.array([300]*a + [1000]*b + [2000]*c), np.array([20.4]*a+ [17.9]*b + [18.7]*c), pltcolours[4], "NH3")
plot_ratefit(np.array([300]*a + [1000]*b + [2000]*c), np.array([23.3]*a+ [22.2]*b + [21.3]*c), pltcolours[5], "H2O")
plt.xlabel('Temperature (K)')
plt.ylabel('logK')
plt.legend()
plt.title("HO2 Reaction: 3-Parameter Fit")
plt.show()
# %%
