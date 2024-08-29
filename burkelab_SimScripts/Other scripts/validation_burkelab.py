
#Validation script
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# bklabct.print_stack_trace_on_segfault()

# from burkelab_SimScripts.ChemicalMechanismCalculationScripts.masterFitter import *

def plot_fixed_T():

    # f, ax = plt.subplots(1, 3, figsize=(4, 3)) 
    # plt.subplots_adjust(wspace=0.2)
    f, ax = plt.subplots(3, 1, figsize=(3, 5)) 
    plt.subplots_adjust(hspace=0.2)
    mpl.rc('font',family='Times New Roman')
    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams['font.size'] = 5
    mpl.rcParams['xtick.labelsize'] = 3
    mpl.rcParams['ytick.labelsize'] = 3
    from matplotlib.legend_handler import HandlerTuple
    plt.rcParams['axes.labelsize'] = 5
    # mpl.rcParams['xtick.major.width'] = 0.5  # Width of major ticks on x-axis
    # mpl.rcParams['ytick.major.width'] = 0.5  # Width of major ticks on y-axis
    # mpl.rcParams['xtick.minor.width'] = 0.5  # Width of minor ticks on x-axis
    # mpl.rcParams['ytick.minor.width'] = 0.5  # Width of minor ticks on y-axis
    # mpl.rcParams['xtick.major.size'] = 2.5  # Length of major ticks on x-axis
    # mpl.rcParams['ytick.major.size'] = 2.5  # Length of major ticks on y-axis
    # mpl.rcParams['xtick.minor.size'] = 1.5  # Length of minor ticks on x-axis
    # mpl.rcParams['ytick.minor.size'] = 1.5  # Length of minor ticks on y-axis

    # PLOG PLOT
    T_ls=[1000]
    P_ls=np.logspace(-5,11,num=50)
    # reaction='H2O2 (+M) <=> 2 OH (+M)'
    reaction='H + O2 (+M) <=> HO2 (+M)'
    # colliders = ["Ar","H2O","CO2", "N2", "H2O2"]
    colliders=["Ar","HE","N2", "H2", "CO2", "NH3", "H2O"]
    colours=[(np.random.rand(), np.random.rand(), np.random.rand()) for _ in range(len(colliders))]
    fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_M_PLOG.yaml'
    gas = ct.Solution(fname)
    for j, X in enumerate(colliders):
        k_TP = []
        for i,P in enumerate(P_ls):
            # temp = []
            # for P in P_ls:
            gas.TPX = T_ls[0], P*ct.one_atm, {X:1.0}
            k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
        if j == 0:
            ax[0].loglog(P_ls,k_TP, linestyle="-",label="M-only",color=colours[j])
        else:
            ax[0].loglog(P_ls,k_TP, linestyle="-",label=None,color=colours[j])

    # fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_purePLOG.yaml'
    # gas = ct.Solution(fname)
    # for j, X in enumerate(colliders):
    #     k_TP = []
    #     for i,P in enumerate(P_ls):
    #         # temp = []
    #         # for P in P_ls:
    #         gas.TPX = T_ls[0], P*ct.one_atm, {X:1.0}
    #         k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
    #     if j == 0:
    #         ax[0].loglog(P_ls,k_TP, linestyle=":",label="pure",color=colours[j])
    #     else:
    #         ax[0].loglog(P_ls,k_TP, linestyle=":",label=None,color=colours[j])
    ax[0].legend()
    # ax[0].set_title("Reaction 25: PLOG (T=1000)")
    ax[0].set_title("Reaction 16: PLOG (T=1000)")
    # plt.savefig('burkelab_SimScripts/rate_constant_plots/'+'Rxn25_PLOG.png', dpi=1000, bbox_inches='tight')


    # TROE PLOT
    # reaction='H2O2 (+M) <=> 2 OH (+M)'
    reaction='H + O2 (+M) <=> HO2 (+M)'
    # colliders = ["Ar","H2O","CO2", "N2", "H2O2"]
    colliders=["Ar","HE","N2", "H2", "CO2", "NH3", "H2O"]
    colours=[(np.random.rand(), np.random.rand(), np.random.rand()) for _ in range(len(colliders))]
    fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_M_Troe.yaml'
    gas = ct.Solution(fname)
    for j, X in enumerate(colliders):
        k_TP = []
        for i,P in enumerate(P_ls):
            # temp = []
            # for P in P_ls:
            gas.TPX = T_ls[0], P*ct.one_atm, {X:1.0}
            k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
        if j == 0:
            ax[1].loglog(P_ls,k_TP, linestyle="-",label="M-only",color=colours[j])
        else:
            ax[1].loglog(P_ls,k_TP, linestyle="-",label=None,color=colours[j])

    # fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_pureTroe.yaml'
    # gas = ct.Solution(fname)
    # for j, X in enumerate(colliders):
    #     k_TP = []
    #     for i,P in enumerate(P_ls):
    #         # temp = []
    #         # for P in P_ls:
    #         gas.TPX = T_ls[0], P*ct.one_atm, {X:1.0}
    #         k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
    #     if j == 0:
    #         ax[1].loglog(P_ls,k_TP, linestyle=":",label="pure",color=colours[j])
    #     else:
    #         ax[1].loglog(P_ls,k_TP, linestyle=":",label=None,color=colours[j])
    ax[1].legend()
    ax[1].set_xlabel("Pressure")
    ax[1].set_title("Reaction 16: Troe (T=1000)")
    # ax[0,1].plot.savefig('burkelab_SimScripts/rate_constant_plots/'+'Rxn25_Troe.png', dpi=1000, bbox_inches='tight')


    # CHEBYSHEV PLOT
    # reaction='H2O2 (+M) <=> 2 OH (+M)'
    reaction='H + O2 (+M) <=> HO2 (+M)'
    # reaction='H2O2 <=> 2 OH'
    # colliders = ["Ar","H2O","CO2", "N2", "H2O2"]
    colliders=["Ar","HE","N2", "H2", "CO2", "NH3", "H2O"]
    colours=[(np.random.rand(), np.random.rand(), np.random.rand()) for _ in range(len(colliders))]
    fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_M_Chebyshev.yaml'
    gas = ct.Solution(fname)
    for j, X in enumerate(colliders):
        k_TP = []
        for i,P in enumerate(P_ls):
            # temp = []
            # for P in P_ls:
            gas.TPX = T_ls[0], P*ct.one_atm, {X:1.0}
            k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
        if j == 0:
            ax[2].loglog(P_ls,k_TP, linestyle="-",label="M-only",color=colours[j])
        else:
            ax[2].loglog(P_ls,k_TP, linestyle="-",label=None,color=colours[j])

    # fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_pureChebyshev.yaml'
    # gas = ct.Solution(fname)
    # for j, X in enumerate(colliders):
    #     k_TP = []
    #     for i,P in enumerate(P_ls):
    #         # temp = []
    #         # for P in P_ls:
    #         gas.TPX = T_ls[0], P*ct.one_atm, {X:1.0}
    #         k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
    #     if j == 0:
    #         ax[2].loglog(P_ls,k_TP, linestyle=":",label="pure",color=colours[j])
    #     else:
    #         ax[2].loglog(P_ls,k_TP, linestyle=":",label=None,color=colours[j])
    ax[2].legend()
    ax[2].set_title("Reaction 16: Chebyshev (T=1000)")
    plt.savefig('burkelab_SimScripts/rate_constant_plots/'+'Rxn25_Plog_Troe_Cheb_fixedT.png', dpi=1000, bbox_inches='tight')


def plot_fixed_P():

    f, ax = plt.subplots(1, 3, figsize=(7, 3)) 
    plt.subplots_adjust(wspace=0.2)
    mpl.rc('font',family='Times New Roman')
    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams['font.size'] = 5
    mpl.rcParams['xtick.labelsize'] = 3
    mpl.rcParams['ytick.labelsize'] = 3
    from matplotlib.legend_handler import HandlerTuple
    plt.rcParams['axes.labelsize'] = 5
    # mpl.rcParams['xtick.major.width'] = 0.5  # Width of major ticks on x-axis
    # mpl.rcParams['ytick.major.width'] = 0.5  # Width of major ticks on y-axis
    # mpl.rcParams['xtick.minor.width'] = 0.5  # Width of minor ticks on x-axis
    # mpl.rcParams['ytick.minor.width'] = 0.5  # Width of minor ticks on y-axis
    # mpl.rcParams['xtick.major.size'] = 2.5  # Length of major ticks on x-axis
    # mpl.rcParams['ytick.major.size'] = 2.5  # Length of major ticks on y-axis
    # mpl.rcParams['xtick.minor.size'] = 1.5  # Length of minor ticks on x-axis
    # mpl.rcParams['ytick.minor.size'] = 1.5  # Length of minor ticks on y-axis

    # PLOG PLOT
    T_ls=np.linspace(200,2000)
    P_ls=[100]
    P_ls=np.logspace(-5,11,num=50)
    # reaction='H2O2 (+M) <=> 2 OH (+M)'
    reaction='H + O2 (+M) <=> HO2 (+M)'
    # colliders = ["Ar","H2O","CO2", "N2", "H2O2"]
    colliders=["Ar","HE","N2", "H2", "CO2", "NH3", "H2O"]
    colours=[(np.random.rand(), np.random.rand(), np.random.rand()) for _ in range(len(colliders))]
    fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_M_PLOG.yaml'
    gas = ct.Solution(fname)
    for j, X in enumerate(colliders):
        k_TP = []
        for i,T in enumerate(T_ls):
            # temp = []
            # for P in P_ls:
            gas.TPX = T, P_ls[0]*ct.one_atm, {X:1.0}
            k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
        if j == 0:
            ax[0].semilogy(T_ls,k_TP, linestyle="-",label="M-only",color=colours[j])
        else:
            ax[0].semilogy(T_ls,k_TP, linestyle="-",label=None,color=colours[j])

    # fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_purePLOG.yaml'
    # gas = ct.Solution(fname)
    # for j, X in enumerate(colliders):
    #     k_TP = []
    #     for i,T in enumerate(T_ls):
    #         # temp = []
    #         # for P in P_ls:
    #         gas.TPX = T, P_ls[0]*ct.one_atm, {X:1.0}
    #         k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
    #     if j == 0:
    #         ax[0].semilogy(T_ls,k_TP, linestyle=":",label="all i",color=colours[j])
    #     else:
    #         ax[0].semilogy(T_ls,k_TP, linestyle=":",label=None,color=colours[j])
    ax[0].legend()
    ax[0].set_title("Reaction 16: PLOG (P=100)")
    # plt.savefig('burkelab_SimScripts/rate_constant_plots/'+'Rxn25_PLOG.png', dpi=1000, bbox_inches='tight')


    # TROE PLOT
    # reaction='H2O2 (+M) <=> 2 OH (+M)'
    reaction='H + O2 (+M) <=> HO2 (+M)'
    # colliders = ["Ar","H2O","CO2", "N2", "H2O2"]
    colliders=["Ar","HE","N2", "H2", "CO2", "NH3", "H2O"]
    colours=[(np.random.rand(), np.random.rand(), np.random.rand()) for _ in range(len(colliders))]
    fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_M_Troe.yaml'
    gas = ct.Solution(fname)
    for j, X in enumerate(colliders):
        k_TP = []
        for i,T in enumerate(T_ls):
            # temp = []
            # for P in P_ls:
            gas.TPX = T, P_ls[0]*ct.one_atm, {X:1.0}
            k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
        if j == 0:
            ax[1].semilogy(T_ls,k_TP, linestyle="-",label="M-only",color=colours[j])
        else:
            ax[1].semilogy(T_ls,k_TP, linestyle="-",label=None,color=colours[j])

    # fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_pureTroe.yaml'
    # gas = ct.Solution(fname)
    # for j, X in enumerate(colliders):
    #     k_TP = []
    #     for i,T in enumerate(T_ls):
    #         # temp = []
    #         # for P in P_ls:
    #         gas.TPX = T, P_ls[0]*ct.one_atm, {X:1.0}
    #         k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
    #     if j == 0:
    #         ax[1].semilogy(T_ls,k_TP, linestyle=":",label="all i",color=colours[j])
    #     else:
    #         ax[1].semilogy(T_ls,k_TP, linestyle=":",label=None,color=colours[j])
    ax[1].legend()
    ax[1].set_xlabel("Temperature")
    ax[1].set_title("Reaction 16: Troe (P=100)")
    # ax[0,1].plot.savefig('burkelab_SimScripts/rate_constant_plots/'+'Rxn25_Troe.png', dpi=1000, bbox_inches='tight')


    # CHEBYSHEV PLOT
    # reaction='H2O2 (+M) <=> 2 OH (+M)'
    reaction='H + O2 (+M) <=> HO2 (+M)'
    # reaction='H2O2 <=> 2 OH'
    # colliders = ["Ar","H2O","CO2", "N2", "H2O2"]
    colliders=["Ar","HE","N2", "H2", "CO2", "NH3", "H2O"]
    colours=[(np.random.rand(), np.random.rand(), np.random.rand()) for _ in range(len(colliders))]
    fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_M_Chebyshev.yaml'
    gas = ct.Solution(fname)
    for j, X in enumerate(colliders):
        k_TP = []
        for i,T in enumerate(T_ls):
            # temp = []
            # for P in P_ls:
            gas.TPX = T, P_ls[0]*ct.one_atm, {X:1.0}
            k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
        if j == 0:
            ax[2].semilogy(T_ls,k_TP, linestyle="-",label="M-only",color=colours[j])
        else:
            ax[2].semilogy(T_ls,k_TP, linestyle="-",label=None,color=colours[j])

    # fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_pureChebyshev.yaml'
    # gas = ct.Solution(fname)
    # for j, X in enumerate(colliders):
    #     k_TP = []
    #     for i,T in enumerate(T_ls):
    #         # temp = []
    #         # for P in P_ls:
    #         gas.TPX = T, P_ls[0]*ct.one_atm, {X:1.0}
    #         k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
    #     if j == 0:
    #         ax[2].semilogy(T_ls,k_TP, linestyle=":",label="all i",color=colours[j])
    #     else:
    #         ax[2].semilogy(T_ls,k_TP, linestyle=":",label=None,color=colours[j])
    ax[2].legend()
    ax[2].set_title("Reaction 16: Chebyshev (P=100)")
    plt.savefig('burkelab_SimScripts/rate_constant_plots/'+'Rxn25_Plog_Troe_Cheb_fixedP.png', dpi=1000, bbox_inches='tight')





# plot_fixed_P()
plot_fixed_T()

# python 'burkelab_SimScripts\validation_burkelab.py'


















# file = 'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\LMRtest_purePLOG.yaml'
# # reactions = ['H + OH (+M) <=> H2O (+M)']
# reactions = ['H + O2 (+M) <=> HO2 (+M)']
# # reactions=[1]
# gas = bklabct.Solution(file)
# #Temp = np.linspace(750,2500,50)
# Temp=[1000]
# Pres = np.logspace(-12,8,100)
# # Pres = [bklabct.one_atm*10**(-6)] # units: Pa
# # Pres=[1e-6,1e-5,1e-4,0.001,0.01,0.1, 1, 10, 100, 1000, 10000, 100000, 1e6, 1e7, 1e8]
# ARlist =[0.9,0.8,0.7,0.6,0.5]
# H2Olist =[0.1,0.2,0.3,0.4,0.5]

# # for x in range(len(ARlist)):

# plt.figure()
# plt.title('H + O2 (+M) <=> HO2 (+M)')
# # labels=["LMRR-PLOG","LMRR-Troe"]
# for i, R in enumerate(reactions):
#     k_list=[]
#     for j, P in enumerate(Pres):
#         temp_list = []
#         # for k,T in enumerate(Temp):
#             # gas.TPX = T,P,{'AR':ARlist[x],'H2O':H2Olist[x]}
#             # gas.TPX = T,P,{'O2':ARlist[x],'H2O':H2Olist[x]}
#             # gas.TPX = T,P,{'O2':1.0}
#         gas.TPX = Temp[0],P*101325,{'AR':1.0}
#             # print(gas.TPX)
#         rc = gas.forward_rate_constants[gas.reaction_equations().index(R)]
#         k_list.append(rc)  
#     # print("%.6e" % k_list[j][0])
#     plt.loglog(Pres,k_list, label="LMRR-PLOG")
    
#     # for j,P in enumerate(Pres):
#     #     plt.plot(Temp,k_list[j],label=str(P)+' atm') 
#     #     #plt.semilogy(Temp,k_list[j],label=str(P)+' atm')    

# # file = 'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\kineticsfromscratch_LMRtest_Troe.yaml'
# # # reactions = ['H + O2 (+M) <=> HO2 (+M)']
# # gas = bklabct.Solution(file)
# # Temp=[1000]
# # Pres = np.logspace(-12,8,50)
# # for i, R in enumerate(reactions):
# #     k_list=[]
# #     for j, P in enumerate(Pres):
# #         temp_list = []
# #         # for k,T in enumerate(Temp):
# #             # gas.TPX = T,P,{'AR':ARlist[x],'H2O':H2Olist[x]}
# #             # gas.TPX = T,P,{'O2':ARlist[x],'H2O':H2Olist[x]}
# #             # gas.TPX = T,P,{'O2':1.0}
# #         gas.TPX = Temp[0],P*101325,{'AR':1.0}
# #             # print(gas.TPX)
# #         rc = gas.forward_rate_constants[gas.reaction_equations().index(R)]
# #         k_list.append(rc)  
# #     # print("%.6e" % k_list[j][0])
# #     plt.loglog(Pres,k_list, label="LMRR-Troe",linestyle="none",marker="o",markersize=3.5,fillstyle='none',color='k')

# # plt.legend()
# # plt.xlabel('Pressure [atm]')
# # plt.ylabel('k')

# # # plt.savefig("reaction1",bbox_inches="tight")
# # plt.show()
    
# # # 1.0810363605840141e-13, 2.1620727211680282e-13, 3.243109081752042e-13, 4.324145442336056e-13,

# # # print()

# # # # %%
# # # T=1000 #K
# # # P=101325 #Pa
# # # def k (dict, T):
# # #     if 'eig0' in dict:
# # #         A = dict['eig0']['A']
# # #         b = dict['eig0']['b']
# # #         Ea = dict['eig0']['Ea']
# # #     if 'plog' in dict:
# # #         A = dict['plog']['A']
# # #         b = dict['plog']['b']
# # #         Ea = dict['plog']['Ea']
# # #     R = 1.987 # cal/molK
# # #     return A*T**(b)*np.exp(-Ea/R/T)
# # # AR = {'X':0.5, 'eig0':{'A': 2.20621e-02, 'b': 4.74036e-01, 'Ea': -1.13148e+02}}
# # # H2O = {'X':0.5, 'eig0':{'A': 1.04529e-01, 'b': 5.50787e-01, 'Ea': -2.32675e+02}}

# # # eig0_mix = AR['X'] * k(AR,T) + H2O['X'] * k(H2O,T)

# # # def Peff(eig0,eig0_mix,P,species):
# # #     Peff = np.exp(np.log(P)+np.log(eig0_mix)-np.log(eig0))/101325
# # #     label = "Peff_"+species+" [atm] = "
# # #     print(label, Peff)
# # #     return Peff
# # # Peff_AR = Peff(k(AR,T),eig0_mix,P,"AR")
# # # Peff_H2O = Peff(k(H2O,T),eig0_mix,P,"H2O")
# # # #%%

# # # def interpolate(p, p1, p2, y1, y2):
# # #     """Interpolates to find y at a specific p using points (p1, y1) and (p2, y2)."""
# # #     return y1 + (y2 - y1) * ((p - p1) / (p2 - p1))

# # # # Given data points
# # # data = {
# # #     1.000e-04: {"A": 5.30514e+12, "b": -2.80725, "Ea": 499.267},
# # #     1.000e-03: {"A": 5.25581e+13, "b": -2.80630, "Ea": 499.946}
# # # }

# # # # Target pressure
# # # p_target = 10**(-3.5)

# # # # Interpolate A, b, Ea
# # # A_interpolated = interpolate(p_target, 1.000e-04, 1.000e-03, data[1.000e-04]["A"], data[1.000e-03]["A"])
# # # b_interpolated = interpolate(p_target, 1.000e-04, 1.000e-03, data[1.000e-04]["b"], data[1.000e-03]["b"])
# # # Ea_interpolated = interpolate(p_target, 1.000e-04, 1.000e-03, data[1.000e-04]["Ea"], data[1.000e-03]["Ea"])


# # # AR = {'X':0.5, 'plog':{'A': 2.20621e-02, 'b': 4.74036e-01, 'Ea': -1.13148e+02}}
# # # H2O = {'X':0.5, 'plog':{'A': 1.04529e-01, 'b': 5.50787e-01, 'Ea': -2.32675e+02}}
