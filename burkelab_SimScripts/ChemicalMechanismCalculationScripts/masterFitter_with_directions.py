"""
Class that allows for fitting of rate constants at various temperatures and pressures (k(T,P))
"""
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.polynomial.chebyshev import chebval
import scipy.optimize
from scipy.optimize import curve_fit
from io import StringIO
import matplotlib as mpl
import re

class masterFitter:
    def __init__(self, T_ls, P_ls, reactions, directions, fname,n_P=7, n_T=7, M_only=False):
        self.T_ls = T_ls
        self.P_ls = P_ls
        self.n_P=n_P
        self.n_T=n_T
        # self.collider = collider
        # self.reaction=reaction
        self.reactions=reactions
        self.directions=directions
        self.fname = fname
        self.P_min = P_ls[0]
        self.P_max = P_ls[-1]
        self.T_min = T_ls[0]
        self.T_max = T_ls[-1]
        # self.colours = colours
        self.M_only=M_only

    def check_if_unimolecular(self,reaction):
        reactants, _ = reaction.split('<=>')
        reactants = reactants.replace('(+M)', '').strip()
        if "H2O2" in reactants:
                return False
        reactant_list = [r.strip() for r in reactants.split('+')]
        num_reactants = 0
        for reactant in reactant_list:
            # Extract numeric prefixes (e.g., '2' in '2 NH2')
            match = re.match(r'(\d*)\s*(\S+)', reactant)
            if match:
                num_molecules = int(match.group(1)) if match.group(1) else 1
                num_reactants += num_molecules
        if num_reactants == 1:
            return True
        else:
            return False

    def get_YAML_kTP(self,reaction,direction,collider):
        gas = ct.Solution(self.fname)
        k_TP = []
        for T in self.T_ls:
            temp = []
            for P in self.P_ls:
                gas.TPX = T, P*ct.one_atm, {collider:1.0}
                if direction=="f":
                    temp.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
                elif direction=="r":
                    temp.append(gas.reverse_rate_constants[gas.reaction_equations().index(reaction)])
            k_TP.append(np.multiply(1000,temp)) # [cc/mol*s]
        return np.array(k_TP)
    
    def first_cheby_poly(self, x, n):
        '''Generate n-th order Chebyshev ploynominals of first kind.'''
        if n == 0: return 1
        elif n == 1: return x
        result = 2. * x * self.first_cheby_poly(x, 1) - self.first_cheby_poly(x, 0)
        m = 0
        while n - m > 2:
            result = 2. * x * result - self.first_cheby_poly(x, m+1)
            m += 1
        # print(result)
        return result
    def reduced_P(self,P):
        '''Calculate the reduced pressure.'''
        P_tilde = 2. * np.log10(P) - np.log10(self.P_min) - np.log10(self.P_max)
        P_tilde /= (np.log10(self.P_max) - np.log10(self.P_min))
        return P_tilde
    def reduced_T(self,T):
        '''Calculate the reduced temperature.'''
        T_tilde = 2. * T ** (-1) - self.T_min ** (-1) - self.T_max ** (-1)
        T_tilde /= (self.T_max ** (-1) - self.T_min ** (-1))
        return T_tilde
    def cheby_poly(self,reaction,direction,collider):
        '''Fit the Chebyshev polynominals to rate constants.
            Input rate constants vector k should be arranged based on pressure.'''
        k_TP = self.get_YAML_kTP(reaction,direction,collider)
        # cheb_mat = np.zeros((len(k_TP)*len(k_TP[0]), self.n_T * self.n_P))
        cheb_mat = np.zeros((len(k_TP.flatten()), self.n_T * self.n_P))
        # print(cheb_mat)
        # for n, P in enumerate(self.P_ls):       # !! assume that at each presssure, we have the same temperateure range
        #     for m, T in enumerate(self.T_ls):
        for m, T in enumerate(self.T_ls):
            for n, P in enumerate(self.P_ls):
                # for j in range(self.n_T):
                for i in range(self.n_T):
                    for j in range(self.n_P):
                    # for i in range(self.n_P):
                        T_tilde = self.reduced_T(T)
                        P_tilde = self.reduced_P(P)
                        T_cheb = self.first_cheby_poly(T_tilde, i)
                        P_cheb = self.first_cheby_poly(P_tilde, j)
                        # cheb_mat[n*len(self.T_ls)+m, i*self.n_P+j] = P_cheb * T_cheb
                        cheb_mat[m * len(self.P_ls) + n, i * self.n_P + j] = T_cheb * P_cheb
        coef = np.linalg.lstsq(cheb_mat, np.log10(k_TP.flatten()),rcond=None)[0].reshape((self.n_T, self.n_P))
        return coef
    def get_cheb_table(self,reaction,direction,collider):
        coef = self.cheby_poly(reaction,direction,collider)
        output2 = StringIO()
        output2.write(f"    temperature-range: [{self.T_min:.1f}, {self.T_max:.1f}]\n")
        output2.write(f"    pressure-range: [{self.P_min:.3e} atm, {self.P_max:.3e} atm]\n")
        output2.write(f"    data:  # {collider}\n")
        for i in range(len(coef)):
            string = ""
            for j in range(len(coef[0])):
                string+=f'{coef[i,j]:.4e}'
                if j<len(coef[0])-1:
                    string+=", "
            output2.write(f"    - [{string}]\n")
        return output2
    
    def get_PLOG_table(self,reaction,direction,collider):
        def arrhenius(T, A, n, Ea):
            return np.log(A) + n*np.log(T)+ (-Ea/(1.987*T))
        gas = ct.Solution(self.fname)
        output2 = StringIO()
        output2.write(f"    rate-constants:  # {collider}\n")
        for P in self.P_ls:
            k_T = []
            for T in self.T_ls:
                gas.TPX = T, P*ct.one_atm, {collider:1}
                if direction == "f":
                    temp = gas.forward_rate_constants[gas.reaction_equations().index(reaction)]
                if direction == "r":
                    temp = gas.reverse_rate_constants[gas.reaction_equations().index(reaction)]
                k_T.append(np.multiply(1000,temp))
            popt, pcov = curve_fit(arrhenius, self.T_ls, np.log(k_T),maxfev = 2000)
            output2.write(f"    - {{P: {P:.3e} atm, A: {popt[0]:.5e}, b: {popt[1]:.5e}, Ea: {popt[2]:.5e}}}\n")
        return output2
    
    def get_Troe_table(self,reaction,direction,collider): #PROBLEM: incorrect treatment of units for unimolecular vs bimolecular reactions
        def f(X,a0,n0,ea0,ai,ni,eai,Fcent):
            N=0.75-np.multiply(1.27,np.log10(Fcent))
            c=-0.4-np.multiply(0.67,np.log10(Fcent))
            d=0.14
            Rcal=1.987
            k0=np.multiply(np.multiply(a0,np.power(X[0,:],n0)),np.exp(np.divide(-ea0,np.multiply(Rcal,X[0,:]))))
            ki=np.multiply(np.multiply(ai,np.power(X[0,:],ni)),np.exp(np.divide(-eai,np.multiply(Rcal,X[0,:]))))
            Rjoule=8.3145
            M=np.divide(X[1,:],np.multiply(Rjoule,X[0,:]))
            M=np.divide(M,1000000.0)
            #print(M,k0,ki)
            logps=np.log10(k0)+np.log10(M)-np.log10(ki)
            den=logps+c
            den=np.divide(den,N-np.multiply(d,den))
            den=np.power(den,2)
            den=den+1.0
            logF=np.divide(np.log10(Fcent),den)
            logk_fit=np.log10(k0)+np.log10(M)+np.log10(ki)+logF-np.log10(ki+np.multiply(k0,M))
            return logk_fit
        if self.check_if_unimolecular(reaction)==False:
            multiplier=1000
        else:
            multiplier=1
            # print(reaction)
        Xvec=[]
        for i,val in enumerate(self.T_ls):
            for j,p in enumerate(self.P_ls):
                Xvec.append([val,p])
        Xvec=np.array(Xvec)
        Xvec=Xvec.T
        gas=ct.Solution(self.fname)
        k_list=np.zeros(np.shape(Xvec)[1])
        for i,temp in enumerate(Xvec[0]):
            gas.TPX=temp,Xvec[1,i],{collider:1.0}
            if direction=="f":
                k_list[i]=gas.forward_rate_constants[gas.reaction_equations().index(reaction)]
            if direction=="r":
                k_list[i]=gas.reverse_rate_constants[gas.reaction_equations().index(reaction)]
            k_list[i]=np.multiply(multiplier,k_list[i])  # [cc/mol*s]
        logk_list=np.log10(k_list)
        gas.TPX=500,0.001*ct.one_atm,{collider:1.0}
        M=np.divide(gas.P,np.multiply(8.3145,gas.T))
        M=np.divide(M,1000000.0)
        
        k0_g=np.divide(gas.forward_rate_constants[gas.reaction_equations().index(reaction)]*multiplier,M)
        gas.TPX=500,1000*ct.one_atm,{collider:1.0}
        ki_g=gas.forward_rate_constants[gas.reaction_equations().index(reaction)]*multiplier
        guess = [k0_g,0.0,0,ki_g,0.0,0,1]

        # bounds = (
        #         [1e-100, -np.inf, -np.inf, 1e-100, -np.inf, -np.inf, 1e-100],  # Lower bounds
        #         [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 1]         # Upper bounds
        #     )

        # popt, pcov = curve_fit(f,Xvec,logk_list,p0=guess,maxfev=100000, bounds=bounds)
        popt, pcov = curve_fit(f,Xvec,logk_list,p0=guess,maxfev=1000000)
        # popt, pcov = curve_fit(f,Xvec,logk_list,maxfev=100000)
        a0,n0,ea0,ai,ni,eai=popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]
        output2 = StringIO()
        output2.write(f"    low-P-rate-constant: {{A: {a0:.3e}, b: {n0:.3e}, Ea: {ea0:.3e}}}\n")
        output2.write(f"    high-P-rate-constant: {{A: {ai:.3e}, b: {ni:.3e}, Ea: {eai:.3e}}}\n")
        output2.write(f"    Troe: {{A: {popt[6]:.3e}, T3: 1.0e-30, T1: 1.0e+30}}\n")
        return output2
    
    def final_yaml(self,foutName,fit_fxn): # returns PLOG in LMRR YAML format
        output = StringIO()
        output.write(open('burkelab_SimScripts\\ChemicalMechanismCalculationScripts\\preamble.txt').read())
        output.write("\n")
        for reaction in self.reactions.keys():
            # print(self.reactions.keys())
            rxn=reaction.replace(" (+M)","")
            output.write(f"- equation: {rxn}\n")
            # output.write("  type: pressure-dependent-Arrhenius\n")
            output.write("  type: LMR_R\n")
            # if reaction == "2 NH2 (+M) <=> N2H4 (+M)" or reaction == "H + OH (+M) <=> H2O (+M)":
            #     output.write("  units: {length: m, quantity: kmol, activation-energy: cal/mol}\n")
            # output.write("  units: {length: m, quantity: kmol, activation-energy: cal/mol}\n")
            output.write("  collider-list:\n")
            # print(reaction)
            for i, collider in enumerate(self.reactions[reaction].keys()):
                # print(self.reactions[reaction].keys())
                # print(collider)
                if self.M_only == True:
                    if i == 0:
                        output.write(f"  - collider: 'M' # {collider} is reference collider\n")
                        output.write(f"    eps: "+self.reactions[reaction][collider]+"\n")
                        # output.write(self.get_PLOG_table(reaction,collider).getvalue())
                        keylist = list(self.reactions.keys())
                        direction = keylist.index(reaction)
                        output.write(fit_fxn(reaction,direction,collider).getvalue())
                    else:
                        output.write(f"  - collider: '{collider}'\n")
                        output.write(f"    eps: "+self.reactions[reaction][collider]+"\n")
                else:
                    if i == 0:
                        output.write(f"  - collider: 'M' # {collider} is reference collider\n")
                    else:
                        output.write(f"  - collider: '{collider}'\n")
                    output.write(f"    eps: "+self.reactions[reaction][collider]+"\n")
                    keylist = list(self.reactions.keys())
                    direction = keylist.index(reaction)
                    output.write(fit_fxn(reaction,direction,collider).getvalue())
            output.write("\n")
        fout = open(f"C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\{foutName}.yaml", "w")
        fout.write(output.getvalue())
        fout.close()
    

    def Troe(self,foutName): # returns PLOG in LMRR YAML format
        self.final_yaml(foutName,self.get_Troe_table)
    
    def PLOG(self,foutName): # returns PLOG in LMRR YAML format
        self.final_yaml(foutName,self.get_PLOG_table)

    def cheb2D(self,foutName): # returns Chebyshev in LMRR YAML format
        self.final_yaml(foutName,self.get_cheb_table)

###################################################
def makeplot(nom_liste,nom_fig):
    colours=[(np.random.rand(), np.random.rand(), np.random.rand()) for _ in range(30)]
    # PLOTTING ACROSS P AT FIXED T 
    titles=["Reaction_13","Reaction_16","Reaction_25","Reaction_72","Reaction_90","Reaction_167"]
    indices=[(0,0),(0,1),(1,0),(1,1),(2,0),(2,1)]
    f, ax = plt.subplots(3, 2, figsize=(10, 10)) 
    # titles=["Reaction_13","Reaction_16","Reaction_25","Reaction_90"]
    # indices=[(0,0),(0,1),(1,0),(1,1)]
    # f, ax = plt.subplots(2, 2, figsize=(10, 10)) 
    plt.subplots_adjust(wspace=0.2)
    mpl.rc('font',family='Times New Roman')
    mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl.rcParams['font.size'] = 5
    mpl.rcParams['xtick.labelsize'] = 5
    mpl.rcParams['ytick.labelsize'] = 5
    from matplotlib.legend_handler import HandlerTuple
    plt.rcParams['axes.labelsize'] = 5
    def get_kTP(fname,P_ls,T_ls,X,reaction,type,linestyle,marker,j,zorder,idx,reverse=False):
        gas = ct.Solution(fname)
        k_TP = []
        for i,P in enumerate(P_ls):
            # temp = []
            # for P in P_ls:
            gas.TPX = T_ls[0], P*ct.one_atm, {X:1.0}
            if reverse==True:
                k_TP.append(gas.reverse_rate_constants[gas.reaction_equations().index(reaction)])
            else:
                k_TP.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
        ax[idx].loglog(P_ls,k_TP, linestyle=linestyle,linewidth=1,markersize=1.5,markeredgewidth=0.6,marker=marker,fillstyle="none",label=f'{X}: {type}',color=colours[j],zorder=zorder)
    T_ls=[1000]
    P_ls=np.logspace(-4,6,num=60)
    nom_PLOG=nom_liste[0]
    nom_troe=nom_liste[1]
    nom_cheb=nom_liste[2]
    # nom_cheb="LMRtest_cheb_M"
    for j,reaction in enumerate(reactions.keys()):
        title=titles[j]
        colliders = reactions[reaction].keys()
        idx=indices[j]
        for j, X in enumerate(colliders):
            fname=f'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\{nom_PLOG}.yaml'
            get_kTP(fname,P_ls,T_ls,X,reaction.replace(" (+M)",""),"PLOG","-","none",j,1,idx)
            fname=f'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\{nom_troe}.yaml'
            get_kTP(fname,P_ls,T_ls,X,reaction.replace(" (+M)",""),"Troe",":","none",j,10,idx)
            fname=f'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\LMRtests\\{nom_cheb}.yaml'
            get_kTP(fname,P_ls,T_ls,X,reaction.replace(" (+M)",""),"Cheb","none","s",j,10,idx)
        ax[idx].set_title(f"{title}: {reaction} (T=1000K)")
        ax[idx].legend()
    plt.savefig('burkelab_SimScripts/rate_constant_plots/'+nom_fig, dpi=1000, bbox_inches='tight')
######################################################################3333


# Note: some of the original rxns already have a PLOG table, which have more limited pressure ranges than the range being explored here (a limitation)
# Note: not all of the Troe entries in yaml sandbox must have efficiencies specified for all colliders (a limitation)

# INPUTS
path='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\'
T_list=np.linspace(200,2000,30)
P_list=np.logspace(-6,12,num=30)

reactions={'H + OH (+M) <=> H2O (+M)':{"N2":"{A: 1, b: 0, Ea: 0}",
                             "AR":"{A: 2.20621e-02, b: 4.74036e-01, Ea: -1.13148e+02}",
                             "H2O":"{A: 1.04529e-01, b: 5.50787e-01, Ea: -2.32675e+02}"},
           'H + O2 (+M) <=> HO2 (+M)':{"AR":"{A: 1, b: 0, Ea: 0}",
                             "HE":"{A: 3.37601e-01, b: 1.82568e-01, Ea: 3.62408e+01}",
                             "N2":"{A: 1.24932e+02, b: -5.93263e-01, Ea: 5.40921e+02}",
                             "H2":"{A: 3.13717e+04, b: -1.25419e+00, Ea: 1.12924e+03}",
                             "CO2":"{A: 1.62413e+08, b: -2.27622e+00, Ea: 1.97023e+03}",
                             "NH3":"{A: 4.97750e+00, b: 1.64855e-01, Ea: -2.80351e+02}",
                             "H2O":"{A: 3.69146e+01, b: -7.12902e-02, Ea: 3.19087e+01}"},
           'H2O2 (+M) <=> 2 OH (+M)':{"AR":"{A: 1, b: 0, Ea: 0}",
                             "N2":"{A: 1.14813e+00, b: 4.60090e-02, Ea: -2.92413e+00}",
                             "CO2":"{A: 8.98839e+01, b: -4.27974e-01, Ea: 2.41392e+02}",
                             "H2O2":"{A: 6.45295e-01, b: 4.26266e-01, Ea: 4.28932e+01}",
                             "H2O":"{A: 1.36377e+00, b: 3.06592e-01, Ea: 2.10079e+02}"},
           'NH3 <=> H + NH2':{"AR":"{A: 1, b: 0, Ea: 0}",
                             "N2":"{A: 4.49281e+00, b: -9.46265e-02, Ea: -1.10071e+02}",
                             "O2":"{A: 1.15210e-01, b: 3.41234e-01, Ea: -3.89210e+02}",
                             "CO2":"{A: 9.19583e+00, b: 6.10696e-02, Ea: 9.01088e+01}",
                             "NH3":"{A: 1.49004e+01, b: 6.06535e-02, Ea: 2.47652e+02}",
                             "H2O":"{A: 1.14560e+01, b: 1.27501e-01, Ea: 3.13959e+02}"},
           '2 NH2 (+M) <=> N2H4 (+M)':{"AR":"{A: 1, b: 0, Ea: 0}",
                             "N2":"{A: 1.46848e+01, b: -3.00962e-01, Ea: 1.65146e+02}",
                             "O2":"{A: 1.53608e+00, b: -3.90857e-02, Ea: 4.43752e+00}",
                             "NH3":"{A: 1.64196e+01, b: -7.29636e-02, Ea: 3.66099e+02}",
                             "H2O":"{A: 2.17658e+01, b: -1.14715e-01, Ea: 4.19216e+02}"},
           'HNO <=> H + NO':{"AR":"{A: 1, b: 0, Ea: 0}",
                             "N2":"{A: 7.56771e+00, b: -1.96878e-01, Ea: 2.53168e+02}",
                             "H2O":"{A: 3.35451e+00, b: 1.33352e-01, Ea: 6.64633e+01}"}
}

# mF = masterFitter(T_list,P_list,reactions,path+"sandbox_substituted.yaml",n_P=7,n_T=7,M_only=True)
# mF.Troe("LMRtest_Troe_M_sub"); mF.PLOG("LMRtest_PLOG_M_sub"); mF.cheb2D("LMRtest_cheb_M_sub")
# makeplot(["LMRtest_PLOG_M_sub","LMRtest_Troe_M_sub","LMRtest_cheb_M_sub"],f'alzSub_Plog_Troe_Cheb_fixedT.png')

reactions={'H2O + M <=> H + OH + M':{"AR":"{A: 1, b: 0, Ea: 0}"},
           'H + O2 (+M) <=> HO2 (+M)':{"AR":"{A: 1, b: 0, Ea: 0}"},
           'H2O2 (+M) <=> 2 OH (+M)':{"AR":"{A: 1, b: 0, Ea: 0}"},
           'H + NH2 (+M) <=> NH3 (+M)':{"AR":"{A: 1, b: 0, Ea: 0}"},
           '2 NH2 (+M) <=> N2H4 (+M)':{"AR":"{A: 1, b: 0, Ea: 0}"},
           'H + NO (+M) <=> HNO (+M)':{"AR":"{A: 1, b: 0, Ea: 0}"}
}

directions=['r','f','f','r','f','r']

mF = masterFitter(T_list,P_list,reactions,directions,path+"sandbox_fullyoriginal.yaml",n_P=7,n_T=7,M_only=True)
mF.Troe("LMRtest_Troe_M_orig")
mF.PLOG("LMRtest_PLOG_M_orig")
mF.cheb2D("LMRtest_cheb_M_orig")
makeplot(["LMRtest_PLOG_M_orig","LMRtest_Troe_M_orig","LMRtest_cheb_M_orig"],f'alzOrig_Plog_Troe_Cheb_fixedT.png')

# python "burkelab_SimScripts\\ChemicalMechanismCalculationScripts\\masterFitter.py" > fitsPLOG.txt





