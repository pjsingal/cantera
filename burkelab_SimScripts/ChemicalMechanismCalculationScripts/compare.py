import sys, os
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import sys, os
import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class chebFit_2D:
    def __init__(self, T_ls, P_ls, n_P, n_T, X, reaction, fname, colours):
        self.T_ls = T_ls
        self.P_ls = P_ls
        self.n_P = n_P
        self.n_T = n_T
        self.X = X
        self.reaction = reaction
        self.fname = fname
        self.P_min = P_ls[0]
        self.P_max = P_ls[-1]
        self.T_min = T_ls[0]
        self.T_max = T_ls[-1]
        self.colours = colours

    def get_YAML_kTP(self):
        gas = ct.Solution(self.fname)
        k_TP = []
        for T in self.T_ls:
            temp = []
            for P in self.P_ls:
                gas.TPX = T, P * ct.one_atm, self.X
                temp.append(gas.forward_rate_constants[gas.reaction_equations().index(self.reaction)])
            k_TP.append(temp)
        return np.array(k_TP)

    def get_cheby_kTP(self, prnt=False):
        '''Compute the rate constant at a given T and P using the Chebyshev coefficients.'''
        alpha = self.cheby_poly()
        k_TP = []
        TP=[]
        for T in self.T_ls:
            temp = []
            t = []
            for P in self.P_ls:
                T_tilde = self.reduced_T(T)
                P_tilde = self.reduced_P(P)
                logk = 0.0
                for i in range(self.n_T):
                    for j in range(self.n_P):
                        T_cheb = self.first_cheby_poly(T_tilde, i)
                        P_cheb = self.first_cheby_poly(P_tilde, j)
                        logk += np.log10(alpha[i, j]) * T_cheb * P_cheb
                temp.append(10 ** logk)
                t.append(T)
                # if prnt:
                #     print(f"T: {T}, P: {P}, T_tilde: {T_tilde}, P_tilde: {P_tilde}, logk: {logk}")
            k_TP.append(temp)
            TP.append(t)
        k_TP = np.array(k_TP)
        # if prnt:
        #     print('\nChebyshev k(TP)')
        #     print(k_TP)
        print(np.array(k_TP).flatten())
        print(np.array(TP))
        return np.array(k_TP)

    def reduced_T(self, T):
        return (2 * T - self.T_min - self.T_max) / (self.T_max - self.T_min)

    def reduced_P(self, P):
        return (2 * P - self.P_min - self.P_max) / (self.P_max - self.P_min)

    def first_cheby_poly(self, x, n):
        if n == 0:
            return 1
        elif n == 1:
            return x
        else:
            return 2 * x * self.first_cheby_poly(x, n - 1) - self.first_cheby_poly(x, n - 2)

    def cheby_poly(self):
        '''Fit the Chebyshev polynomials to rate constants.
           Input rate constants vector k should be arranged based on pressure.'''
        k_TP = self.get_YAML_kTP()
        cheb_mat = np.zeros((len(k_TP.flatten()), self.n_T * self.n_P))
        for n, P in enumerate(self.P_ls):
            for m, T in enumerate(self.T_ls):
                T_tilde = self.reduced_T(T)
                P_tilde = self.reduced_P(P)
                for i in range(self.n_T):
                    for j in range(self.n_P):
                        T_cheb = self.first_cheby_poly(T_tilde, i)
                        P_cheb = self.first_cheby_poly(P_tilde, j)
                        cheb_mat[n*len(self.T_ls)+m, i*self.n_P+j] = P_cheb * T_cheb
        coef = np.linalg.lstsq(cheb_mat, np.log10(k_TP.flatten()), rcond=None)[0].reshape((self.n_T, self.n_P))
        return coef

# Parameters for testing
T_list = np.linspace(200, 3000, 10)
P_list = np.linspace(9.869232667160128e-03, 100, 10)
n_T = 5
n_P = 5
X = {'Ar': 1.0}
reaction = 'H2O2 (+M) <=> 2 OH (+M)'
fname = 'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\sandbox.yaml'
colours = [(np.random.rand(), np.random.rand(), np.random.rand()) for _ in range(len(P_list))]

# Create chebFit_2D object
cheb_fit = chebFit_2D(T_list, P_list, n_P, n_T, X, reaction, fname, colours)

# Get k(T,P) data from both methods
yaml_kTP = cheb_fit.get_YAML_kTP()
cheby_kTP = cheb_fit.get_cheby_kTP(prnt=True)

# Compare the results
df_yaml = pd.DataFrame(yaml_kTP, index=T_list, columns=P_list)
df_cheby = pd.DataFrame(cheby_kTP, index=T_list, columns=P_list)

print("YAML k(T,P) DataFrame:\n", df_yaml)
print("Chebyshev k(T,P) DataFrame:\n", df_cheby)

comparison = df_yaml.compare(df_cheby)
print("Comparison of YAML and Chebyshev k(T,P) DataFrames:\n", comparison)
