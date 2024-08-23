
import sys, os
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy.polynomial.chebyshev import chebval

class chebFit_2D:
    def __init__(self, T_ls, P_ls, n_P, n_T, collider, reaction, fname,colours):
        self.T_ls = T_ls
        self.P_ls = P_ls
        self.n_P=n_P
        self.n_T=n_T
        self.collider = collider
        self.reaction=reaction
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
                gas.TPX = T, P*ct.one_atm, {self.collider:1.0}
                temp.append(gas.forward_rate_constants[gas.reaction_equations().index(self.reaction)])
            k_TP.append(temp)
        return np.array(k_TP)
    
    # def getKfromTroe(self):
    #     gas = ct.Solution(self.fname)
    #     def getK(Temp,Pres,X) :
    #         gas.TPX = Temp,Pres,X
    #         k = gas.forward_rate_constants[gas.reaction_equations().index(self.reaction)]
    #         return k
    #     # plt.figure()
    #     k_P=[]
    #     for P in self.P_ls:
    #         k_T = []
    #         for T in self.T_ls:
    #             k_T.append(getK(T,P,self.X))
    #         k_P.append(k_T)
    #         # plt.plot(T_list,np.log(k_T),label=str(P/101325)+" atm")
    #     # plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    #     # plt.title("Troe Rates for "+self.reaction)
    #     # plt.xlabel("Temperature")
    #     # plt.ylabel("logK [units still undetermined]")
    #     # plt.show()
    #     return k_P
    
    # def test_YAML_kTP(self, T_test, P_test,prnt=False):
    #     gas = ct.Solution(self.fname)
    #     k_TP = []
    #     for T in T_test:
    #         temp = []
    #         for P in P_test:
    #             gas.TPX = T, P*ct.one_atm, self.X
    #             temp.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
    #         k_TP.append(temp)
    #     if prnt == True:
    #         print('\nYAML k(TP)')
    #         print(k_TP)
    #     return np.array(k_TP)

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
    
    # def first_cheby_poly(self, x, n):
    #     '''Generate n-th order Chebyshev polynomials of the first kind.'''
    #     coeffs = np.zeros(n + 1)
    #     coeffs[n] = 1
    #     return chebval(x, coeffs)

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


    def cheby_poly(self):
        '''Fit the Chebyshev polynominals to rate constants.
            Input rate constants vector k should be arranged based on pressure.'''
        k_TP = self.get_YAML_kTP()
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


    # def cheby_poly(self):
    #     '''Fit the Chebyshev polynomials to rate constants.
    #     Input rate constants vector k should be arranged based on pressure.'''
    #     k_TP = self.get_YAML_kTP()
    #     cheb_mat = np.zeros((len(k_TP.flatten()), self.n_T * self.n_P))
    #     # for n, P in enumerate(self.P_ls):
    #     #     for m, T in enumerate(self.T_ls):
    #     for m, T in enumerate(self.T_ls):
    #         for n, P in enumerate(self.P_ls):
    #             T_tilde = self.reduced_T(T)
    #             P_tilde = self.reduced_P(P)
    #             T_cheb_vals = chebval(T_tilde, np.eye(self.n_T))
    #             P_cheb_vals = chebval(P_tilde, np.eye(self.n_P))
    #             for i in range(self.n_T):
    #                 for j in range(self.n_P):
    #                     # cheb_mat[n * len(self.T_ls) + m, i * self.n_P + j] = T_cheb_vals[i] * P_cheb_vals[j]
    #                     cheb_mat[m * len(self.P_ls) + n, i * self.n_P + j] = T_cheb_vals[i] * P_cheb_vals[j]
    #     coef = np.linalg.lstsq(cheb_mat, np.log10(k_TP.flatten()), rcond=None)[0].reshape((self.n_T, self.n_P))
    #     return coef

    def get_cheby_kTP(self, prnt=False):
        '''Compute the rate constant at a given T and P using the Chebyshev coefficients.'''
        alpha = self.cheby_poly()
        k_TP = []
        for T in self.T_ls:
            temp = []
            for P in self.P_ls:
                T_tilde = self.reduced_T(T)
                P_tilde = self.reduced_P(P)
                T_cheb_vals = chebval(T_tilde, np.eye(self.n_T))
                P_cheb_vals = chebval(P_tilde, np.eye(self.n_P))
                logk = np.sum(alpha * T_cheb_vals[:, None] * P_cheb_vals[None, :])
                temp.append(10 ** logk)
            k_TP.append(temp)
        k_TP = np.array(k_TP)
        if prnt:
            print('\nChebyshev k(TP)')
            print(k_TP)
        return np.array(k_TP)
    
    # def get_cheby_kTP(self, prnt=False):
    #     '''Compute the rate constant at a given T and P using the Chebyshev coefficients.'''
    #     alpha = self.cheby_poly()
    #     k_TP = []
    #     for T in self.T_ls:
    #         temp = []
    #         for P in self.P_ls:
    #             T_tilde = self.reduced_T(T)
    #             P_tilde = self.reduced_P(P)
    #             logk = np.polynomial.chebyshev.chebval(T_tilde, alpha) @ np.polynomial.chebyshev.chebval(P_tilde, alpha.T)
    #             temp.append(10 ** logk)
    #         k_TP.append(temp)
    #     k_TP = np.array(k_TP)
    #     if prnt:
    #         print('\nChebyshev k(TP)')
    #         print(k_TP)
    #     return np.array(k_TP)

# def calc_polynomial(self):
#     #calculate rate constants helper function
#     alpha = self.cheby_poly()
#     T_tilde = self.reduced_T(T)
#     values = np.polynomial.chebyshev.chebval(T_tilde,alpha)

#     alpha = self.cheby_poly()
#     k_TP = []
#     for T in self.T_ls:
#         temp = []
#         for P in self.P_ls:
#             T_tilde = self.reduced_T(T)
#             P_tilde = self.reduced_P(P)
#             T_cheb_vals = chebval(T_tilde, np.eye(self.n_T))
#             P_cheb_vals = chebval(P_tilde, np.eye(self.n_P))
#             logk = np.sum(alpha * T_cheb_vals[:, None] * P_cheb_vals[None, :])
#             temp.append(10 ** logk)
#         k_TP.append(temp)
#     k_TP = np.array(k_TP)
#     if prnt:
#         print('\nChebyshev k(TP)')
#         print(k_TP)
#     return np.array(k_TP)
#     return values

    # def get_cheby_kTP(self, prnt=False):
    #     '''Compute the rate constant at a given T and P using the Chebyshev coefficients.'''
    #     alpha = self.cheby_poly()
    #     k_TP = []
    #     for T in self.T_ls:
    #         temp = []
    #         for P in self.P_ls:
    #             T_tilde = self.reduced_T(T)
    #             P_tilde = self.reduced_P(P)
    #             logk = 0.0
    #             for i in range(self.n_T):
    #                 for j in range(self.n_P):
    #                     T_cheb = self.first_cheby_poly(T_tilde, i)
    #                     P_cheb = self.first_cheby_poly(P_tilde, j)
    #                     logk += alpha[i, j] * T_cheb * P_cheb
    #                     # logk += np.log10(alpha[i, j]) * T_cheb * P_cheb
    #             temp.append(10 ** logk)
    #             # temp.append(logk)
    #         k_TP.append(temp)
    #         # for P in self.P_ls:
    #         #     T_tilde = self.reduced_T(T)
    #         #     P_tilde = self.reduced_P(P)
    #         #     T_cheb_vals = np.polynomial.chebyshev.chebval(T_tilde, np.identity(self.n_T))
    #         #     P_cheb_vals = np.polynomial.chebyshev.chebval(P_tilde, np.identity(self.n_P))
    #         #     logk = np.sum(alpha * T_cheb_vals[:, None] * P_cheb_vals[None, :])
    #         #     temp.append(10 ** logk)
    #         # k_TP.append(temp)
    #     k_TP = np.array(k_TP)
    #     if prnt == True:
    #         print('\nChebyshev k(TP)')
    #         print(k_TP)
    #     return np.array(k_TP)
    
    # def test_cheby_kTP(self, T_test, P_test, prnt=False):
    #     '''Compute the rate constant at a given T and P using the Chebyshev coefficients.'''
    #     alpha = self.cheby_poly()
    #     k_TP = []
    #     for T in T_test:
    #         temp = []
    #         for P in P_test:
    #             T_tilde = self.reduced_T(T)
    #             P_tilde = self.reduced_P(P)
    #             logk = 0.0
    #             for i in range(self.n_T):
    #                 for j in range(self.n_P):
    #                     T_cheb = self.first_cheby_poly(T_tilde, i)
    #                     P_cheb = self.first_cheby_poly(P_tilde, j)
    #                     logk += alpha[i, j] * T_cheb * P_cheb
    #             temp.append(10 ** logk)
    #         k_TP.append(temp)
    #     if prnt == True:
    #         print('\nChebyshev k(TP)')
    #         print(k_TP)
    #     return np.array(k_TP)
    
    def plot_cheby_kTP(self):
        k_TP = self.get_cheby_kTP()
        for j in range(len(self.P_ls)):
            k_T = []
            for row in k_TP:
                k_T.append(row[j])
            # label = "cheb "+ f'{self.P_ls[j]:.2e}'+" atm"
            if j ==0:
                plt.semilogy(self.T_ls,k_T,label="Chebyshev",linestyle='none',marker="x",markersize=5,color=self.colours[j])
            else:
                plt.semilogy(self.T_ls,k_T,label=None,linestyle='none',marker="x",markersize=5,color=self.colours[j])
    
    def plot_Troe_kTP(self, linestyle):
        k_TP = self.get_YAML_kTP()
        for j in range(len(self.P_ls)):
            k_T = []
            for row in k_TP:
                k_T.append(row[j])
            # label = "troe "+ f'{self.P_ls[j]:.2e}'+" atm"
            if j ==0:
                plt.semilogy(self.T_ls,k_T,label="Troe",linestyle=linestyle,color=self.colours[j])
            else:
                plt.semilogy(self.T_ls,k_T,label=None,linestyle=linestyle,color=self.colours[j])
    
    def print_cheby_YAML_format(self):
        coef = self.cheby_poly()
        # formatted_list = ', '.join([f'{x:.4f}' for x in coef)
        print(("- equation: %s")%(self.reaction))
        print("  type: Chebyshev")
        print(("  temperature-range: [%.1f, %.1f]")%(self.T_min,self.T_max))
        print(("  pressure-range: [%.3e atm, %.3e atm]")%(self.P_min,self.P_max))
        print("  data:  # "+self.collider)
        for i in range(len(coef)):
            string = ""
            for j in range(len(coef[0])):
                string+=f'{coef[i,j]:.4e}'
                if j<len(coef[0])-1:
                    string+=", "
            print("  - ["+string+"]")
        print("\n\n\n")

