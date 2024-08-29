import os 
import io
import numpy as np
import sys
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
    

def first_cheby_poly(x, n):
    '''Generate n-th order Chebyshev ploynominals of first kind.'''
    if n == 0: return 1
    elif n == 1: return x
    result = 2. * x * first_cheby_poly(x, 1) - first_cheby_poly(x, 0)
    m = 0
    while n - m > 2:
        result = 2. * x * result - first_cheby_poly(x, m+1)
        m += 1
    return result

def reduced_T( T, T_min, T_max):
    '''Calculate the reduced temperature.'''
    T = np.array(T)
    T_tilde = 2. * T ** (-1.0) - T_min ** (-1.0) - T_max ** (-1.0)
    T_tilde /= (T_max ** (-1.0) - T_min ** (-1.0))
    return T_tilde
   
def calc_polynomial(T,alpha):
    #calculate rate constants helper function
    T_reduced_list = reduced_T(T,200,2400)
    # T_reduced_list = reduced_T(T,T[0],T[-1])
    values = np.polynomial.chebyshev.chebval(T_reduced_list,alpha)
    return values

def fit_cheby_poly_1d(n_T, k, T_ls):
    #T needs to be a lsit 
    '''Fit the Chebyshev polynominals to rate constants.
       Input rate constants vector k should be arranged based on pressure.'''
    cheb_mat = np.zeros((len(k), n_T))
    for m, T in enumerate(T_ls):
        T_min = T_ls[0]
        T_max = T_ls[-1]
        for i in range(n_T):
            T_tilde = reduced_T(T, T_min, T_max)
            T_cheb = first_cheby_poly(T_tilde, i)

            cheb_mat[m,i] =  T_cheb
            #log_k = np.log10(np.array(k))
    coef,b,c,d = np.linalg.lstsq(cheb_mat,k,rcond=-1)
    
    return coef

def run_cantera_calculate_rate_constant(T,cti_file):
    gas = ct.Solution(cti_file)
    cantera_k = []
    for Temp in T:
        gas.TPX = Temp,101325,{'Ar':1}
        cantera_k.append(gas.forward_rate_constants[0])
    return np.array(cantera_k)*1000

def run_fitter(n_T,T_ls,cti_file):
    k_from_troe = run_cantera_calculate_rate_constant(T_ls,cti_file)
    alpha = fit_cheby_poly_1d(n_T,k_from_troe,T_ls)
    k_calculated = calc_polynomial(T_ls,alpha)
    # plt.figure()
    # plt.semilogy(T_ls,k_calculated, label='k from Cheb')
    # plt.semilogy(T_ls,k_from_troe, label='k from Troe')
    # plt.show()
    print(alpha)
    
# T_ls = temperature list (in increasing order)
# n_T number of rows in data matrix. This number is the highest order of T_i minus one that the matrix will use
# (e.g., n_T=3 means the bottom row of "data" represents a second-order polynomial)

# k = number of columns in data matrix


# gas = ct.Solution(file)

# k_list=[]
# for i, T in enumerate(T_ls):
#     gas.TPX = T_ls[i],101325,{'AR':1.0}
#     rc = gas.forward_rate_constants[gas.reaction_equations().index('H + O2 (+M) <=> HO2 (+M)')]
#     k_list.append(rc)
# k_list = run_cantera_calculate_rate_constant(T_ls,file)

n_T=20
T_ls=np.linspace(200,2000)
file = 'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\sandbox.yaml'
run_fitter(n_T,T_ls,file)