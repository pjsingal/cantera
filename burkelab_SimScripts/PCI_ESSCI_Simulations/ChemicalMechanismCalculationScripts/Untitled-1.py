#%%
import cantera as ct
import numpy as np

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

def reduced_T(T, T_min, T_max):
    '''Calculate the reduced temperature.'''
    T_tilde = 2. * T ** (-1) - T_min ** (-1) - T_max ** (-1)
    T_tilde /= (T_max ** (-1) - T_min ** (-1))
    return T_tilde

def reduced_P(P, P_min, P_max):
    '''Calculate the reduced pressure.'''
    P_tilde = 2. * np.log(P) - np.log(P_min) - np.log(P_max)
    P_tilde /= (np.log(P_max) - np.log(P_min))
    return P_tilde

def cheby_poly(n_T, n_P, k, T_ls, P_ls, P_min=0.01, P_max=100, T_min=200, T_max=3000):
    '''Fit the Chebyshev polynominals to rate constants.
        Input rate constants vector k should be arranged based on pressure.'''
    # modified for abstraction reactions
    if P_ls == ['--']:
        P_ls = [1.0]

    cheb_mat = np.zeros((len(k), n_T * n_P))
    for n, P in enumerate(P_ls):       # !! assume that at each presssure, we have the same temperateure range
        for m, T in enumerate(T_ls):
            for i in range(n_T):
                T_tilde = reduced_T(T, T_min, T_max)
                T_cheb = first_cheby_poly(T_tilde, i)
                for j in range(n_P):
                    P_tilde = reduced_P(P, P_min, P_max)
                    P_cheb = first_cheby_poly(P_tilde, j)
                    cheb_mat[n*len(T_ls)+m, i*n_P+j] = P_cheb * T_cheb
    coef = np.linalg.lstsq(cheb_mat, np.log10(np.array(k)))[0]
    return coef

model = 'C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\sandbox.yaml'
gas = ct.Solution(model)

reaction = 'H + O2 (+M) <=> HO2 (+M)'

T_list = np.linspace(200,4000)
# P_list = np.linspace(200,4000)
P_list=[1]
k = []
for T in T_list:
    temp = []
    for P in P_list:
        gas.TPX = T, P*ct.one_atm, 'Ar:1'
        temp.append(gas.forward_rate_constants[gas.reaction_equations().index(reaction)])
    k.append(temp)

k = np.array(k)
coef = cheby_poly(7, 1, k, T_list, P_list)
# %%
print(coef)
# %%
