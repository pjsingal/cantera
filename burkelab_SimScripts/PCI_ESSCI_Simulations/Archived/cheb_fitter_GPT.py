import sys
sys.path.append("C:/Users/pjsin/Documents/cantera/build/python")
import cantera as ct
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def troe_rate_constants(low_rate, high_rate, troe_params, T, P):
    """
    Calculate the rate constants for a Troe reaction.
    """
    k0 = low_rate[0] * (T / 300)**low_rate[1] * np.exp(-low_rate[2] / (ct.gas_constant * T))
    k_inf = high_rate[0] * (T / 300)**high_rate[1] * np.exp(-high_rate[2] / (ct.gas_constant * T))
    Pr = k0 * P / k_inf
    F_cent = (1 - troe_params[0]) * np.exp(-T / troe_params[1]) + \
             troe_params[0] * np.exp(-T / troe_params[2]) + \
             np.exp(-troe_params[3] / T)
    logF_cent = np.log10(F_cent)
    c = -0.4 - 0.67 * logF_cent
    n = 0.75 - 1.27 * logF_cent
    f = (logF_cent / (n - 0.14 * logF_cent))**2
    logF = logF_cent / (1 + f)
    F = 10**logF
    k = k_inf * (Pr / (1 + Pr)) * F
    return k

def chebyshev_fit(rate_constants, T_min, T_max, P_min, P_max, n_T, n_P):
    """
    Fit the rate constants to a Chebyshev polynomial.
    """
    T_samples = np.linspace(T_min, T_max, n_T)
    P_samples = np.linspace(P_min, P_max, n_P)
    T_mid = (T_max + T_min) / 2
    P_mid = (P_max + P_min) / 2
    T_scale = (T_max - T_min) / 2
    P_scale = (P_max - P_min) / 2
    
    def scaled_T(T):
        return (T - T_mid) / T_scale

    def scaled_P(P):
        return (P - P_mid) / P_scale

    T_mesh, P_mesh = np.meshgrid(T_samples, P_samples)
    k_mesh = np.log(rate_constants(T_mesh, P_mesh))

    coeffs = np.polynomial.chebyshev.chebfit(np.ravel(scaled_T(T_mesh)), np.ravel(scaled_P(P_mesh)), k_mesh.flatten(), [n_T-1, n_P-1])
    return coeffs

# Define the Troe reaction parameters
low_rate = [6.366e+20, -1.72, 524.8 * 4184]  # A, b, Ea in J/mol
high_rate = [4.7e+12, 0.44, 0.0]  # A, b, Ea in J/mol
troe_params = [0.5, 1.0e-30, 1.0e+30, 0]  # A, T3, T1, T2

# Define the temperature and pressure range for fitting
T_min = 300
T_max = 2000
P_min = ct.one_atm
P_max = 100 * ct.one_atm

# Generate the rate constants
def rate_constants(T, P):
    return troe_rate_constants(low_rate, high_rate, troe_params, T, P)

# Fit the Chebyshev polynomial
n_T = 6  # Number of terms in temperature
n_P = 4  # Number of terms in pressure
coeffs = chebyshev_fit(rate_constants, T_min, T_max, P_min, P_max, n_T, n_P)

# Display the Chebyshev coefficients
print("Chebyshev coefficients:")
print(coeffs)

# Optionally, plot the fitted rate constants
T_fit = np.linspace(T_min, T_max, 100)
P_fit = np.linspace(P_min, P_max, 100)
T_mesh_fit, P_mesh_fit = np.meshgrid(T_fit, P_fit)
k_fit = np.polynomial.chebyshev.chebval2d((T_fit - (T_min + T_max) / 2) / ((T_max - T_min) / 2),
                                          (P_fit - (P_min + P_max) / 2) / ((P_max - P_min) / 2),
                                          coeffs)
plt.contourf(T_mesh_fit, P_mesh_fit, k_fit)
plt.colorbar(label='Log(rate constant)')
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (Pa)')
plt.title('Chebyshev Polynomial Fit')
plt.show()
