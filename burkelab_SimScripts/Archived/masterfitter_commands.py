# Master Fitter

from masterFitter import *

T_list=np.linspace(200,2000,50)
# P_list=np.logspace(-6,12,num=140)
P_list=np.logspace(-6,12,num=10)
fname='C:\\Users\\pjsin\\Documents\\cantera\\test\\data\\sandbox.yaml'
# reaction='H2O2 (+M) <=> 2 OH (+M)'

reactions={'H + OH <=> H2O':{"N2":"{A: 1, b: 0, Ea: 0}",
                             "Ar":"{A: 2.20621e-02, b: 4.74036e-01, Ea: -1.13148e+02}",
                             "H2O":"{A: 1.04529e-01, b: 5.50787e-01, Ea: -2.32675e+02}"},
           'H + O2 <=> HO2':{"Ar":"{A: 1, b: 0, Ea: 0}",
                             "HE":"{A: 3.37601e-01, b: 1.82568e-01, Ea: 3.62408e+01}",
                             "N2":"{A: 1.24932e+02, b: -5.93263e-01, Ea: 5.40921e+02}",
                             "H2":"{A: 3.13717e+04, b: -1.25419e+00, Ea: 1.12924e+03}",
                             "CO2":"{A: 1.62413e+08, b: -2.27622e+00, Ea: 1.97023e+03}",
                             "NH3":"{A: 4.97750e+00, b: 1.64855e-01, Ea: -2.80351e+02}",
                             "H2O":"{A: 3.69146e+01, b: -7.12902e-02, Ea: 3.19087e+01}"},
           'H2O2 <=> 2 OH':{"Ar":"{A: 1, b: 0, Ea: 0}",
                             "N2":"{A: 1.14813e+00, b: 4.60090e-02, Ea: -2.92413e+00}",
                             "CO2":"{A: 8.98839e+01, b: -4.27974e-01, Ea: 2.41392e+02}",
                             "H2O2":"{A: 6.45295e-01, b: 4.26266e-01, Ea: 4.28932e+01}",
                             "H2O":"{A: 1.36377e+00, b: 3.06592e-01, Ea: 2.10079e+02}"},
           'NH3 <=> H + NH2':{"Ar":"{A: 1, b: 0, Ea: 0}",
                             "N2":"{A: 4.49281e+00, b: -9.46265e-02, Ea: -1.10071e+02}",
                             "O2":"{A: 1.15210e-01, b: 3.41234e-01, Ea: -3.89210e+02}",
                             "CO2":"{A: 9.19583e+00, b: 6.10696e-02, Ea: 9.01088e+01}",
                             "NH3":"{A: 1.49004e+01, b: 6.06535e-02, Ea: 2.47652e+02}",
                             "H2O":"{A: 1.14560e+01, b: 1.27501e-01, Ea: 3.13959e+02}"},
           '2 NH2 (+M) <=> N2H4 (+M)':{"Ar":"{A: 1, b: 0, Ea: 0}",
                             "N2":"{A: 1.46848e+01, b: -3.00962e-01, Ea: 1.65146e+02}",
                             "O2":"{A: 1.53608e+00, b: -3.90857e-02, Ea: 4.43752e+00}",
                             "NH3":"{A: 1.64196e+01, b: -7.29636e-02, Ea: 3.66099e+02}",
                             "H2O":"{A: 2.17658e+01, b: -1.14715e-01, Ea: 4.19216e+02}"},
           'HNO <=> H + NO':{"Ar":"{A: 1, b: 0, Ea: 0}",
                             "N2":"{A: 7.56771e+00, b: -1.96878e-01, Ea: 2.53168e+02}",
                             "H2O":"{A: 3.35451e+00, b: 1.33352e-01, Ea: 6.64633e+01}"}
}
# reactions=reaction_dict.keys()
# colliders=reaction_dict.values()
n_T=7
n_P=7
#consider adding more rows/cols to cheb mat

colours=[(np.random.rand(), np.random.rand(), np.random.rand()) for _ in range(len(P_list))]

# colliders = ["Ar","H2O","CO2", "N2", "H2O2"]

# for reaction in reactions.keys():
#     # for collider in reactions[reaction]:
#     #     cheb_fit = chebFit_2D(T_list,P_list,n_P, n_T, collider, reaction, fname,colours)
#     #     cheb_fit.print_cheby_YAML_format()  
#     for collider in reactions[reaction]:
plog_fit = PLOGFit(T_list,P_list,reactions,fname)
plog_fit.print_PLOG_YAML_format()       

# for j,reaction in enumerate(reactions):
#     # for i in range(len(colliders[j])):
#     #     cheb_fit = chebFit_2D(T_list,P_list,n_P, n_T, colliders[j][i], reaction, fname,colours)
#     #     cheb_fit.print_cheby_YAML_format()

#     for i,c in enumerate(colliders):
#         plog_fit = PLOGFit(T_list,P_list,c,reaction,fname)
#         plog_fit.print_PLOG_YAML_format()

# python "burkelab_SimScripts\\ChemicalMechanismCalculationScripts\\masterfitter_commands.py"


