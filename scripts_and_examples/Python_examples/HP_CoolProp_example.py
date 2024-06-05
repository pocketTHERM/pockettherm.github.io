# python script that serves multiple purposes:
#  1) in the first instance, this script serves as an example
#     of how the 'thermo_props.py' and 'hp_simulator.py' modules
#     can be used outside of the pocketTHERM environment to 
#     model and simulate heat pump systems
#  2) secondly, this script demonstrates how to use the same 
#     'hp_simulator.py' module, but rely on CoolProp for 
#     carrying out thermodynamic property calculations, rather
#     than using the Peng-Robinson model
#  3) conduct a direct comparison between cycle predictions 
#     obtained using the Peng-Robinson model and using CoolProp
#
# Note:
#   in order for this code to work, the CoolProp package needs
#   to be imported within the 'thermo_props.py' module. 
#
# Martin T. White, University of Sussex, 05/06/2024

import sys
sys.path.append('../python_modules/')
import thermo_props
import hp_simulator as hp
import matplotlib.pyplot as plt
import numpy as np

# heat exchanger discretisation:
n_hxc = [2,2,2,5,2,2,2,2]

# fluid properties:
fluids = [ "R1233zd(E)", "Pentane",     "MM" ]		# fluid names
Tc     = [  439.6,        469.7,        518.7     ]	# critical temperature, K
Pc     = [  3.6237e6,     3.3675e6,     1.93113e6 ]	# critical pressure, Pa
om     = [  0.3025,       0.2510,       0.4180    ]	# acentric factor
wm     = [  130.4962,     72.1488,      162.3768  ]	# molecular weight, g/mol
cp1    = [  33.3490,      12.9055,      64.6367   ]	# polynomial coefficients: 
cp2    = [  0.2823,       0.3906,       0.6695    ]     #   c_p(T) = cp1 + ... 
cp3    = [ -0.1523e-3,   -0.1036e-3,   -0.2895e-3 ]     #     cp2*T + cp3*T**2

# heat source and heat sink parameters:
cld = [323,1,4200]  # inlet temp. K, flow rate kg/s, specific heat J/kg K
hot = [328,1,4200]  # inlet temp. K, flow rate kg/s, specific heat J/kg K

# fixed ORC parameters:
Tevap  = 313	# evaporation temperature, K
dtsh   = 5	# compressore superheating, K
dtsc   = 5      # expansion valve subcooling, K
dt_ihe = 25	# internal heat exchanger temperature difference, K
Tco    = 318	# heat-source outlet temperature, K

# isentropic efficiency (compressor):
eta = [0.7]

# number of points in array:
n = 25

# pressure ratio array:
Prat = np.linspace(2,10,n)

# initialise arrays for data storage:
COP_PR = np.zeros([3,n])
COP_CP = np.zeros([3,n])

# run parametric study:
for i in range(len(fluids)):

    # setup Peng-Robinson fluid:
    cp = [ cp1[i], cp2[i], cp3[i] ] 
    PR = thermo_props.pr_fluid(fluids[i],Tc[i],Pc[i],om[i],cp,300,0.01,wm[i])

    # setup CoolProp fluid:
    CP = thermo_props.coolprop_fluid(fluids[i],'HEOS')

    for j in range(len(Prat)):

        # create array of inputs:
        x = [Tevap,dtsh,Prat[j],dtsc,dt_ihe,Tco]
    
        # run model using Peng-Robinson model:
        [props,cycle_out,_,pp,*_] = hp.simulate_HP(PR,x,eta,hot,cld,n_hxc,0)
        COP_PR[i,j] = cycle_out[0]

        # run model using CoolProp:
        try:
            [props,cycle_out,_,pp,*_] = hp.simulate_HP(CP,x,eta,hot,cld,n_hxc,0)
            COP_CP[i,j] = cycle_out[0]
        except:
            COP_CP[i,j] = np.nan

plt.plot(Prat,COP_CP[0,:],'k-',label='R1233zd(E) - HEOS')
plt.plot(Prat,COP_CP[1,:],'r-',label='n-pentane - HEOS')
plt.plot(Prat,COP_CP[2,:],'b-',label='MM - HEOS')
plt.plot(Prat,COP_PR[0,:],'k--',label='R1233zd(E) - PR')
plt.plot(Prat,COP_PR[1,:],'r--',label='n-pentane - PR')
plt.plot(Prat,COP_PR[2,:],'b--',label='MM - PR')
plt.legend()
plt.xlabel('Pressure ratio')
plt.ylabel('Coefficient of performance')
plt.show()

# End of file
