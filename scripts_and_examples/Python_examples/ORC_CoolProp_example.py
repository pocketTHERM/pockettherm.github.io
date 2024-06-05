# python script that serves multiple purposes:
#  1) in the first instance, this script serves as an example
#     of how the 'thermo_props.py' and 'orc_simulator.py' modules
#     can be used outside of the pocketTHERM environment to 
#     model and simulate ORC power systems
#  2) secondly, this script demonstrates how to use the same 
#     'orc_simulator.py' module, but rely on CoolProp for 
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
import orc_simulator as orc
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
hot = [523,1,1000]  # inlet temp. K, flow rate kg/s, specific heat J/kg K
cld = [288,1,4200]  # inlet temp. K, flow rate kg/s, specific heat J/kg K

# fixed ORC parameters:
Tcond = 313	# condensation temperature, K
dtsc  = 2	# pump subcooling, K
x3    = 1.2     # expander inlet parameter
eff_r = 0.7     # recuperator effectivness
Tho   = 423	# heat-source outlet temperature, K

# isentropic efficiencies (pump and turbine):
eta = [0.7,0.8]

# number of points in array:
n = 25

# maxium reduced pressure:
Pred_max = 0.9

# initialise arrays for data storage:
Pev_PR = np.zeros([3,n])
eta_PR = np.zeros([3,n])
Pev_CP = np.zeros([3,n])
eta_CP = np.zeros([3,n])

# run parametric study:
for i in range(len(fluids)):

    # setup Peng-Robinson fluid:
    cp = [ cp1[i], cp2[i], cp3[i] ] 
    PR = thermo_props.pr_fluid(fluids[i],Tc[i],Pc[i],om[i],cp,300,0.01,wm[i])

    # setup CoolProp fluid:
    CP = thermo_props.coolprop_fluid(fluids[i],'HEOS')

    # setup array of pressure ratios:
    PR.saturation_t(Tcond,0)
    Pcond = PR.fluid.p()
    Prat_max = Pred_max * PR.pc / Pcond
    Prat = np.linspace(2,Prat_max,n)

    for j in range(len(Prat)):

        # create array of inputs:
        x = [Tcond,dtsc,Prat[j],x3,eff_r,Tho]
    
        # run model using Peng-Robinson model:
        [props,cycle_out,_,pp,*_] = orc.simulate_ORC(PR,x,eta,hot,cld,n_hxc,0)
        Pev_PR[i,j] = props[1,1]/1e5
        if min(pp[0]) > 0:
            eta_PR[i,j] = cycle_out[1]
        else:
            eta_PR[i,j] = np.nan

        # run model using CoolProp:
        [props,cycle_out,_,pp,*_] = orc.simulate_ORC(CP,x,eta,hot,cld,n_hxc,0)
        Pev_CP[i,j] = props[1,1]/1e5
        if min(pp[0]) > 0:
            eta_CP[i,j] = cycle_out[1]
        else:
            eta_CP[i,j] = np.nan

plt.plot(Pev_CP[0,:],eta_CP[0,:],'k-',label='R1233zd(E) - HEOS')
plt.plot(Pev_CP[1,:],eta_CP[1,:],'r-',label='n-pentane - HEOS')
plt.plot(Pev_CP[2,:],eta_CP[2,:],'b-',label='MM - HEOS')
plt.plot(Pev_PR[0,:],eta_PR[0,:],'k--',label='R1233zd(E) - PR')
plt.plot(Pev_PR[1,:],eta_PR[1,:],'r--',label='n-pentane - PR')
plt.plot(Pev_PR[2,:],eta_PR[2,:],'b--',label='MM - PR')
plt.xlabel('Evaporation pressure [bar]')
plt.ylabel('Cycle thermal efficiency [%]')
plt.legend()
plt.show()

# End of file
