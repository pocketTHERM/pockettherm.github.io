from js import document
from pyodide.ffi import create_proxy 
import thermo_props
import orc_simulator
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# initialise figures:
fig, ax = plt.subplots()
fig.set_size_inches(4, 4)
ax.set_position([0.175,0.125,0.80,0.85])
display(fig,target="plt-hot-out")
display(fig,target="plt-cold-mass")
plt.close(fig)

# initialise fluid:
fluid = thermo_props.pr_fluid("n-pentane",469.7,3.3675e6,0.2510,
                               [12.9055,0.3906,-0.1036e-3],300,0.01,72.1488)

# load saturation curve:
df   = pd.read_csv(r'./n-pentane.csv')
ssat = df['s_sat']
tsat = df['T_sat']

# fontsize for in-figure text:
fs = 7

# saturation limits:
smin = min(ssat)
smax = max(ssat)
tmin = min(tsat)
tmax = max(tsat)
ds   = smax - smin
dt   = tmax - tmin
smin = smin - 0.1*ds
smax = smax + 0.1*ds
tmin = tmin - 0.1*dt
tmax = tmax + 0.1*dt

# heat exchanger discretisation (not user defined):
n_hxc = [2,2,2,5,2,2,2,2]

# inputs:
hot = [423,10,1000]            # heat-source conditions
cld = [288,1000,4200]         # heat-sink conditions
eta = [0.7,0.8]               # component efficiencies
x   = [303,2,8,1.5,0.8,423]  # cycle variables
   
def _heat_source_outlet(*args, **kwargs):
    
    # get inputs:
    x[5] = float(document.getElementById("inpt-hot-out").value)

    # simulate cycle:
    [props,cycle_out,UA,pp,th,tc] = orc_simulator.simulate_ORC(fluid,x,eta,hot,cld,n_hxc,0)

    # print cycle outputs:
    eta_str = 'Thermal efficiency = ' + "{:.2f}".format(cycle_out[1]) + ' %'
    m_str   = 'Cycle mass-flow rate = ' + "{:.2f}".format(cycle_out[8]) + ' kg/s'
    Wn_str  = 'Net power generated = ' + "{:.2f}".format(cycle_out[0]) + ' kW'
    
    # clear T-s plot:   
    document.getElementById("plt-hot-out").innerHTML = ""
    
    # update T-s plot:
    fig, ax = orc_simulator.plot_cycle_ts(fluid,props,th,tc,df)
#    ax.set_ylim((tmin,tmax))
    ax.set_xlim(-1200,600)
    ax.set_ylim(273,473)
    ax.text(-1150,473-8,eta_str,color='k',size=fs)
    ax.text(-1150,473-16,m_str,color='k',size=fs)
    ax.text(-1150,473-24,Wn_str,color='k',size=fs)    
    ax.set_position([0.175,0.125,0.80,0.85])
    fig.set_size_inches(4, 4)
    display(fig,target="plt-hot-out")
    plt.close(fig)
    
def _heat_sink_mass(*args, **kwargs):
    
    x[5] = 373
    cld[1] = float(document.getElementById("inpt-cold-mass").value)
    
    # simulate cycle:
    [props,cycle_out,UA,pp,th,tc] = orc_simulator.simulate_ORC(fluid,x,eta,hot,cld,n_hxc,0)
    
    # clear T-s plot:   
    document.getElementById("plt-cold-mass").innerHTML = ""
    
    # update T-s plot:
    fig, ax = orc_simulator.plot_cycle_ts(fluid,props,th,tc,df)
#    ax.set_ylim((tmin,tmax))
    ax.set_xlim(-1200,600)
    ax.set_ylim(273,473)
    ax.set_position([0.175,0.125,0.80,0.85])
    fig.set_size_inches(4, 4)
    display(fig,target="plt-cold-mass")
    plt.close(fig)
    
    
# run function on click:    
heat_source_outlet = create_proxy(_heat_source_outlet)
document.getElementById(
    "btn-hot-out").addEventListener("click", heat_source_outlet)

heat_sink_mass = create_proxy(_heat_sink_mass)
document.getElementById(
    "btn-cold-mass").addEventListener("click", heat_sink_mass)
