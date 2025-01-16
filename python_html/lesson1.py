from pyscript import document
from pyscript import display
from js import console
from pyodide.ffi import create_proxy 
import csv
import matplotlib.pyplot as plt
import numpy as np
import thermo_props

# initialise figures:
fig, ax = plt.subplots()
fig.set_size_inches(4, 4)
ax.set_position([0.175,0.125,0.80,0.85])
display(fig,target="plt-psat-single")
display(fig,target="plt-tcond-pr")

plt.tight_layout()
#plt.close(fig)

# fontsize for in-figure text:
fs = 7

# initialise fluid:
fluid = thermo_props.pr_fluid("n-pentane",469.7,3.3675e6,0.2510,
                               [12.9055,0.3906,-0.1036e-3],300,0.01,72.1488)

# read csv file containing saturation curve:
with open('n-pentane.csv', mode='r') as file:
    csv_reader = csv.DictReader(file, quoting=csv.QUOTE_NONNUMERIC)

    # initialize an empty dictionary for columns
    df = {header: [] for header in csv_reader.fieldnames}

    # populate the dictionary
    for row in csv_reader:
        for key in row:
            df[key].append(row[key]) 

ssat = df['s_sat']
tsat = df['T_sat']

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

# update saturation curve:
ssat[:] = [x - smin for x in ssat]

def _p_saturation_single(*args, **kwargs):
    
    # get pressure:
    p = float(document.getElementById("inpt-psat-single").value)*1e5
    
    # saturated liquid properties:
    fluid.saturation_p(p,0)
    sl = fluid.fluid.smass()
    Tl = fluid.fluid.T()
    dl = fluid.fluid.rhomass()
    hl = fluid.fluid.hmass()
    
    # saturated vapour properties:
    fluid.saturation_p(p,1)
    sv = fluid.fluid.smass()
    Tv = fluid.fluid.T()
    dv = fluid.fluid.rhomass()
    hv = fluid.fluid.hmass()
    
    # clear T-s plot:   
    document.getElementById("plt-psat-single").innerHTML = ""
    
    # plot figure
    fig, ax = plt.subplots()
    ax.plot(ssat,tsat,'k-')
    ax.plot(sl-smin,Tl,'bo')
    ax.plot(sv-smin,Tv,'ro')
    ax.set_xlabel('Entropy, s [J/(kg K)]')
    ax.set_ylabel('Temperature, T [K]')
    ax.set_xlim((0,smax-smin))
    ax.set_ylim((tmin,tmax))
    
    # strings for plot:
    T_sat_str  = 'Saturation temperature = ' + "{:.2f}".format(Tl) + ' K'
    p_sat_str  = 'Saturation pressure = ' + "{:.2f}".format(p/1e5) + ' bar'
    d_satl_str = 'Liquid density = ' + "{:.2f}".format(dl) + ' kg/m3'
    h_satl_str = 'Liquid enthalpy = ' + "{:.2f}".format(hl/1e3) + ' kJ/kg'
    s_satl_str = 'Liquid entropy = ' + "{:.2f}".format(sl-smin) + ' J/(kg K)'
    d_satv_str = 'Vapour density = ' + "{:.2f}".format(dv) + ' kg/m3'
    h_satv_str = 'Vapour enthalpy = ' + "{:.2f}".format(hv/1e3) + ' kJ/kg'
    s_satv_str = 'Vapour entropy = ' + "{:.2f}".format(sv-smin) + ' J/(kg K)'
    
    # add strings to plot:
    ax.text(50,tmax-10,T_sat_str,color='k',size=fs)
    ax.text(50,tmax-20,p_sat_str,color='k',size=fs)
    ax.text(50,tmax-30,d_satl_str,color='b',size=fs)
    ax.text(50,tmax-40,h_satl_str,color='b',size=fs)
    ax.text(50,tmax-50,s_satl_str,color='b',size=fs)
    ax.text(50,tmax-60,d_satv_str,color='r',size=fs)
    ax.text(50,tmax-70,h_satv_str,color='r',size=fs)
    ax.text(50,tmax-80,s_satv_str,color='r',size=fs)
    
    # set position of figure:
    ax.set_position([0.175,0.125,0.80,0.85])
    fig.set_size_inches(4, 4)
    display(fig,target="plt-psat-single")
    #plt.close(fig)

def _evaporation_condensation(*args, **kwargs):
    
    # get inputs:
    Tcond = float(document.getElementById("inpt-tcond").value)
    pr    = float(document.getElementById("inpt-pr").value)
    
    # condensation conditions:
    fluid.saturation_t(Tcond,0)
    p1 = fluid.fluid.p()
    s1l = fluid.fluid.smass()
    fluid.saturation_t(Tcond,1)
    s1v = fluid.fluid.smass()
    
    # evaporation conditions:
    p2 = p1*pr
    fluid.saturation_p(p2,0)
    T2 = fluid.fluid.T()
    s2l = fluid.fluid.smass()
    fluid.saturation_p(p2,1)
    s2v = fluid.fluid.smass()
    
    # update T-s plot:   
    document.getElementById("plt-tcond-pr").innerHTML = ""
    
    #
    s1av = 0.5*(s1l+s1v)
    s2av = 0.5*(s2l+s2v)
    
    # plot figure
    fig, ax = plt.subplots()
    ax.plot(ssat,tsat,'k-',linewidth=1)
    ax.plot(np.array([s1l,s1v])-smin,np.array([Tcond,Tcond]),'go-',linewidth=1,markersize=3)
    ax.plot(np.array([s2l,s2v])-smin,np.array([T2,T2]),'go-',linewidth=1,markersize=3)
    ax.set_xlabel('Entropy, s [J/(kg K)]')
    ax.set_ylabel('Temperature, T [K]')
    
    ax.text(s1av-smin,Tcond-10,'Condensation',color='g',ha='center')
    ax.text(s2av-smin,T2-10,'Evaporation',color='g',ha='center')
    
    Tc_str  = 'Condensation temperature = ' + "{:.2f}".format(Tcond) + ' K'
    Pc_str  = 'Condensation pressure = ' + "{:.2f}".format(p1/1e5) + ' bar'
    Te_str  = 'Evaporation temperature = ' + "{:.2f}".format(T2) + ' K'
    Pe_str  = 'Evaporation pressure = ' + "{:.2f}".format(p2/1e5) + ' bar'
    
    ax.text(50,tmax-10,Tc_str,color='k',size=fs)
    ax.text(50,tmax-20,Pc_str,color='k',size=fs)
    ax.text(50,tmax-30,Te_str,color='k',size=fs)
    ax.text(50,tmax-40,Pe_str,color='k',size=fs)
    
    ax.set_xlim((0,smax-smin))
    ax.set_ylim((tmin,tmax))
    ax.set_position([0.175,0.125,0.80,0.85])
    fig.set_size_inches(4, 4)
    display(fig,target="plt-tcond-pr")
    #plt.close(fig)

# run function on click:    
p_saturation_single = create_proxy(_p_saturation_single)
document.getElementById(
    "btn-psat-single").addEventListener("click", p_saturation_single)

evaporation_condensation = create_proxy(_evaporation_condensation)
document.getElementById(
    "btn-tcond-pr").addEventListener("click", evaporation_condensation)
