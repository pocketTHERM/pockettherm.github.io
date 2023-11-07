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
display(fig,target="plt-sub-super")
display(fig,target="plt-4-variables")
plt.close(fig)

# fontsize for in-figure text:
fs = 7

# initialise fluid:
fluid = thermo_props.pr_fluid("n-pentane",469.7,3.3675e6,0.2510,
                               [12.9055,0.3906,-0.1036e-3],300,0.01,72.1488)

# load saturation curve:
df   = pd.read_csv(r'./n-pentane.csv')
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
    
def _subcooled_superheated(*args, **kwargs):
    
    # pressure
    p = 5e5
    
    # get input temperature:
    T = float(document.getElementById("inpt-sub-super").value)
    
    # saturation temperature:
    fluid.saturation_p(p,0)
    Tsat = fluid.fluid.T()
    
    # construct liquid subcooling:
    Tl = np.linspace(273,Tsat,50)
    sl = np.zeros(len(Tl))
    for i in range(len(Tl)):
        fluid.tpflash(Tl[i],p,0)
        sl[i] = fluid.fluid.smass()
        
    # construct vapour subcooling:
    Tv = np.linspace(Tsat,473,50)
    sv = np.zeros(len(Tv))
    for i  in range(len(Tv)):
        fluid.tpflash(Tv[i],p,1)
        sv[i] = fluid.fluid.smass()
        
    # point properties:
    if T < Tsat:
        region = 0
    else: 
        region = 1
    fluid.tpflash(T,p,region)
    s = fluid.fluid.smass()
    h = fluid.fluid.hmass()
    d = fluid.fluid.rhomass()
    
    # clear T-s plot:   
    document.getElementById("plt-sub-super").innerHTML = ""
    
    # plot:
    fig, ax = plt.subplots()
    ax.plot(ssat-smin,tsat,'k-',linewidth=1)
    ax.plot(sl-smin,Tl,'g-',linewidth=1)
    ax.plot(sv-smin,Tv,'g-',linewidth=1)
    ax.plot(np.array([sl[49],sv[0]])-smin,np.array([Tsat,Tsat]),'g-',linewidth=1)
    ax.plot(s-smin,T,'go',markersize=4)
    ax.set_xlabel('Entropy, s [J/(kg K)]')
    ax.set_ylabel('Temperature, T [K]')
    T_str = 'Temperature = ' + "{:.2f}".format(T) + ' K'
    P_str = 'Pressure = ' + "{:.2f}".format(p/1e5) + ' bar'
    h_str = 'Enthalpy = ' + "{:.2f}".format(h/1e3) + ' kJ/kg'
    s_str = 'Entropy = ' + "{:.2f}".format(s-smin) + ' J/(kg K)'
    d_str = 'Density = ' + "{:.2f}".format(d) + ' kg/m3'
    ax.text(70,tmax-10,T_str,color='k',size=fs)
    ax.text(70,tmax-20,P_str,color='k',size=fs)
    ax.text(70,tmax-30,h_str,color='k',size=fs)
    ax.text(70,tmax-40,s_str,color='k',size=fs)
    ax.text(70,tmax-50,d_str,color='k',size=fs)
    ax.set_xlim((0,smax-smin+200))
    ax.set_ylim((tmin,tmax))
    ax.set_position([0.175,0.125,0.80,0.85])
    fig.set_size_inches(4, 4)
    display(fig,target="plt-sub-super")
    plt.close(fig)

    
def _four_variables(*args, **kwargs):
    
    # get inputs:
    Tcond   = float(document.getElementById("inpt-tcond").value)
    subcool = float(document.getElementById("inpt-subcool").value)
    pr      = float(document.getElementById("inpt-pr").value)
    supheat = float(document.getElementById("inpt-superheat").value)
        
    # condensation conditions:
    fluid.saturation_t(Tcond,0)
    p1 = fluid.fluid.p()
    s1l = fluid.fluid.smass()
    fluid.saturation_t(Tcond,1)
    s1v = fluid.fluid.smass()
    
    # pump inlet conditions:
    T1 = Tcond - subcool
    fluid.tpflash(T1,p1,0)
    s1 = fluid.fluid.smass()
    
    # evaporation conditions:
    p2 = p1*pr
    fluid.saturation_p(p2,0)
    T2l = fluid.fluid.T()
    s2l = fluid.fluid.smass()
    p3  = p2
    fluid.saturation_p(p3,1)
    T3v = fluid.fluid.T()
    s3v = fluid.fluid.smass()
    
    # turbine inlet conditions:
    T3 = T3v + supheat
    fluid.tpflash(T3,p3,1)
    s3 = fluid.fluid.smass()
    
    # clear T-s plot:   
    document.getElementById("plt-4-variables").innerHTML = ""
    
    # update plot:
    fig, ax = plt.subplots()
    ax.plot(ssat-smin,tsat,'k-',linewidth=1)
    ax.plot(np.array([s1,s1l,s1v])-smin,np.array([T1,Tcond,Tcond]),'go-',linewidth=1,markersize=4)
    ax.plot(np.array([s2l,s3v,s3])-smin,np.array([T2l,T3v,T3]),'go-',linewidth=1,markersize=4)
    ax.set_xlabel('Entropy, s [J/(kg K)]')
    ax.set_ylabel('Temperature, T [K]')
    T1_str = 'Pump inlet temperature = ' + "{:.2f}".format(T1) + ' K'
    P1_str = 'Pump inlet pressure = ' + "{:.2f}".format(p1/1e5) + ' bar'
    T3_str = 'Expander inlet temperature = ' + "{:.2f}".format(T3) + ' K'
    P3_str = 'Expander inlet pressure = ' + "{:.2f}".format(p3/1e5) + ' bar'
    ax.text(50,tmax-10,T1_str,color='k',size=fs)
    ax.text(50,tmax-20,P1_str,color='k',size=fs)
    ax.text(50,tmax-30,T3_str,color='k',size=fs)
    ax.text(50,tmax-40,P3_str,color='k',size=fs)
    ax.set_xlim((0,smax-smin))
    ax.set_ylim((tmin,tmax))
    ax.set_position([0.175,0.125,0.80,0.85])
    fig.set_size_inches(4, 4)
    display(fig,target="plt-4-variables")
    plt.close(fig)

# run function on click:    
subcooled_superheated = create_proxy(_subcooled_superheated)
document.getElementById(
    "btn-sub-super").addEventListener("click", subcooled_superheated)


four_variables = create_proxy(_four_variables)
document.getElementById(
    "btn-4-variables").addEventListener("click", four_variables)

