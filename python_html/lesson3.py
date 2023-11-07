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
display(fig,target="plt-ideal-cycle")
display(fig,target="plt-real-cycle")
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

def _ideal_cycle(*args, **kwargs):
    
    # get inputs:
    Tcond   = float(document.getElementById("inpt-tcond-ideal").value)
    subcool = float(document.getElementById("inpt-subcool-ideal").value)
    pr      = float(document.getElementById("inpt-pr-ideal").value)
    supheat = float(document.getElementById("inpt-superheat-ideal").value)
        
    # condensation conditions:
    fluid.saturation_t(Tcond,0)
    p1 = fluid.fluid.p()
    T1l = fluid.fluid.T()
    s1l = fluid.fluid.smass()
    fluid.saturation_t(Tcond,1)
    T4v = fluid.fluid.T()
    s4v = fluid.fluid.smass()
    
    # pump inlet conditions:
    T1 = Tcond - subcool
    fluid.tpflash(T1,p1,0)
    s1 = fluid.fluid.smass()
    h1 = fluid.fluid.hmass()
    
    # isentropic pump outlet conditions:
    p2 = p1*pr
    fluid.psflash_mass(p2,s1,0)
    T2 = fluid.fluid.T()
    s2 = s1
    h2 = fluid.fluid.hmass()
    
    # evaporation conditions:
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
    h3 = fluid.fluid.hmass()
    
    # isentropic turbine outlet conditions:
    p4 = p1
    fluid.psflash_mass(p4,s3,1)
    s4 = fluid.fluid.smass()
    h4 = fluid.fluid.hmass()
    T4 = fluid.fluid.T()
    
    # efficiency:
    wp = h2 - h1
    wt = h3 - h4
    wn = wt - wp
    qh = h3 - h2
    eta = 100*(wn/qh)
    
    # clear T-s plot:   
    document.getElementById("plt-ideal-cycle").innerHTML = ""
    
    # cycle points:
    s_ts = np.array([s1,s2,s2l,s3v,s3,s4,s4v,s1l,s1])
    t_ts = np.array([T1,T2,T2l,T3v,T3,T4,T4v,T1l,T1])
    
    # strings for plot:
    eta_str = 'Thermal efficiency = ' + "{:.2f}".format(eta) + ' %'
    wn_str  = 'Net specific work = ' + "{:.2f}".format(wn/1000) + ' kJ/kg'
    qh_str  = 'Specific thermal input = ' + "{:.2f}".format(qh/1000) + ' kJ/kg'

    # update plot:
    fig, ax = plt.subplots()
    ax.plot(ssat-smin,tsat,'k-',linewidth=1)
    ax.plot(np.array([s1,s2,s3,s4])-smin,np.array([T1,T2,T3,T4]),'go',markersize=3)
    ax.fill(s_ts-smin,t_ts,facecolor='g', alpha=0.5, linewidth=1,)
    ax.plot(s_ts-smin,t_ts,'g-',linewidth=1)
    ax.plot(np.array([s1,s3])-smin,np.array([423,423]),'ro-',linewidth=1,markersize=3)
    ax.plot(np.array([s1,s3])-smin,np.array([293,293]),'bo-',linewidth=1,markersize=3)
    ax.set_xlabel('Entropy, s [J/(kg K)]')
    ax.set_ylabel('Temperature, T [K]')
    ax.set_xlim((0,smax-smin))
    ax.set_ylim((tmin,tmax))
    ax.text(50,tmax-8,eta_str,color='k',size=fs)
    ax.text(50,tmax-16,wn_str,color='k',size=fs)
    ax.text(50,tmax-24,qh_str,color='k',size=fs)
    ax.set_position([0.175,0.125,0.80,0.85])
    fig.set_size_inches(4, 4)
    display(fig,target="plt-ideal-cycle")
    plt.close(fig)

    
def _real_cycle(*args, **kwargs):
    
    # get inputs:
    etap = float(document.getElementById("inpt-etap").value)/100
    etae = float(document.getElementById("inpt-etae").value)/100
    
    # fixed inputs:
    Tcond   = 313
    subcool = 2
    pr      = 10
    supheat = 10
        
    # condensation conditions:
    fluid.saturation_t(Tcond,0)
    p1 = fluid.fluid.p()
    T1l = fluid.fluid.T()
    s1l = fluid.fluid.smass()
    fluid.saturation_t(Tcond,1)
    T4v = fluid.fluid.T()
    s4v = fluid.fluid.smass()
    
    # pump inlet conditions:
    T1 = Tcond - subcool
    fluid.tpflash(T1,p1,0)
    s1 = fluid.fluid.smass()
    h1 = fluid.fluid.hmass()
    
    # isentropic pump outlet conditions:
    p2 = p1*pr
    fluid.psflash_mass(p2,s1,0)
    h2s = fluid.fluid.hmass()
    
    # pump outlet conditions:
    h2 = h1 + (h2s - h1)/etap
    fluid.phflash_mass(p2,h2,0)
    T2 = fluid.fluid.T()
    s2 = fluid.fluid.smass()
    
    # evaporation conditions:
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
    h3 = fluid.fluid.hmass()
    
    # isentropic turbine outlet conditions:
    p4 = p1
    fluid.psflash_mass(p4,s3,1)
    h4s = fluid.fluid.hmass()
    
    # turbine outlet conditions:
    h4 = h3 - etae*(h3 - h4s)
    fluid.phflash_mass(p4,h4,1)
    s4 = fluid.fluid.smass()
    T4 = fluid.fluid.T()
    
    # efficiency:
    wp = h2 - h1
    wt = h3 - h4
    wn = wt - wp
    qh = h3 - h2
    eta = 100*(wn/qh)
    
    # clear T-s plot:   
    document.getElementById("plt-real-cycle").innerHTML = ""
    
    # cycle points:
    s_ts = np.array([s1,s2,s2l,s3v,s3,s4,s4v,s1l,s1])
    t_ts = np.array([T1,T2,T2l,T3v,T3,T4,T4v,T1l,T1])
    
    # strings for plot:
    eta_str = 'Thermal efficiency = ' + "{:.2f}".format(eta) + ' %'
    wn_str  = 'Net specific work = ' + "{:.2f}".format(wn/1000) + ' kJ/kg'
    qh_str  = 'Specific thermal input = ' + "{:.2f}".format(qh/1000) + ' kJ/kg'

    # update plot:
    fig, ax = plt.subplots()
    ax.plot(ssat-smin,tsat,'k-',linewidth=1)
    ax.plot(np.array([s1,s2,s3,s4])-smin,np.array([T1,T2,T3,T4]),'go',markersize=3)
    ax.fill(s_ts-smin,t_ts,facecolor='g', alpha=0.5, linewidth=1)
    ax.plot(s_ts-smin,t_ts,'g-',linewidth=1)
    ax.plot(np.array([s2,s3])-smin,np.array([423,423]),'ro-',linewidth=1,markersize=3)
    ax.plot(np.array([s1,s4])-smin,np.array([293,293]),'bo-',linewidth=1,markersize=3)
    ax.set_xlabel('Entropy, s [J/(kg K)]')
    ax.set_ylabel('Temperature, T [K]')
    ax.set_xlim((0,smax-smin))
    ax.set_ylim((tmin,tmax))
    ax.text(50,tmax-8,eta_str,color='k',size=fs)
    ax.text(50,tmax-16,wn_str,color='k',size=fs)
    ax.text(50,tmax-24,qh_str,color='k',size=fs)
    
    # inset axes showing pump compression process:
    axins = ax.inset_axes([0.051, 0.45, 0.3, 0.3])
    x1, x2, y1, y2 = s1-20-smin, s1+20-smin, T1-2, T2+2
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    axins.set_xticklabels([])
    axins.set_yticklabels([])
    axins.plot(ssat-smin,tsat,'k-',linewidth=1)
    axins.plot(np.array([s1,s2,s3,s4])-smin,np.array([T1,T2,T3,T4]),'go')
    axins.fill(s_ts-smin,t_ts,facecolor='g', alpha=0.5, linewidth=1)
    axins.plot(s_ts-smin,t_ts,'g-',linewidth=1)
    ax.indicate_inset_zoom(axins, edgecolor="black")

    # show plot:
    ax.set_position([0.175,0.125,0.80,0.85])
    fig.set_size_inches(4, 4)
    display(fig,target="plt-real-cycle")
    plt.close(fig)
    
    
# run function on click:    
ideal_cycle = create_proxy(_ideal_cycle)
document.getElementById(
    "btn-ideal-cycle").addEventListener("click", ideal_cycle)

real_cycle = create_proxy(_real_cycle)
document.getElementById(
    "btn-real-cycle").addEventListener("click", real_cycle)

