from pyscript import document
from pyscript import display
from js import console
from pyodide.ffi import create_proxy 
import csv
import matplotlib.pyplot as plt
import numpy as np
import thermo_props

# initialise figures:
fig, ax = plt.subplots(1,2)
fig.set_size_inches(8, 4)
ax[0].set_position([0.09,0.125,0.39,0.85])
ax[1].set_position([0.59,0.125,0.39,0.85])
display(fig,target="plt-recup-eff")
#plt.close(fig)

# ([0.175,0.125,0.80,0.85])

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

# update saturation curve:
ssat[:] = [x - smin for x in ssat]

   
def _recuperated_cycle(*args, **kwargs):
    
    # get inputs:
    eff_r = float(document.getElementById("inpt-recup-eff").value)/100
    
    # fixed inputs:
    Tcond   = 313
    subcool = 0
    pr      = 10
    supheat = 10
    etap    = 0.7
    etae    = 0.8
        
    # condensation conditions:
    fluid.saturation_t(Tcond,0)
    p1 = fluid.fluid.p()
    T1l = fluid.fluid.T()
    s1l = fluid.fluid.smass()
    fluid.saturation_t(Tcond,1)
    T4v = fluid.fluid.T()
    s4v = fluid.fluid.smass()
    h4v = fluid.fluid.hmass()
    
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
    
    # recuperator:
    fluid.tpflash(T4,p2,0)
    h2max = fluid.fluid.hmass()
    fluid.tpflash(T2,p4,1)
    h4min = fluid.fluid.hmass()
    dh_max = min([h4-h4min,h2max-h2])
    h4r = h4 - eff_r*dh_max
    h2r = h2 + eff_r*dh_max
    if eff_r < 1e-3:
        T4r = T4
        s4r = s4
        T2r = T2
        s2r = s2
    else:
        fluid.phflash_mass(p2,h2r,0)
        T2r = fluid.fluid.T()
        s2r = fluid.fluid.smass()
        fluid.phflash_mass(p4,h4r,1)
        T4r = fluid.fluid.T()
        s4r = fluid.fluid.smass()

    # efficiency:
    wp = h2 - h1
    wt = h3 - h4
    wn = wt - wp
    qh = h3 - h2r
    eta = 100*(wn/qh)
    
    # clear T-s plot:   
    document.getElementById("plt-recup-eff").innerHTML = ""
    
    # cycle points:
    s_ts = np.array([s1,s2,s2r,s2l,s3v,s3,s4,s4r,s4v,s1l,s1])
    t_ts = np.array([T1,T2,T2r,T2l,T3v,T3,T4,T4r,T4v,T1l,T1])
    
    # recuperator points:
    s_rec = np.array([s2,s2r,s4,s4r,s2])
    t_rec = np.array([T2,T2r,T4,T4r,T2])
    
    # strings for plot:
    eta_str = 'Thermal efficiency = ' + "{:.2f}".format(eta) + ' %'
    wn_str  = 'Net specific work = ' + "{:.2f}".format(wn/1000) + ' kJ/kg'
    qh_str  = 'Specific thermal input = ' + "{:.2f}".format(qh/1000) + ' kJ/kg'

    # update plot:
    fig, ax = plt.subplots(1,2)
    ax[0].plot(ssat,tsat,'k-',linewidth=1)
    ax[0].plot(np.array([s1,s2,s2r,s3,s4,s4r])-smin,np.array(
        [T1,T2,T2r,T3,T4,T4r]),'go',markersize=3)
    ax[0].fill(s_ts-smin,t_ts,facecolor='g', alpha=0.2, linewidth=1)
    ax[0].fill(s_rec-smin,t_rec,facecolor='g', alpha=0.5, linewidth=2)
    ax[0].plot(s_ts-smin,t_ts,'g-',linewidth=1)
    ax[0].plot(np.array([s2,s2r])-smin,np.array([T2,T2r]),'g-',linewidth=2)
    ax[0].plot(np.array([s4,s4r])-smin,np.array([T4,T4r]),'g-',linewidth=2)
    ax[0].plot(np.array([s2r,s3])-smin,np.array([423,423]),'ro-',markersize=3,linewidth=1)
    ax[0].plot(np.array([s1,s4r])-smin,np.array([293,293]),'bo-',markersize=3,linewidth=1)
    ax[0].set_xlabel('Entropy, s [J/(kg K)]')
    ax[0].set_ylabel('Temperature, T [K]')
    ax[0].set_xlim((0,smax-smin))
    ax[0].set_ylim((tmin,tmax))
    ax[0].text(50,tmax-8,eta_str,color='k',size=fs)
    ax[0].text(50,tmax-16,wn_str,color='k',size=fs)
    ax[0].text(50,tmax-24,qh_str,color='k',size=fs)
    
    ax[1].plot(np.array([0,1]),np.array([T4r,T4]),'r-',linewidth=1)
    ax[1].plot(np.array([0,1]),np.array([T2,T2r]),'b-',linewidth=1)
    ax[1].set_xlabel('Normalised heat exchanger location [-]')
    ax[1].set_ylabel('Temperature, T [K]')
    ax[1].set_xlim((0,1))
    ax[1].set_ylim((T2-10,T4+10))
    ax[1].text(0.67,T4+2,'Expander outlet',color='r',size=fs)
    ax[1].text(0.02,T2-3,'Pump outlet',color='b',size=fs)

    # show plot:
    ax[0].set_position([0.09,0.125,0.39,0.85])
    ax[1].set_position([0.59,0.125,0.39,0.85])
    display(fig,target="plt-recup-eff")
    #plt.close(fig)
    
    
# run function on click:    
recuperated_cycle = create_proxy(_recuperated_cycle)
document.getElementById(
    "btn-recup-eff").addEventListener("click", recuperated_cycle)

