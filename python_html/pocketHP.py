from js import document
from js import console
from pyodide.ffi import create_proxy
import hp_simulator
import thermo_props
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# heat exchanger discretisation:
n_hxc = [2, 2, 2, 5, 2, 2, 2, 2]
# similar to user-defined functions but w/o name
func = lambda x:0.1 ** 2 

# initialise the graph using matplotlib:
fig, ax = plt.subplots()
fig.set_size_inches(4, 4)
ax.set_position([0.175, 0.125, 0.80, 0.85])
display(fig,target="single_study_plot")
plt.tight_layout()
plt.close(fig)

def setup_fluid(fluid_name):
    
    # avilable fluids:
    fluids = [ "R1233zd(E)", "n-pentane",  "MM" ]
    Tc     = [  439.6,        469.7,        518.7     ]
    Pc     = [  3.6237e6,     3.3675e6,     1.93113e6 ]
    om     = [  0.3025,       0.2510,       0.4180    ]
    wm     = [  130.4962,     72.1488,      162.3768  ]
    cp1    = [  33.3490,      12.9055,      64.6367   ]
    cp2    = [  0.2823,       0.3906,       0.6695    ]
    cp3    = [ -0.1523e-3,   -0.1036e-3,   -0.2895e-3 ] 

    # set index and read csv file containing saturation curve:
    if fluid_name == "R1233zd(E)":
        i  = 0
        df = pd.read_csv(r'./R1233zd(E).csv')
    elif fluid_name == "n-pentane":
        i  = 1
        df = pd.read_csv(r'./n-pentane.csv')
    elif fluid_name == "MM":
        i  = 2
        df = pd.read_csv(r'./MM.csv')
        
    # setup fluid:
    fluid = thermo_props.pr_fluid(
        fluids[i],Tc[i],Pc[i],om[i],[cp1[i],cp2[i],cp3[i]],300,0.01,wm[i])
        
    # return fluid and datafile:
    return fluid, df


# ------------------------------------------------------------------------ #
def single_Cycle(*args, **kwargs):

    name = document.getElementById('working-fluid-selecter').value
    [fluid, df] = setup_fluid(name)

    # get heat-source conditions:
    cld = [0,0,0]
    cld[0] = float(document.getElementById("heat-source-temperature").value)
    cld[1] = float(document.getElementById("heat-source-mass-flow").value)
    cld[2] = float(document.getElementById("heat-source-heat-capacity").value)
    
    # get heat-sink conditions:
    hot = [0,0,0]
    hot[0] = float(document.getElementById("heat-sink-temperature").value)
    hot[1] = float(document.getElementById("heat-sink-mass-flow").value)
    hot[2] = float(document.getElementById("heat-sink-heat-capacity").value)

    # get component efficiencies: (fixed)
    eta = [0]
    eta[0] = float(document.getElementById("compressor-efficiency").value)/100
    
    # cycle variables
    x = [0, 0, 0, 0, 0, 0]
    x[0] = float(document.getElementById("evaporation-temperature").value)
    x[1] = float(document.getElementById("compressor-superheat").value)
    x[2] = float(document.getElementById("pressure-ratio").value)
    x[3] = float(document.getElementById("expansion-value-subcooling").value)
    x[4] = float(document.getElementById("internal-heat-exchange-temp").value)
    x[5] = float(document.getElementById("heat-source-outlet").value)

    [props, cycle_out, UA, pp, th, tc] = hp_simulator.simulate_HP(fluid, x, eta, hot, cld, n_hxc, 0)

    # output: hp_simulator --> cycle_performance
    # change #
    document.getElementById("coefficient-of-performance").innerHTML = "{:.3f}".format(cycle_out[0]) 
    document.getElementById("compressor-work").innerHTML = "{:.3f}".format(cycle_out[1])
    document.getElementById("heat-input").innerHTML = "{:.3f}".format(cycle_out[2])
    document.getElementById("heat-output").innerHTML = "{:.3f}".format(cycle_out[3])
    document.getElementById("heat-recuperated").innerHTML = "{:.3f}".format(cycle_out[4])
    document.getElementById("mass-flow-rate").innerHTML = "{:.3f}".format(cycle_out[5])

    # heat output
    document.getElementById("ua-hot").innerHTML = "{:.3f}".format(sum(UA[0])/1000)
    document.getElementById("ua-cold").innerHTML = "{:.3f}".format(sum(UA[1])/1000)
    document.getElementById("ua-recup").innerHTML = "{:.3f}".format(UA[2]/1000)
    document.getElementById("pp-hot").innerHTML = "{:.3f}".format(min(pp[0]))
    document.getElementById("pp-cold").innerHTML = "{:.3f}".format(min(pp[1]))
    document.getElementById("pp-recup").innerHTML = "{:.3f}".format(pp[2])
    
    # change text color if pinch point is negative:
    if min(pp[0]) < 1:
        document.getElementById("pp-hot").style.color = '#d00'
    else:
        document.getElementById("pp-hot").style.color = '#000000'
    if min(pp[1]) < 1:
        document.getElementById("pp-cold").style.color = '#d00'
    else:
        document.getElementById("pp-cold").style.color = '#000000'
    if pp[2] < 1:
        document.getElementById("pp-recup").style.color = '#d00'
    else:
        document.getElementById("pp-recup").style.color = '#000000'
    
    # reset the plot:
    document.getElementById("single_study_plot").innerHTML = ""

    # update plot:
    singleParam = document.getElementById("single-study-select").value
    if singleParam == "temperature-entropy":    
        fig, ax = hp_simulator.plot_cycle_ts(fluid,props,th,tc, df)
    elif singleParam == "pressure-enthalpy":
        fig, ax = hp_simulator.plot_cycle_ph(fluid,props,df)
    elif singleParam == "pressure-volume":
        fig, ax = hp_simulator.plot_cycle_pv(fluid,props,df)
    ax.set_position([0.175, 0.125, 0.80, 0.85])
    fig.set_size_inches(4, 4)
    display(fig, target = "single_study_plot")
    plt.close(fig)

    # print Temperature [K]: 
    document.getElementById("pl-1").innerHTML ="{:.2f}".format(props[0][0])
    document.getElementById("po-1").innerHTML ="{:.2f}".format(props[0][1])
    document.getElementById("roh-1").innerHTML ="{:.2f}".format(props[0][2])
    document.getElementById("slh-1").innerHTML ="{:.2f}".format(props[0][3])
    document.getElementById("svh-1").innerHTML ="{:.2f}".format(props[0][4])
    document.getElementById("ei-1").innerHTML ="{:.2f}".format(props[0][5])
    document.getElementById("eo-1").innerHTML ="{:.2f}".format(props[0][6])
    document.getElementById("rol-1").innerHTML ="{:.2f}".format(props[0][7])
    document.getElementById("svl-1").innerHTML ="{:.2f}".format(props[0][8])

    # print Pressure [kPa]: 
    document.getElementById("pl-2").innerHTML ="{:.2f}".format(props[1][0]/1e3)
    document.getElementById("po-2").innerHTML ="{:.2f}".format(props[1][1]/1e3)
    document.getElementById("roh-2").innerHTML ="{:.2f}".format(props[1][2]/1e3)
    document.getElementById("slh-2").innerHTML ="{:.2f}".format(props[1][3]/1e3)
    document.getElementById("svh-2").innerHTML ="{:.2f}".format(props[1][4]/1e3)
    document.getElementById("ei-2").innerHTML ="{:.2f}".format(props[1][5]/1e3)
    document.getElementById("eo-2").innerHTML ="{:.2f}".format(props[1][6]/1e3)
    document.getElementById("rol-2").innerHTML ="{:.2f}".format(props[1][7]/1e3)
    document.getElementById("svl-2").innerHTML ="{:.2f}".format(props[1][8]/1e3)

    # print Enthalpy [kJ/kg]: 
    document.getElementById("pl-3").innerHTML ="{:.2f}".format(props[2][0]/1e3)
    document.getElementById("po-3").innerHTML ="{:.2f}".format(props[2][1]/1e3)
    document.getElementById("roh-3").innerHTML ="{:.2f}".format(props[2][2]/1e3)
    document.getElementById("slh-3").innerHTML ="{:.2f}".format(props[2][3]/1e3)
    document.getElementById("svh-3").innerHTML ="{:.2f}".format(props[2][4]/1e3)
    document.getElementById("ei-3").innerHTML ="{:.2f}".format(props[2][5]/1e3)
    document.getElementById("eo-3").innerHTML ="{:.2f}".format(props[2][6]/1e3)
    document.getElementById("rol-3").innerHTML ="{:.2f}".format(props[2][7]/1e3)
    document.getElementById("svl-3").innerHTML ="{:.2f}".format(props[2][8]/1e3)

    # print Entropy [kJ/(kg K)]: 
    document.getElementById("pl-4").innerHTML ="{:.3f}".format(props[3][0]/1e3)
    document.getElementById("po-4").innerHTML ="{:.3f}".format(props[3][1]/1e3)
    document.getElementById("roh-4").innerHTML ="{:.3f}".format(props[3][2]/1e3)
    document.getElementById("slh-4").innerHTML ="{:.3f}".format(props[3][3]/1e3)
    document.getElementById("svh-4").innerHTML ="{:.3f}".format(props[3][4]/1e3)
    document.getElementById("ei-4").innerHTML ="{:.3f}".format(props[3][5]/1e3)
    document.getElementById("eo-4").innerHTML ="{:.3f}".format(props[3][6]/1e3)
    document.getElementById("rol-4").innerHTML ="{:.3f}".format(props[3][7]/1e3)
    document.getElementById("svl-4").innerHTML ="{:.3f}".format(props[3][8]/1e3)

    # print Density [kg/m^3]: 
    document.getElementById("pl-5").innerHTML ="{:.2f}".format(props[4][0])
    document.getElementById("po-5").innerHTML ="{:.2f}".format(props[4][1])
    document.getElementById("roh-5").innerHTML ="{:.2f}".format(props[4][2])
    document.getElementById("slh-5").innerHTML ="{:.2f}".format(props[4][3])
    document.getElementById("svh-5").innerHTML ="{:.2f}".format(props[4][4])
    document.getElementById("ei-5").innerHTML ="{:.2f}".format(props[4][5])
    document.getElementById("eo-5").innerHTML ="{:.2f}".format(props[4][6])
    document.getElementById("rol-5").innerHTML ="{:.2f}".format(props[4][7])
    document.getElementById("svl-5").innerHTML ="{:.2f}".format(props[4][8])
    
def paramCycle(*args, **kwargs):

    name = document.getElementById('working-fluid-selecter').value
    [fluid, df] = setup_fluid(name)

    # get heat-source conditions:
    cld = [0,0,0]
    cld[0] = float(document.getElementById("heat-source-temperature").value)
    cld[1] = float(document.getElementById("heat-source-mass-flow").value)
    cld[2] = float(document.getElementById("heat-source-heat-capacity").value)
    
    # get heat-sink conditions:
    hot = [0,0,0]
    hot[0] = float(document.getElementById("heat-sink-temperature").value)
    hot[1] = float(document.getElementById("heat-sink-mass-flow").value)
    hot[2] = float(document.getElementById("heat-sink-heat-capacity").value)

    # get component efficiencies:
    eta = [0,0]
    eta[0] = float(document.getElementById("compressor-efficiency").value)/100

    # get cycle variables: fix
    x = [0, 0, 0, 0, 0, 0]
    x[0] = float(document.getElementById("evaporation-temperature").value)
    x[1] = float(document.getElementById("compressor-superheat").value)
    x[2] = float(document.getElementById("pressure-ratio").value)
    x[3] = float(document.getElementById("expansion-value-subcooling").value)
    x[4] = float(document.getElementById("internal-heat-exchange-temp").value)
    x[5] = float(document.getElementById("heat-source-outlet").value)

    # extract parametric study:
    param_var = document.getElementById("param-study-selecter").value
    param_min = float(document.getElementById("min-Param").value)
    param_max = float(document.getElementById("max-Param").value)
    x_param   = np.linspace(param_min,param_max,10)

    # set index of parametric variable: for output
    if param_var == "Evaporation-temp":
        param_indx = 0
    elif param_var == "Compressor-superheat":
        param_indx = 1
    elif param_var == "Pressure":
        param_indx = 2
    elif param_var == "Expansion-subcooling":
        param_indx = 3
    elif param_var == "Internal-heat":
        param_indx = 4
    elif param_var == "Heat-source-outlet":
        param_indx = 5
        
    # display name on first column:
    document.getElementById("param-name-hp").innerHTML = param_var

     # run parametric study:
    for i in range(len(x_param)):
        
        # run for current set of inputs:
        x[param_indx] = x_param[i]
        [props,cycle_out,UA,pp,th,tc] = hp_simulator.simulate_HP(fluid,x,eta,hot,cld,n_hxc,0)
        console.log(props)

        # print cycle outputs:
        document.getElementById("var-"+str(i+1)).innerHTML = "{:.3f}".format(x_param[i])
        document.getElementById("cop-"+str(i+1)).innerHTML = "{:.3f}".format(cycle_out[0])
        document.getElementById("Wc-"+str(i+1)).innerHTML = "{:.3f}".format(cycle_out[1])
        document.getElementById("Qh-"+str(i+1)).innerHTML = "{:.3f}".format(cycle_out[2])
        document.getElementById("Qc-"+str(i+1)).innerHTML = "{:.3f}".format(cycle_out[3])
        document.getElementById("Qr-"+str(i+1)).innerHTML = "{:.3f}".format(cycle_out[4])
        document.getElementById("mdot-"+str(i+1)).innerHTML = "{:.3f}".format(cycle_out[5])
        document.getElementById("uah-"+str(i+1)).innerHTML = "{:.3f}".format(sum(UA[0])/1000)
        document.getElementById("uac-"+str(i+1)).innerHTML = "{:.3f}".format(sum(UA[1])/1000)
        document.getElementById("uar-"+str(i+1)).innerHTML = "{:.3f}".format(UA[2]/1000)
        document.getElementById("pph-"+str(i+1)).innerHTML = "{:.3f}".format(min(pp[0]))
        document.getElementById("ppc-"+str(i+1)).innerHTML = "{:.3f}".format(min(pp[1]))
        document.getElementById("ppr-"+str(i+1)).innerHTML = "{:.3f}".format(pp[2])
        
        # change text color if pinch point is negative:
        if min(pp[0]) < 1:
            document.getElementById("pph-"+str(i+1)).style.color = '#d00'
        else:
            document.getElementById("pph-"+str(i+1)).style.color = '#000000'
        if min(pp[1]) < 1:
            document.getElementById("ppc-"+str(i+1)).style.color = '#d00'
        else:
            document.getElementById("ppc-"+str(i+1)).style.color = '#000000'
        if pp[2] < 1:
            document.getElementById("ppr-"+str(i+1)).style.color = '#d00'
        else:
            document.getElementById("ppr-"+str(i+1)).style.color = '#000000'
            

# run function on click:    
singleCycle = create_proxy(single_Cycle)
document.getElementById("single-cycle-hp").addEventListener("click", singleCycle)

# run function on click:    
param_Cycle = create_proxy(paramCycle)
document.getElementById("param-cycle-hp").addEventListener("click", param_Cycle)