import thermo_props
import matplotlib.pyplot as plt
import csv
import pandas as pd

# fluid properties:
fluids = [ "R1233zd(E)", "n-pentane",  "MM" ]
Tc     = [  439.6,        469.7,        518.7     ]
Pc     = [  3.6237e6,     3.3675e6,     1.93113e6 ]
om     = [  0.3025,       0.2510,       0.4180    ]
wm     = [  130.4962,     72.1488,      162.3768  ]
cp1    = [  33.3490,      12.9055,      64.6367   ]
cp2    = [  0.2823,       0.3906,       0.6695    ]
cp3    = [ -0.1523e-3,   -0.1036e-3,   -0.2895e-3 ]  
pmin   = [  47.7939e3,    24.2887e3,    1.3428e3  ]
pmax   = [  0.98,         0.98,         0.97      ]

for i in range(len(fluids)):
    
    # setup fluid:
    fluid = thermo_props.pr_fluid(
        fluids[i],Tc[i],Pc[i],om[i],[cp1[i],cp2[i],cp3[i]],300,0.01,wm[i])
    
    # construct saturation curves:
    fluid.saturation_curve(500,pmin[i],pmax[i])
   
    # write saturation curve to file:
    file_name = fluids[i] + ".csv"
    with open(file_name, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["s_sat","T_sat","p_sat","h_sat","d_sat"])
        for i in range(len(fluid.ssat)): writer.writerow([fluid.ssat[i],fluid.Tsat[i],fluid.psat[i],fluid.hsat[i],fluid.dsat[i]])
