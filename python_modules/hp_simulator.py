import numpy as np
import matplotlib.pyplot as plt

################################################################################
def simulate_HP(flu,x,eta,hot,cld,n_hxc,ts_flg,*args):
    
    # inputs:
    #   x    = [ t1   , sh    , pr  , sc , dt_ihe , tco ]
    #   eta  = [ eta_c ]
    #   hot  = [ thi  , mh    , cph ]
    #   cld  = [ tci  , mc    , cpc ]
    #   args = [ plt_axes ]
    
    # simulate cycle:
    props = HP_state_points(flu,x[0],x[1],x[2],x[3],x[4],eta[0])
       
    # size heat addition heat exchanger:
    [mdot,UAc,ppc,ppc_min] = addition_hxc(
        flu,props[:,6:9],cld[0],x[5],cld[1],cld[2],n_hxc[0:2])

    # size heat rejection heat exchanger:
    [tho,UAh,pph,pph_min] = rejection_hxc(
        flu,props[:,1:5],mdot,hot[0],hot[1],hot[2],n_hxc[3:6])

    # size recuperator:
    if (x[4] > 0):
        [UAr,ppr] = UA_sizing_recup(
            flu,mdot,props[:,4:6],0,props[:,8:10],1,n_hxc[6],props[:,3])
    else:
        UAr = 0
        ppr = 999

    # evaluate cycle performance:
    cycle_outpts = cycle_performance(flu,mdot,props)

    # temperature profiles of source and sink:
    tc = [cld[0],x[5]]
    th = [hot[0],tho]

    # plot cycle t-s diagram:
    if (ts_flg == 1):
        if (len(args) == 0):
            plot_cycle_ts(flu,props,th,tc)
        else:
            plot_cycle_ts(
                flu,props,th,tc,args[0],args[1])

    # combine UA and pinch points:
    UA = [UAh,UAc,UAr]
    pp = [pph,ppc,ppr]

    # return outputs:
    return props, cycle_outpts, UA, pp, th, tc
    
################################################################################
def HP_state_points(flu,t_evap,dt_sh,pr,dt_sc,dt_ihe,etac):

    # initiate props array:
    props = np.zeros([5,10])

    # get low pressure saturation conditions:
    flu.saturation_t(t_evap,2)
    p1  = flu.psat
    h1_ = flu.hmass_v
    s1_ = flu.smass_v
    d1_ = flu.rhomass_v
    p4  = p1
    h4_ = flu.hmass_l
    s4_ = flu.smass_l
    d4_ = flu.rhomass_l
    
    # store low pressure saturation conditions:
    props[:,7] = [t_evap,p1,h1_,s1_,d1_]

    # set initial guess for low pressure pressure-enthalpy flash:
    satl_lp = [t_evap,p4,h4_,s4_,d4_]
    satv_lp = props[:,7]    

    # properties at inlet to compression process:
    t1 = t_evap + dt_sh + dt_ihe
    flu.tpflash(t1,p1,1)
    h1 = flu.fluid.hmass()
    s1 = flu.fluid.smass()
    d1 = flu.fluid.rhomass()

    # store properties at inlet to compression process: 
    props[:,0] = [t1,p1,h1,s1,d1]
    props[:,9] = props[:,0]

    # calculate high pressure of system:
    p2 = p1*pr

    # get high pressure saturation conditions:
    if (p2 < flu.pc): # then subcritical cycle

        # saturated liquid and vapour properties:
        flu.saturation_p(p2,2)
        t2_ = flu.tsat
        h2_ = flu.hmass_v
        s2_ = flu.smass_v
        d2_ = flu.rhomass_v
        t3_ = t2_
        h3_ = flu.hmass_l 
        s3_ = flu.smass_l
        d3_ = flu.rhomass_l 

        # initial guess for high pressure pressure-enthalpy flash:
        satv_hp = [t2_,p2,h2_,s2_,d2_]
        satl_hp = [t3_,p2,h3_,s3_,d3_]

    else: # then supercritical cycle
        
        # initial guess for high pressure pressure-enthalpy flash:
        satv_hp = None
        satl_hp = None
        
    # run compression process:
    props[:,1] = compression(flu,props[:,0],1,p2,etac,satv_hp,satl_hp) 

    # properties at outlet of condenser:
    p3 = p2
    if (p3 < flu.pc): # then subcritical cycle
        
        # expander inlet conditions:
        if (dt_sc <= 0.001):

            # saturated liquid:
            t3 = t2_
            h3 = h3_
            s3 = s3_
            d3 = d3_
        
        else:
            
            # subcooled liquid:
            t3 = t3_ - dt_sc
            flu.tpflash(t3,p3,0)
            h3 = flu.fluid.hmass()
            s3 = flu.fluid.smass()
            d3 = flu.fluid.rhomass()

        # store saturation properties:
        props[:,2] = [t2_,p2,h2_,s2_,d2_]
        props[:,3] = [t3_,p3,h3_,s3_,d3_]
        props[:,4] = [t3, p3,h3 ,s3 ,d3 ]

    else: # then supercritical cycle
        
        t3 = flu.tc - dt_sc
        flu.tpflash(t3,p3,2)
        h3 = flu.fluid.hmass()
        s3 = flu.fluid.smass()
        d3 = flu.fluid.rhomass()
        
        # store properties:
        props[:,2] = props[:,1]
        props[:,3] = props[:,2] 
        props[:,4] = [t3, p3,h3 ,s3 ,d3 ]
        
    if dt_ihe > 0:
        
        # recuperator energy balance (low pressure):
        t1r = t1 - dt_ihe
        flu.tpflash(t1r,p1,1)
        h1r = flu.fluid.hmass()
        s1r = flu.fluid.smass()
        d1r = flu.fluid.rhomass()

        # recuperator energy balance (high pressure):
        h3r = h3 - (h1 - h1r)
        flu.phflash_mass(p3,h3r,0,props[:,4],t1)
        t3r = flu.fluid.T()
        s3r = flu.fluid.smass()
        d3r = flu.fluid.rhomass()

        # store recuperator outlet properties:
        props[:,5] = [t3r,p3,h3r,s3r,d3r]
        props[:,8] = [t1r,p1,h1r,s1r,d1r]
        
    else:
        props[:,5] = props[:,4]
        props[:,8] = props[:,0]

    # run expansion stage:
    props[:,6] = expansion(flu,props[:,5],0,p4,0,satl_lp,satv_lp)

    # return array of thermodynamic properties:
    return props

################################################################################
# run compression process: 
def compression(flu,props_in,region,p_out,eta,satv,satl):

    # unpack compressor inlet conditions:
    h_in = props_in[2]
    s_in = props_in[3]

    # unpack saturated conditions at outlet pressure:
    [tsat,dum,hsatl,ssatl,dsatl] = satl
    [tsat,dum,hsatv,ssatv,dsatv] = satv
    
    # calculate outlet enthalpy after isentropic compression:
    if (s_in > ssatv):
        flu.psflash_mass(p_out,s_in,region,satv,satv[0]+50)
        h_out_s = flu.fluid.hmass()
    else:
        x_s = (s_in - ssatl)/(ssatv - ssatl)
        h_out_s = hsatl + x_s*(hsatv - hsatl)

    # calculate oulet enthalpy after real compression: 
    h_out = h_in + (h_out_s - h_in)/eta

    # calculate thermodynamic properties at outlet of compressor:
    if (h_out > hsatv):
        flu.phflash_mass(p_out,h_out,region,satv,satv[0]+50)
        t_out = flu.fluid.T()
        s_out = flu.fluid.smass()
        d_out = flu.fluid.rhomass()
    else:
        x_s   = (h_out - hsatl)/(hsatv - hsatl)
        t_out = tsat
        s_out = ssatl + x_s*(ssatv - ssatl)
        #v_out = vsatl + x_s*(vsatv - vsatl)
        d_out = 1/((1/dsatl) + x_s*(1/dsatv - 1/dsatl))

    # return properties:
    return np.array([t_out, p_out, h_out, s_out, d_out])

################################################################################
# run expansion process:
def expansion(flu,props_in,region,p_out,eta,satl,satv):

    # unpack expander inlet conditions:
    T_in = props_in[0]
    h_in = props_in[2]
    s_in = props_in[3]

    # unpack saturated conditions at outlet pressure:
    [tsat,dum,hsatl,ssatl,dsatl] = satl
    [tsat,dum,hsatv,ssatv,dsatv] = satv

    # calculate outlet enthalpy after isentropic expansion:
    if (s_in > ssatv):
        flu.psflash_mass(p_out,s_in,region,satv,T_in)
        h_out_s = flu.fluid.hmass()
    else:
        x_s = (s_in - ssatl)/(ssatv - ssatl)
        h_out_s = hsatl + x_s*(hsatv - hsatl)

    # calculate outlet enthalpy after real expansion:
    h_out = h_in - eta*(h_in - h_out_s)

    # calculate thermodynamic properties at outlet of expander:
    if (h_out > hsatv):
        flu.phflash_mass(p_out,h_out,region,satv,T_in)
        t_out = flu.fluid.T()
        s_out = flu.fluid.smass()
        d_out = flu.fluid.rhomass()
    else:
        x_s   = (h_out - hsatl)/(hsatv - hsatl)
        t_out = tsat
        s_out = ssatl + x_s*(ssatv - ssatl)
        #v_out = vsatl + x_s*(vsatv - vsatl)
        d_out = 1/((1/dsatl) + x_s*(1/dsatv - 1/dsatl))
 
    # return properties:
    return np.array([t_out, p_out, h_out, s_out, d_out])

################################################################################
# heat-addition heat exchanger function:
def addition_hxc(flu,props,tci,tco,mc,cpc,n):

    # working fluid mass-flow rate:
    mdot = (mc*cpc*(tci - tco))/(props[2,2] - props[2,0])

    # determine heat-source intermediate temperature:
    tc_ev_ou = tci - mdot*(props[2,2] - props[2,1])/(mc*cpc)

    # heat exchanger UA and pinch point calculation via discretisation:
    [UAev,pp_ev] = UA_sizing(flu,mdot,props[:,0:2],tco,mc,cpc,n[0],0)
    if (tc_ev_ou < tci):
        [UAsh,pp_sh] = UA_sizing(flu,mdot,props[:,1:3],tc_ev_ou,mc,cpc,n[1],1)
    else:
        UAsh  = 0
        pp_sh = 999

    # store values in arrays:
    UA = np.array([UAev, UAsh])
    pp = np.array([pp_ev,pp_sh])
    
    # minimum pinch point:
    pp_min = min(pp)

    # return outputs:
    return mdot, UA, pp, pp_min

################################################################################
# heat-rejection heat exchanger function:
def rejection_hxc(flu,props,mdot,thi,mh,cph,n):

    # working fluid pressure:
    p_orc = props[1,0]

    # heat-sink outlet temperature:
    tho = thi + mdot*(props[2,0] - props[2,3])/(mh*cph)

    if (p_orc < flu.pc): # then subcritical cycle:

        # determine heat-source intermediate temperatures:
        th_ev_in = thi + mdot*(props[2,1] - props[2,3])/(mh*cph)
        th_ev_ou = thi + mdot*(props[2,2] - props[2,3])/(mh*cph)

        # heat exchanger UA and pinch point calculation via discretisation:
        [UApc,pp_pc] = UA_sizing(flu,mdot,props[:,0:2],tho,mh,cph,n[0],1)
        [UAco,pp_co] = UA_sizing(flu,mdot,props[:,1:3],th_ev_in,mh,cph,n[1],0)
        if (th_ev_ou > thi):
            [UAsc,pp_sc] = UA_sizing(flu,mdot,props[:,2:4],th_ev_ou,mh,cph,n[2],1)
        else:
            UAsh  = 0
            pp_sh = 999

        # store values in arrays:
        UA = np.array([UApc, UAco, UAsc])
        pp = np.array([pp_pc,pp_co,pp_sc])

        # minimum pinch point:
        pp_min = min(pp)

    else: # then supercritical cycle
        
        # size supercritical heat exchanger:
        [UA,pp] = UA_sizing(flu,mdot,props[:,2:4],tho,mh,cph,n[3],1)

        # conver to numpy array for consistency with subcritical:
        UA = np.array([UA])
        pp = np.array([pp])

        # minimum pinch point:
        pp_min = pp

    # return outputs:
    return tho, UA, pp, pp_min

################################################################################
# UA sizing function when one stream has a fixed cp:
def UA_sizing(flu,mdot,props,t_stream_out,m_stream,cp_stream,n,region):

    # ensure n is an interger:
    n = int(n)

    # create temperature array for organic fluid:
    t_orc = np.linspace(props[0,0],props[0,1],n)

    # initialise stream temperature array:
    t_str = np.zeros(n)

    # determine enthalpies of organic fluid:
    if (props[0,0] == props[0,1]): # two-phase region
        h_orc = np.linspace(props[2,0],props[2,1],n)
    else:
        h_orc      = np.zeros(n)
        h_orc[0]   = props[2,0]
        h_orc[n-1] = props[2,1]
        for i in range(1,n-1,1):
            flu.tpflash(t_orc[i],props[1,0],region)
            h_orc[i] = flu.fluid.hmass()
    
    # determine stream temperatures:
    t_str[:] = t_stream_out + mdot*(h_orc - h_orc[0])/(m_stream*cp_stream)

    # compute pinch points:
    pp = np.zeros(n)
    if (props[2,0] < props[2,1]): # working fluid increasing in temperature
        pp[:] = t_str - t_orc
    else: # working fluid decreasing in temperature
        pp[:] = t_orc - t_str

    # get maximum and minimum pinch points:
    pp_min = min(pp)
    pp_max = max(pp)

    # complete UA sizing:
    if (np.sign(pp_min) == np.sign(pp_max)):
        UA = 0
        for i in range(1,n,1):
            dtA  = t_orc[i]   - t_str[i]
            dtB  = t_orc[i-1] - t_str[i-1]
            lmtd = (dtB - dtA)/np.log(dtB/dtA)
            Q    = mdot*(h_orc[i] - h_orc[i-1]) 
            UA   = UA + abs(Q/lmtd)
    else:
        UA = 999

    # return UA and minimum pinch point:
    return UA, pp_min

################################################################################
# UA sizing function for recuperator:
def UA_sizing_recup(
    flu,mdot,props_hot,region_hot,props_cld,region_cld,n,sat_hot):

    # ensure n is an interger:
    n = int(n)

    # create temperature array for cold fluid:
    t_cld = np.linspace(props_cld[0,0],props_cld[0,1],n)

    # determine enthalpies for cold fluid:
    h_cld      = np.zeros(n)
    h_cld[0]   = props_cld[2,0]
    h_cld[n-1] = props_cld[2,1]
    for i in range(1,n-1,1):
        flu.tpflash(t_cld[i],props_cld[1,0],region_cld)
        h_cld[i] = flu.fluid.hmass()

    # create enthalpy array for hot fluid:
    h_hot    = np.zeros(n)
    h_hot[:] = props_hot[2,1] + (h_cld - h_cld[0])

    # determine temperatures for hot fluid:
    t_hot      = np.zeros(n)
    t_hot[0]   = props_hot[0,1]
    t_hot[n-1] = props_hot[0,0]
    for i in range(1,n-1,1):
        flu.phflash_mass(props_hot[1,0],h_hot[i],region_hot,sat_hot,t_hot[i-1])
        t_hot[i] = flu.fluid.hmass()

    # compute pinch points:
    pp    = np.zeros(n)
    pp[:] = abs(t_hot - t_cld)

    # complete UA sizing:
    UA = 0
    for i in range(1,n,1):
        dtA  = t_hot[i]   - t_cld[i]
        dtB  = t_hot[i-1] - t_cld[i-1]
        lmtd = (dtB - dtA)/np.log(dtB/dtA)
        Q    = mdot*(h_hot[i] - h_hot[i-1])
        UA   = UA + abs(Q/lmtd)

    # return UA and minimum pinch point:
    return UA, min(pp)

################################################################################
# evaluate cycle performance:
def cycle_performance(flu,mdot,props):

    # enthalpies within cycle:
    h_orc = props[2,:]

    # cycle specific peformance:
    wc = h_orc[1] - h_orc[0]
    qc = h_orc[8] - h_orc[6]
    qh = h_orc[1] - h_orc[4]
    qr = h_orc[9] - h_orc[8]
    
    # coefficient of performance:
    cop = qh/wc

    # energy exchanges [kW]:
    Wc = mdot*wc/1000
    Qh = mdot*qh/1000
    Qc = mdot*qc/1000
    Qr = mdot*qr/1000

    # return values:
    return np.array([cop, Wc, Qh, Qc, Qr, mdot])

################################################################################
# plot t-s diagram:
def plot_cycle_ts(flu,props,th,tc,*args):

    # construct saturation curve:
    if ((len(args) == 0) or ((len(args) == 1) and (args[0] is None))): 
        flu.saturation_curve(100)
        ssat = flu.ssat
        tsat = flu.Tsat
    else:
        df   = args[0]
        ssat = df['s_sat']
        tsat = df['T_sat']

    # state points for t-s diagram:
    t = props[0,:]
    p = props[1,:]
    h = props[2,:]   
    s = props[3,:]

    # construct high pressure isobar:
    if (p[1] < flu.pc): # then subcritical cycle
        
        # construct precooling region:
        [t_pc,s_pc] = construct_isobar(flu,props[:,1],props[:,2],1,10)
        
        # construct subcooling region:
        if (t[5] < t[3]):
            [t_sc,s_sc] = construct_isobar(flu,props[:,3],props[:,5],0,5)
        else:
            t_sc = t[5]
            s_sc = s[5]

        # combine:
        t_hp = np.concatenate((t_pc,t_sc),axis=None)
        s_hp = np.concatenate((s_pc,s_sc),axis=None)

    else: # supercritical cycle
        [t_hp,s_hp] = construct_isobar(flu,props[:,1],props[:,5],2,50)
 
    # construct superheating region:
    [t_ph,s_ph] = construct_isobar(flu,props[:,7],props[:,9],1,3) 

    # combine isoarbars:
    t_ts = np.concatenate((t[0],t_hp,t[6],t_ph),axis=None)
    s_ts = np.concatenate((s[0],s_hp,s[6],s_ph),axis=None)

    # determine axis limits and scaling:
    s_cyc_min = s[5]
    s_cyc_max = s[1]
    t_cyc_min = min((min(tc),t[6]))
    t_cyc_max = max((max(th),t[1]))   
    ds = s_cyc_max - s_cyc_min
    dt = t_cyc_max - t_cyc_min
    smin = s_cyc_min - 0.1*ds
    smax = s_cyc_max + 0.1*ds
    tmin = t_cyc_min - 0.1*dt
    tmax = t_cyc_max + 0.1*dt
    
    fig, ax = plt.subplots()

    # plot ts diagram:
    ax.plot(ssat-smin,tsat,'k-',linewidth=1)
    ax.fill(s_ts-smin,t_ts,facecolor='green', alpha=0.5, linewidth=1)
    ax.plot(s_ts-smin,t_ts,'g-',linewidth=1)
    ax.plot(s-smin,t,'go',markersize=3)
    ax.plot([s[4]-smin,s[1]-smin],th,'ro-',markersize=3,linewidth=1)
    ax.plot([s[8]-smin,s[6]-smin],tc,'bo-',markersize=3,linewidth=1)
    ax.set_xlabel('Entropy, s [J/(kg K)]')
    ax.set_ylabel('Temperature, T [K]')
    ax.set_xlim((0,smax-smin))
    ax.set_ylim((tmin,tmax))
    
    # range of axis:
    ds = smax-smin
    dt = tmax-tmin
        
    # label cycle points:
    txt_index = np.array([  0,  1 ,  5 ,  6  ]) 
    txt_labs  = [ '1','2','3','4' ]   
    s_shift   = np.array([ 0.03,  0.03, -0.05, -0.05])*ds
    t_shift   = np.array([-0.00,  0.00,  0.00,  0.00])*dt    
    for i in range(len(txt_index)):
        j = txt_index[i]
        ax.text(s[j]+s_shift[i]-smin,t[j]+t_shift[i],txt_labs[i],color='k',size=10)

    return fig, ax

################################################################################
# plot p-h diagram:
def plot_cycle_ph(flu,props,*args):

    # construct saturation curve:
    if ((len(args) == 0) or ((len(args) == 1) and (args[0] is None))): 
        flu.saturation_curve(100)
        psat = flu.psat
        hsat = flu.hsat
    else:
        df   = args[0]
        psat = df['p_sat']
        hsat = df['h_sat']

    # state points for p-h diagram:
    t = props[0,:]
    p = props[1,:]
    h = props[2,:]   
    s = props[3,:]

    # determine axis limites and scaling:
    h_cyc_min = h[6]
    h_cyc_max = h[1]
    p_cyc_min = p[0]
    p_cyc_max = p[1]  
    dh = h_cyc_max - h_cyc_min
    dp = p_cyc_max - p_cyc_min
    hmin = h_cyc_min - 0.1*dh
    hmax = h_cyc_max + 0.1*dh
    pmin = p_cyc_min - 0.1*dp
    pmax = p_cyc_max + 0.1*dp

    fig, ax = plt.subplots()

    # plot ph diagram:
    ax.plot(hsat/1e3,psat/1e5,'k-',linewidth=1)
    ax.fill(h/1e3,p/1e5,facecolor='green', alpha=0.5, linewidth=1)
    ax.plot(h/1e3,p/1e5,'g-',linewidth=1)
    ax.plot(h/1e3,p/1e5,'go',markersize=3)
    ax.set_xlabel('Enthalpy, h [kJ/kg]')
    ax.set_ylabel('Pressure, P [bar]')
    ax.set_xlim((hmin/1e3,hmax/1e3))
    ax.set_ylim((pmin/1e5,pmax/1e5))
    
    # range of axis:
    dp = pmax-pmin
    dh = hmax-hmin
    
    # label cycle points:
    txt_index = np.array([  0,  1  , 5 ,  6  ]) 
    txt_labs  = [ '1','2','3','4' ]
    h_shift   = np.array([ 0.03, 0.03, -0.05, -0.05])*dh
    p_shift   = np.array([ 0.00, 0.00, -0.00, -0.00])*dp    
    for i in range(len(txt_index)):
        j = txt_index[i]
        ax.text((h[j]+h_shift[i])/1e3,(p[j]+p_shift[i])/1e5,txt_labs[i],color='k',size=10)

    return fig, ax

################################################################################
# plot p-v diagram:
def plot_cycle_pv(flu,props,*args):

    # construct saturation curve:
    if ((len(args) == 0) or ((len(args) == 1) and (args[0] is None))): 
        flu.saturation_curve(100)
        psat = flu.psat
        vsat = 1/flu.dsat
    else:
        df   = args[0]
        psat = df['p_sat']
        dsat = df['d_sat']
        vsat = 1/dsat

    # state points for p-v diagram:
    t = props[0,:]
    p = props[1,:]
    h = props[2,:]   
    s = props[3,:]
    v = 1/props[4,:]
    
    # create compression:
    p_com = np.linspace(p[0],p[1],50)
    grad  = (t[1]-t[0])/(p[1]-p[0])
    incp  = t[0] - grad*p[0]
    t_com = grad*p_com + incp
    v_com = np.zeros(len(p_com))
    v_com[0] = v[0]
    for i in range(1,len(p_com)-1,1):
        flu.tpflash(t_com[i],p_com[i],1)
        v_com[i] = 1/flu.fluid.rhomass()
    v_com[len(p_com)-1] = v[1]
    
    # create expansion:
    p_exp = np.linspace(p[5],p[6],20)
    h_exp = h[5]
    v_exp = np.zeros(len(p_exp))
    v_exp[0] = v[5]
    for i in range(1,len(p_exp)-1):
        flu.saturation_p(p_exp[i],2)
        tsat = flu.tsat
        hl = flu.hmass_l
        vl = 1/flu.rhomass_l
        hv = flu.hmass_v
        vv = 1/flu.rhomass_v
        if h_exp < hl:
            flu.phflash_mass(p_exp[i],h_exp,0,[tsat,p_exp[i],hl])
            v_exp[i] = 1/flu.fluid.rhomass()
        else:
            q = (h_exp - hl)/(hv - hl)
            v_exp[i] = vl + q*(vv - vl)
    v_exp[len(p_exp)-1] = v[6]
    
    # combine isoarbars:
    p_ts = np.concatenate((p_com,p[1:6],p_exp,p[6:9],p[0]),axis=None)
    v_ts = np.concatenate((v_com,v[1:6],v_exp,v[6:9],v[0]),axis=None)
    
    # determine axis limites and scaling:
    v_cyc_min = v[5]
    v_cyc_max = v[0]
    p_cyc_min = p[0]
    p_cyc_max = p[1]  
    dv = v_cyc_max - v_cyc_min
    dp = p_cyc_max - p_cyc_min
    vmin = v_cyc_min - 0.1*dv
    vmax = v_cyc_max + 0.1*dv
    pmin = p_cyc_min - 0.1*dp
    pmax = p_cyc_max + 0.1*dp

    fig, ax = plt.subplots()
    
    # plot pv diagram:
    ax.plot(vsat,psat/1e5,'k-',linewidth=1)
    ax.fill(v_ts,p_ts/1e5,facecolor='green', alpha=0.5, linewidth=1)
    ax.plot(v_ts,p_ts/1e5,'g-',linewidth=1)
    ax.plot(v,p/1e5,'go',markersize=3)
    ax.set_xlabel('Specific volume, v [$\mathregular{m^{3}}$/kg]')
    ax.set_ylabel('Pressure, P [bar]')
    ax.set_xlim((vmin,vmax))
    ax.set_ylim((pmin/1e5,pmax/1e5))
    
    # range of axis:
    dp = pmax-pmin
    dv = vmax-vmin
    
    # label cycle points:
    txt_index = np.array([  0,  1 , 5 ,  6  ]) 
    txt_labs  = [ '1','2','3','4']
    v_shift   = np.array([ 0.03, 0.03, -0.05, -0.05])*dv
    p_shift   = np.array([-0.00, 0.00, 0.00, -0.02])*dp    
    for i in range(len(txt_index)):
        j = txt_index[i]
        ax.text((v[j]+v_shift[i]),(p[j]+p_shift[i])/1e5,txt_labs[i],color='k',size=10)

    return fig, ax
    
################################################################################
def construct_isobar(flu,props_in,props_out,region,n):

    # initialise arrayss:
    t = np.linspace(props_in[0],props_out[0],n)
    s = np.zeros(n)
    
    # set inlet and outlet properties:
    t[0]   = props_in[0]
    s[0]   = props_in[3]
    t[n-1] = props_out[0]
    s[n-1] = props_out[3]
    
    # construct intermediate points:
    for i in range(1,n-1,1):
        flu.tpflash(t[i],props_in[1],region)
        s[i] = flu.fluid.smass()

    # return properties:
    return t, s

# =========================================================================
# END OF FILE
# =========================================================================
