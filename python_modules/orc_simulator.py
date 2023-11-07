import numpy as np
import matplotlib.pyplot as plt

################################################################################
# function to complete cycle analysis:
def simulate_ORC(flu,x,eta,hot,cld,n_hxc,ts_flg,*args):

    # inputs:
    #   x    = [ t1   , sc    , pr  , x3 , eff, tho ]
    #   eta  = [ eta_p, eta_t ]
    #   hot  = [ thi  , mh    , cph ]
    #   cld  = [ tci  , mc    , cpc ]
    #   args = [ plt_axes ]

    # simulate cycle:
    props = cycle_state_points(
        flu,x[0],x[1],x[2],x[3],hot[0],eta[0],eta[1],x[4])
    
    # size heat addition heat exchanger:
    [mdot,UAh,pph,pph_min] = addition_hxc(
        flu,props[:,2:6],hot[0],x[5],hot[1],hot[2],n_hxc[0:4])

    # size heat rejection heat exchanger:
    [tco,UAc,ppc,ppc_min] = rejection_hxc(
        flu,props[:,7:11],mdot,cld[0],cld[1],cld[2],n_hxc[4:7])

    # size recuperator:
    if ((x[4] > 0) and props[0,6] > x[0]):
        [UAr,ppr] = UA_sizing_recup(
            flu,mdot,props[:,6:8],1,props[:,1:3],0,n_hxc[7],props[:,8])
    else:
        UAr = 0
        ppr = 999

    # evaluate cycle performance:
    cycle_outpts = cycle_performance(
        flu,mdot,props,hot[0],cld[0],hot[1],hot[2])

    # plot cycle t-s diagram:
    if (ts_flg == 1):
        if (len(args) == 0):
            plot_cycle_ts(flu,props,[x[5],hot[0]],[tco,cld[0]])
        else:
            plot_cycle_ts(
                flu,props,[x[5],hot[0]],[tco,cld[0]],args[0],args[1])

    # combine UA and pinch points:
    UA = [UAh,UAc,UAr]
    pp = [pph,ppc,ppr]
    
    th = [x[5],hot[0]]
    tc = [tco,cld[0]]

    # return outputs:
    return props, cycle_outpts, UA, pp, th, tc

################################################################################
# function to analyse thermodynamic cycle:
def cycle_state_points(flu,t_cond,dt_sc,pr,exp_in,thi,etap,etat,effr):

    # initiate props array:
    props = np.zeros([5,11])

    # get low pressure saturation conditions:
    flu.saturation_t(t_cond,2)
    p1  = flu.psat
    h1_ = flu.hmass_l
    s1_ = flu.smass_l
    d1_ = flu.rhomass_l
    p4  = p1
    h4_ = flu.hmass_v 
    s4_ = flu.smass_v
    d4_ = flu.rhomass_v 
    
    # store low pressure saturation conditions:
    props[:,9] = [t_cond,p1,h1_,s1_,d1_]
    props[:,8] = [t_cond,p4,h4_,s4_,d4_]

    # set initial guess for low pressure pressure-enthalpy flash:
    satl_lp = props[:,9]
    satv_lp = props[:,8]    

    # properties at inlet to compression process:
    t1 = t_cond - dt_sc
    flu.tpflash(t1,p1,0)
    h1 = flu.fluid.hmass()
    s1 = flu.fluid.smass()
    d1 = flu.fluid.rhomass()

    # store properties at inlet to compression process: 
    props[:,0]  = [t1,p1,h1,s1,d1]
    props[:,10] = props[:,0]

    # calculate high pressure of system:
    p2 = p1*pr

    # get high pressure saturation conditions:
    if (p2 < flu.pc): # then subcritical cycle

        # saturated liquid and vapour properties:
        flu.saturation_p(p2,2)
        t2_ = flu.tsat
        h2_ = flu.hmass_l
        s2_ = flu.smass_l
        d2_ = flu.rhomass_l
        t3_ = t2_
        h3_ = flu.hmass_v 
        s3_ = flu.smass_v
        d3_ = flu.rhomass_v 

        # initial guess for high pressure pressure-enthalpy flash:
        satl_hp = [t2_,p2,h2_,s2_,d2_]
        satv_hp = [t2_,p2,h3_,s3_,d3_]

    else: # then supercritical cycle
        
        # initial guess for high pressure pressure-enthalpy flash:
        satl_hp = None
        satv_hp = None

    # run compression process:
    props[:,1] = compression(flu,props[:,0],0,p2,etap,satl_hp) 

    # properties at inlet to expansion process:
    p3 = p2
    if (p3 < flu.pc): # then subcritical cycle
        
        # expander inlet conditions:
        if (exp_in <= 1):

            # two-phase expansion:
            t3 = t3_
            h3 = h2_ + exp_in*(h3_ - h2_)
            s3 = s2_ + exp_in*(s3_ - s2_)
            #v3 = v2_ + exp_in*(v3_ - v2_)
            d3 = 1/((1/d2_) + exp_in*(1/d3_ - 1/d2_))

            # store saturation properties:
            props[:,3] = [t2_,p2,h2_,s2_,d2_]
            props[:,4] = [t3, p3,h3 ,s3 ,d3 ]

        else:

            # dry vapour expansion:
            t3 = t3_ + (exp_in - 1)*(thi - t3_)
            flu.tpflash(t3,p3,1)
            h3 = flu.fluid.hmass()
            s3 = flu.fluid.smass()
            d3 = flu.fluid.rhomass()

            # store saturation properties:
            props[:,3] = [t2_,p2,h2_,s2_,d2_]
            props[:,4] = [t3_,p3,h3_,s3_,d3_]

    else: # then supercritical cycle

        # dry vapour expansion:
        t3 = flu.tc + (exp_in - 1)*(thi - flu.tc)
        flu.tpflash(t3,p3,1)
        h3 = flu.fluid.hmass()
        s3 = flu.fluid.smass()
        d3 = flu.fluid.rhomass()

        # store saturation properties:
        props[:,3] = [t3,p3,h3,s3,d3]
        props[:,4] = props[:,3]

    # store expansion inlet properties:
    props[:,5] = [t3,p3,h3,s3,d3]

    # run expansion stage:
    props[:,6] = expansion(flu,props[:,5],1,p4,etat,satl_lp,satv_lp)

    # run recuperation:
    if (effr > 0 and props[0,6] > t_cond):
        [props[:,7],props[:,2]] = recuperation(
            flu,props[:,6],1,props[:,1],0,effr,satv_lp,satl_hp)
    else:
        props[:,2] = props[:,1]
        props[:,7] = props[:,6]    
        if (effr > 0):
            props[:,8] = props[:,7]
            
    # check if expander outlet is two-phase:
    if props[2,6] < props[2,8]:
        props[:,8] = props[:,6]

    # return array of thermodynamic properties:
    return props

################################################################################
# run compression process: 
def compression(flu,props_in,region,p_out,eta,sat):

    # unpack compressor inlet conditions:
    h_in = props_in[2]
    s_in = props_in[3]
    
    # calculate outlet enthalpy after isentropic compression:
    flu.psflash_mass(p_out,s_in,region,sat)
    h_out_s = flu.fluid.hmass()

    # calculate oulet enthalpy after real compression: 
    h_out = h_in + (h_out_s - h_in)/eta

    # calculate thermodynamic properties at outlet of compressor:
    flu.phflash_mass(p_out,h_out,region,sat)
    t_out = flu.fluid.T()
    s_out = flu.fluid.smass()
    d_out = flu.fluid.rhomass()

    # return properties:
    return np.array([t_out, p_out, h_out, s_out, d_out])

################################################################################
# run expansion process:
def expansion(flu,props_in,region,p_out,eta,satl,satv):

    # unpack expander inlet conditions:
    h_in = props_in[2]
    s_in = props_in[3]

    # unpack saturated conditions at outlet pressure:
    [tsat,dum,hsatl,ssatl,dsatl] = satl
    [tsat,dum,hsatv,ssatv,dsatv] = satv

    # calculate outlet enthalpy after isentropic expansion:
    if (s_in > ssatv):
        flu.psflash_mass(p_out,s_in,region,satv)
        h_out_s = flu.fluid.hmass()
    else:
        x_s = (s_in - ssatl)/(ssatv - ssatl)
        h_out_s = hsatl + x_s*(hsatv - hsatl)

    # calculate outlet enthalpy after real expansion:
    h_out = h_in - eta*(h_in - h_out_s)

    # calculate thermodynamic properties at outlet of expander:
    if (h_out > hsatv):
        flu.phflash_mass(p_out,h_out,region,satv)
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
# recuperator function:
def recuperation(
    flu,props_hot,region_hot,props_cld,region_cld,effr,sat_hot,sat_cld):

    # unpack fluid properties (hot side):
    th_in = props_hot[0]
    ph_in = props_hot[1]
    hh_in = props_hot[2]

    # unpack fluid properties (cold side):
    tc_in = props_cld[0]
    pc_in = props_cld[1]
    hc_in = props_cld[2]

    # calculate enthalpies for maximum heat exchange:
    flu.tpflash(tc_in,ph_in,region_hot)
    hh_out_min = flu.fluid.hmass()
    flu.tpflash(th_in,pc_in,region_cld)
    hc_out_max = flu.fluid.hmass()
    
    # adjust hot side in case of two-phase conditions:
    if hh_out_min < sat_hot[2]: hh_out_min = sat_hot[2]
    
    # calculate maximum enthalpy changes:
    dhh_max = hh_in - hh_out_min
    dhc_max = hc_out_max - hc_in

    # calculate maximum possible enthalpy change:
    dh_max = min([dhh_max, dhc_max])
    
    # compute real enthalpy change:
    dh = effr*dh_max

    # calculate outlet enthalpies:
    hh_out = hh_in - dh
    hc_out = hc_in + dh

    # calculate outlet properties:
    flu.phflash_mass(ph_in,hh_out,region_hot,sat_hot)
    th_out = flu.fluid.T()
    sh_out = flu.fluid.smass()
    dh_out = flu.fluid.rhomass()
    flu.phflash_mass(pc_in,hc_out,region_cld,sat_cld)
    tc_out = flu.fluid.T()
    sc_out = flu.fluid.smass()
    dc_out = flu.fluid.rhomass()
 
    # store properties in array:
    props_hot_out = np.array([th_out, ph_in, hh_out, sh_out, dh_out])
    props_cld_out = np.array([tc_out, pc_in, hc_out, sc_out, dc_out])

    # return hot and cold side properties:
    return props_hot_out, props_cld_out

################################################################################
# heat-addition heat exchanger function:
def addition_hxc(flu,props,thi,tho,mh,cph,n):

    # orc pressure:
    p_orc = props[1,0]

    # orc mass-flow rate:
    mdot = (mh*cph*(thi - tho))/(props[2,3] - props[2,0])

    if (p_orc < flu.pc): # then subcritical cycle:

        # determine heat-source intermediate temperatures:
        th_ev_in = thi - mdot*(props[2,3] - props[2,2])/(mh*cph)
        th_ev_ou = thi - mdot*(props[2,3] - props[2,1])/(mh*cph)

        # heat exchanger UA and pinch point calculation via discretisation:
        [UAph,pp_ph] = UA_sizing(flu,mdot,props[:,0:2],tho,mh,cph,n[0],0)
        [UAev,pp_ev] = UA_sizing(flu,mdot,props[:,1:3],th_ev_ou,mh,cph,n[1],0)
        if (th_ev_in < thi):
            [UAsh,pp_sh] = UA_sizing(
                flu,mdot,props[:,2:4],th_ev_in,mh,cph,n[2],1)
        else:
            UAsh  = 0
            pp_sh = 999

        # store values in arrays:
        UA = np.array([UAph, UAev, UAsh])
        pp = np.array([pp_ph,pp_ev,pp_sh])

        # minimum pinch point:
        pp_min = min(pp)

    else: # then supercritical cycle
        
        # size supercritical heat exchanger:
        [UA,pp] = UA_sizing(flu,mdot,props[:,0:2],tho,mh,cph,n[3],1)

        # conver to numpy array for consistency with subcritical:
        UA = np.array([UA])
        pp = np.array([pp])

        # minimum pinch point:
        pp_min = pp

    # return outputs:
    return mdot, UA, pp, pp_min

################################################################################
# heat-rejection heat exchanger function:
def rejection_hxc(flu,props,mdot,tci,mc,cpc,n):

    # orc pressure:
    p_orc = props[0,1]

    # determine heat-sink temperatures:
    tco      = tci + mdot*(props[2,0] - props[2,3])/(mc*cpc)    
    tc_co_in = tci + mdot*(props[2,1] - props[2,3])/(mc*cpc)
    tc_co_ou = tci + mdot*(props[2,2] - props[2,3])/(mc*cpc)

    # heat exchanger UA and pinch point calculation via discretisation:
    if (tc_co_in < tco):
        [UApc,pp_pc] = UA_sizing(flu,mdot,props[:,0:2],tco,mc,cpc,n[0],1)
    else:
        UApc  = 0
        pp_pc = 999
    [UAco,pp_co] = UA_sizing(flu,mdot,props[:,1:3],tc_co_in,mc,cpc,n[1],0)
    if (tc_co_ou > tci):
        [UAsc,pp_sc] = UA_sizing(flu,mdot,props[:,2:4],tc_co_ou,mc,cpc,n[2],0)
    else:
        UAsc  = 0
        pp_sc = 999

    # store values in arrays:
    UA = np.array([UApc, UAco, UAsc])
    pp = np.array([pp_pc,pp_co,pp_sc])

    # minimum pinch point:
    pp_min = min(pp)

    # return outputs:
    return tco, UA, pp, pp_min

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
        flu.phflash_mass(props_hot[1,0],h_hot[i],region_hot,sat_hot)
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
def cycle_performance(flu,mdot,props,thi,tci,mh,cph):

    # enthalpies within cycle:
    h_orc = props[2,:]

    # cycle specific peformance:
    wp = h_orc[1] - h_orc[0]
    wt = h_orc[5] - h_orc[6]
    qh = h_orc[5] - h_orc[2]
    qr = h_orc[2] - h_orc[1]
    qc = h_orc[7] - h_orc[10]
    wn = wt - wp
    
    # cycle thermal efficiency [%]:
    eta_th = 100*(wn/qh)

    # energy exchanges [kW]:
    Wp = mdot*wp/1000
    Wt = mdot*wt/1000
    Wn = mdot*wn/1000
    Qh = mdot*qh/1000
    Qc = mdot*qc/1000
    Qr = mdot*qr/1000

    # heat-source exergy [kW]:
    ex = ((mh*cph)*((thi - tci) - tci*np.log(thi/tci)))/1000

    # second-law efficiency [%]:
    eta_ex = 100*(Wn/ex)

    # return values:
    return np.array([Wn, eta_th, eta_ex, Wp, Wt, Qh, Qc, Qr, mdot])

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
        
        # construct preheating region:
        [t_ph,_,s_ph,_] = construct_isobar(flu,props[:,1],props[:,3],0,10)
        
        # construct superheating region:
        if (t[4] < t[5]):
            [t_sh,_,s_sh,_] = construct_isobar(flu,props[:,4],props[:,5],1,5)
        else:
            t_sh = t[5]
            s_sh = s[5]

        # combine: between t_hp and s_hp 
        t_hp = np.concatenate((t_ph,t_sh),axis=None)
        s_hp = np.concatenate((s_ph,s_sh),axis=None)

    else: # supercritical cycle
        [t_hp,_,s_hp,_] = construct_isobar(flu,props[:,1],props[:,5],2,50)
 
    # construct precooling region:
    if (t[6] > t[8]):
        [t_pc,_,s_pc,_] = construct_isobar(flu,props[:,6],props[:,8],1,5)
    else:
        t_pc = t[6]
        s_pc = s[6]

    # construct subcooling region:
    if (t[9] > t[10]): 
        [t_sc,_,s_sc,_] = construct_isobar(flu,props[:,9],props[:,10],0,3) 
    else:
        t_sc = t[10]
        s_sc = s[10]

    # combine isoarbars:
    t_ts = np.concatenate((t_hp,t_pc,t_sc),axis=None)
    s_ts = np.concatenate((s_hp,s_pc,s_sc),axis=None)

    # determine axis limites and scaling:
    s_cyc_min = s[0]
    s_cyc_max = s[5]
    t_cyc_min = min((min(tc),t[0]))
    t_cyc_max = max((max(th),t[5]))   
    ds = s_cyc_max - s_cyc_min
    dt = t_cyc_max - t_cyc_min
    smin = s_cyc_min - 0.1*ds
    smax = s_cyc_max + 0.1*ds
    tmin = t_cyc_min - 0.1*dt
    tmax = t_cyc_max + 0.1*dt
    
    fig, ax = plt.subplots()

    # plot ts diagram:
    ax.plot(ssat,tsat,'k-',linewidth=1)
    ax.fill(s_ts,t_ts,facecolor='green', alpha=0.5, linewidth=1)
    ax.plot(s_ts,t_ts,'g-',linewidth=1)
    ax.plot(s,t,'go',markersize=3)
    ax.plot([s[2],s[5]],th,'ro-',markersize=3,linewidth=1)
    ax.plot([s[7],s[10]],tc,'bo-',markersize=3,linewidth=1)
    ax.set_xlabel('Entropy, s [J/(kg K)]')
    ax.set_ylabel('Temperature, T [K]')
    ax.set_xlim((smin,smax))
    ax.set_ylim((tmin,tmax))
    
    # range of axis:
    ds = smax-smin
    dt = tmax-tmin
        
    # label cycle points:
    txt_index = np.array([  0,  1 , 2 , 5 ,  6 ,  7  ]) 
    txt_labs  = [ '1','2','2r','3','4','4r' ]
    if s[2] - s[1] < 1:
        txt_labs[2] = ''
        txt_labs[5] = ''    
    s_shift   = np.array([-0.02, -0.02, -0.03, 0.01, 0.01, 0.01])*ds
    t_shift   = np.array([-0.05,  0.01,  0.01, 0.00, 0.00, -0.02])*dt    
    for i in range(len(txt_index)):
        j = txt_index[i]
        ax.text(s[j]+s_shift[i],t[j]+t_shift[i],txt_labs[i],color='k',size=10)

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
    h_cyc_min = h[0]
    h_cyc_max = h[5]
    p_cyc_min = p[0]
    p_cyc_max = p[5]  
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
    txt_index = np.array([  0,  1 , 2 , 5 ,  6 ,  7  ]) 
    txt_labs  = [ '1','2','2r','3','4','4r' ]
    h_shift   = np.array([-0.03, -0.03, -0.00, 0.00, 0.00, 0.00])*dh
    p_shift   = np.array([ 0.02, 0.02,  0.02, 0.02, -0.05, -0.05])*dp    
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
    
    # create expansion:
    p_exp = np.linspace(p[5],p[6],50)
    grad  = (t[6]-t[5])/(p[6]-p[5])
    incp  = t[5] - grad*p[5]
    t_exp = grad*p_exp + incp
    v_exp = np.zeros(len(p_exp))
    v_exp[0] = v[5]
    for i in range(1,len(p_exp)-1,1):
        flu.tpflash(t_exp[i],p_exp[i],1)
        v_exp[i] = 1/flu.fluid.rhomass()
    v_exp[len(p_exp)-1] = v[6]  
    
    # combine isoarbars:
    p_ts = np.concatenate((p[0:6],p_exp,p[6:10]),axis=None)
    v_ts = np.concatenate((v[0:6],v_exp,v[6:10]),axis=None)
    
    # determine axis limites and scaling:
    v_cyc_min = v[1]
    v_cyc_max = v[6]
    p_cyc_min = p[0]
    p_cyc_max = p[5]  
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
    v_shift   = np.array([-0.03,-0.03, 0.01, 0.01])*dv
    p_shift   = np.array([-0.04, 0.01, 0.01, 0.01])*dp    
    for i in range(len(txt_index)):
        j = txt_index[i]
        ax.text((v[j]+v_shift[i]),(p[j]+p_shift[i])/1e5,txt_labs[i],color='k',size=10)

    return fig, ax

################################################################################
def construct_isobar(flu,props_in,props_out,region,n):

    # initialise arrays:
    t = np.linspace(props_in[0],props_out[0],n)
    s = np.zeros(n)
    h = np.zeros(n)
    d = np.zeros(n)
    
    # set inlet and outlet properties:
    t[0]   = props_in[0]
    h[0]   = props_in[2]
    s[0]   = props_in[3]
    d[0]   = props_in[4]
    t[n-1] = props_out[0]
    h[n-1] = props_out[2]
    s[n-1] = props_out[3]
    d[n-1] = props_out[4]
    
    # construct intermediate points:
    for i in range(1,n-1,1):
        flu.tpflash(t[i],props_in[1],region)
        h[i] = flu.fluid.hmass()
        s[i] = flu.fluid.smass()
        d[i] = flu.fluid.rhomass()

    # return properties:
    return t, h, s, d

# =========================================================================
# END OF FILE
# =========================================================================
