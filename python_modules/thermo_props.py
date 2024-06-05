# thermo_props.py 
# 
# python module for thermodynamic property calculation. this module 
# contains a class that can be used to predict fluid properties using an
# in-house Peng-Robinson equation of state. a separate class is available
# that interfaces with the CoolProp package using the same syntax to allow
# easy transition between different models in external scripts.
# 
#   example:
#
#     # using CoolProp get properties at given temp. & press.:
#     fluid = thermo_props.coolprop_fluid(fluid,backend)
#     fluid.tpflash(t,p)
#     fluid.fluid.hmass()
#     fluid.fluid.smass()
#
#     # using Peng-Robinson get properties at given temp. & press.:
#     fluid = thermo_props.pr_fluid(name,tc,pc,om,cp,tr,vr,wm)
#     fluid.tpflash(t,p)
#     fluid.fluid.hmass()
#     fluid.fluid.smass()
#
# =========================================================================
#
# contents:
#   class:  pr_fluid
#       class for calling a fluid defined via the internal Peng-Robinson
#       equation of state.
#             departure(self,t,p,z)
#             tvflash_mol(self,t,v)
#             tpflash(self,t,p,region)
#             phflash_mol(self,p,h,region,*args)
#             psflash_mol(self,p,s,region,*args)
#             hsflash_mol(self,h,s,region,guess1,guess2,*args)
#             psflash_mass(self,p,s_mass,region,*args)
#             phflash_mass(self,p,h_mass,region,*args)
#             hsflash_mass(self,h_mass,s_mass,region,guess1,guess2,*args)
#             saturation_t(self,t,q)
#             saturation_p(self,p,q)
#             saturation_curve(self,n,*args)
#             soundspeed_mol(self,t,v)
# 
#   class:  coolprop_fluid
#       class for calling a fluid defined in CoolProp. can be setup to
#       call using either the CoolProp or the REFPROP backend.
#             __init__(self,name,backend,*args)
#             tdflash(self,t,d,*args)
#             tpflash(self,t,p,*args) 
#             pdflash(self,p,d,*args)
#             psflash_mol(self,p,s,*args)
#             phflash_mol(self,p,h,*args)
#             hsflash_mol(self,h,s,*args)  
#             pdflash_mass(self,p,d,*args)
#             psflash_mass(self,p,s,*args)
#             phflash_mass(self,p,h,*args)
#             tsflash_mass(self,t,s,*args)
#             hsflash_mass(self,h,s,*args)
#             duflash_mass(self,d,u,*args)
#             hsflash_mass_backup(self,h,s,p_guess1,p_guess2)
#             saturation_t(self, t, q)
#             saturation_p(self, p, q)
#             saturation_curve(self, n *args)
# =========================================================================

import numpy as np
#import CoolProp.CoolProp as CoolProp

# =========================================================================
# define generic fluid class for use with Peng-Robinson fluid. Required so
#  the Peng-Robinson fluid class behaves in the same way as the CoolProp.
class fluid:
    
    ##########################################################################
    # initialisation:
    def __init__(self):
        self.temperature   = np.nan
        self.pressure      = np.nan
        self.enthalpy_mol  = np.nan
        self.entropy_mol   = np.nan
        self.volume_mol    = np.nan
        self.enthalpy_mass = np.nan
        self.entropy_mass  = np.nan
        self.density_mass  = np.nan
        self.sound_speed   = np.nan
        
    ##########################################################################
    # return functions:
    def p(self): return self.pressure
    def T(self): return self.temperature
    def hmolar(self): return self.enthalpy_mol
    def smolar(self): return self.entropy_mol
    def vmolar(self): return self.volume_mol
    def hmass(self): return self.enthalpy_mass
    def smass(self): return self.entropy_mass
    def rhomass(self): return self.density_mass
    def speed_sound(self): return self.sound_speed

# ===========================================================================
# define class for Peng-Robinson fluid:
class pr_fluid:

    ##########################################################################
    # initialisation of class:    
    def __init__(self,name,tc,pc,omega,cp,tr,vr,wm):
        
        # universal gas constant:
        r = 8.3144621         # [J/(mol K)]
        
        # determine equation paramters:
        a0 = (0.45724*(r**2)*(tc**2))/pc
        b  = 0.0778*(r*tc)/pc
        n  = 0.37464 + 1.54226*omega - 0.26992*(omega**2)
    
        # reference paraemters:
        ar = a0*((1 + n*(1 - (tr/tc)**(1/2)))**2)
        pr = (r*tr/(vr-b)) - (ar/(vr**2 + 2*b*vr - b**2))
        zr = pr*vr/(r*tr)
    
        # define fluid structure:
        self.name = name
        self.r  = r
        self.a0 = a0
        self.b  = b
        self.n  = n
        self.tc = tc
        self.pc = pc
        self.tr = tr
        self.vr = vr
        self.pr = pr
        self.cp = cp
        self.om = omega
        self.wm = wm/1000
        
        # saturation properties (single call for both phases):     
        self.tsat      = np.nan
        self.psat      = np.nan
        self.smass_l   = np.nan
        self.hmass_l   = np.nan
        self.rhomass_l = np.nan
        self.smass_v   = np.nan
        self.hmass_v   = np.nan
        self.rhomass_v = np.nan

        # initiate fluid substructure:
        self.fluid = fluid()
        
        # calulate departure function of reference parameters:
        (self.hr_dep,self.sr_dep) = self.departure(tr,pr,zr)
        
    ##########################################################################
    # departure function for enthalpy and entropy (molar):
    def departure(self,t,p,z):

        # reduced temperature and pressure:   
        tr = t/self.tc
        pr = p/self.pc

        # parameters for calculation of departure enthalpy and entropy:
        a  = ((1 + self.n*(1 - (t/self.tc)**(1/2)))**2)
        B  = 0.07780*pr/tr
        zf = (z + 2.414*B)/(z - 0.414*B)
        
        # departure enthalpy and entropy [molar basis]:
        h_dep = self.r*self.tc*(tr*(z-1) \
                  - 2.078*(1+self.n)*(a**0.5)*np.log(zf))
        s_dep = self.r*(np.log(z-B) \
                  - 2.078*self.n*((1+self.n)/(tr**0.5) - self.n) *np.log(zf))

        # return calculated departure values:
        return h_dep, s_dep

    ##########################################################################
    # temperature-volume flash:
    def tvflash_mol(self,t,v):

        # solve cubic equation of state:
        a = self.a0*((1 + self.n*(1 - (t/self.tc)**0.5))**2)
        p = (self.r*t)/(v - self.b) - a/(v**2 + 2*self.b*v - self.b**2)

        # compute compressibility factor:
        z = p*v/(self.r*t)
   
        # enthalpy and entropy departure functions:
        (h_dep,s_dep) = self.departure(t,p,z)
   
        # enthalpy [J/mol]:
        h_ideal = self.cp[0]*(t - self.tr) \
                    + (self.cp[1]/2)*(t**2 - self.tr**2) \
                    + (self.cp[2]/3)*(t**3 - self.tr**3)
        h = h_dep - self.hr_dep + h_ideal
   
        # entropy [J/(mol K)]:
        s_ideal = self.cp[0]*np.log(t/self.tr) \
                    + (self.cp[1])*(t - self.tr) \
                    + (self.cp[2]/2)*(t**2 - self.tr**2) \
                    - self.r*np.log(p/self.pr);
        s = s_dep - self.sr_dep + s_ideal
        
        # determine speed of sound:
        a = self.soundspeed_mol(t,np.real(v))

        # set to real:
        self.fluid.temperature   = t
        self.fluid.pressure      = p
        self.fluid.enthalpy_mol  = np.real(h)
        self.fluid.entropy_mol   = np.real(s)
        self.fluid.volume_mol    = np.real(v)
        self.fluid.enthalpy_mass = self.fluid.enthalpy_mol/self.wm
        self.fluid.entropy_mass  = self.fluid.entropy_mol/self.wm
        self.fluid.density_mass  = self.wm/self.fluid.volume_mol
        self.fluid.sound_speed   = a/np.sqrt(self.wm)

    ##########################################################################
    # temperature-pressure flash:
    def tpflash(self,t,p,region):
        
        # solve cubic equation of state:
        a   = self.a0*((1 + self.n*(1 - (t/self.tc)**0.5))**2)
        c0  = p
        c1  = p*self.b - self.r*t
        c2  = -(3*p*self.b**2 + 2*self.b*self.r*t - a)
        c3  = p*self.b**3 + self.r*t*self.b**2 - a*self.b
        rts = np.roots([c0, c1, c2, c3])
        
        # fluid specific volume [m^3/mol]:
        if region == 0:
            v = rts[2]
        elif region == 1:
            v = rts[0]
        elif region == 2:
            v = rts[np.isreal(rts) & (rts > self.b)]

        # compute compressibility factor:
        z = p*v/(self.r*t)
    
        # enthalpy and entropy departure functions:
        (h_dep,s_dep) = self.departure(t,p,z)
    
        # enthalpy [J/mol]:
        h_ideal = self.cp[0]*(t - self.tr) \
                    + (self.cp[1]/2)*(t**2 - self.tr**2) \
                    + (self.cp[2]/3)*(t**3 - self.tr**3)
        h = h_dep - self.hr_dep + h_ideal
    
        # entropy [J/(mol K)]:
        s_ideal = self.cp[0]*np.log(t/self.tr) \
                    + (self.cp[1])*(t - self.tr) \
                    + (self.cp[2]/2)*(t**2 - self.tr**2) \
                    - self.r*np.log(p/self.pr);
        s = s_dep - self.sr_dep + s_ideal
        
        # determine speed of sound:
        a = self.soundspeed_mol(t,np.real(v))

        # set to real:
        self.fluid.temperature   = t
        self.fluid.pressure      = p
        self.fluid.enthalpy_mol  = np.real(h)
        self.fluid.entropy_mol   = np.real(s)
        self.fluid.volume_mol    = np.real(v)
        self.fluid.enthalpy_mass = self.fluid.enthalpy_mol/self.wm
        self.fluid.entropy_mass  = self.fluid.entropy_mol/self.wm
        self.fluid.density_mass  = self.wm/self.fluid.volume_mol
        self.fluid.sound_speed   = a/np.sqrt(self.wm)

    ########################################################################
    # pressure-enthalpy flash (molar input):
    def phflash_mol(self,p,h,region,*args):
 
        if ((len(args) == 0) or ((len(args) == 1) and (args[0] is None))): 
            # no initial guess provided
            if (region == 0):
                # use saturated liquid conditions:
                self.saturation_p(p,0)
                t1 = self.fluid.T()
                h1 = self.fluid.hmolar()

                # set second temperature to a value below initial guess:
                t2 = t1 - 25

            elif (region == 1):
                # use saturated vapour conditions:
                self.saturation_p(p,1)
                t1 = self.fluid.T()
                h1 = self.fluid.hmolar()

                # set second temperature to a value above initial guess:
                t2 = t1 + 25

            elif (region == 2):
                # determine conditions at critical temperature:
                t1 = self.tc
                self.tpflash(t1,p,region)
                h1 = self.fluid.hmolar()

                # set second temperature to a value above initial guess:
                t2 = t1 + 25
        else: 
            # initial guess provided
            t1 = args[0][0]
            h1 = args[0][2] * self.wm 

	    # second temperature:
            t2 = args[1]

        # compute properties at next temperature guess:
        self.tpflash(t2,p,region)
        h2 = self.fluid.hmolar()
 
        # use secant method to determine temperature: 
        while (abs(h1 - h2) > 1):
            m  = (t2 - t1)/(h2 - h1)
            c  = t1 - m*h1
            t1 = t2
            t2 = m*h + c
            h1 = h2
            self.tpflash(t2,p,region)
            h2 = self.fluid.hmolar()    

        # compute final properties:
        t = t2;
        self.tpflash(t,p,region)
        
    ###########################################################################
    # pressure-entropy flash (molar input):
    def psflash_mol(self,p,s,region,*args):

        if ((len(args) == 0) or ((len(args) == 1) and (args[0] is None))):
            # no initial guess provided
            if (region == 0):
                # use saturated liquid conditions:
                self.saturation_p(p,0)
                t1 = self.fluid.T()
                s1 = self.fluid.smolar()

                # set second temperature to a value below initial guess:
                t2 = t1 - 25

            elif (region == 1):
                # use saturated vapour conditions:
                self.saturation_p(p,1)
                t1 = self.fluid.T()
                s1 = self.fluid.smolar()

                # set second temperature to a value above initial guess:
                t2 = t1 + 25

            elif (region == 2):
                # determine conditions at critical temperature:
                t1 = self.tc
                self.tpflash(t1,p,region)
                s1 = self.fluid.smolar()

                # set second temperature to a value above initial guess:
                t2 = t1 + 25
        else:
            # initial guess provided
            t1 = args[0][0]
            s1 = args[0][3] * self.wm

            # second temperature:
            t2 = args[1]

        # compute properties at next temperature guess:
        self.tpflash(t2,p,region)
        s2 = self.fluid.smolar() 

        # use secant method to solve for temperature: 
        while (abs(s1 - s2) > 1e-5):
            m  = (t2 - t1)/(s2 - s1)
            c  = t1 - m*s1
            t1 = t2
            t2 = m*s + c
            s1 = s2
            self.tpflash(t2,p,region)
            s2 = self.fluid.smolar()

        # compute final properties:
        t = t2
        self.tpflash(t,p,region)

    ##########################################################################
    # enthalpy-entropy flash (molar input):
    def hsflash_mol(self,h,s,region,guess1,guess2,*args):

        # initial guesses for pressures:
        pa = guess1
        pb = guess2

        # compute enthalpy for initial guesses:
        self.psflash_mol(pa,s,region)
        ha = self.fluid.hmolar()
        self.psflash_mol(pb,s,region)
        hb = self.fluid.hmolar()

        # use secant method to solve for pressure:
        while (abs(pa - pb) > 10):
            m  = (pb - pa)/(hb - ha)
            c  = pa - m*ha
            pa = pb
            ha = hb
            pb = m*h + c
            self.psflash_mol(pb,s,region)
            hb = self.fluid.hmolar()

        # compute final properties:
        p = pb
        self.psflash_mol(p,s,region)   
        
    #######################################################################
    # conversion functions from molar to mass basis:
    def psflash_mass(self,p,s_mass,region,*args):
        s_mol = s_mass*self.wm
        self.psflash_mol(p,s_mol,region,*args)
    def phflash_mass(self,p,h_mass,region,*args):
        h_mol = h_mass*self.wm
        self.phflash_mol(p,h_mol,region,*args)
    def hsflash_mass(self,h_mass,s_mass,region,guess1,guess2,*args):
        h_mol = h_mass*self.wm
        s_mol = s_mass*self.wm
        self.hsflash_mol(h_mol,s_mol,region,guess1,guess2,*args)
        
    #######################################################################
    # saturation properties at given temperature:
    def saturation_t(self,t,q):

        # determine temperature dependent factor:
        a = self.a0*((1 + self.n*(1 - (t/self.tc)**0.5))**2)

        # iteration to determine saturation pressure:
        x1 = 0
        x2 = 100
        dp = 100;
        while (abs(x1/x2 - 1) > 1e-5):
        
            x1 = x2
        
            # determine difference in fugacity at pressure 1:
            p1 = x1
 
            # find roots of cubic equation (volumes):
            c0  = p1
            c1  = p1*self.b - self.r*t
            c2  = -(3*p1*self.b**2 + 2*self.b*self.r*t - a)
            c3  = p1*self.b**3 + self.r*t*self.b**2 - a*self.b
            rts = np.roots([c0, c1, c2, c3])
        
            # determine fugacities:
            A   = a*p1/((self.r*t)**2)
            B   = self.b*p1/(self.r*t)
            z   = (p1*rts)/(self.r*t)
            phi = np.exp(z - 1 - np.log(z - B) - (A/(2*np.sqrt(2)*B)) \
                    *np.log((z + (1 + np.sqrt(2))*B)/(z + (1 - np.sqrt(2))*B)))

            # determine objective function:
            f1 = phi[0] - phi[2]

            # determine difference in fugacity at pressure 2:
            p2 = p1 + dp
        
            # find roots of cubic equation (volumes):
            c0  = p2
            c1  = p2*self.b - self.r*t
            c2  = -(3*p2*self.b**2 + 2*self.b*self.r*t - a)
            c3  = p2*self.b**3 + self.r*t*self.b**2 - a*self.b
            rts = np.roots([c0, c1, c2, c3])
        
            # determine fugacities:
            A   = a*p2/((self.r*t)**2)
            B   = self.b*p2/(self.r*t)
            z   = (p2*rts)/(self.r*t)
            phi = np.exp(z - 1 - np.log(z - B) - (A/(2*np.sqrt(2)*B)) \
                    *np.log((z + (1 + np.sqrt(2))*B)/(z + (1 - np.sqrt(2))*B)))
  
            # determine objective function:
            f2 = phi[0] - phi[2]
        
            # next guess for saturation pressure:
            df_dx = (f2 - f1)/(p2 - p1)
            x2 = x1 - f1/df_dx
       
        # final pressure and solve equation: 
        p   = x2
        c0  = p
        c1  = p*self.b - self.r*t
        c2  = -(3*p*self.b**2 + 2*self.b*self.r*t - a)
        c3  = p*self.b**3 + self.r*t*self.b**2 - a*self.b
        rts = np.roots([c0, c1, c2, c3])
           
        # specific volumes:
        vl = rts[2]
        vv = rts[0]
 
        # ideal terms:
        h_ideal = self.cp[0]*(t - self.tr) \
                    + (self.cp[1]/2)*(t**2 - self.tr**2) \
                    + (self.cp[2]/3)*(t**3 - self.tr**3)

        s_ideal = self.cp[0]*np.log(t/self.tr) \
                + (self.cp[1])*(t - self.tr) \
                + (self.cp[2]/2)*(t**2 - self.tr**2) \
                - self.r*np.log(p/self.pr)

        # compressibility factors:
        zl = p*vl/(self.r*t)
        zv = p*vv/(self.r*t)
    
        # enthalpy and entropy departure functions:
        (hl_dep,sl_dep) = self.departure(t,p,zl)
        (hv_dep,sv_dep) = self.departure(t,p,zv)
    
        # enthalpy [J/mol]:
        hl = hl_dep - self.hr_dep + h_ideal
        hv = hv_dep - self.hr_dep + h_ideal    

        # entropy [J/(mol K)]:
        sl = sl_dep - self.sr_dep + s_ideal
        sv = sv_dep - self.sr_dep + s_ideal
        
        # set values to real in case complex:
        p  = np.real(p)
        hl = np.real(hl)
        sl = np.real(sl)
        vl = np.real(vl)
        hv = np.real(hv)
        sv = np.real(sv)
        vv = np.real(vv)

        # return saturation pressure, enthalpy, entropy and specific volume:
        if (q == 0): # return liquid properties
            self.fluid.temperature   = t
            self.fluid.pressure      = p
            self.fluid.enthalpy_mol  = hl
            self.fluid.entropy_mol   = sl
            self.fluid.volume_mol    = vl
            self.fluid.enthalpy_mass = hl/self.wm
            self.fluid.entropy_mass  = sl/self.wm
            self.fluid.density_mass  = self.wm/vl
        elif (q == 1): # return vapour properties
            self.fluid.temperature   = t
            self.fluid.pressure      = p 
            self.fluid.enthalpy_mol  = hv
            self.fluid.entropy_mol   = sv
            self.fluid.volume_mol    = vv
            self.fluid.enthalpy_mass = hv/self.wm
            self.fluid.entropy_mass  = sv/self.wm
            self.fluid.density_mass  = self.wm/vv
        elif (q == 2): # return liquid and vapour properties:
            self.tsat = t
            self.psat = p
            #self.hmolar_l  = hl
            #self.smolar_l  = sl
            #self.vmolar_l  = vl
            self.hmass_l   = hl/self.wm
            self.smass_l   = sl/self.wm
            self.rhomass_l = self.wm/vl
            #self.hmolar_v  = hv
            #self.smolar_v  = sv
            #self.vmolar_v  = vv
            self.hmass_v   = hv/self.wm
            self.smass_v   = sv/self.wm
            self.rhomass_v = self.wm/vv
    
    ###########################################################################
    # saturation properties at specified pressure:
    def saturation_p(self,p,q):    

        # initial guesses for temperature:
        t1 = self.tr
        t2 = self.tc

        # determine pressures at initial guesses:
        self.saturation_t(t1,q)
        if q != 2:
            p1 = self.fluid.p()
        else:
            p1 = self.psat
        self.saturation_t(t2,q)
        if q != 2:
            p2 = self.fluid.p()
        else:
            p2 = self.psat

        # iteration (secant method):
        while (abs(p1 - p2) > 10):
            m  = (t2 - t1)/(p2 - p1)
            c  = t1 - m*p1
            t1 = t2
            t2 = m*p + c
            p1 = p2
            self.saturation_t(t2,q)
            if q != 2:
                p2 = self.fluid.p()
            else:
                p2 = self.psat
    
        # converged temperature:
        t = t2
        self.saturation_t(t,q)
        
    #########################################################################
    # construct saturation curve
    def saturation_curve(self,n,*args):

        # set minimum pressure:
        if (len(args) == 1):
            pmin = args[0]
        else:
            pmin = 1e3

        # define pressure array:
        psat = np.linspace(pmin,0.98*self.pc,n)

        # initialise arrays:
        Tsat  = []
        ssatl = []
        ssatv = []
        hsatl = []
        hsatv = []
        dsatl = []
        dsatv = []

        # compute saturation properties:
        for i in range(n):

            # saturated liquid properties:
            self.saturation_p(psat[i],2)
            Tsat.append(self.tsat)
            ssatl.append(self.smass_l)
            hsatl.append(self.hmass_l)
            dsatl.append(self.rhomass_l)
            ssatv.append(self.smass_v)
            hsatv.append(self.hmass_v)
            dsatv.append(self.rhomass_v)

        # combine liquid and vapour properties:
        self.Tsat = np.concatenate([Tsat,np.flip(Tsat)])
        self.psat = np.concatenate([psat,np.flip(psat)])
        self.hsat = np.concatenate([hsatl,np.flip(hsatv)])
        self.ssat = np.concatenate([ssatl,np.flip(ssatv)])
        self.dsat = np.concatenate([dsatl,np.flip(dsatv)])
        
    ##########################################################################
    # calculation of speed of sound:
    def soundspeed_mol(self,t,v):

        # ideal specific heat capacities:
        cp_id = self.cp[0] + self.cp[1]*t + self.cp[2]*t**2
        cv_id = cp_id - self.r

        # first derivative of pressure wrt temperature @ fixed volume:
        dpdt_A = self.r/(-self.b+v)
        dpdt_B = self.a0*self.n*(1+self.n*(1-np.sqrt(t/self.tc)))
        dpdt_C = (np.sqrt(self.tc*t)*(v**2 + 2*self.b*v - self.b**2))
        dpdt   = dpdt_A + dpdt_B/dpdt_C 

        # first derivative of pressure wrt volume @ fixed temperature:
        dpdv_A = -self.r*t/((-self.b+v)**2)
        dpdv_B = self.a0*(2*(self.b+v))*((1 + self.n*(1-np.sqrt(t/self.tc)))**2)
        dpdv_C = (v**2 + 2*self.b*v - self.b**2)**2
        dpdv   = dpdv_A + dpdv_B/dpdv_C   
 
        # first derivative of volume wrt temperature @ fixed pressure:
        dvdt = -dpdt/dpdv

        # real specific volume term:
        d2alpha_dT2 = self.n*(1+self.n)*((self.tc/t)**0.5)/(2*self.tc*t)
        intgral = np.log((v-self.b*(np.sqrt(2)-1))/(v+self.b*(np.sqrt(2)+1)))
        cv_re = -self.a0*t*d2alpha_dT2*intgral/(2*np.sqrt(2)*self.b)
    
        # specific heat capacity at constant volume:    
        cv = cv_id + cv_re

        # specific heat capacity at constant pressure:
        cp = cv + t*dpdt*dvdt

        # speed of sound:
        a = np.sqrt(-(v**2)*(cp/cv)*dpdv)

        # return speed of sound:
        return a
    
# =========================================================================
# class for a CoolProp fluid using various backends
class coolprop_fluid:    
    
    #######################################################################
    # initialisation of class: 
    def __init__(self,name,backend):

        # create abstract state:
        self.fluid = CoolProp.AbstractState(backend,name)

        # append fluid parameters:
        self.pc = self.fluid.p_critical()
        self.tc = self.fluid.T_critical()
        self.wm = self.fluid.molar_mass()*1000
        
        # saturation properties (single call for both phases):     
        self.tsat      = np.nan
        self.psat      = np.nan
        self.smass_l   = np.nan
        self.hmass_l   = np.nan
        self.rhomass_l = np.nan
        self.smass_v   = np.nan
        self.hmass_v   = np.nan
        self.rhomass_v = np.nan
        
    #######################################################################
    # direct CoolProp calls:
    def tdflash(self,t,d,*args):
        self.fluid.update(CoolProp.DmassT_INPUTS,d,t)
    def tpflash(self,t,p,*args): 
        self.fluid.update(CoolProp.PT_INPUTS,p,t)
    def pdflash(self,p,d,*args):
        self.fluid.update(CoolProp.DmassP_INPUTS,d,p)
    def psflash_mol(self,p,s,*args):
        self.fluid.update(CoolProp.PSmolar_INPUTS,p,s)
    def phflash_mol(self,p,h,*args):
        self.fluid.update(CoolProp.HmolarP_INPUTS,h,p)
    def hsflash_mol(self,h,s,*args):
        self.fluid.update(CoolProp.HmolarSmolar_INPUTS,h,s)    
    def pdflash_mass(self,p,d,*args):
        self.fluid.update(CoolProp.DmassP_INPUTS,d,p)
    def psflash_mass(self,p,s,*args):
        self.fluid.update(CoolProp.PSmass_INPUTS,p,s)
    def phflash_mass(self,p,h,*args):
        self.fluid.update(CoolProp.HmassP_INPUTS,h,p)
    def hsflash_mass(self,h,s,*args):
        self.fluid.update(CoolProp.HmassSmass_INPUTS,h,s)
            
    ##########################################################################
    def saturation_t(self,t,q,*args):
        if (q != 2): # return specified phase
            self.fluid.update(CoolProp.QT_INPUTS,q,t)
        elif (q == 2): # return liquid and vapour properties
            
            # saturated liquid properties:
            self.fluid.update(CoolProp.QT_INPUTS,0,t)
            #self.hmolar_l  = self.fluid.hmolar()
            #self.smolar_l  = self.fluid.smolar()
            self.tsat      = t
            self.psat      = self.fluid.p()
            self.hmass_l   = self.fluid.hmass()
            self.smass_l   = self.fluid.smass()
            self.rhomass_l = self.fluid.rhomass()

            # saturated vapour properties:
            self.fluid.update(CoolProp.QT_INPUTS,1,t)    
            #self.hmolar_v  = self.fluid.hmolar()
            #self.smolar_v  = self.fluid.smolar()
            self.hmass_v   = self.fluid.hmass()
            self.smass_v   = self.fluid.smass()
            self.rhomass_v = self.fluid.rhomass()
            
    ##########################################################################
    def saturation_p(self,p,q,*args): 
        if (q != 2): # return specified phase
            self.fluid.update(CoolProp.PQ_INPUTS,p,q)
        elif (q == 2): # return liquid and vapour properties
            
            # saturated liquid properties:
            self.fluid.update(CoolProp.PQ_INPUTS,p,0)
            #self.hmolar_l  = self.fluid.hmolar()
            #self.smolar_l  = self.fluid.smolar()
            self.tsat      = self.fluid.T()
            self.psat      = p
            self.hmass_l   = self.fluid.hmass()
            self.smass_l   = self.fluid.smass()
            self.rhomass_l = self.fluid.rhomass()

            # saturated vapour properties:
            self.fluid.update(CoolProp.PQ_INPUTS,p,1)
            #self.hmolar_v  = self.fluid.hmolar()
            #self.smolar_v  = self.fluid.smolar()
            self.hmass_v   = self.fluid.hmass()
            self.smass_v   = self.fluid.smass()
            self.rhomass_v = self.fluid.rhomass()
            
    ###########################################################################
    # construct saturation curve
    def saturation_curve(self,n,*args):

        # set minimum pressure:
        if (len(args) == 1):
            pmin = args[0]
        else:
            pmin = 1e3 

        # define pressure array:
        psat  = np.linspace(pmin,pmax,n)
        
        # initialise arrays:
        Tsatl = []
        Tsatv = []
        ssatl = [] 
        ssatv = []
        hsatl = []
        hsatv = []

        # compute saturation properties:
        for i in range(n):

            # saturated liquid properties:
            self.saturation_p(psat[i],0)
            ssatl.append(self.fluid.smass())
            hsatl.append(self.fluid.hmass())
            Tsatl.append(self.fluid.T())

            # saturated vapour properties:
            self.saturation_p(psat[i],1)
            ssatv.append(self.fluid.smass())
            hsatv.append(self.fluid.hmass())            
            Tsatv.append(self.fluid.T())

        # combine liquid and vapour properties:
        self.Tsat = np.concatenate([Tsatl,np.flip(Tsatv)])
        self.psat = np.concatenate([psat ,np.flip(psat)])
        self.hsat = np.concatenate([hsatl,np.flip(hsatv)])
        self.ssat = np.concatenate([ssatl,np.flip(ssatv)])

# =========================================================================
# END OF FILE
# =========================================================================
