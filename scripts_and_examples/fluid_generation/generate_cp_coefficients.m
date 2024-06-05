clear variables; clc;

% import CoolProp:
CoolProp = py.importlib.import_module('CoolProp.CoolProp');

% fluids:
fluids = { 'R1233zdE','pentane','MM' };

% temperature values used to construct polynomial:
T = linspace(273,573,20);

% minimum pressure:
p_min = 1e3;

% generate coefficients:
for i = 1:numel(fluids)

    % setup REFPROP:
    fluid = CoolProp.AbstractState('REFPROP',fluids{i});

    % fluid parameters:
    Tc(i) = fluid.T_critical(); %#ok<*SAGROW> 
    Pc(i) = fluid.p_critical();
    om(i) = fluid.acentric_factor();
    wm(i) = fluid.molar_mass()*1000;

    % saturation temperature at 1 kPa:
    fluid.update(CoolProp.PQ_INPUTS,1e3,0)
    Tmin = fluid.T();
    if Tmin < 273
        fluid.update(CoolProp.QT_INPUTS,0,273)
        pm(i) = fluid.p();
    else
        pm(i) = pmin;
    end

    % calculate cp value at each temperature value:
    cp = zeros(1,numel(T));
    for j = 1:numel(T)
        try
            fluid.update(CoolProp.PT_INPUTS,1e-99,T(j));
            cp(j) = fluid.molar_mass()*fluid.cpmass();
        catch
            cp(j) = nan;
        end
    end
    
    % fit second order polynomial:
    warning off
    [p,~] = polyfit(T(~isnan(cp)),cp(~isnan(cp)),2);
    warning on    
    
    % reconstruct polynomial vector:
    cp_coeffs(i,:) = [ p(3) p(2) p(1) ];

end