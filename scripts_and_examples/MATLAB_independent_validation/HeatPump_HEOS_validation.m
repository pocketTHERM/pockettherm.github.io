% MATLAB script to carry out a comparison between the predictions obtained
% from the pocketTHERM online calculator to an independent MATLAB model
% that calculates fluid properties using the Helmholtz equation of state
% available within CoolProp.

clear variables; clc;

% import pocketORC results:
R1233zdE = readtable('pocketTHERM_results.xlsx','Sheet','HP_R1233zd(E)');
pentane  = readtable('pocketTHERM_results.xlsx','Sheet','HP_n-pentane');

% import CoolProp module:
CoolProp = py.importlib.import_module('CoolProp.CoolProp');
backend = 'REFPROP';    % 'REFPROP' or 'HEOS' (HEOS = CoolProp)

% fluids:
fluids = { 'R1233zdE','Pentane'};

% ranges for parametric study (min. and max. pressure ratios):
pr_min = [ 2  2  ];
pr_max = [ 10 10 ];

% resolution of CoolProp simulations:
n = 100;

% heat-source conditions:
tci = 323;                  % heat-source inlet temperature, K
tco = 318;                  % heat-source outlet temperature, K
mc  = 1;                    % heat-source mass-flow rate, kg/s
cpc = 4200;                 % heat-source specific-heat capacity, J/kg K

% heat-sink conditions:
thi = 328;                  % heat-sink inlet temperature, K
mh  = 1;                    % heat-sink mass-flow rate, kg/s
cph = 4200;                 % heat-source specific heat capacity, J/kg K

% component performance:
etac = 0.7;                 % compressor isentropic efficiency

% fixed ORC parameters:
tevap  = 313;               % evaporation temperature, K
dtsh   = 5;                 % compresor superheating, K
dtsc   = 5;                 % expansion valve subcooling, K
dt_ihe = 15;                % internal heat exchanger temperature diff, K

% run parametric study:
COP = zeros(numel(fluids),n); Prat = COP; Qh_kW = COP;
for i = 1:numel(fluids)

    % setup CoolProp with Helmholtz equation of state:
    fluid = CoolProp.AbstractState(backend,fluids{i});

    % setup array of pressure ratios:
    pr_array = linspace(pr_min(i),pr_max(i),n);

    for j = 1:numel(pr_array)

        % evaporator outlet conditions:
        fluid.update(CoolProp.QT_INPUTS,0,tevap)
        p1 = fluid.p;
        t1 = tevap + dtsh;
        fluid.update(CoolProp.PT_INPUTS,p1,t1)
        h1 = fluid.hmass();
        s1 = fluid.smass();
  
        % compressor inlet conditions:
        t2 = t1 + dt_ihe;
        p2 = p1;
        fluid.update(CoolProp.PT_INPUTS,p2,t2)
        h2 = fluid.hmass();
        s2 = fluid.smass();
    
        % compressor outlet conditions:
        p3 = p2*pr_array(j);
        fluid.update(CoolProp.PSmass_INPUTS,p3,s2)
        h3s = fluid.hmass();
        h3 = h2 + (h3s - h2)/etac;
        fluid.update(CoolProp.HmassP_INPUTS,h3,p3)
        t3 = fluid.T();
        s3 = fluid.smass();

        % condenser outlet conditions:
        p4 = p3;
        fluid.update(CoolProp.PQ_INPUTS,p4,1)
        tcond = fluid.T();
        t4 = tcond - dtsc;
        fluid.update(CoolProp.PT_INPUTS,p4,t4)
        h4 = fluid.hmass();
        s4 = fluid.smass();

        % expansion valve inlet conditions:
        p5 = p4;
        h5 = h4 - (h2 - h1);

        % expansion valuve outlet conditions:
        p6 = p1;
        h6 = h5;

        % heat pump mass-flow rate:
        m  = (mc*cpc*(tci - tco))/(h1 - h6);

        % compressor work:
        Wc = m*(h3 - h2);
        Qh = m*(h3 - h4);

        % coefficient of performance:
        COP(i,j) = Qh/Wc;

        % heat generated [kW]:
        Qh_kW(i,j) = Qh/1000;

        % evaporation pressure [bar]:
        Prat(i,j) = pr_array(j);

    end
end

% plot results - net power:
figure
hold on
plot(Prat(1,:),Qh_kW(1,:),'k-')
plot(Prat(2,:),Qh_kW(2,:),'r-')
plot(R1233zdE.PressureRatio,R1233zdE.Heater_kW_,...
    'ko','markerfacecolor','k')
plot(pentane.PressureRatio,pentane.Heater_kW_,...
    'ro','markerfacecolor','r')
hold off
xlabel('Pressure ratio [-]','fontsize',14)
ylabel('Heat generated [kW]','fontsize',14)
set(gca,'fontsize',14,'box','on','xgrid','on','ygrid','on')
legend('R1233zd(E): Helmholtz','n-pentane: Helmholtz',...
    'R1233zd(E): Peng-Robinson','n-pentane: Peng-Robinson',...
    'location','southeast','fontsize',12)
    
% plot results - cycle thermal efficiency:
figure
hold on
plot(Prat(1,:),COP(1,:),'k-')
plot(Prat(2,:),COP(2,:),'r-')
plot(R1233zdE.PressureRatio,R1233zdE.CoefficientOfPerformance___,...
    'ko','markerfacecolor','k')
plot(pentane.PressureRatio,pentane.CoefficientOfPerformance___,...
    'ro','markerfacecolor','r')
hold off
xlabel('Pressure ratio [-]','fontsize',14)
ylabel('Coefficient of performance (heating) [-]')
set(gca,'fontsize',14,'box','on','xgrid','on','ygrid','on')
legend('R1233zd(E): Helmholtz','n-pentane: Helmholtz',...
    'R1233zd(E): Peng-Robinson','n-pentane: Peng-Robinson',...
    'MM - Peng-Robinson','location','northeast','fontsize',14)

% END OF FILE