% MATLAB script to carry out a comparison between the predictions obtained
% from the pocketTHERM online calculator to an independent MATLAB model
% that calculates fluid properties using the Helmholtz equation of state
% available within CoolProp.

clear variables; clc;

% import pocketORC results:
R1233zdE = readtable('pocketTHERM_results.xlsx','Sheet','ORC_R1233zd(E)');
pentane  = readtable('pocketTHERM_results.xlsx','Sheet','ORC_n-pentane');
MM       = readtable('pocketTHERM_results.xlsx','Sheet','ORC_MM');

% import CoolProp module:
CoolProp = py.importlib.import_module('CoolProp.CoolProp');
backend = 'REFPROP';    % 'REFPROP' or 'HEOS' (HEOS = CoolProp)

% fluids:
fluids = { 'R1233zdE','Pentane','MM' };

% ranges for parametric study (min. and max. pressure ratios):
pr_min = [ 2  2  2  ];
pr_max = [ 15 20 40 ];

% resolution of CoolProp simulations:
n = 100;

% heat-source conditions:
thi = 473;                  % heat-source inlet temperature, K
tho = 423;                  % heat-source outlet temperature, K
mh  = 1;                    % heat-source mass-flow rate, kg/s
cph = 1000;                 % heat-source specific-heat capacity, J/kg K

% heat-sink conditions:
tci = 288;                  % heat-sink inlet temperature, K
mc  = 1;                    % heat-sink mass-flow rate, kg/s
cpc = 4200;                 % heat-source specific heat capacity, J/kg K

% component performance:
etap = 0.7;                 % pump isentropic efficiency
etae = 0.8;                 % expander isentropic efficiency
effr = 0.7;                 % recuperator effectiveness

% fixed ORC parameters:
tcond = 313;                % condensation temperature, K
dtsc  = 2;                  % pump subcooling, K
x3    = 1.2;                % expander inlet parameter

% run parametric study:
Wnet = zeros(numel(fluids),n); eta_th = Wnet; Pevap = Wnet;
for i = 1:numel(fluids)

    % setup CoolProp with Helmholtz equation of state:
    fluid = CoolProp.AbstractState(backend,fluids{i});

    % setup array of pressure ratios:
    pr_array = linspace(pr_min(i),pr_max(i),n);

    for j = 1:numel(pr_array)

        % pump inlet conditions:
        fluid.update(CoolProp.QT_INPUTS,0,tcond)
        p1 = fluid.p();
        t1 = tcond - dtsc;
        fluid.update(CoolProp.PT_INPUTS,p1,t1)
        h1 = fluid.hmass();
        s1 = fluid.smass();
    
        % pump outlet conditions:
        p2 = p1*pr_array(j);
        fluid.update(CoolProp.PSmass_INPUTS,p2,s1)
        h2s = fluid.hmass();
        h2 = h1 + (h2s - h1)/etap;
        fluid.update(CoolProp.HmassP_INPUTS,h2,p2)
        t2 = fluid.T();
        s2 = fluid.smass();

        % expander inlet conditions:
        p3 = p2;
        fluid.update(CoolProp.PQ_INPUTS,p3,1)
        tevap = fluid.T();
        t3 = tevap + (x3-1)*(thi - tevap);
        fluid.update(CoolProp.PT_INPUTS,p3,t3)
        h3 = fluid.hmass();
        s3 = fluid.smass();

        % expander outlet conditions:
        p4 = p1;
        fluid.update(CoolProp.PSmass_INPUTS,p4,s3)
        h4s = fluid.hmass();
        h4 = h3 - etae*(h3 - h4s);
        fluid.update(CoolProp.HmassP_INPUTS,h4,p4)
        t4 = fluid.T();
        s4 = fluid.smass();

        % recuperator:
        fluid.update(CoolProp.PT_INPUTS,p2,t4)
        h2max = fluid.hmass();
        if t2 > tcond
            fluid.update(CoolProp.PT_INPUTS,p4,t2)
        else
            fluid.update(CoolProp.PQ_INPUTS,p4,1)
        end
        h4min = fluid.hmass();
        dhmax = min([h2max-h2, h4-h4min]);
        h2r = h2 + effr*dhmax;
        h4r = h4 - effr*dhmax;

        % ORC mass-flow rate:
        m  = (mh*cph*(thi - tho))/(h3 - h2r);

        % thermal and mechanical powers:
        Wp = m*(h2 - h1);
        We = m*(h3 - h4);
        Qh = m*(h3 - h2r);
        
        % net power [kW]:
        Wnet(i,j) = (We - Wp)/1000;

        % cycle thermal efficiency efficiency [%]:
        eta_th(i,j) = 100*((We-Wp)/Qh);

        % evaporation pressure [bar]:
        Pevap(i,j) = p2/1e5;

    end
end

% plot results - net power:
figure
hold on
plot(Pevap(1,:),Wnet(1,:),'k-')
plot(Pevap(2,:),Wnet(2,:),'r-')
plot(Pevap(3,:),Wnet(3,:),'b-')
plot(R1233zdE.EvaporationPressure_bar_,R1233zdE.Net_kW_,...
    'ko','markerfacecolor','k')
plot(pentane.EvaporationPressure_bar_,pentane.Net_kW_,...
    'ro','markerfacecolor','r')
plot(MM.EvaporationPressure_bar_,MM.Net_kW_,...
    'bo','markerfacecolor','b')
hold off
xlabel('Evaporation pressure [bar]','fontsize',14)
ylabel('Net power [kW]','fontsize',14)
set(gca,'fontsize',14,'box','on','xgrid','on','ygrid','on')
legend('R1233zd(E): Helmholtz','n-pentane: Helmholtz','MM: Helmholtz',...
    'R1233zd(E): Peng-Robinson','n-pentane: Peng-Robinson',...
    'MM: Peng-Robinson','location','southeast','fontsize',14)
    
% plot results - cycle thermal efficiency:
figure
hold on
plot(Pevap(1,:),eta_th(1,:),'k-')
plot(Pevap(2,:),eta_th(2,:),'r-')
plot(Pevap(3,:),eta_th(3,:),'b-')
plot(R1233zdE.EvaporationPressure_bar_,R1233zdE.Efficiency_cycle____,...
    'ko','markerfacecolor','k')
plot(pentane.EvaporationPressure_bar_,pentane.Efficiency_cycle____,...
    'ro','markerfacecolor','r')
plot(MM.EvaporationPressure_bar_,MM.Efficiency_cycle____,'bo',...
    'markerfacecolor','b')
hold off
xlabel('Evaporation pressure [bar]')
ylabel('Cycle thermal efficiency [%]')
set(gca,'fontsize',14,'box','on','xgrid','on','ygrid','on')
legend('R1233zd(E): Helmholtz','n-pentane: Helmholtz','MM: Helmholtz',...
    'R1233zd(E): Peng-Robinson','n-pentane: Peng-Robinson',...
    'MM: Peng-Robinson','location','southeast','fontsize',14)

% END OF FILE