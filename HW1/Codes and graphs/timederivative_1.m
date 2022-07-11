function dPdt=timederivative(P)

%Define variables for the membrane properties
Cm=100; %pF
N=100;
E=-70; %mV
Gisingle=100; %nS

%Define variables for the occupancy probability
kappa=1.38e-23; %J/K
h=6.626e-34; %J/s
T=293; %K
R=8.314; %J/(K*mol)
F=96485; %C/mol
z=1;
%Assumptions
G1=6e4; %J
G2=1e3; %J
G3=6e4; %J

k12=  kappa*T/h * exp(-1/(R*T) * (G1 +      (z*1e-3*F*P(4)/4)));
k21=  kappa*T/h * exp(-1/(R*T) * (G1 - G2 - (z*1e-3*F*P(4)/4)));
k23=  kappa*T/h * exp(-1/(R*T) * (G3 - G2 + (z*1e-3*F*P(4)/4)));
k32=  kappa*T/h * exp(-1/(R*T) * (G3 -      (z*1e-3*F*P(4)/4)));

dPdt = [k21 * P(2) - k12 * P(1);...%dP1/dt
        k12 * P(1) + k32 * P(3) - (k21+k23) * P(2);...%dP2/dt
        k23 * P(2) - k32 * P(3);...%dP3/dt
        -P(3) * N * Gisingle *(P(4)-E)/Cm];%%dV/dt 