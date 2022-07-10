figure(1)
clf

% Define constants
Gi = 10; % nS
Ei = 55; % mV
CA = 1e6;  % uF/cm^2
A = 1e-4; %cm^2
C=CA*A;%pF
Iapp=-100; %pA

%Case i 
t_points=[0:0.01:35];
tl=length(t_points);
Iapp_m=zeros(1,tl);
pulselogic=t_points<=30 & t_points>=5;
Iapp_m(pulselogic)=-100;
Vmi=NaN(1,tl);
for i=1:tl
    Vmi(i)=(-Iapp_m(i)+Ei*Gi)/Gi;
end

%Case ii
%Between time slice 5ms to 30ms
dViidt = @(Vii) -Iapp/C  ;

% Define span of time to evalute over (in miliseconds)
% Evaluates from 0 to 5 ms at intervals of 0.01 ms
t_span1 = [5:0.01:30];

% Define initial conditions to ODE
V0 = 55; % mV

% Run ODE23 to get V over time
[tii,Vii] = ode23(@(tii,Vii) dViidt(Vii), t_span1, V0);
Vii=Vii.';

t1=[0:0.01:4.99];
tl1=length(t1);
Vprior=Ei*ones(1,tl1);
t2=[30.01:0.01:35];
tl2=length(t2);
Vlatter=Ei*ones(1,tl2);
Viif=horzcat(Vprior,Vii,Vlatter);


%Case iii
%-------------------
% Define ODE as a function of membrane potential, V
dVdt = @(V1)  (Gi*(Ei-V1))/C;

% Define span of time to evalute over (in miliseconds)
% Evaluates from 0 to 5 ms at intervals of 0.01 ms
t_span1 = [0:0.01:5];

% Define initial conditions to ODE
V0 = 55; % mV

% Run ODE23 to get V over time
[t1,V1] = ode23(@(t1,V1) dVdt(V1), t_span1, V0);
V1=V1.';
t1=t1.';

%---------------------------

% Define ODE as a function of membrane potential, V
dVdt = @(V2)  (Gi*(Ei-V2)-Iapp)/C;

% Define span of time to evalute over (in miliseconds)
% Evaluates from 0 to 35 ms at intervals of 0.01 ms
t_span2 = [5:0.01:30];

% Define initial conditions to ODE
V0 = V1(end); % mV

% Run ODE23 to get V over time
[t2,V2] = ode23(@(t2,V2) dVdt(V2), t_span2, V0);
V2=V2.';
t2=t2.';



%-------------------
% Define ODE as a function of membrane potential, V
dVdt = @(V3)  (Gi*(Ei-V3))/C;

% Define span of time to evalute over (in miliseconds)
% Evaluates from 0 to 35 ms at intervals of 0.01 ms
t_span3 = [30:0.01:35];

% Define initial conditions to ODE
VX = V2(end); % mV

% Run ODE23 to get V over time
[t3,V3] = ode23(@(t3,V3) dVdt(V3), t_span3, VX);
V3=V3.';
t3=t3.';

t=horzcat(t1,t2,t3);
V=horzcat(V1,V2,V3);

%-------------------------------------------------------
% Plot
plot(t_points,Vmi,'--', 'DisplayName', 'Situation i')
title('RC-circuit with current model')
xlabel('Time [ms]')
ylabel('V(t) [mV]')

hold on
plot(t_points,Viif,'--','DisplayName','Situation ii')
plot(t,V,'--','DisplayName','Situation iii')
legend('show','Location','Best')
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
figure(2)
t=[0:0.01:35];
length_t=length(t);
Iapp=zeros(length_t,1);
a=t>=5 & t<=30;
Iapp(a)=-100;

plot(t,Iapp, 'DisplayName', 'Numerical')
title('Current Pulse Over Time')
xlabel('Time [ms]')
ylabel('I(t) [pA]')