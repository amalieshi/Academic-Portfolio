%Hui Shi
%created on 2017-09-17
figure(1)
tspan = [0:.001:0.2];

Voltage_clamp=[-100:5:100];
vcl=length(Voltage_clamp);
P_st=NaN(vcl,3);
I=NaN(vcl,1);

for i=1:vcl
    vc=Voltage_clamp(i);
    inits = [1; 0; 0; vc];
    [t,P] = ode23(@(t,P) timederivative_1(P),tspan,inits);
    P_ss= P(end,1:3) ;%steady state
    P_st(i,:)=P_ss;
    I(i)=P(3)*100*(vc-70)*1e-3;
end

plot(Voltage_clamp,P_st(:,1),'--','DisplayName', 'P1')
title('Occupancy Probability of State 1, 2 & 3 with different voltage clamp')
xlabel('Voltage[mV]')

hold on
plot(Voltage_clamp,P_st(:,2),'--','DisplayName', 'P2')
plot(Voltage_clamp,P_st(:,3),'--','DisplayName', 'P3')
legend('show','Location','Best')

figure(2)
plot(Voltage_clamp,I,'--')
title('I-V Curve')
xlabel('Voltage [mV]')
ylabel('Current [pA]')



