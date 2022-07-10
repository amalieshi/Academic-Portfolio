%Hui Shi
%created on 2017-09-17
figure(1)
tspan = [0:.001:0.2];
inits = [1; 0; 0; 0];

[t,P] = ode23(@(t,P) timederivative_1(P),tspan,inits);

plot(t,P(:,1),'--','DisplayName', 'P1')
title('Occupancy Probability of State 1, 2 & 3')
xlabel('Time [s]')

hold on
plot(t,P(:,2),'--','DisplayName', 'P2')
plot(t,P(:,3),'--','DisplayName', 'P3')
legend('show','Location','Best')




figure(2)
plot(t,P(:,4),'--')
title('Membrane Potential Versus Time')
xlabel('Time [s]')
ylabel('Voltage [mV]')

