%Hui Shi
%created on 2017-09-17
figure(1)
tspan = [0:.001:0.05];
inits = [1; 0; 0; 0];

[t,P] = ode23(@(t,P) timederivative_h1(P),tspan,inits);
[t1,P1] = ode23(@(t,P) timederivative_1(P),tspan,inits);
plot(t,P(:,1),'-g','DisplayName', 'P1_New')
title('Occupancy Probability of State 1, 2 & 3')
xlabel('Time [s]')

hold on
plot(t,P(:,2),'-k','DisplayName', 'P2_New')
plot(t,P(:,3),'-r','DisplayName', 'P3_New')
plot(t,P1(:,1),'og','DisplayName', 'P1')
plot(t,P1(:,2),'ok','DisplayName', 'P2')
plot(t,P1(:,3),'or','DisplayName', 'P3')
legend('show','Location','Best')


