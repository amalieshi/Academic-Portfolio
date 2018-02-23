%% HW 3. Problem 1
y1=@(L) cosh(L) + (1+tanh(L)) * sinh(L);
y2=@(L) cosh(L) + (1+1/tanh(L)) * sinh(L);
x0=0.001;
one_input=fzero(@(L) y1(L)-4,x0);
two_inputs=fzero(@(L) y2(L)-8,x0);
l = linspace(0,2,100); % Distance from soma

disp([num2str(one_input), ' < L < ',num2str(two_inputs)])

voltage1 = NaN(100,1);
voltage2 = voltage1;
for i = 1:100
voltage1(i) = (y1(l(i))).^-1;
voltage2(i) = 2* (y2(l(i))).^-1;
end

voltage11 = (y1(one_input)).^-1;
voltage22 = 2*(y1(two_inputs)).^-1;

figure(1); 
subplot(121)
plot(l,voltage1,'DisplayName','One Input')
xlabel('Electrotonic length L'), ylabel('Membrane potential at the soma [mV]')
vline(one_input);hline(voltage11);
text(one_input-0.5,voltage11,'One Input Steady State')
title('One Synpatic Input Voltage versus L')

subplot(122)
plot(l,voltage2,'DisplayName','Two Inputs')
xlabel('Electrotonic length L'), ylabel('Membrane potential at the soma [mV]')
vline(two_inputs);hline(voltage22);
text(two_inputs-0.85,voltage22,'Two Inputs Steady State')
title('Two Synpatic Inputs Voltage versus L')
saveas(figure(1),'hw2.png')

