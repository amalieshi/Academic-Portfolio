%% HW 4. Problem 1
y1=@(x) 1/(cosh(x)+(tanh(x)*sinh(x))+sinh(x));
y2=@(x) 1/(cosh(x)+(1/tanh(x)*sinh(x))+sinh(x));
x0=0.001;
one_input=fzero(@(x) y1(x)-1/4,x0);
two_inputs=fzero(@(x) 2*y2(x)-1/4,x0);

lvalue = linspace(0,5,100);

voltage1 = NaN(100,1);
voltage2 = voltage1;
for i = 1:100
voltage1(i) = y1(lvalue(i));
voltage2(i) = 2*y2(lvalue(i));
end

figure; 
subplot(121)
plot(lvalue,voltage1-1/4,'DisplayName','One Input')
xlabel('Dimensionless Length'), ylabel('Voltage (mV)')
legend('Show','Location','Best')
hline(0);vline(one_input);
text(one_input,0,'One Input Steady State')
title('One Synpatic Input Voltage versus L')

subplot(122)
plot(lvalue,voltage2-1/4,'DisplayName','Two Inputs')
xlabel('Dimensionless Length'), ylabel('Voltage (mV)')
legend('Show','Location','Best')
hline(0);vline(two_inputs);
text(two_inputs,0,'Two Inputs Steady State')
title('Two Synpatic Inputs Voltage versus L')

