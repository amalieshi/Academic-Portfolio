% Script glm_part1.m
% MATLAB code to visualize data, fit a GLM model of the relation between
% spiking and the rat's position, and visualize this model for the
% Neuroinformatics GLM problem set.
% The code is initialized with an overly simple GLM model construction.
% Please improve it!

% load the rat trajectory and spiking data;
close all;
load('glm_data.mat');

% visualize the raw data
figure;
plot(xN,yN,x_at_spiketimes,y_at_spiketimes,'r.');
axis tight square;
xlabel('x position (m)'); ylabel('y position (m)');

% fit a GLM model to the x and y positions.  Can you think of a better
% model form based on the visualized data?
[b,dev,stats] = glmfit([xN yN],spikes_binned,'poisson');

%visualize your model
% construct a grid of positions to plot the model against...
figure;
[x_new,y_new]=meshgrid(-1:.1:1);
y_new = flipud(y_new);
x_new = fliplr(x_new);

% compute lambda for each point on this grid using the GLM model
lambda = exp(b(1) + b(2)*x_new + b(3)*y_new);
lambda(find(x_new.^2+y_new.^2>1))=nan;

%plot lambda as a function position over this grid
h_mesh = mesh(x_new,y_new,lambda,'AlphaData',0);
get(h_mesh,'AlphaData');
set(h_mesh,'AlphaData',0);
hold on;
plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
