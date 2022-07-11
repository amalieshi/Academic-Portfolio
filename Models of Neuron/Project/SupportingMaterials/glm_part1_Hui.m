% Script glm_part1.m
% MATLAB code to visualize data, fit a GLM model of the relation between
% spiking and the rat's position, and visualize this model for the
% Neuroinformatics GLM problem set.
% The code is initialized with an overly simple GLM model construction.
% Please improve it!

% load the rat trajectory and spiking data;
close all;
load('train.mat');

% fit a GLM model to the x and y positions.  Can you think of a better
% model form based on the visualized data?
[b,dev,stats] = glmfit([xN yN],spikes_binned,'poisson');

%visualize your model
% construct a grid of positions to plot the model against...
figure;
[x_new,y_new]=meshgrid(-1:.1:1);
y_new = flipud(y_new);
x_new = fliplr(x_new);



X = NaN(T-5,8);

% compute lambda for each point on this grid using the GLM model
lambda = exp(B'*T);
lambda(find(x_new.^2+y_new.^2>1))=nan;

% dN1 =0;
% L1 = exp(dN1.*log(lambda)-log(factorial(dN1))-lambda);
% for i = 1:length(spikes_binned)
%   dNi=spikes_binned(i);
%   L=exp(dNi.*log(lambda)-log(factorial(dNi))-lambda); 
%   if i ==1
%     L2=L1.*L; 
%   else 
%     L2=L.*L2;
%   end
%       
% end


%plot lambda as a function position over this grid
h_mesh = mesh(x_new,y_new,lambda,'AlphaData',0);
get(h_mesh,'AlphaData');
set(h_mesh,'AlphaData',0);
hold on;
plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
xlabel('x{k}'); ylabel('y{k}'); zlabel('lambda value')
title('predicted lambda strength based on the two covariates')


