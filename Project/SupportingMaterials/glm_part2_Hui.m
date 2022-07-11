% script glm_part2.m
% MATLAB script to constuct and visualize GLM models of spiking against 
% other covariates individually.  
% The script is initialized to visualize the raw data in the form of 
% occupancy normalized histograms.  Use this to constuct GLM models 
% and compare them to the histograms 

% load the rat trajectory and spiking data;
load('glm_data.mat');

% Find spike bin index for occupancy normalized histograms
ind = [];
for numspikes = 1:max(spikes_binned),
    ind = [ind; find(spikes_binned >= numspikes)];
end;

figure;
warning('off');
set(gcf,'Name','Occupancy Normalized Histograms');

% Histogram of spiking to x-velocity
subplot(2,2,1);
    velocities = -.04:.001:.04;
    bar(velocities,hist(vxN(ind),velocities)./hist(vxN,velocities));
    ylabel('normalized spike counts');
    xlabel('x velocity');
     b=glmfit(exp(vxN),spikes_binned,'poisson');
     hold on;
     plot(velocities,exp(b(1)+b(2).*exp((hist(vxN(ind),velocities)./hist(vxN,velocities)))),'r');
     
% Histogram of spiking to y-velocity 
% fit change covariate to square root
subplot(2,2,2);
    bar(velocities,hist(vyN(ind),velocities)./hist(vyN,velocities));
     ylabel('normalized spike counts');
     xlabel('y velocity');
     pow=0.5;
     b=glmfit((vyN).^(pow),spikes_binned,'poisson');
     hold on;
     plot(velocities,(b(1)+b(2).*sqrt(hist(vyN(ind),velocities)./hist(vyN,velocities))).^(pow),'r');
     
    

% Histogram of spiking to movement speed 
subplot(2,2,3);
    rs = 0:.001:.04;
    bar(rs,hist(r(ind),rs)./hist(r,rs));
    xlabel('speed');
    b=glmfit(r.^(1/2),spikes_binned,'poisson');    
    hold on;
    plot(rs,exp(b(1)+b(2).*sqrt(hist(r(ind),rs)./hist(r,rs)))-0.1,'r');

% Histogram of spiking to movement direction
subplot(2,2,4);
    phis = -pi:.1:pi;
    bar(phis,hist(phi(ind),phis)./hist(phi,phis));
    xlabel('direction');
    pow=1/2;
    b=glmfit((phi).^(pow),spikes_binned,'poisson');    
    hold on;
    plot(phis,(b(1)+b(2).*(hist(phi(ind),phis)./hist(phi,phis))).^(pow),'r');
