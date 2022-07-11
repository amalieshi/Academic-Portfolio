% Script glm_part1_ks.m
% MATLAB code to fit a GLM model of the relation between
% spiking and the rat's position, and draw K-S plot for the
% second Neuroinformatics 2005 GLM problem set.

% load the rat trajectory and spiking data;
load('glm_data.mat');

% visualize the raw data
figure;
title('Visualize the raw data');
plot(xN,yN,x_at_spiketimes,y_at_spiketimes,'r.');


% fit a GLM model to the x and y positions. (ADD ALL YOUR MODEL CANDIDATES
% HERE!!! Two possible models are listed below.

[b_lin,dev_lin,stats_lin] = glmfit([xN yN],spikes_binned,'poisson');

%*******  K-S Plot  *******************
%Note that the ks plot for one model is given below. You should overlay ks
%plots for each of your models.

%graph the K-S plot and confidence intervals for the K-S statistic

%first generate the conditional intensity at each timestep
% ** Adjust the below line according to your choice of model.
% remember to include a column of ones to multiply the default constant GLM parameter beta_0**

 % based on our GLM model with the log "link function"
                       % Use your parameter estimates (b) from glmfit along
                       % with the covariates you used (xN, yN, ...)

lambdaEst = exp(b_lin(1)+b_lin(2)*xN+ b_lin(3)*yN); 
                       
timestep = 1;
lambdaInt = 0;
j=0;
for t=1:length(spikes_binned),
    lambdaInt = lambdaInt + lambdaEst(t)*timestep;
    if (spikes_binned(t)),
        j = j + 1;
        KS(j) = 1-exp(-lambdaInt);
        lambdaInt = 0;
    end;
end;
KSSorted = sort( KS );
N = length( KSSorted);
figure(2);
plot( ([1:N]-.5)/N, KSSorted, 'b', 0:.01:1,0:.01:1, 'g',0:.01:1, [0:.01:1]+1.36/sqrt(N), 'r', 0:.01:1,[0:.01:1]-1.36/sqrt(N), 'r' );
axis( [0 1 0 1] );
xlabel('Uniform CDF');
ylabel('Empirical CDF of Rescaled ISIs');
title('KS Plot with 95% Confidence Intervals');

ks_stat = max(abs(KSSorted - ([1:N]-.5)/N));

