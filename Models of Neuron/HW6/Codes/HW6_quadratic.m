load glm_data.mat
[b,dev,stats] = glmfit([xN.^2 yN.^2],spikes_binned,'poisson');

%visualize your model
% construct a grid of positions to plot the model against...
figure;
[x_new,y_new]=meshgrid(-1:.1:1);
y_new = flipud(y_new);
x_new = fliplr(x_new);

% compute lambda for each point on this grid using the GLM model
lambda = exp(b(1) + b(2)*x_new.^2 + b(3)*y_new.^2);
lambda(find(x_new.^2+y_new.^2>1))=nan;

%plot lambda as a function position over this grid
h_mesh = mesh(x_new,y_new,lambda,'AlphaData',0);
get(h_mesh,'AlphaData');
set(h_mesh,'AlphaData',0);
hold on;
plot3(cos(-pi:1e-2:pi),sin(-pi:1e-2:pi),zeros(size(-pi:1e-2:pi)));
xlabel('x{k}'); ylabel('y{k}'); zlabel('lambda value')
title('predicted lambda strength based on the two covariates')

%Examine the confidence intervals computed for each parameter of two models
%based on the square root of the inverse of the observed Fisher information
figure;
Fisher=stats.se;
CI = sqrt(1./Fisher)
errorbar(b,2.*Fisher);
xlim([0,4]);
% All of the b values are statistically different from zeros

pval = stats.p;
%All of the p values are very small and they are below 0.05 so 
%they are all significant

[yhat,dylo,dyhi] = glmval(b,[[-1:.01:1].^2; [-1:.01:1].^2]','log',stats);
figure;
errorbar(-1:.01:1,yhat,dylo,dyhi);

%% 3. Characterize the relative goodness-of-fit among competing models
elements = [xN.^2 yN.^2];
parameters = numel(stats.se);
AIC = dev + 2*parameters

%% 4. 
lambdaEst = exp(b(1)+b(2)*xN.^2+ b(3)*yN.^2); 
                       
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
figure;
plot( ([1:N]-.5)/N, KSSorted, 'b', 0:.01:1,0:.01:1, 'g',0:.01:1, [0:.01:1]+1.36/sqrt(N), 'r', 0:.01:1,[0:.01:1]-1.36/sqrt(N), 'r' );
axis( [0 1 0 1] );
xlabel('Uniform CDF');
ylabel('Empirical CDF of Rescaled ISIs');
title('KS Plot with 95% Confidence Intervals');

ks_stat = max(abs(KSSorted - ([1:N]-.5)/N))


