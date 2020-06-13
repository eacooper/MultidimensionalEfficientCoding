% This script generates the manifold potentials and displacement fields for
% efficiently encoding 2D stimulus probability distributions. The script
% generates a large number of example distributions from bivariate
% Gaussians, to be used in the analysis characterizing correlations between
% tuning curve properties and stimulus probabilty
clear all; close all;

%% Define Parameters

% stimulus space is centered on 0 and goes from -1 - 1
stim_range   = [-1 1]; 

% number of bins for probability distribution in each dimension, must be odd
nbin = 201;

% sample stimulus space
sx = linspace(stim_range(1),stim_range(2),nbin+1);
sy = linspace(stim_range(2),stim_range(1),nbin+1);

% convert lattice points to bin centers
sbx  = sx(1:end-1) + diff(sx(1:2))/2;  % bin center x
sby  = sy(1:end-1) + diff(sy(1:2))/2;  % bin center y

% generate mesh 
[sbmx, sbmy] = meshgrid(sbx, sby);

% number of simulations
sims_start  = 165;
sims_end    = 500;

% possible sigmas for probability distribution
stdev_min = 0.1;
stdev_max = 0.4;

%% Run simulations

% for each sim
for n = sims_start:sims_end
    
    % stimulus probability
    stdev_px = rand*(stdev_max - stdev_min) + stdev_min;
    stdev_py = rand*(stdev_max - stdev_min) + stdev_min;
    theta = 180*rand*(pi/180);
            
    p2 = @(x,y) exp( -( ( ( (cos(theta).^2)./(2.*stdev_px.^2) + (sin(theta).^2)./(2.*stdev_py.^2) ).*x.^2 ) ...
                + (2.* ( -(sin(2.*theta))./(4.*stdev_px.^2) + (sin(2.*theta))./(4.*stdev_py.^2) ).*x.*y) ...
                + ( ( (sin(theta).^2)./(2.*stdev_px.^2) + (cos(theta).^2)./(2.*stdev_py.^2) ).*y.^2) ) );
  
    % evaluate function
    pD = p2(sbmx,sbmy);
    
    % normalize probability
    pD = pD./sum(pD(:));
    
    % determine manifold potential
    [F,Fx,Fy,Fxi,Fyi,pU,fhand] = compute_uniform_transform(nbin,pD,'circle');
    
    % save results and diagnostic figure
    saveas(fhand,['sims/manifolds/plots/sim' num2str(n) '_sx' num2str(stdev_px,2) '_sy' num2str(stdev_py,2) '_th' num2str(round(theta*180/pi)) '_optimization_result.eps'],'epsc')
    close(fhand);
    save(['sims/manifolds/data/sim' num2str(n) '_sx' num2str(stdev_px,2) '_sy' num2str(stdev_py,2) '_th' num2str(round(theta*180/pi)) '_data.mat'],'sbx','sby','sbmx','sbmy','stdev_px','stdev_py','theta','pD','F','Fx','Fy','Fxi','Fyi','pU');
    
end



