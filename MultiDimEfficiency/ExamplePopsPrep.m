% This script generates the manifold potentials and displacement fields for
% efficiently encoding a set of example 2D stimulus probability distributions
clear all; close all;


%% Define Stimulus Space Parameters

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


%% Generate examples
examples = 1:5;

% for each example stimulus distribution
for p = examples
    
    %% Stimulus probability
    switch p
        
        case 1
            
            % gaussian in s1, uniform in s2
            stdev_px    = 0.25;
            stdev_py    = 10000;
            theta       = 0*(pi/180);
            
            p2 = @(x,y) exp( -( ( ( (cos(theta).^2)./(2.*stdev_px.^2) + (sin(theta).^2)./(2.*stdev_py.^2) ).*x.^2 ) ...
                + (2.* ( -(sin(2.*theta))./(4.*stdev_px.^2) + (sin(2.*theta))./(4.*stdev_py.^2) ).*x.*y) ...
                + ( ( (sin(theta).^2)./(2.*stdev_px.^2) + (cos(theta).^2)./(2.*stdev_py.^2) ).*y.^2) ) );
            
            % evaluate function
            pD = p2(sbmx,sbmy);
            
            % shape of windowing for manifold
            window_type = 'square';
            
        case 2
            
            % isotropic Gaussian in s1 and s2
            stdev_px    = 0.25;
            stdev_py    = 0.25;
            theta       = 0*(pi/180);
            
            p2 = @(x,y) exp( -( ( ( (cos(theta).^2)./(2.*stdev_px.^2) + (sin(theta).^2)./(2.*stdev_py.^2) ).*x.^2 ) ...
                + (2.* ( -(sin(2.*theta))./(4.*stdev_px.^2) + (sin(2.*theta))./(4.*stdev_py.^2) ).*x.*y) ...
                + ( ( (sin(theta).^2)./(2.*stdev_px.^2) + (cos(theta).^2)./(2.*stdev_py.^2) ).*y.^2) ) );
            
            % evaluate function
            pD = p2(sbmx,sbmy);
            
            % shape of windowing for manifold
            window_type = 'circle';
            
        case 3
            
            % weird distribution
            [pD,~] = weird_pdf1(sbmx,sbmy);
            
            % shape of windowing for manifold
            window_type = 'circle';
            
        case 4
           
            % weird distribution
            [pD,~] = weird_pdf2(sbmx,sbmy);
            
            % shape of windowing for manifold
            window_type = 'circle';
            
        case 5
            
            % leptokurtotic generalized normal
            stdev_px    = 0.75;
            stdev_py    = 0.25;
            pwr         = 1.1;
            
            p2 = @(x,y) exp( -( ((abs(x).^pwr)/(stdev_px^pwr)) + ((abs(y).^pwr)/(stdev_py^pwr)) ) );
            
            
            % evaluate function
            pD = p2(sbmx,sbmy);
            
            % shape of windowing for manifold
            window_type = 'circle';
            
        otherwise
            
    end
    
    % normalize probability
    pD = pD./sum(pD(:));
    
    % separable approximation
    [u,s,v] = svd(pD);
    p_s1    = v(:,1);
    p_s2    = u(:,1);
    psep    = p_s2*p_s1';
    
    
    %% Solve for manifolds
    
    % determine manifold potential
    [F,Fx,Fy,Fxi,Fyi,pU,fhand] = compute_uniform_transform(nbin,pD,window_type);
    
    % save diagnostic figure
    saveas(fhand,['panels/Examples/p' num2str(p) '_optimization_result.eps'],'epsc')
    close(fhand);
    
    % determine manifold potential for separable approximation
    [Fsep,Fxsep,Fysep,Fxisep,Fyisep,pUsep,fhandsep] = compute_uniform_transform(nbin,psep,window_type);
    
    % save diagnostic figure
    saveas(fhandsep,['panels/Examples/p' num2str(p) '_optimization_result_sep.eps'],'epsc')
    close(fhandsep);
    
    % save results for loading in next step
    save(['panels/Examples/p' num2str(p) 'data.mat']);
    
end



