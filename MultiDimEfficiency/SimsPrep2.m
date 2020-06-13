% This script examines the correlation between tuning curve properties and
% stimulus probably for the large sample of simulated probability
% distributions generated in SimsPrep
clear all; close all;

% files to load
listing = dir('./sims/manifolds/data/*.mat');
n_sims = numel(listing);

%% Define parameters

% stimulus space
stim_range   = [-1 1];          % stimulus space is centered on 0 and goes from -1 - 1
axis_range   = [-0.55 0.55];    % crop edges to remove boundary artifacts

% possible sigmas for initial tuning curves
stdev_minh = 0.03;
stdev_maxh = 0.07;
    
% number of neurons in each slice sample
n_neurons   = 25;

% neuronal tuning preferences
m   = linspace(stim_range(1), stim_range(2), n_neurons);
    
% flag for whether to plot each simulation
plot_it = 1;

%% Load and analyze each simulation

for n = 1:n_sims
    
    display(num2str(n));
    
    % load data
    load(['./sims/manifolds/data/' listing(n).name]);
    
    % scale probability from 0-1
    pD = pD./max(pD(:));
    
    % define uniform tuning curves
    stdev(n) = rand*(stdev_maxh - stdev_minh) + stdev_minh; % of tuning curves
    
    % tuning curve shape
    h2 = @(x,y) exp(-0.5 * (x.^2 + y.^2)./(stdev(n)^2));
    
    % random 1D slice to take, within axis range
    slices      = find(abs(sbx) <= abs(axis_range(1)));
    slice(n)    = randsample(slices,1);
    
    % which dimension to slice
    dim(n)      = randsample([1 2],1);

    % clear variables
    prob = []; pref = []; gain = []; width = [];
    
    % initiate figure
    figure; hold on;
    subplot(2,2,1); title('Warped Population');
    set(gca,'ytick',[],'xtick',[],'ztick',[]);axis([axis_range axis_range 0 1.2]); view(2);
    box on; axis square; colormap parula;
    subplot(2,2,2); title('1D tuning slices');
   
    
    % start counter
    cnt = 1;
    
    % Warped tuning curves
    for n1 = 1:n_neurons
        
        for n2 = 1:n_neurons
            
            % warped curve
            curve2_warped_n = h2(sbmx + Fx - m(n1),sbmy + Fy - m(n2));
            
            % plot every 4th iteration
            if ~mod(n1,4) && ~mod(n2,4)
                subplot(2,2,1); hold on;
                ha = surf(sbmx,sbmy,curve2_warped_n);
                set(ha,'edgecolor','none');
            end
            
            % Simulated 1-D measurements
            
            % sample differently depending on dimension
            if dim(n) == 1
                
                % x slice (horizontal)
                this_slice = curve2_warped_n(slice(n),:);
                
                % plot horizontal line
                if n1 == 1 && n2 == 1
                    subplot(2,2,1); hold on;
                    plot3(stim_range,[sby(slice(n)) sby(slice(n))],[1.1 1.1],'r-');
                end
                
                % what is the max response (gain) in this slice
                gainHtmp = max(this_slice);
                
                % what is the tuning width in this slice
                xtmp        = sbx(this_slice >= gainHtmp*0.5);
                widthHtmp   = abs(xtmp(1)-xtmp(end));
                
                % what is the preferred stimulus in this slice
                tmp         = find(this_slice == gainHtmp);
                prefHtmp    = tmp(1);
                
                % what is the probability at the max response in this slice
                probHtmp = pD(slice(n), prefHtmp);
                
                % if this neuron responds in this slice & doesn't exceed the axis boundary, store the results
                if max(this_slice) >= 0.2 && ~any(abs(xtmp) > abs(axis_range(1)))
                    
                    % more precise width estimate by linear interpolation
                    curve_x     = sbx(this_slice >= gainHtmp*0.1);
                    curve_y     = this_slice(this_slice >= gainHtmp*0.1);
                    curve_xH    = linspace(min(sbx(this_slice >= gainHtmp*0.1)),max(sbx(this_slice >= gainHtmp*0.1)),50);
                    curve_yH    = interp1(curve_x,curve_y,curve_xH);
                    ztmp        = curve_xH(curve_yH >= gainHtmp*0.5);
                    widthHtmp   = abs(ztmp(1)-ztmp(end));
                    
                    % store variables
                    gain(cnt)  = gainHtmp;
                    pref(cnt)  = prefHtmp;
                    width(cnt) = widthHtmp;
                    prob(cnt)  = probHtmp;
                    
                    % plot response
                    subplot(2,2,2); hold on;
                    plot(this_slice,'k-');
                    
                    cnt = cnt + 1;
                end
                
            elseif dim(n) == 2
                
                % y slice (vertical)
                this_slice = curve2_warped_n(:,slice(n));
                
                % plot vertical line
                if n1 == 1 && n2 == 1
                    subplot(2,2,1); hold on;
                    plot3([sbx(slice(n)) sbx(slice(n))],stim_range,[1.1 1.1],'r-');
                end
                
                % what is the max response (gain) in this slice
                gainVtmp = max(this_slice);
                
                % what is the tuning width in this slice
                ytmp        = sby(this_slice >= gainVtmp*0.5);
                widthVtmp   = abs(ytmp(1)-ytmp(end));
                
                % what is the preferred stimulus in this slice
                tmp         = find(this_slice == gainVtmp);
                prefVtmp    = tmp(1);
                
                % what is the probability at the max response in this slice
                probVtmp    = pD(prefVtmp, slice(n));
                
                % if this neuron responds in this slice & doesn't exceed the axis boundary, store the results
                if gainVtmp >= 0.2 && ~any(abs(ytmp) > abs(axis_range(1)))
                    
                    % more precise width estimate by linear interpolation
                    curve_x     = sby(this_slice >= gainVtmp*0.1);
                    curve_y     = this_slice(this_slice >= gainVtmp*0.1);
                    curve_xH    = linspace(min(curve_x),max(curve_x),50);
                    curve_yH    = interp1(curve_x,curve_y,curve_xH);
                    ztmp        = curve_xH(curve_yH >= gainVtmp*0.5);
                    widthVtmp   = abs(ztmp(1)-ztmp(end));
                    
                    % store variables
                    gain(cnt)  = gainVtmp;
                    pref(cnt)  = prefVtmp;
                    width(cnt) = widthVtmp;
                    prob(cnt)  = probVtmp;
                    
                    % plot response
                    subplot(2,2,2); hold on;
                    plot(this_slice,'k-');
                    
                    cnt = cnt + 1;
                end
                
            end
            
        end
        
    end
    
    % compute correlations if enough valid samples were found
    if cnt > 5
        r_gain = corr(prob',gain');
        r_width = corr(prob',1./width');
    else
        r_gain = NaN;
        r_width = NaN;
    end
    
    %% correlation plots
    
    % gain plot
    subplot(2,2,3); hold on;
    scatter(prob,gain,'ko','filled');
    xlabel('2D Probability'); ylabel('1D Gain'); box on;
    
    % bandwidth plot
    subplot(2,2,4); hold on;
    scatter(prob,1./width,'ko','filled');
    xlabel('2D Probability'); ylabel('1D Bandwidth'); box on;
    
    % save results
    if ~exist('sims/correlations/plots/'); mkdir('sims/correlations/plots/'); end
    if ~exist('sims/correlations/data/'); mkdir('sims/correlations/data/'); end
    
    saveas(gcf,['sims/correlations/plots/sim' num2str(n) '_sx' num2str(stdev_px,2) '_sy' num2str(stdev_py,2) '_th' num2str(round(theta*180/pi)) '_correlations.eps'],'epsc')
    close(gcf);
    save(['sims/correlations/data/sim' num2str(n) '_sx' num2str(stdev_px,2) '_sy' num2str(stdev_py,2) '_th' num2str(round(theta*180/pi)) '_data.mat'],'r_gain','r_width');
    
end

