% This script examines the correlation between tuning curve properties and
% stimulus probably for example neuronal populations
clear all; close all;

%% Select example

Example = 5;

% load data
load(['panels/Examples/p' num2str(Example) 'data.mat']);

% scale probability from 0-1
pD = pD./max(pD(:));

%% Define parameters

% stimulus space
stim_range   = [-1 1];          % stimulus space is centered on 0 and goes from -1 - 1
axis_range   = [-0.55 0.55];    % crop edges to remove boundary artifacts

% neural population
stdev     = 0.05; % of tuning curves
n_neurons = 20; % number of neurons in each dimension
m         = linspace(stim_range(1), stim_range(2), n_neurons); % neuronal tuning preferences
h2        = @(x,y) exp(-0.5 * (x.^2 + y.^2)./(stdev^2)); % tuning curve shape

%% Set up plotting

% colormap
rcmap = parula(128);

figHori = figure; axHori = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]); axis([axis_range axis_range 0 1.2]); box on; axis square; set(gca,'linewidth',3); colormap(rcmap);

figVert = figure; axVert = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]); axis([axis_range axis_range 0 1.2]); box on; axis square; set(gca,'linewidth',3); colormap(rcmap);

figGain = figure; axGain = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([-0.1 1.1 0.05 1.15]); box on; axis square; set(gca,'linewidth',3);

figBandW = figure; axBandW = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([-0.1 1.1 0 50]); box on; axis square; set(gca,'linewidth',3);

figBandW1D = figure; axBandW1D = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([-0.1 1.1 0 50]); box on; axis square; set(gca,'linewidth',3);

fig1Ds1 = figure; ax1Ds1 = axes(); hold on; box on; axis square; set(gca,'linewidth',3);
fig1Ds2 = figure; ax1Ds2 = axes(); hold on; box on; axis square; set(gca,'linewidth',3);


%% Perform sampling

% 1D slices to take
num_slices  = 31;
slices      = find(abs(sbx) <= abs(axis_range(1)));
slices      = floor(linspace(slices(1),slices(end),num_slices));
slicesV     = slices + sign(rand(1,31)).*randi(5,1,num_slices); % vertical
slicesH     = slices + sign(rand(1,31)).*randi(5,1,num_slices); % horizontal


% initialize counters
cntV = 1;
cntH = 1;

% Warped tuning curves
for n1 = 1:n_neurons
    
    for n2 = 1:n_neurons
        
        % slightly randomize neuron peak
        mn1 = m(n1) + 0.5*randn;
        mn2 = m(n2) + 0.5*randn;
        
        % warped curve
        curve2_warped_n = h2(sbmx + Fx - mn1,sbmy + Fy - mn2);
        
        % plot the full curve - will overlay 1D prefs horizontal
        ha = surf(axHori,sbmx,sbmy,curve2_warped_n);
        set(ha,'edgecolor','none');
        
        % plot the full curve - will overlay 1D prefs vertical
        ha = surf(axVert,sbmx,sbmy,curve2_warped_n);
        set(ha,'edgecolor','none');
        
        %% Simulated 1-D measurements
        
        % for each y slice (vertical)
        for ss = 1:numel(slicesV)
            
            % plot vertical line
            if n1 == 1 && n2 == 1
                plot3(axVert,[sbx(slicesV(ss)) sbx(slicesV(ss))],stim_range,[1.1 1.1],'r-');
            end
            
            % get this slice
            this_slice = curve2_warped_n(:,slicesV(ss));
  
            % what is the max response (gain) in this slice
            gainVtmp = max(this_slice);
            
            % what is the tuning width in this slice
            ytmp        = sby(this_slice >= gainVtmp*0.5);
            widthVtmp   = abs(ytmp(1)-ytmp(end));
              
            % what is the preferred stimulus in this slice
            tmp      = find(this_slice == gainVtmp);
            prefVtmp = tmp(1);
            
            % what is the probability at the max response in this slice
            probVtmp = pD(prefVtmp, slicesV(ss));
            
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
                gainV(cntV)  = gainVtmp;
                prefV(cntV)  = prefVtmp;
                widthV(cntV) = widthVtmp;
                probV(cntV)  = probVtmp;
                
                % plot peak and bandwidth
                plot3(axVert,sbx(slicesV(ss)),sby(prefV(cntV)),1.2,'rs');
                plot3(axVert,[sbx(slicesV(ss)) sbx(slicesV(ss))],[ztmp(1) ztmp(end)],[1.2 1.2],'k-','linewidth',3);
                
                cntV = cntV + 1;
            end
            
        end
        
        
        
        % for each x slice (horizontal)
        for ss = 1:numel(slicesH)
            
            % plot horizontal line
            if n1 == 1 && n2 == 1
                plot3(axHori,stim_range,[sby(slicesH(ss)) sby(slicesH(ss))],[1.1 1.1],'k-');
            end
            
            % get this slice
            this_slice = curve2_warped_n(slicesH(ss),:);

            % what is the max response (gain) in this slice
            gainHtmp = max(this_slice);
            
            % what is the tuning width in this slice
            xtmp      = sbx(this_slice >= gainHtmp*0.5);
            widthHtmp = abs(xtmp(1)-xtmp(end));
                
            % what is the preferred stimulus in this slice
            tmp      = find(this_slice == gainHtmp);
            prefHtmp = tmp(1);
            
            % what is the probability at the max response in this slice
            probHtmp = pD(slicesH(ss), prefHtmp);
            
            % if this neuron responds in this slice & doesn't exceed the axis boundary, store the results
            if max(this_slice) >= 0.2 && ~any(abs(xtmp) > abs(axis_range(1)))
                
                % more precise width estimate by linear interpolation
                curve_x   = sbx(this_slice >= gainHtmp*0.1);
                curve_y   = this_slice(this_slice >= gainHtmp*0.1);
                curve_xH  = linspace(min(sbx(this_slice >= gainHtmp*0.1)),max(sbx(this_slice >= gainHtmp*0.1)),50);
                curve_yH  = interp1(curve_x,curve_y,curve_xH);
                ztmp      = curve_xH(curve_yH >= gainHtmp*0.5);
                widthHtmp = abs(ztmp(1)-ztmp(end));
                
                % store variables
                gainH(cntH)  = gainHtmp;
                prefH(cntH)  = prefHtmp;
                widthH(cntH) = widthHtmp;
                probH(cntH)  = probHtmp;
                
                % plot peak and bandwidth
                plot3(axHori,sbx(prefH(cntH)),sby(slicesH(ss)),1.1,'ko');
                plot3(axHori,[ztmp(1) ztmp(end)],[sby(slicesH(ss)) sby(slicesH(ss))],[1.1 1.1],'r-','linewidth',2);
                
                cntH = cntH + 1;
            end
            

        end
        
    end
    
end

display(['num neurons vertical = ' num2str(numel(gainV))]);
display(['num neurons horizontal = ' num2str(numel(gainH))]);

%% correlation plots

% gain plot
scatter(axGain,probV,gainV,'ro','filled');
scatter(axGain,probH,gainH,'ko','filled');

% bandwidth plot
scatter(axBandW,probV,1./widthV,'ro','filled');
scatter(axBandW,probH,1./widthH,'ko','filled');

%% for sanity check, calculate the predictions for marginalized probabilities and 1-D populations

% marginalized 1D probabilities
pDs1 = sum(pD,1);
pDs1 = pDs1/max(pDs1);

pDs2 = sum(pD,2);
pDs2 = pDs2/max(pDs2);

% cumulatives
pDs1C = cumsum(pDs1)/sum(pDs1);
pDs2C = cumsum(pDs2)/sum(pDs2);

% 1D tuning curve shape and spacing
h     = @(x) exp(-0.5 * (x.^2)./(stdev^2));
m1D   = linspace(stim_range(1), stim_range(2), n_neurons);

% measure bandwidth for 1-D population
for n = 1:n_neurons
    
    % nth tuning curve
    curve_n  = h(sbx-m1D(n));
    
    % s1
    
    % compute warped tuning curve
    curve_warped_n = h((pDs1C-0.5)*2 - m1D(n));
    
    % preferred stim
    tmp          = sbx(curve_warped_n == max(curve_warped_n));
    pref1Ds1(n)  = mean(tmp);
    pref1Dcoords = find(abs(sbx-pref1Ds1(n)) == min(abs(sbx-pref1Ds1(n))));
    pref1Dcoords = min(pref1Dcoords);
    
    % fwhm
    tmp = sbx(curve_warped_n >= 0.5);
    if ~isempty(tmp)
        width1Ds1(n) = abs(tmp(1)-tmp(end));
    else
        width1Ds1(n) = NaN;
    end
    
    % probability at preference
    prob1Ds1(n) = pDs1(pref1Dcoords);
    
    % plot warped population
    plot(ax1Ds1, sbx, curve_warped_n);
    
    % s2
    
    % compute warped tuning curve
    curve_warped_n = h((pDs2C-0.5)*2 - m1D(n));
    
    % preferred stim
    tmp          = sby(curve_warped_n == max(curve_warped_n));
    pref1Ds2(n)  = mean(tmp);
    pref1Dcoords = find(abs(sby-pref1Ds2(n)) == min(abs(sby-pref1Ds2(n))));
    pref1Dcoords = min(pref1Dcoords);
    
    % fwhm
    tmp = sby(curve_warped_n >= 0.5);
    if ~isempty(tmp)
        width1Ds2(n) = abs(tmp(1)-tmp(end));
    else
        width1Ds2(n) = NaN;
    end
    
    % probability at preference
    prob1Ds2(n) = pDs2(pref1Dcoords);
    
    % plot warped population
    plot(ax1Ds2, sby, curve_warped_n);
    
end



%% 1D correlation plotplots

% gain plot
scatter(axBandW1D,prob1Ds1,1./width1Ds1,'ro','filled');
scatter(axBandW1D,prob1Ds2,1./width1Ds2,'ko','filled');


%% save plots
if ~exist('panels/CorrelationAnalysis/'); mkdir('panels/CorrelationAnalysis/'); end

saveas(figGain,['panels/CorrelationAnalysis/fig' num2str(Example) '_prob_gain_1D.eps'],'epsc')
saveas(figBandW,['panels/CorrelationAnalysis/fig' num2str(Example) '_prob_width_1D.eps'],'epsc')
saveas(figHori,['panels/CorrelationAnalysis/fig' num2str(Example) 'b.eps'],'epsc')
saveas(figVert,['panels/CorrelationAnalysis/fig' num2str(Example) 'c.eps'],'epsc')
saveas(figBandW1D,['panels/CorrelationAnalysis/fig' num2str(Example) '_prob_width_1Dmarg.eps'],'epsc')

saveas(fig1Ds1,['panels/CorrelationAnalysis/fig' num2str(Example) '_marg_s1.eps'],'epsc')
saveas(fig1Ds2,['panels/CorrelationAnalysis/fig' num2str(Example) '_marg_s2.eps'],'epsc')