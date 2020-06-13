% This script shows how to visualize the allocation of neural encoding resources by calculating the determinant of the Fisher
% information matrix for example 2D populations
clear all; close all;

%% Select example

Example = 5;

% load manifold potential for this panel
load(['panels/Examples/p' num2str(Example) 'data.mat']);

% flag for whether to increase population sampling by a factor of 2 in each dimension
% (increases FIM and creates a better estimate of the stimulus probability based on FIM)
scale_up = 1;

%% Define parameters

% stimulus space
stim_range   = [-1 1];          % stimulus space is centered on 0 and goes from -1 - 1
axis_range   = [-0.65 0.65];    % crop edges to remove boundary artifacts

% tuning curve parameters
stdev        = 0.05; % of tuning curves
spacing      = 0.2; % spacing between neuron centers

if scale_up
    spacing = spacing/2;
end

% tuning curve shape
h2 = @(x,y) exp(-0.5 * (x.^2 + y.^2)./(stdev^2));

%% Set up plotting

pcmap = cool(128);
rcmap = parula(128);
fmap  = copper(128);

% set up figures
figFish     = figure; axFish = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(fmap);
figFishSep  = figure; axFishSep = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(fmap);
figProb     = figure; axProb = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(pcmap);
figProbSep  = figure; axProbSep = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(pcmap);
figProbEst  = figure; axProbEst = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(pcmap);
figPop      = figure; axPop = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]); axis([axis_range axis_range 0 1.2]); colormap(rcmap); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figPopSep   = figure; axPopSep = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]); axis([axis_range axis_range 0 1.2]); colormap(rcmap); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figPopCum   = figure; axPopCum = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]); axis([axis_range axis_range 0 1.2]); colormap(rcmap); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figFishCum  = figure; axFishCum = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]); box on;axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(fmap);
figDnstySample = figure; axDnstySample = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(pcmap);


%% Probability plots

% plot stimulus probability
ha = surf(axProb,sbmx,sbmy,pD);
set(ha,'edgecolor','none');
set(ha,'FaceAlpha',1)

% plot separable approximation of stimulus probability
ha = surf(axProbSep,sbmx,sbmy,psep);
set(ha,'edgecolor','none');
set(ha,'FaceAlpha',1)

% plot marginals of stimulus probability
figMarg = figure; hold on;
plot(sbx,p_s1,'k');
plot(sby,p_s2,'r');

%% Neural population

% generate hexagonal lattice for neurons
[ms1,ms2] = generate_hexagonal_lattice(spacing,stim_range);


% initialize fisher matrix for full and separable approximation populations
fisher      = zeros(201,201,2,2);
fisher_sep  = zeros(201,201,2,2);

% Plot the warped tuning curves
for n = 1:numel(ms1)
    
    % full optimization
    curve2_warped_n = h2(sbmx + Fx - ms1(n),sbmy + Fy - ms2(n));
    
    % calculate fisher information
    fisher_n = single_neuron_FIM(curve2_warped_n);
    
    % add the FIM from each neuron
    fisher   = fisher + fisher_n;
    
    % plot tuning curves
    ha = surf(axPop,sbmx,sbmy,curve2_warped_n);
    set(ha,'edgecolor','none');
    
    % also plot for separable approximation
    curve2_warped_n_sep = h2(sbmx + Fxsep - ms1(n),sbmy + Fysep - ms2(n));
    
    % calculate fisher information
    fisher_n_sep = single_neuron_FIM(curve2_warped_n_sep);
    
    % add the FIM from each neuron
    fisher_sep   = fisher_sep + fisher_n_sep;
    
    % plot tuning curves
    ha = surf(axPopSep,sbmx,sbmy,curve2_warped_n_sep);
    set(ha,'edgecolor','none');
    
end

%% Fisher Information

% compute the determinant of the FIM
for s1 = 1:size(fisher,2)
    
    for s2 = 1:size(fisher,1)
        
        fisherD(s1,s2)      = det(squeeze(fisher(s1,s2,:,:)));
        fisherD_sep(s1,s2)  = det(squeeze(fisher_sep(s1,s2,:,:)));
        
    end
    
end

% determine color axis range (max within window without border artifacts)
fisher_max = max(fisherD(abs(sbmx(1:end,1:end))<=abs(axis_range(1)) & abs(sbmy(1:end,1:end))<=abs(axis_range(1))));
fisher_sep_max = max(fisherD_sep(abs(sbmx(1:end,1:end))<=abs(axis_range(1)) & abs(sbmy(1:end,1:end))<=abs(axis_range(1))));

display(['FIM min = 0 ; FIM max = ' num2str(fisher_max)]);


% plot determinant of FIM
ha = surf(axFish,sbmx,sbmy,fisherD);
set(ha,'edgecolor','none');
set(ha,'FaceAlpha',1);
caxis(axFish,[0 fisher_max]);

ha = surf(axFishSep,sbmx,sbmy,fisherD_sep);
set(ha,'edgecolor','none');
set(ha,'FaceAlpha',1);
caxis(axFishSep,[0 fisher_sep_max]);

% plot recovered estimate of stimulus probability
ha = surf(axProbEst,sbmx,sbmy,sqrt(fisherD));
set(ha,'edgecolor','none');
set(ha,'FaceAlpha',1)
caxis(axProbEst,[0 sqrt(fisher_max)]);

%% Density mapping

% to confirm density mapping, just select a bunch of points with uniform density in stimulus space,
% warp, and plot new density (should be proportionate to probability)

points_x = randsample(1:201,5000,true);
points_y = randsample(1:201,5000,true);

points_x_warped = sbmx(points_x,points_y) + Fxi(points_x,points_y);
points_y_warped = sbmy(points_x,points_y) + Fyi(points_x,points_y);

[N,Xedges,Yedges] = histcounts2(points_x_warped(:),points_y_warped(:),51);

N = N'; % hist counts treats x as vertical

ha = surf(axDnstySample,Xedges(1:end-1),Yedges(1:end-1),N);
set(ha,'edgecolor','none');
set(ha,'FaceAlpha',1)

%% 1-D approximation

% perform 1-D separable calculation to see how close a match this is for optimizing density/FI

% 1-D normalized coords
s1D = (sbx - min(sbx))/range(sbx);

% 1-D cumulatives
p_s1C = cumsum(p_s1)/sum(p_s1);
p_s2C = cumsum(p_s2)/sum(p_s2);

% show these cumulatives
figMargCum = figure; hold on;
plot(s1D,p_s1C,'k');
plot(s1D,p_s2C,'r');

% make 2D
pDs1C = repmat(p_s1C',201,1);
pDs2C = repmat(p_s2C,1,201);

% scale to align with stimulus space (-1 to 1)
pDs1C_scaled = ((2*pDs1C)-1);
pDs2C_scaled = ((2*pDs2C)-1);

% 2D cumulative
pDs12C = p_s2C*p_s1C';

% initialize fisher info matrix
fisher_cum = zeros(201,201,2,2);

% Plot the warped tuning curves based on separable calculation
for n = 1:numel(ms1)
    
    % compute warped tuning curve
    curve2_warped_n = h2(pDs1C_scaled - ms1(n),pDs2C_scaled - ms2(n));
    
    % plot
    ha = surf(axPopCum,sbmx,sbmy,curve2_warped_n);
    set(ha,'edgecolor','none');
    
    % calculate fisher information
    fisher_n = single_neuron_FIM(curve2_warped_n);
    
    % add the FIM from each neuron
    fisher_cum   = fisher_cum + fisher_n;
    
end

% compute the determinant of the FIM
for s1 = 1:size(fisher_cum,2)
    
    for s2 = 1:size(fisher_cum,1)
        
        fisherD_cum(s1,s2) = det(squeeze(fisher_cum(s1,s2,:,:)));
        
    end
    
end

% plot determinant of FIM for calculation using cumulative
ha = surf(axFishCum,sbmx(1:end,1:end),sbmy(1:end,1:end),fisherD_cum);
set(ha,'edgecolor','none');
set(ha,'FaceAlpha',1)
caxis(axFishCum,[0 ...
    max(fisherD_cum(abs(sbmx(1:end,1:end))<=abs(axis_range(1)) & abs(sbmy(1:end,1:end))<=abs(axis_range(1))))]);

%% save panels

if ~scale_up
    saveas(figProb,['panels/Fisher/fig' num2str(Example) '_prob.eps'],'epsc')
    saveas(figProbSep,['panels/Fisher/fig' num2str(Example) '_probsep.eps'],'epsc')
    saveas(figProbEst,['panels/Fisher/fig' num2str(Example) '_probest.eps'],'epsc')
    
    saveas(figFish,['panels/Fisher/fig' num2str(Example) '_FIMfull_maxval' num2str(fisher_max,3) '.eps'],'epsc')
    saveas(figFishSep,['panels/Fisher/fig' num2str(Example) '_FIMsep.eps'],'epsc')
    saveas(figFishCum,['panels/Fisher/fig' num2str(Example) '_FIMcum.eps'],'epsc')
    
    saveas(figMarg,['panels/Fisher/fig' num2str(Example) '_1Dmarginals.eps'],'epsc')
    saveas(figMargCum,['panels/Fisher/fig' num2str(Example) '_1Dcumulatives.eps'],'epsc')
    
    saveas(figPop,['panels/Fisher/fig' num2str(Example) '_neuronsfull.eps'],'epsc')
    saveas(figPopSep,['panels/Fisher/fig' num2str(Example) '_neuronssep.eps'],'epsc')
    saveas(figPopCum,['panels/Fisher/fig' num2str(Example) '_neuronscum.eps'],'epsc')
    
    saveas(figDnstySample,['panels/Fisher/fig' num2str(Example) '_densityfull.eps'],'epsc')
    
else
    saveas(figProb,['panels/Fisher/fig' num2str(Example) '_prob_scaledup.eps'],'epsc')
    saveas(figProbSep,['panels/Fisher/fig' num2str(Example) '_probsep_scaledup.eps'],'epsc')
    saveas(figProbEst,['panels/Fisher/fig' num2str(Example) '_probest_scaledup.eps'],'epsc')
    
    saveas(figFish,['panels/Fisher/fig' num2str(Example) '_FIMfull_scaledup_maxval' num2str(fisher_max,3) '.eps'],'epsc')
    saveas(figFishSep,['panels/Fisher/fig' num2str(Example) '_FIMsep_scaledup.eps'],'epsc')
    saveas(figFishCum,['panels/Fisher/fig' num2str(Example) '_FIMcum_scaledup.eps'],'epsc')
    
    saveas(figMarg,['panels/Fisher/fig' num2str(Example) '_1Dmarginals_scaledup.eps'],'epsc')
    saveas(figMargCum,['panels/Fisher/fig' num2str(Example) '_1Dcumulatives_scaledup.eps'],'epsc')
    
    saveas(figPop,['panels/Fisher/fig' num2str(Example) '_neuronsfull_scaledup.eps'],'epsc')
    saveas(figPopSep,['panels/Fisher/fig' num2str(Example) '_neuronssep_scaledup.eps'],'epsc')
    saveas(figPopCum,['panels/Fisher/fig' num2str(Example) '_neuronscum_scaledup.eps'],'epsc')
    
    saveas(figDnstySample,['panels/Fisher/fig' num2str(Example) '_densityfull_scaledup.eps'],'epsc')
    
end












