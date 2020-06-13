% This script plots the results for a specified example distribution from
% ExampleFigPrep script
clear all; close all;

%% Select example
Example = 1;

load(['panels/Examples/p' num2str(Example) 'data.mat']);


%% Define parameters

% stimulus space
stim_range   = [-1 1];          % stimulus space is centered on 0 and goes from -1 - 1
axis_range   = [-0.65 0.65];    % crop edges to remove boundary artifacts

% neural population
stdev        = 0.05;        % of tuning curves
spacing      = 0.2;         % spacing between neuron centers in s1

% tuning curve shape
h2 = @(x,y) exp(-0.5 * (x.^2 + y.^2)./(stdev^2));

%% Set up plotting

% colormap
pcmap = cool(128);
rcmap = parula(128);

% set up figures
figProb = figure; axProb = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]); box on;
axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(pcmap);

figProbSep = figure; axProbSep = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]); box on;
axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(pcmap);

figDnsty = figure; axDnsty = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]); box on;
axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(pcmap);

figMani = figure; axMani = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]);
box on; axis square; set(gca,'linewidth',3,'DefaultLineLineWidth', 1.5); view(axMani,-45,25); colormap gray;

figUni = figure; axUni = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]); colormap(rcmap);
axis([axis_range axis_range]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); view(axUni,-11,84);

figWarpTilt = figure; axWarpTilt = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]); axis([axis_range axis_range 0 1.2]);
box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); view(axWarpTilt,-5,86); colormap(rcmap);

figWarp = figure; axWarp = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]);axis([axis_range axis_range 0 1.2]);
box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(rcmap);

figWarpSep = figure; axWarpSep = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]); axis([axis_range axis_range 0 1.2]);
box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(rcmap);

figArrow = figure; axArrow = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]);
box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);

figSlice1 = figure; axSlice1 = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]);
axis([axis_range(1) axis_range(2) 0 3 0 1.2]); box on; axis square; set(gca,'linewidth',1, 'DefaultLineLineWidth', 1.5); view(axSlice1,-1,1); colormap(rcmap);

figSlice2 = figure; axSlice2 = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]);
axis([axis_range(1) axis_range(2) 0 3 0 1.2]); box on; axis square; set(gca,'linewidth',1, 'DefaultLineLineWidth', 1.5); view(axSlice2,-1,1); colormap(rcmap);



%% Stimulus plots

% plot probability of stimulus
ha = surf(axProb,sbmx,sbmy,pD);
set(ha,'edgecolor','none');
set(ha,'FaceAlpha',1)

% plot optimal density (proporational to probability)
ha = surf(axDnsty,sbmx,sbmy,pD);
set(ha,'edgecolor','none');
set(ha,'FaceAlpha',1)

% plot separable approximation of probability
ha = surf(axProbSep,sbmx,sbmy,psep);
set(ha,'edgecolor','none');
set(ha,'FaceAlpha',1)

% Plot the manifold potential
ha = surf(axMani,sbmx,sbmy,F);
set(ha,'edgealpha', 0.1);

%% Neuronal populations

% generate hexagonal lattice for neurons
[ms1,ms2] = generate_hexagonal_lattice(spacing,stim_range);

% uniform population
for n = 1:numel(ms1)
    
    curve2_n  = h2(sbmx-ms1(n),sbmy-ms2(n));
    
    % plot population in tilted view
    ha = surf(axUni,sbmx,sbmy,curve2_n);
    set(ha,'edgecolor','none');
    
end

% create down-sampled lattice to plot displacement arrows
downs   = 20;
tbmx    = sbmx(1:downs:end,1:downs:end);
tbmy    = sbmy(1:downs:end,1:downs:end);

tFxi    = Fxi(1:downs:end,1:downs:end);
tFyi    = Fyi(1:downs:end,1:downs:end);

X       = [tbmx(:) tbmy(:)]';
V       = [tFxi(:) tFyi(:)]';

% This is the INVERSE of the displacement field, which is more intuitive for arrows
figure(figArrow)
hq = equiver(X,V,'k-','linewidth',2);

% coordinates for two horizontal slices through the population
slices = [ceil(nbin*0.5) ceil(nbin*0.68)];

% Plot the warped tuning curves in 2D and 1D slices
for n = 1:numel(ms1)
    
    curve2_warped_n = h2(sbmx + Fx - ms1(n),sbmy + Fy - ms2(n));
    
    % plot population tilted
    ha = surf(axWarpTilt,sbmx,sbmy,curve2_warped_n);
    set(ha,'edgecolor','none');
    
    % plot flat
    ha = surf(axWarp,sbmx,sbmy,curve2_warped_n);
    set(ha,'edgecolor','none');
    
    % plot 2 horizontal 1D slices
    
    % if this neuron responds to stimuli in this slice, plot it's 1-D tuning curve
    if max(curve2_warped_n(slices(1),:)) >= 0.2
        ha = surf(axSlice1,[sbmx(1,:) ; sbmx(1,:)],[ones(size(sbmx(1,:))) ; 2*ones(size(sbmx(1,:)))],[curve2_warped_n(slices(1),:) ; curve2_warped_n(slices(1),:)]);
        set(ha,'edgecolor','none');
    end
    
    % same for second slice
    if max(curve2_warped_n(slices(2),:)) >= 0.2
        ha = surf(axSlice2,[sbmx(1,:) ; sbmx(1,:)],[ones(size(sbmx(1,:))) ; 2*ones(size(sbmx(1,:)))],[curve2_warped_n(slices(2),:) ; curve2_warped_n(slices(2),:)]);
        set(ha,'edgecolor','none');
    end
    
    % also plot tuning curves for separable approximation
    curve2_warped_n_sep = h2(sbmx + Fxsep - ms1(n),sbmy + Fysep - ms2(n));
    
    ha = surf(axWarpSep,sbmx,sbmy,curve2_warped_n_sep);
    set(ha,'edgecolor','none');
    
end



% save panels
saveas(figProb,['panels/Examples/fig' num2str(Example) 'a.eps'],'epsc')
saveas(figDnsty,['panels/Examples/fig' num2str(Example) 'c.eps'],'epsc')
saveas(figProbSep,['panels/Examples/fig' num2str(Example) 'a2.eps'],'epsc')
saveas(figMani,['panels/Examples/fig' num2str(Example) 'f.eps'],'epsc')
saveas(figUni,['panels/Examples/fig' num2str(Example) 'b.eps'],'epsc')
saveas(figArrow,['panels/Examples/fig' num2str(Example) 'e.eps'],'epsc')
saveas(figWarpTilt,['panels/Examples/fig' num2str(Example) 'd.eps'],'epsc')
saveas(figWarp,['panels/Examples/fig' num2str(Example) 'd2.eps'],'epsc')
saveas(figSlice1,['panels/Examples/fig' num2str(Example) 'g.eps'],'epsc')
saveas(figSlice2,['panels/Examples/fig' num2str(Example) 'h.eps'],'epsc')
saveas(figWarpSep,['panels/Examples/fig' num2str(Example) 'd2_sep.eps'],'epsc')







