% This script shows an example of warping a uniform 1D neural population to an efficient population
% following the method described in:
%
% Efficient sensory encoding and Bayesian inference with heterogeneous neural populations.  
% Ganguli D, Simoncelli EP.  
% Neural computation. 2014;26(10):2103?2134

clear all; close all;

%% Define Parameters

% stimulus space
stim_range = [0 1];     % stimulus space is centered on 0.5 and goes from 0-1
axis_range = [0.2 0.8]; % crop edges to remove boundary artifacts

% number of bins for probability distribution
nbin = 4001;

% neural population
stdev      = 0.05/2;    % of tuning curves
spacing    = 0.2/2;     % spacing between tuning curves


%% Set up plotting

% colormaps
rcmap = parula(128);
pcmap = flipud(autumn(128));
dcmap = cool(128);

% figure panels
figProb     = figure; axProb = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 1.1]); box on; axis square; set(gca,'linewidth',3,'DefaultLineLineWidth', 1.5);
figUni      = figure; axUni = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 1.1]); box on; axis square; set(gca,'linewidth',3,'DefaultLineLineWidth', 1.5);
figCml      = figure; axCml = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 1]); box on; axis square; set(gca,'linewidth',3,'DefaultLineLineWidth', 1.5);
figWarp     = figure; axWarp = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 1.1]); box on; axis square; set(gca,'linewidth',3,'DefaultLineLineWidth', 1.5);
figDnsty    = figure; axDnsty = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 1.1]); box on; axis square; set(gca,'linewidth',3,'DefaultLineLineWidth', 1.5);


%% Define space and warping

% sample stimulus space
s  = linspace(stim_range(1), stim_range(2), nbin+1);
s  = s(1:end-1) + diff(s(1:2))/2; 

% spacing of neurons
m = stim_range(1):spacing:stim_range(2);

% stimulus probability
p = @(x) exp(-0.5 * ((x-0.5).^2)./(0.1^2));

% cumulative of stimulus probability
pc = @(x) 0.5*(1 + erf((x-0.5)/(sqrt(2)*0.1)));

% plot stimulus probability, cumulative, and optimal density
plot(axProb, s, p(s),'k-');    % probability
plot(axCml, s, pc(s),'k-');   % cumulative
plot(axDnsty, s, p(s),'k-');    % density (= probability)

%% Warp population

% tuning curve shape
h = @(x) exp(-0.5 * (x.^2)./(stdev^2));

% for each neuron
for n = 1:numel(m)
    
    % nth tuning curve
    curve_n  = h(s-m(n));
    
    % compute warped tuning curve
    curve_warped_n = h(pc(s) - m(n));
    
    % plot uniform population
    plot(axUni, s, curve_n,'k-');

    % plot warped population
    plot(axWarp, s, curve_warped_n,'k-');

end

%% Save plots
saveas(figProb,'panels/Summary/fig1a.eps','epsc')
saveas(figUni,'panels/Summary/fig1b.eps','epsc')
saveas(figCml,'panels/Summary/fig1c.eps','epsc')
saveas(figWarp,'panels/Summary/fig1d.eps','epsc')
saveas(figDnsty,'panels/Summary/fig1e.eps','epsc')
