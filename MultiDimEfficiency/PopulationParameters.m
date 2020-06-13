% This script illustrates the parameterization of a heterogeneous neuronal
% population in terms of density and gain, as described in
%
% Efficient sensory encoding and Bayesian inference with heterogeneous neural populations.  
% Ganguli D, Simoncelli EP.  
% Neural computation. 2014;26(10):2103?2134
close all; clear all;

%%  Define Parameters

% stimulus space
stim_range  = [-12, 12];    % edges of stimulus space
axis_range  = [-8 8];       % crop edges to remove boundary artifacts

% number of bins for probability distribution
nbin = 4001;

% neural population
stdev       = 1.0;
n_neurons   = 13;

%% Set up plotting

figDispGood = figure; axDispGood = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range -1.2 1.2]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figDispBad  = figure; axDispBad = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range -1.2 1.2]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figGainGood = figure; axGainGood = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 1.2]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); 
figGainBad  = figure; axGainBad = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 1.2]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figUni      = figure; axUni = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 1.2]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figWarpGood1 = figure; axWarpGood1 = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 1.2]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figWarpBad1 = figure; axWarpBad1 = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 1.2]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figWarpGood2 = figure; axWarpGood2 = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 1.2]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figWarpBad2 = figure; axWarpBad2 = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 1.2]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figFishUni = figure; axFishUni = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 3]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figFishGood = figure; axFishGood = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 3]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);
figFishBad = figure; axFishBad = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range 0 10]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5);


%% Define space and warping

% Define axis and neuron positions. 
s   = linspace(stim_range(1), stim_range(2), nbin);
m   = linspace(stim_range(1), stim_range(2), n_neurons);

% tuning curve shape and derivative
h = @(x) exp(-0.5 * (x.^2)./(stdev^2));
h_p = @(x) -x.*exp(-0.5 * (x.^2)./(stdev^2));

% GOOD displacement field -- slowly varying
f_good = @(x) cos(x./2);
df_good = @(x) -0.5.*sin(x./2);

% GOOD gain function -- slowly varying
g_good = @(x) 0.4*cos(x./4) + 0.6;
g_p_good = @(x) -0.1.*sin(x./4);

% BAD displacement field -- steep gradients
f_bad = @(x) cos(1.5.*x);
df_bad = @(x) -1.5*sin(1.5.*x);

% BAD gain function -- steep gradients
g_bad = @(x) 0.4.* cos(5.*x) + 0.6;
g_p_bad = @(x) -2.*sin(5.*x);

% plot the displacement fields
plot(axDispGood, s, f_good(s),'k-');
plot(axDispGood,axis_range,[0 0],'k:');
plot(axDispBad, s, f_bad(s),'k-');
plot(axDispBad,axis_range,[0 0],'k:');

% plot the gain functions
plot(axGainGood, s, g_good(s),'k-');
plot(axGainBad, s, g_bad(s),'k-');


%% Compute Fisher information and plot neuronal populations

% Initialize Fisher Info
fisher              = zeros(size(s));
fisher_scaled_good  = zeros(size(s));
fisher_scaled_bad   = zeros(size(s));

for n = 1:n_neurons
    
    % for the uniform population, nth tuning curve
    curve_n  = h(s-m(n));
    deriv_n  = h_p(s-m(n));
    
    % corresponding fisher information
    fisher_n = deriv_n.^2 ./ curve_n;
    fisher   = fisher + fisher_n;
    
    % warped curves for GOOD and BAD displacement fields
    curve_warped_n_good = double(h(s + f_good(s) - m(n)));
    curve_warped_n_bad  = double(h(s + f_bad(s) - m(n)));
    
    % warped & scaled curves for GOOD and BAD displacement and gain
    curve_scaled_n_good = double(h(s + f_good(s) - m(n)).*g_good(s));
    curve_scaled_n_bad = double(h(s + f_bad(s) - m(n)).*g_bad(s));
    
    % Fisher information for warped & scales curves
    
    % good
    deriv_scaled_n_good = (1 + df_good(s)) .* h_p(s + f_good(s) - m(n)) .*g_good(s);
    deriv_scaled_n_good = deriv_scaled_n_good + h(s + f_good(s)-m(n)) .* g_p_good(s);
    
    fisher_scaled_n_good = deriv_scaled_n_good.^2 ./ curve_scaled_n_good;
    fisher_scaled_good   = fisher_scaled_good + fisher_scaled_n_good;
    
    % bad
    deriv_scaled_n_bad = (1 + df_bad(s)) .* h_p(s + f_bad(s) - m(n)) .*g_bad(s);
    deriv_scaled_n_bad = deriv_scaled_n_bad + h(s + f_bad(s)-m(n)) .* g_p_bad(s);
    
    fisher_scaled_n_bad = deriv_scaled_n_bad.^2 ./ curve_scaled_n_bad;
    fisher_scaled_bad   = fisher_scaled_bad + fisher_scaled_n_bad;
    
    % plot uniform population
    plot(axUni, s, curve_n);
    
    % plot both sets of warped curves curves
    plot(axWarpGood1, s, curve_warped_n_good);
    plot(axWarpBad1, s, curve_warped_n_bad);
    
    % plot both sets of warped & scaled curves
    plot(axWarpGood2, s, curve_scaled_n_good);
    plot(axWarpBad2, s, curve_scaled_n_bad);
    
end

% plot both sets of fisher information for the uniform populations
I_conv = mean(fisher(round(nbin*0.05):round(nbin*0.95)));
plot(axFishUni, s, fisher);
plot(axFishUni, s, ones(size(s)).*I_conv);

% plot fisher information for the warped & scaled populations
plot(axFishGood, s, fisher_scaled_good);
plot(axFishGood, s, ones(size(s)) .* I_conv.*(double((1+df_good(s))).^2) .*g_good(s));

plot(axFishBad, s, fisher_scaled_bad);
plot(axFishBad, s, ones(size(s)) .* I_conv.*(double((1+df_bad(s))).^2) .*g_bad(s));

%% Save plots
saveas(figDispGood,'panels/Parameters/fig1a.eps','epsc')
saveas(figDispBad,'panels/Parameters/fig1d.eps','epsc')
saveas(figWarpGood1,'panels/Parameters/fig1b.eps','epsc')
saveas(figWarpGood2,'panels/Parameters/fig1c.eps','epsc')
saveas(figWarpBad1,'panels/Parameters/fig1e.eps','epsc')
saveas(figWarpBad2,'panels/Parameters/fig1f.eps','epsc')
saveas(figGainGood,'panels/Parameters/fig1g.eps','epsc')
saveas(figGainBad,'panels/Parameters/fig1h.eps','epsc')
saveas(figUni,'panels/Parameters/fig1i.eps','epsc')
saveas(figFishUni,'panels/Parameters/fig1j.eps','epsc')
saveas(figFishGood,'panels/Parameters/fig1k.eps','epsc')
saveas(figFishBad,'panels/Parameters/fig1l.eps','epsc')









