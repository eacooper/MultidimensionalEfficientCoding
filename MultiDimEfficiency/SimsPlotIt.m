% This script examines the correlation between tuning curve properties and
% stimulus probably for the large sample of simulated probability
% distributions generated in SimsPrep
clear all; close all;

% files to load
listing = dir('./sims/correlations/data/*.mat');
n_sims = numel(listing);

%% Load  each simulation

for n = 1:n_sims
    
    display(num2str(n));
    
    % load data
    load(['./sims/correlations/data/' listing(n).name]);
    
    % grab correlations
    r_gain_all(n) = r_gain;
    r_width_all(n) = r_width;

end

fig1 = figure; hold on;
h1 = histogram(r_gain_all,linspace(-1,1,21));
set(gca,'ytick',[]); set(gca,'xtick',[0],'xticklabel',[]); axis([-1 1 0 n_sims]); box on; axis square; set(gca,'linewidth',1);
set(h1,'FaceColor',[0.5 0.5 0.5]);
set(h1,'FaceAlpha',1)

fig2 = figure; hold on;
h2 = histogram(r_width_all,linspace(-1,1,21));
set(gca,'ytick',[]); set(gca,'xtick',[0],'xticklabel',[]); axis([-1 1 0 n_sims]); box on; axis square; set(gca,'linewidth',1);
set(h2,'FaceColor',[0.5 0.5 0.5]);
set(h2,'FaceAlpha',1)

if ~exist('panels/CorrelationSims/'); mkdir('panels/CorrelationSims/'); end
saveas(fig1,['panels/CorrelationSims/figprob_gain_1D.eps'],'epsc')
saveas(fig2,['panels/CorrelationSims/figprob_width_1D.eps'],'epsc')

% stats
[p_gain,h_gain,stats_gain] = signrank(r_gain_all);
[p_width,h_width,stats_width] = signrank(r_width_all);

m_gain = median(r_gain_all);
m_width = median(r_width_all);

display(['median gain = ' num2str(m_gain,4)]);
display(['median width = ' num2str(m_width,4)]);

display(['p gain = ' num2str(p_gain,4)]);
display(['p width = ' num2str(p_width,4)]);

