% This script shows an example of warping a uniform 2D neural population to an efficient population
clear all; close all;

%% Define Parameters

% stimulus space
stim_range   = [-1 1];      % stimulus space is centered on 0 and goes from -1 - 1
axis_range   = [-0.6 0.6];  % crop edges to remove boundary artifacts

% neural population
% (stimulus space is 2x G&S space, so we multiply parameters from SummaryFig1D by 2 to match)
stdev        = 0.025*2;     % of tuning curves
spacing      = 0.1*2;       % spacing between tuning curves

% tuning curve shape
h2 = @(x,y) exp(-0.5 * (x.^2 + y.^2)./(stdev^2));

% number of bins for probability distribution in each dimension, must be odd
nbin = 201;


%% Set up plotting

% colormaps
rcmap = parula(128);
dcmap = cool(128);
pcmap = cool(128);

% figure panels
figProb     = figure; axProb = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(pcmap);
figDnsty    = figure; axDnsty = axes(); hold on; set(gca,'ytick',[],'xtick',[]); axis([axis_range axis_range]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(dcmap);
figUni      = figure; axUni = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]); axis([axis_range axis_range]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(rcmap);
figWarp2D   = figure; axWarp2D = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]); axis([axis_range axis_range]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); colormap(rcmap);
figArrow    = figure; axArrow = axes(); hold on; set(gca,'ytick',[],'xtick',[]);box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); view(axArrow,[-45 25]); axis equal; 
figMani     = figure; axMani = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); view(axMani,-45,25); colormap gray;
figWarp1D   = figure; axWarp1D = axes(); hold on; set(gca,'ytick',[],'xtick',[],'ztick',[]); axis([axis_range(1) axis_range(2) 0 3 0 1.2]); box on; axis square; set(gca,'linewidth',3, 'DefaultLineLineWidth', 1.5); view(axWarp1D,-1,1); colormap(rcmap);


%% Define space and warping

% sample stimulus space
sx = linspace(stim_range(1),stim_range(2),nbin+1);
sy = linspace(stim_range(2),stim_range(1),nbin+1);

% convert lattice points to bin centers
sbx  = sx(1:end-1) + diff(sx(1:2))/2;  % bin center x
sby  = sy(1:end-1) + diff(sy(1:2))/2;  % bin center y

% generate mesh 
[sbmx, sbmy] = meshgrid(sbx, sby);

% evaluate stimulus probability function over this lattice
[pD,~]  = weird_pdf2(sbmx,sbmy);
pD      = pD./sum(pD(:));

% plot stimulus probability
hs = surf(axProb,sbmx,sbmy,pD);
set(hs,'edgecolor','none','FaceAlpha',1)

% plot optimal density (proportional to probability)
hs = surf(axDnsty,sbmx,sbmy,pD);
set(hs,'edgecolor','none','FaceAlpha',1)

% determine manifold potential and displacement field
[F,Fx,Fy,Fxi,Fyi,pU,fapprox] = compute_uniform_transform(nbin,pD,'circle');

% Plot the manifold potential
hs = surf(axMani,sbmx,sbmy,F);
set(hs,'edgealpha', 0.1);

% save diagnostic figure from numerical approximation
saveas(fapprox,'./panels/Summary/optimization_result.eps','epsc')
close(fapprox);

%% Warp population

% generate hexagonal lattice for neurons
[ms1,ms2] = generate_hexagonal_lattice(spacing,stim_range);

% uniform population
for n = 1:numel(ms1)
    
    curve2_n  = h2(sbmx-ms1(n),sbmy-ms2(n));
    
    % plot it
    hs = surf(axUni,sbmx,sbmy,curve2_n);
    set(hs,'edgecolor','none');

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
hq = equiver(X,V,'k-');
axis([-0.6 0.6 -0.6 0.6 0 0.001]);

% warped tuning curves
for n = 1:numel(ms1)
    
    curve2_warped_n = h2(sbmx + Fx - ms1(n),sbmy + Fy - ms2(n));
    
    % in 2D
    hs = surf(axWarp2D,sbmx,sbmy,curve2_warped_n);
    set(hs,'edgecolor','none');
    
    % in 1D
    hs = surf(axWarp1D,[sbmx(1,:) ; sbmx(1,:)],[ones(size(sbmx(1,:))) ; 2*ones(size(sbmx(1,:)))],[curve2_warped_n(ceil(nbin/2),:) ; curve2_warped_n(ceil(nbin/2),:)]);
    set(hs,'edgecolor','none');
    
end

%% Save plots
if ~exist('panels/Summary/'); mkdir('panels/Summary/'); end

saveas(figProb,'panels/Summary/fig2a.eps','epsc')
saveas(figUni,'panels/Summary/fig2b.eps','epsc')
saveas(figDnsty,'panels/Summary/fig2c.eps','epsc')
saveas(figWarp2D,'panels/Summary/fig2d.eps','epsc')
saveas(figArrow,'panels/Summary/fig2e.eps','epsc')
saveas(figMani,'panels/Summary/fig2f.eps','epsc')
saveas(figWarp1D,'panels/Summary/fig2g.eps','epsc')
