function [F,Fx,Fy,Fxi,Fyi,pU,fhand] = compute_uniform_transform(nbin,pD,win_type)
% FUNCTION [F,Fx,Fy,Fxi,Fyi,pU,fhand] = compute_uniform_transform(nbin,pD,win_type)
%
% This function runs the numerical method for estimation the transformation
% from non-uniform to uniform probabilities for 2-D probability
% distributions describe in:
%
% Uniform Transformation of Non-Separable Probability Distributions  
% Eric Kee
% axXiv preprint598arXiv:160901982. 2016
%
% PARAMETERS
%       nbin:       number of bins used to sample each dimension of the probability distribution
%       pD:         non-uniform probability distribution
%       win_type:   shape of the window used to match the uniform and
%                   non-uniform distributions at the edges of the domain
%                   ('circle' or 'square')
%
% OUTPUT
%       F:          manifold potential
%       Fx:         gradient field in X
%       Fy:         gradient fiel in Y
%       Fxi:        inverse gradient field in X
%       Fyi:        inverse gradient field in Y
%       pU:         windowed uniform probability distribution
%       fhand:      figure handle for diagnostic figure showing result

%% Set up probability distributions

% Build grid for data distribution
beg = -1;
ned =  1;

ex  = linspace(beg, ned, nbin+1);    % bin edge x
ey  = linspace(ned, beg, nbin+1);    % bin edge y

bx  = ex(1:end-1) + diff(ex(1:2))/2;  % bin center x
by  = ey(1:end-1) + diff(ey(1:2))/2;  % bin center y

% make 2D mesh with center coordinates of each bin
[bmx, bmy] = meshgrid(bx, by);

% spacing between elements, used for normalization
h   = ex(2)-ex(1);

% Make uniform density of same resolution/support as pD and compute change in density
[pU,p] = compute_density_change(win_type,bmx,bmy,h,pD);

% just put these variables into a single structure, useful for debugging
meta = estruct(p, pU, pD, h, bmx, bmy);


%% Solve for manifold potential

% number of filter taps for first derivative (must be odd)
ntaps = 5;

% solve for the manifold potential
F    = potential(p, bx, by, meta, 1000, ntaps);

% get handle of final diagnostic figure
fhand = gcf;

% Compute gradient field & inverse
Fx 	= deriv(F, 'x', h, ntaps, 'replicate');
Fy	= deriv(F, 'y', h, ntaps, 'replicate');

Fx 	= double(Fx);
Fy 	= double(Fy);

[Fxi, Fyi] = invertfield(Fx, Fy, bmx, bmy);

% indices to use (invalid points are outside the window)
M = double(pU > 0);

% set invalid indices to zero
Fx = M.*Fx;
Fy = M.*Fy;

Fxi = M.*Fxi;
Fyi = M.*Fyi;

