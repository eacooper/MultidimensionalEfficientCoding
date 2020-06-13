function [pU,p] = compute_density_change(win_type,bmx,bmy,h,pD)
% FUNCTION [pU,p] = compute_density_change(win_type,bmx,bmy,h,pD)
%
% This function makes a uniform density distribution of same
% resolution/support as pD and computes change in density between the two
%
%
% PARAMETERS
%       nbin:       number of bins used to sample probability distribution
%       pD:         non-uniform probability distribution
%       win_type:   shape of the window used to match the uniform and
%                   non-uniform distributions at the edges of the domain
%                   ('circle' or 'square')
%       bmy,bmx     2D mesh with center coordinates of each bin
%       h           spacing between elements, used for normalization
%
% OUTPUT
%       pU:         windowed uniform probability distribution
%       p:          normalized change in density

switch win_type
    
    case 'circle'
        
        % Make uniform density of same resolution/support as pD, but with a circular window after which pU = 0
        npad   = 1;
        d      = bmx(npad+1:end-npad+1, npad+1:end-npad+1);
        d      = range(d(:));
        pU     = circledensity(d, bmx, bmy, 'strict');
        
        % make sure pD == zero wherever pU == 0
        border     = pU==0;
        pD(border) = 0;
        
        % renormalize so that windowed pD is probability density
        pD = pD/sumn(pD*h^2);
        
        % Compute normalized change in density
        p       = (pD - pU) ./ pU;
        
        % replace nans from 0/0 with 0
        p(isnan(p)) = 0;
        
        
    case 'square'
        
        % Make uniform density of same resolution/support as pD, but with a square window after which pU = 0
        npad   = 4;
        pU                      = ones(size(pD));
        pU(1:npad,:)            = 0;
        pU(end-npad+1:end,:)    = 0;
        pU(:,1:npad)            = 0;
        pU(:,end-npad+1:end)    = 0;
        pU                      = pU/sumn(pU*h^2);
        
        % make sure pD == zero wherever pU == 0
        border     = pU==0;
        pD(border) = 0;
        
        % renormalize so that windowed pD is probability density
        pD = pD/sumn(pD*h^2);
        
        % Compute normalized change in density
        p  = (pD - pU) ./ pU;
        
        % set all nans to zeros
        p(isnan(p)) = 0;
        
        
    otherwise
        error('invalid window type');
        
end