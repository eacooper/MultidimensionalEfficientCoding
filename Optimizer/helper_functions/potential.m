function F = potential(p, bx, by, meta, iter, ntaps)
% FUNCTION F = potential(p, bx, by, iter, meta)
%   Computes a scalar potential function F that
%   embodies the gradient field that serves as a map between stimulus space and sensory space.
%   The potential function
%   is constrained by the normalized change in
%   density, from stimulus space to sensory space
%     p = (pD - pU) / pU
%
% PARAMETERS
%       p : A MxM matrix
%      bx : Bin x-positions
%      by : Bin y-positions
%    iter : number of iterations to run in optimization
%    meta : Metadata for error function
%

% Compute number of pyramid levels (course to fine scales for optimization)
nbin = size(p,1);
nl   = floor(log2(nbin/30));

% For each level
lvls = nl:-1:0;

for l = 1:numel(lvls)
    
    % Compute size of pyramid level
    nbinl = nbin/2^lvls(l); % 2^lvl(l) represents the degree of downsampling
    nbinl = ceil(nbinl);
    
    % Make low-resolution version of p
    bxl = linspace(bx(1), bx(end), nbinl);
    byl = linspace(by(1), by(end), nbinl);
    hl  = bxl(2) - bxl(1);
    plo = imresize(p, nbinl*[1 1], 'bicubic');
    
    % Make low-resolution version of pD
    [bmxl, bmyl] = meshgrid(bxl,byl);                        % mesh of low-resolution support
    pDlo        = imresize(meta.pD, nbinl*[1 1], 'bicubic'); % low-resultion probability distribution
    pDlo        = pDlo / sumn(pDlo*hl^2);                    % normalized
    
    % Make metadata for debug minimization
    metalo      = meta;
    metalo.p    = plo;
    metalo.pD   = pDlo;
    metalo.h    = hl;
    metalo.bmx  = bmxl;
    metalo.bmy  = bmyl;
    
    % add a scalar value for probability contained within each bin in a uniform probability distribution (no windowing)
    mE          = 1/numel(bmxl);
    metalo.pU   = mE/hl^2; 
    
    % Initialize potential function
    if l==1
        F = zeros(nbinl, nbinl);
    else
        F = imresize(F, nbinl*[1 1], 'bicubic');
    end
    
    % Solve current level
    f = packfun(F, ntaps); % vectorize the potential function
    f = minimize(f, @dEdF, iter, plo, hl, ntaps, metalo);
    F = unpackfun(f, size(plo), ntaps);
end
end