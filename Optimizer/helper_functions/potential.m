function F = potential(p, bx, by, meta, iter, ntaps)
% FUNCTION F = potential(p, bx, by, meta)
%   Computes a scalar potential function F that
%   embodies the gradient field that maps Euclidean
%   space into data space. The potential function
%   is constrained by the normalized change in 
%   density, from Euclidean to data space
%     p = (pD - pU) / pU
%
% PARAMETERS
%       p : A MxM matrix, p = (pD-pU)/pU
%      bx : Bin x-positions
%      by : Bin y-positions
%    meta : Metadata for error function
%
  % Compute number of pyramid levels
  nbin = size(p,1);
  nl   = floor(log2(nbin/30));

  % For each level
  lvls = nl:-1:0;
  for l = 1:numel(lvls)
    % Compute size of pyramid level
    nbinl = nbin/2^lvls(l);
    nbinl = ceil(nbinl);
    
    % Make low-resolution version of p
    bxl = linspace(bx(1), bx(end), nbinl);
    byl = linspace(by(1), by(end), nbinl);
    hl  = bxl(2) - bxl(1);
    plo = imresize(p, nbinl*[1 1], 'bicubic');

    % Make metadata for debug minimization
    [bmxl bmyl] = meshgrid(bxl,byl);
    mE          = 1/numel(bmxl);
    pDlo        = imresize(meta.pD, nbinl*[1 1], 'bicubic');
    pDlo        = pDlo / sumn(pDlo*hl^2);
    metalo      = meta;
    metalo.p    = plo;
    metalo.pU   = mE/hl^2;
    metalo.pD   = pDlo;
    metalo.h    = hl;
    metalo.bmx  = bmxl;
    metalo.bmy  = bmyl;

    % Initialize potential function
    if l==1
      F = zeros(nbinl, nbinl);
    else
      F = imresize(F, nbinl*[1 1], 'bicubic');
    end

    % Solve current level  
    f = packfun(F, ntaps);
    f = minimize(f, @dEdF, iter, plo, hl, ntaps, metalo);
    F = unpackfun(f, size(plo), ntaps);
  end
end