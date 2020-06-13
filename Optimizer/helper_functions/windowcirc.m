function p = windowcirc(p, d, bmx, bmy)
% FUNCTION p = windowrho(p, d)
%   This tapers the values at the boundary of p
%   gradually to zero 
%
  % Compute radial distances
  midx = mean(bmx(1,:));
  midy = mean(bmy(:,1));
  r    = sqrt( [bmx-midx].^2 + [bmy-midy].^2);

  % Compute max distance
  mmax = maxn(abs(bmx-midx));

  % Make sigmoidal falloff
  fmax = 0.999999;
  fmin = 0.01;
  f  = sigmoidfun(d/2, mmax, fmax, fmin);
  m  = f(r);
  m(m>fmax) = 1;

  % Window the function
  p = p.*m;
end