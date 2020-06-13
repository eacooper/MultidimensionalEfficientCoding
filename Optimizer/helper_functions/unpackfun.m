function F = unpackfun(f, sz, ntaps)
% FUNCTION F = unpackfun(f, sz, ntaps)
%   Creates a potential function F from a vector that was returned by the 
%   nonlinear minimizer.
%  
% PARAMETERS
%   f  :  A vector, Nx1
%   sz :  The size of the potential function
%
% RETURNS
%   F  :  A potential function
%
  % Unpack f 
  b     = borderwid(ntaps); 
  szlow = sz-2*b;
  F     = reshape(f, szlow);
  F     = padarray(F, [b b], 'replicate');
end