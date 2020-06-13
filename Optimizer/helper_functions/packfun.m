function f = packfun(F, ntaps)
% FUNCTION f = packfun(F, ntaps)
%   Creates a vector v from F that is passed to the nonlinear 
%   minimizer to start an optimization.  User unpackfun(..)
%   to transform f back into F.
%  
% PARAMETERS
%   F    :  A potential function
%  ntaps :  Number of filter taps in the first derivative filter
%
% RETURNS
%   f  :  Vector of potential values
%  
  % Convert central region into a vector
  b = borderwid(ntaps);
  s = size(F,1);
  t = b+1:s-b;
  f = F(t,t);
  f = f(:);
end