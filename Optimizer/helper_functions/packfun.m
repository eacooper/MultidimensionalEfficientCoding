function f = packfun(F, ntaps)
% FUNCTION f = packfun(F, ntaps)
%   Creates a vector v from F that is passed to the nonlinear 
%   minimizer to start an optimization.  Use unpackfun(..)
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
  b = borderwid(ntaps); % compute width of padding border, which is the amount of the edges of F that have been replicated
  s = size(F,1);        % size of F along one dimension
  t = b+1:s-b;          % indices that exclude padded region
  f = F(t,t);           % values of F excluding padded region
  f = f(:);             % vectorize these values
end