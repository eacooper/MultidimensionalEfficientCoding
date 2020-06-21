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
  b     = borderwid(ntaps);                 % compute width of padding border, which is the amount of the edges of F that have been replicated
  szlow = sz-2*b;                           % size of F along one dimension without padding
  F     = reshape(f, szlow);                % reshape f into F   
  F     = padarray(F, [b b], 'replicate');  % add the replicate padding in
end