function n = borderlen(fsz, ntaps)
% FUNCTION n = borderlen(fsz, ntaps)
%   Computes the number of samples along the border
%   of a potential function fsz that is being processed
%   with a first derivative filter having number 
%   of taps, 'ntaps'.
% 
% PARAMETERS
%   fsz   : Size of potential function 
%   ntaps : Number of taps in first derivative function
%
% RETURNS
%   n  :  Number of samples along the border    
%
%
  b   = borderwid(ntaps);
  sz1 = fsz - 2*b;
  sz0 = fsz - 2*(b+1);
  n   = prod(sz1) - prod(sz0);
end