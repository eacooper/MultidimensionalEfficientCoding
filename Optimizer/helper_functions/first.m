function y = first(x, n)
% FUNCTION y = first(x, n)
%   Finds the first n elements.
%   The code is:
%     y = find(x, n, 'first');
%   
%   n is optional!
%
% PARAMETERS
%   x  : An array
%   n  : (optional) The number of first items to return.
%                   If not specified, n = 1.
%
% RETURNS
%   y  : The index of the first n elements
%

  if ~exist('n','var')
    n = 1;
  end
  y = find(x, n, 'first');

end