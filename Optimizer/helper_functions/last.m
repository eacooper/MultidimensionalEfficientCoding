function y = last(x, n)
% FUNCTION y = last(x, n)
%   Finds the last n elements.
%   The code is:
%     y = find(x, n, 'last');
%   
%   n is optional!
%
% PARAMETERS
%   x  : An array
%   n  : (optional) The number of last items to return.
%                   If not specified, n = 1.
%
% RETURNS
%   y  : The index of the last n elements
%

  if ~exist('n','var')
    n = 1;
  end
  y = find(x, n, 'last');

end