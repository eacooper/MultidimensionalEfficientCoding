function S = estruct(varargin)
% FUNCTION S = estruct(varargin)
%   Returns a structure S containing each of the varaibles passed as parameters.
%   This simplifies matlab's typical struct constructor, which requires variable
%   name/value pairs as arguments.  estruct uses the same varaible names as passed
%   to it.
% 
% PARAMETERS
%   Any set of variables.  Note: you cannot call estruct with literals as parameters.
%     e.g.  estruct(t)   is allowed, but not   estruct('t')
%
% RETURNS
%   S  :  A structure containing the variables in varargin     
%
% REVISION
%   2010.04.25 : Eric Kee
%
  for k = 1:numel(varargin)
    S.(inputname(k)) = varargin{k};
  end

end