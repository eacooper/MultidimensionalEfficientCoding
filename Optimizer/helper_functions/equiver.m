function q = equiver(p, n, varargin)
% FUNCTION q = equiver(p, n, varargin)
%   Allows you to call quiver with a matrix of points and vectors
%   rather than having to specify each component individually.
%
% PARAMETERS
%   p  : Points  2xN or 3xN
%   n  : Vectors 2xN or 3xN
%
%   varargin : Other parameters that you would pass to quiver()
%
% RETURNS
%   q  : Handle of quivergroup object
% 
%

  d = size(p,1);
  if d == 2
    q = quiver(p(1,:), p(2,:), n(1,:), n(2,:), varargin{:});
  elseif d==3
    q = quiver3(p(1,:), p(2,:), p(3,:), n(1,:), n(2,:), n(3,:), varargin{:});
  else
    error('Quiver points or normals have wrong dimensions.');
  end

end