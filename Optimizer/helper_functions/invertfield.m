function [Fxi Fyi] = invertfield(Fx, Fy, xm, ym)
% FUNCTION [Fxi Fyi] = invertfield()
%   Computes the inverse of the field, (Fx, Fy)
%
% PARAMETERS
%   Fx  :  x-component of vector field, MxN
%   Fy  :  y-component of vector field, MxN
%   xm  :  Mesh of x-coordinates, MxN
%   ym  :  Mesh of y-coordinates, MxN
%
% RETURNS
%   Fxi :  x-component of inverse field at points (xm,ym)
%   Fyi :  y-component of inverse field at points (xm,ym)
%
  xFx = xm(:) + Fx(:);
  yFy = ym(:) + Fy(:);
  Qx  = scatteredInterpolant(xFx, yFy, -Fx(:));
  Qy  = scatteredInterpolant(xFx, yFy, -Fy(:));
  Fxi = Qx(xm,ym);
  Fyi = Qy(xm,ym);
end