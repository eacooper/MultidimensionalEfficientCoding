function p = circledensity(d, bmx, bmy, opt)
% FUNCTION p = circledensity(d, bmx, bmy, opt)
%   Returns PDF for a pillbox function having
%   diameter d over domain bmx, bmy.
%
% PARAMETERS
%   d    :  Diameter of circle
%   bmx  :  x-coordinates of lattice
%   bmy  :  y-coordinates of lattice
%
% RETURNS
%     p  :  A PDF, size(p)==size(bmx)
%
  h = bmx(1,2)-bmx(1,1);
  p = fspecial('disk', d/h/2);
  g = size(bmx) - size(p);
  p = padarray(p, g/2);

  if exist('opt') && isequal(opt, 'strict')
    p = double(p==maxn(p));
  end
  p = p/sumn(p*h^2);
end