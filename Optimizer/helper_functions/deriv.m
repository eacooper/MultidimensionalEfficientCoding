function y = deriv(x, dtype, h, ntaps, boundopt)
% Function y = deriv(x, dtype, h, ntaps, boundopt)
%   Computes the derivative of a 2D function, x.
%   
% PARAMETERS
%          x : A 2D function, MxN
%      dtype : Kind of derivative, may be:
%                'x', 'y', 'xx', 'yy', or 'xy'
%          h : Lattice spacing, e.g. 1
%      ntaps : Number of filter taps, 3, 5, or 9
%   boundopt : Boundary options for imfilter
%                          0 : pad with zero
%                'replicate' : replicate
%
% RETURNS
%      y :  The derivative of x
%
  % Pick filter
  switch ntaps
   case 3
    p = [0.229879  0.540242  0.229879];
    d = [0.425287  0.000000 -0.425287];
   case 5
    p = [0.030320  0.249724  0.439911  0.249724  0.030320];
    d = [0.104550  0.292315  0.000000 -0.292315 -0.104550];
   case 9
    p = [0.000721  0.015486  0.090341  0.234494  0.317916   0.234494   0.090341   0.015486   0.000721];
    d = [0.003059  0.035187  0.118739  0.143928  0.000000  -0.143928  -0.118739  -0.035187  -0.003059];
  end

  % First derivative
  if dtype(1) == 'x'
    y = imfilter( imfilter(x, d, 'conv', boundopt), p',         'conv', boundopt) / h;
  else
    y = imfilter( imfilter(x, p, 'conv', boundopt), flipud(d'), 'conv', boundopt) / h;
  end

  % Second derivative
  n = numel(dtype);
  if n==2 && dtype(2) == 'x'
    y = imfilter( imfilter(y, d, 'conv', boundopt), p',         'conv', boundopt) / h;
  elseif n==2
    y = imfilter( imfilter(y, p, 'conv', boundopt), flipud(d'), 'conv', boundopt) / h;
  end
end