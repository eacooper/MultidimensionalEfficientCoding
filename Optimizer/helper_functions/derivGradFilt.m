function [u v j i] = derivGradFilt(I, dtype, h, ntaps, boundopt)
% FUNCTION [u v j i] = derivGradFilt(I, dtype, h, ntaps, boundopt)
%   Computes a separable filter that computes
%   a special gradient at the edge of the domain.
%
% PARAMETERS
%      I    :  Array containing 1s and 0s. I should have a value 
%              of 1 at all positions where some function F is
%              set by the same variable e.g. F(1,1) = x_1 = F(1,2)
%              means I(1,1:2) = 1.
%   dtype   :  The kind of derivative to compute. This is 
%                'x' 'y' 'xx' 'yy' 'xy' 'yx'
%            
%      h    :  The lattice spacing that your derivatives
%              would be computed with
% 
%   ntaps   :  Number of taps to use in the derivative filters
%
%  boundopt :  What to do at boundary,
%                          0 : pad with zero
%                'replicate' : replicate
%
% RETURNS
%   u   :  Second 1D filter to apply via convolution, Nx1
%   v   :  First 1D filter to apply via convolution,  1xM
%   j,i :  Apply the filter to I(j,i); all other positions 
%
  % Compute SVD to make separable version
  dd = deriv_new(I, dtype, h, ntaps, boundopt);  
  [u s v] = svd(dd,0); 
  u = sqrt(s(1))*u(:,1); 
  v = sqrt(s(1))*v(:,1)';  

  % Compute where the filter is non-zero
  uf = any(dd ~= 0, 2);  % Intentionally use hard ~= 0.  That can only be true for positions
  vf = any(dd ~= 0, 1);  % of dd where the convolution kernal overlapped all exact-zeros in I

  % Save range of indices at which to apply the filter
  j = [first(uf) last(uf)];
  i = [first(vf) last(vf)];

  % Keep only non-zero portion
  u  = u(j(1):j(2));
  v  = v(i(1):i(2));

  % Pre-flip filter for convolution
  u = flipud(u);
  v = fliplr(v);
end