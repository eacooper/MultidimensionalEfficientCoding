function n = num2ndDerivTaps(ntaps)
% FUNCTION n = num2ndDerivTaps(ntaps)
%   Returns number of taps in a 2nd derivative
%   filter when using computing that second deriatives
%   with 2 first derivative filters
%
% PARAMETERS
%   ntaps : Number of taps in first deriative filter
% 
% RETURNS
%       n : Number of taps in second derivative filter
%
  % Compute width of filter
  filt = deriv(1, 'xx', 1, ntaps, 'full');  
  n    = size(filt,1);
end