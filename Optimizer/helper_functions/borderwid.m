function b = borderwid(ntaps)
% FUNCTION b = borderwid(ntaps)
%   Computes the width of the padding border that 
%   should be added to the potential function if 
%   ntaps are used to compute first derivatives.
%
% PARAMETERS
%   ntaps : Number of taps in first derivative filter
%  
% RETURNS
%       b : Number of border samples
%
  % Compute the border with. This may be confusing. We need to allow
  % the partial derivatives to be non-zero at the edge of F. The boundary
  % derivatives will be zero because we are adding that boundary to F.
  % That boundary is created by duplicating the edge values of F.
  % With that as context, for the partial at the edge of F to be non-zero,
  % the derviative filter at the edge must overlap some non-replicated 
  % values in F. If we were to replicate a full filter-width of samples
  % at the edge of F, then the derivative filter at the edge would have
  % to give zero. The correct solution is to replicate one less than the
  % filter width of samples. 
  %   That explanation is horrible, so you are going to have to think. Sorry.
  b = floor(ntaps/2);% - 1; 
end