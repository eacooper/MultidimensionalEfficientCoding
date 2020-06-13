function [U nu] = packit(U)
% HELPER FUNCTION 
%   Packs a cell array of vectors U into a matrix 
%   when those vectors may have varying length
%   (pads with zeros). Assumes vectors are homogeneous
%   in their orientations as rows or columns, and 
%   concatenates them appropriately
%
  if isempty(U) 
    return;
  end
  nu   = cellfun(@(x)numel(x), U);
  maxn = max(nu);
  for k = 1:numel(nu)
    tt         = U{k};
    tt(maxn+1) = 0;
    U{k}       = tt(1:maxn);
  end
  sz = size(U{1});
  if sz(1) > sz(2)
    U = cell2mat(U(:)');
  else
    U = cell2mat(U(:));
  end
end