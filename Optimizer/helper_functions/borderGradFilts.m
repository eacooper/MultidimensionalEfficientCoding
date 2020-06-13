function [U V D n] = borderGradFilts(dtype, h, ntaps, fsz)
% FUNCTION [U V D n] = borderGradFilts(dtype, h, ntaps, fsz)
%   Computes the filters that are needed to compute the 
%   error gradient along the boundary of the domain.
%   
%   NOTE:  Tested for 5 and 9 taps; 3 taps is not defined
%          for the proposed method.
%
% PARAMETERS
%   dtype  :  What kind of convolution are you computing 
%             gradient of? 'xx'  'yy' or  'xy'?
%
%     fsz  :  The size of the matrix to which these 
%             filters are going to be applied. This way
%             we can pre-compute the indices that are filtered over
%
% RETURNS
%   U,V    :  Matrices specifying separable filters
%               U(:,k) and V(k,:) specify a separable filter
%     D    :  A matrix, each column specifies
%               nu     :  Number of valid elements in a column of U
%                          e.g. u1 = U(1:nu(1),1);
%               nv     :  Valid elements in a row of V
%               bj, ej :  Begin and end row indices to pass filter over
%               bi, ei :  Begin and end col indices to pass filter over
%               bj, ej :  Begin and end row indices to save result to
%               bi, ei :  Begin and end col indices to save result to 
%     n    :  The total number of edge samples that are computed
%             by the filter bank
%
  % Intialize lists of filter data
  k = 0;  U = {};  V = {};  D = {};

  % Corners
  types = {'ul' 'ur' 'lr' 'll'};
  for k = 1:numel(types)
    type       = types{k};
    [u v inds] = cornerGradFilt(type, dtype, h, ntaps, fsz);
    U{end+1}   = u;
    V{end+1}   = v;
    D{end+1}   = inds; 
  end

  % Edges
  types = {'t' 'r' 'b' 'l'};
  for k = 1:numel(types)
    type       = types{k};
    [u v inds] = edgeGradFilts(type, dtype, h, ntaps, fsz); 
    U = [U u];
    V = [V v];
    D = [D inds];
  end

  % Pack U and V into matrices for processing speed
  [U nu] = packit(U);
  [V nv] = packit(V);
  D      = [nu; nv; packit(D)];

  % Compute the length of the border for reference
  n = 0;
  for k = 1:size(D,2)
    aj = D( 7, k);
    bj = D( 8, k);
    ai = D( 9, k);
    bi = D(10, k);
    nj = numel(aj:bj);
    ni = numel(ai:bi);
    n  = n + nj*ni;
  end
end
