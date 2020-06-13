function [u v inds] = cornerGradFilt(ctype, dtype, h, ntaps, fsz)
% FUNCTION [u v nu nv inds] = cornerGradFilt(ctype, dtype, h, ntaps, fsz)
%   Computes the filter that should be applied to the corner of a function
%   to compute the partial w.r.t a derivative filter.
%
%  WARNING: STILL HARD CODED FOR ONLY 5 TAP FILTERS
% 
% PARAMETERS
%   ctype : Corner, 'ul' 'ur', 'lr, 'll'
%   dtype : Type of derivative, 'xx' 'yy', 'yx', 'x', 'y'
%   h     : Lattice spacing
%   ntaps : Number of filter taps
%   fsz   : Size of 2D function F that filter will be applied to
%
% RETURNS
%    u,  v : Filters for rows and cols respectively
%   nu, nv : Number of elements in filters (convienence)
%    inds  : The indices that filter should be applied over
%              inds = [jbeg jend ibeg iend]
%
  % Setup
  isz = fsz;        % Size of indicator array
  b   = borderwid(ntaps);   

  % Compute indices of 1's in indicator matrix
  switch ctype
   case 'ul'
    j = 1:b+1;   % b+1 is the unknown, which is repeated onto 
    i = 1:b+1;   %  ..the border elements 1:b
   case 'ur'
    j = 1:b+1;
    i = isz(2)-b:isz(2);
   case 'lr'
    j = isz(1)-b:isz(1);
    i = isz(2)-b:isz(2);
   case 'll'
    j = isz(1)-b:isz(1);
    i = 1:b+1;
  end

  % Compute the filter
  I         = zeros(isz);   % Indicator array marks where the variable appears
  I(j,i)    = 1;            % ...
  [u v j i] = derivGradFilt(I, dtype, h, ntaps, 'replicate');
  nu        = numel(u);  
  nv        = numel(v);  

  % Compute indices where result should be saved
  switch ctype
   case 'ul'
    jo = (b+1)*[1 1];
    io = (b+1)*[1 1];
   case 'ur'
    jo =      (b+1)*[1 1];
    io = (fsz(2)-b)*[1 1];
   case 'lr'
    jo = (fsz(1)-b)*[1 1];
    io = (fsz(2)-b)*[1 1];
   case 'll'
    jo = (fsz(1)-b)*[1 1];
    io =      (b+1)*[1 1];
  end

  % Returnd the indices to filter
  inds   = [j i jo io]';
end