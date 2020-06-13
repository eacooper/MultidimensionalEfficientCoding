function [U V inds] = edgeGradFilts(etype, dtype, h, ntaps, fsz)
% FUNCTION [U V inds] = edgeGradFilts(etype, dtype, h, ntaps, fsz)
%   Computes the distinct filters along a particular edge.
%
%   WARNING: Assumes 2nd derivatives are being computed by
%            applying two first deriative operations, which
%            makes the 2nd deriative filter larger than ntaps.
%
%   NOTE   : Has been manually tested for 5 and 9 taps, and
%            3 taps are degenerate given how we impose the
%            boundary constraints.
%
% PARAMETERS
%   etype  :  Edge type, a character: l t r b
%   dtype  :  Kind of derivative to apply, e.g. 'xx'
%       h  :  Lattice spacing
%   ntaps  :  Number taps in SECOND derivative filter
%     fsz  :  Size of the matrix that filters will be applied to 
%
% RETURNS
%     U,V  :  Matrix of u and v filters
%    inds  :  Indicies that should be convolved over, and indices
%             where the results should be stored.
%
  % Compute width of filter
  fwid2 = floor( num2ndDerivTaps(ntaps)/2 );

  % Length of edge
  if (etype=='l') || (etype=='r')
    edgelen = fsz(1);
  else
    edgelen = fsz(2);
  end
  
  % Build indices along edge that we need to compute filters for
  minsz     = 2*fwid2 + 1;   
  b         = borderwid(ntaps);   % ntaps, not ntaps2
  c         = b+2;                % b+1 is the corner sample on edge; b+2 is the first non-corner sample
  midpoint  = fwid2+1;            % Midpointx along edge is where filter is not cropped by a corner    
  leftInds  = c:midpoint-1;
  rightInds = edgelen - fliplr(leftInds) + 1;
  edgeInds  = [leftInds midpoint rightInds];

  % For each edge index
  for j = 1:numel(edgeInds)
    i  = edgeInds(j);
    I  = zeros(fsz);

    % Set indicator array according to which edge
    switch etype
     case 't'
      jj      = 1:b+1;
      I(jj,i) = 1;
     case 'r'
      jj      = fsz-b:fsz;
      I(i,jj) = 1;
     case 'b'
      jj      = fsz-b:fsz;
      I(jj,i) = 1;
     case 'l'
      jj      = 1:b+1;
      I(i,jj) = 1;
    end

    % Compute the filter
    [u v uj vi] = derivGradFilt(I, dtype, h, ntaps, 'replicate');  
    U{j} = u;
    V{j} = v;

    % Build input and output indices
    %   Input indices:  indices into the function that will be convolved
    %     in_inds = [jbegin jend ibegin iend]
    %   Output indices: indices where the convolution result should be stored
    %     out_inds = [jbegin jend ibegin iend]
    %
    %   inds = [in_inds out_inds]
    % 
    % Shorthand variables

    hei = fsz(1);           % Shorthand for height of matrix that filters will be applied to
    wid = fsz(2);           % Shorthand for width ...
    switch etype
     % Top edge
     case 't'
      if (i == midpoint)  
        cent    = [b+1 b+1 midpoint wid-midpoint+1];  
        inds{j} = [uj  1 wid  cent]';
      else
        inds{j} = [uj vi  b+1 b+1  i i ]';
      end

     % Right edge
     case 'r'
      if (i==midpoint)
        cent = [midpoint hei-midpoint+1 wid-b wid-b];  % Indices for central part of edge
        inds{j} = [1 hei vi cent]';
      else
        inds{j} = [uj vi  i i wid-b wid-b ]';
      end

     % Bottom edge
     case 'b'
      if (i==midpoint)
        cent = [hei-b hei-b midpoint wid-midpoint+1];  % Indices for central part of edge
        inds{j} = [uj  1 wid  cent]';
      else
        inds{j} = [uj  vi  hei-b  hei-b  i  i ]';
      end

     % Left edge
     case 'l'
      if (i==midpoint)
        cent = [midpoint hei-midpoint+1 b+1 b+1];  % Indices for central part of edge
        inds{j} = [1  hei  vi  cent]';
      else
        inds{j} = [uj  vi  i  i  b+1  b+1 ]';
      end
    end
  end
end