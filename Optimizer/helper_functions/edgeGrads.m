function [dFs ids] = edgeGrads(F, dtype, h, ntaps)
% FUNCTION [dFs ids] = edgeGrads(F, dtype, h, ntap)
%   Computes the gradient of function f around its edge, assuming that its
%   edge value is repeated, and its gradient on the interior is given
%   by a derivative filter of kind, 'kind'
%
% PARAMETERS
%   F     :  Function to compute gradient along edge
%   dtype :  Derivative filter type, 'xx' 'yy' or 'xy'
%         
% RETURNS
%   dFs   :  A vector of the edge samples in F
%   ids   :  A vector of linear indices
%              e.g. F(ids) corresponds to dFs
%
  persistent Uxx Vxx Dxx Uyy Vyy Dyy Uxy Vxy Dxy n;
  if isempty(n) || (n ~= borderlen(size(F),ntaps))
    [Uxy Vxy Dxy n] = borderGradFilts('xy', h, ntaps, size(F));
    [Uxx Vxx Dxx n] = borderGradFilts('xx', h, ntaps, size(F));
    [Uyy Vyy Dyy n] = borderGradFilts('yy', h, ntaps, size(F));
  end

  % Make space for each border element
  n   = borderlen(size(F), ntaps);
  dFs = zeros(1,n);    
  ids = zeros(1,n);

  % Pick filter bank
  switch dtype 
   case 'xx'
    U = Uxx;
    V = Vxx;
    D = Dxx;
   case 'yy'
    U = Uyy;
    V = Vyy;
    D = Dyy;
   case 'xy'
    U = Uxy;
    V = Vxy;
    D = Dxy;
  end

  % Compute partial derivatives
  c = 0;
  for k = 1:size(D,2)
    % Get filter
    nu = D(1,k);
    nv = D(2,k);
    u  = U(1:nu,k);
    v  = V(k,1:nv);

    % Compute the partials
    j    = D(3,k):D(4,k);
    i    = D(5,k):D(6,k);
    vals = conv2( conv2(F(j,i), v, 'valid'), u, 'valid');

    % Construct linear indices where gradients should be stored
    j    = D(7,k):D(8,k);
    i    = D(9,k):D(10,k);
    inds = j + (i-1)*size(F,1);   % This is manual sub2ind to avoid error when size(j)~=size(i)
                 
    % Store gradient results
    m      = numel(vals);
    j      = c+1 : c+m;
    dFs(j) = vals;
    ids(j) = inds;
    c      = c + m;
  end
end