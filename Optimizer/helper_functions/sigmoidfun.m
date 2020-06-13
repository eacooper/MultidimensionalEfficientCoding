function f = sigmoidfun(a,b,l,h)
% FUNCTION f = sigmoidfun(a, b, l, h)
%   Returns a 1D sigmoid function such that
%     f(a) == 0.01 
%     f(b) == 0.99
%   and
%     f([a b]) == [0.01 0.99]
% 
% PARAMETERS
%   a : The lower threshold of the sigmoid
%   b : The upper threshold of the sigmoid
%   l : (Optional) f(a) == l
%   h : (Optional) f(b) == h
%
% RETURNS
%   f : A function pointer to the sigmoid
%         f(a) = 0.01
%         f(b) = 0.09
%

  % Defaults
  if ~exist('l', 'var')
    l = 0.01;
  end
  if ~exist('h', 'var')
    h = 0.99;
  end

  % Solve Mp = B
  M = [a 1; b 1];
  B = [H(l); H(h)];
  p = M\B;

  % Return f
  f = @(x) reshape(1./(1+exp(-[x(:) ones(numel(x),1)]*p)), size(x));

end


function y = H(x)
  y = -log((1-x)/x);
end















