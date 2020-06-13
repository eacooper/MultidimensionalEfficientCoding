function [p,psep] = weird_pdf2(x,y)
%
% generate a non-separable 2-D probability distribution

stdev_x = 0.15;
stdev_y = 0.6;
mu_y    = -0.65;

p = abs((y-mu_y)*1.5).*exp(-(((x.^2)./(stdev_x.*abs((y-mu_y).^1.5)+0.0001)) + (((y-mu_y).^2)./stdev_y)));

p(y < mu_y) = 0;

% compute separable approximation
[u,s,v] = svd(p);
psep = u(:,1)*v(:,1)';