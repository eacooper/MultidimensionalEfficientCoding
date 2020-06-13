function [p,psep] = weird_pdf1(x,y)
%
% generate a non-separable 2-D probability distribution

stdev_x = 0.15;
stdev_y = 0.15;
p       = exp(-(((x.^2)./(stdev_x.*abs(y*1.5)+0.005)) + ((y.^2)./stdev_y)));

% compute separable approximation
[u,s,v] = svd(p);
psep    = u(:,1)*v(:,1)';