function [err dEdF pest] = dEdF(f, p, h, ntaps, meta)
% FUNCTION [err dEdF pest] = dEdF(f, p, h, ntaps, meta)
%   Computes an error, and gradient of that error.
%   Variable f is a potential function that embodies a gradient field
%   which warps a uniform probability density into a desired
%   non-uniform density. The transform to that non-uniform density 
%   embodied by p, the normalized change in density at each point.
%   p can be computed from f. This function does that, 
%   and then computes an error between p and the estimate from f.
%   
% PARAMATERS
%   f    :  A potential function, MxN vectorized
%   p    :  Normalized change in density, MxN
%   h    :  Lattice spacing
%  ntaps :  Number of taps to specify to deriv() function
%  meta  :  meta structure used for visualization and debugging
%
% RETURNS:
%   err  :  the error between p and the estimate of p derived from f

  % Get potential function from solution vector
  F = unpackfun(f, size(p), ntaps);

  % Derivatives
  Fx   = deriv_new(F,  'x', h, ntaps, 'replicate');
  Fy   = deriv_new(F,  'y', h, ntaps, 'replicate');
  Fxx  = deriv_new(Fx, 'x', h, ntaps, 'replicate');
  Fyy  = deriv_new(Fy, 'y', h, ntaps, 'replicate');
  Fxy  = deriv_new(Fx, 'y', h, ntaps, 'replicate');

  % Residual
  pest = Fxx.*Fyy - Fxy.^2 + Fxx + Fyy; % Eq 17 from Kee paper
  R    = pest - p;

  % Error
  err = 0.5*sum(R(:).^2);

  % Terms for computing partials in central area (Eqs 18-23 in Kee paper)
  FxxR = Fxx.*R;
  FyyR = Fyy.*R;
  FxyR = Fxy.*R;

  T1 = deriv_new(FxxR, 'yy', h, ntaps, 0);
  T2 = deriv_new(FyyR, 'xx', h, ntaps, 0);
  T3 = deriv_new(FxyR, 'xy', h, ntaps, 0);
  T4 = deriv_new(   R, 'xx', h, ntaps, 0);
  T5 = deriv_new(   R, 'yy', h, ntaps, 0); 

  % Partials in central area 
  dEdF = T1 + T2 - 2*T3 + T4 + T5;

  % Partials along edges 
  [E1 ids] = edgeGrads(FxxR, 'yy', h, ntaps);
  [E2 ids] = edgeGrads(FyyR, 'xx', h, ntaps);
  [E3 ids] = edgeGrads(FxyR, 'xy', h, ntaps);  
  [E4 ids] = edgeGrads(   R, 'xx', h, ntaps);  
  [E5 ids] = edgeGrads(   R, 'yy', h, ntaps);  
  
  % partials along edges
  dEdF(ids) = E1 + E2 - 2*E3 + E4 + E5; 

  % Plot to show how it's going
  if (2*rand-0.005 < 0)
      
    pD      = meta.pU*p    + meta.pU;
    pDest   = meta.pU*pest + meta.pU;
    
    vdir = [1 1 1];

    figure(21); hold on; set(gcf, 'pos', [881 5 560 800]);
    subplot(321);
    surf(meta.bmx, meta.bmy, meta.p, 'edgealpha', 0.1); grid on; axis xy tight square; box on;
    title('Actual Rho..');
    view(vdir);

    subplot(322); 
    surf(meta.bmx, meta.bmy, pest, 'edgealpha', 0.1); grid on; axis xy tight square; box on;
    title('Estimated Rho..');
    view(vdir);

    subplot(323); 
    imagesc(meta.bmx(1,:), meta.bmy(:,1)', R); grid on; axis xy tight square; box on;
    title('Residual..'); colorbar;
    
    b  = borderwid(ntaps);
    vj = b+1:size(dEdF,1)-b;
    vi = b+1:size(dEdF,2)-b;

    subplot(324); 
    imagesc(meta.bmx(1,:), meta.bmy(:,1), dEdF); grid on; axis xy tight equal; box on; grid on;
    title('Gradient..'); colorbar;

    subplot(325); hold off;
    mid = round(size(meta.bmx,1)/2);
    plot(meta.bmx(mid,:), p(mid,:)); grid on; axis xy tight; hold on;
    plot(meta.bmx(mid,:), pest(mid,:));
    title('Change in density..');

    subplot(326); hold off;
    surf(meta.bmx, meta.bmy, F, 'edgealpha', 0.1); grid on; axis xy tight square; box on;
    title('Estimated Potential..');
    view(vdir);

    drawnow;
  end

  % Return only the unknowns that we are solving for 
  dEdF = packfun(dEdF, ntaps);

end