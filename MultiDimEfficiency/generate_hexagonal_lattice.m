function [ms1,ms2] = generate_hexagonal_lattice(spacing,stim_range)
% generate hexagonal lattice for neurons

m1          = 0:spacing:stim_range(2);
m1          = [-m1(end:-1:2) m1];
m2          = 0:sqrt(0.75)*spacing:stim_range(2);
m2          = [-m2(end:-1:2) m2];
[ms1,ms2]   = meshgrid(m1,m2);

% shift spacing on odd rows
for r = 1:size(ms1,2)
    if mod(r,2)
        ms1(r,:) = ms1(r,:)+(spacing/2);
    end
end

% crop to circle
ms1     = ms1(:);
ms2     = ms2(:);
keep    = abs(ms1) <= 1.5 & abs(ms2) <= 1.5;
ms1     = ms1(keep);
ms2     = ms2(keep);
