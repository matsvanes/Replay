function [pval z] = circ_rtest_multivar(alpha, w, d, dim)
% Adaptation of circ_rtest from Matlab's circular statistics toolbox (see
% original documentation below). Adapted to work for multivariate data.
%
% Mats van Es, 2022
%
% [pval, z] = circ_rtest(alpha,w)
%   Computes Rayleigh test for non-uniformity of circular data.
%   H0: the population is uniformly distributed around the circle
%   HA: the populatoin is not distributed uniformly around the circle
%   Assumption: the distribution has maximally one mode and the data is
%   sampled from a von Mises distribution!
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%
%   Output:
%     pval  p-value of Rayleigh's test
%     z     value of the z-statistic
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 3 || isempty(d)
  % per default do not apply correct for binned data
  d = 0;
end

if nargin < 4 || isempty(dim)
  dim = 1;
end

reshapeflag=0;
if length(size(alpha))>2 || dim>2 % reshape alpha to be 2-dimensional, with the first dimension to compute the statistic over
  origdim=dim;
  s=size(alpha);
  alpha = reshape(permute(alpha, [dim, setdiff(1:length(s), dim)]), s(dim), []);
  dim=1;
  reshapeflag=1;
end

if nargin < 2 || isempty(w)
  w=[];
  n=size(alpha,dim);
else
  n=sum(w);
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
    error('Input dimensions do not match');
  end
end

r =  circ_r(alpha,w(:),d, dim);

% compute Rayleigh's R (equ. 27.1)
R = n*r;

% compute Rayleigh's z (equ. 27.2)
z = R.^2 / n;

% compute p value using approxation in Zar, p. 617
pval = exp(sqrt(1+4*n+4*(n^2-R.^2))-(1+2*n));

if reshapeflag
  pval = reshape(pval, s(setdiff(1:length(s), origdim)));
  z = reshape(z, s(setdiff(1:length(s), origdim)));
end



% outdated version:
% compute the p value using an approximation from Fisher, p. 70
% pval = exp(-z);
% if n < 50
%   pval = pval * (1 + (2*z - z^2) / (4*n) - ...
%    (24*z - 132*z^2 + 76*z^3 - 9*z^4) / (288*n^2));
% end

