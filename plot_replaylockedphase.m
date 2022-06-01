function [N, EDGES] = plot_replaylockedphase(topphase, nbins, doplot)

if nargin<2 || isempty(nbins), nbins=12; end
if nargin<3 || isempty(doplot), doplot=false; end


[N, EDGES] = histcounts(topphase, nbins);

if doplot
  plot(2*pi/nbins:2*pi/nbins:2*pi, N)
end
