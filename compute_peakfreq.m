function [pow, peakfreq] = compute_peakfreq(D, freq, taper, smoothing, doplot)

if ~exist('freq','var'), freq=10; end
if ~exist('taper','var'), taper='hanning'; end
if ~exist('smoothing', 'var'), smoothing=5; end

% get alpha power and phase time course
d = spm2fieldtrip(D);
d.hdr.label=d.label;

cfg         = [];
cfg.length  = 10; % 10 second windows (0.1 Hz resolution)
cfg.overlap = 0.5; % half overlap
d = ft_redefinetrial(cfg, d); 


cfg=[];
cfg.method = 'mtmfft';
cfg.taper = taper;
if strcmp(taper, 'dpss')
  cfg.tapsmofrq = smoothing;
end
cfg.foi = freq;
cfg.pad = ceil(d.time{1}(end));
cfg.output = 'pow';
pow = ft_freqanalysis(cfg, d);

sel = [1:4,11:16,21,22,25,26];
x=mean(pow.powspctrm(sel,:));
if doplot
  figure; plot(pow.freq, x)
  hold on
  [val, idx] = max(x(find(pow.freq==6):find(pow.freq==14)));
  peakfreq = pow.freq(find(pow.freq==6)-1+idx);
  vline(peakfreq)
end
end

