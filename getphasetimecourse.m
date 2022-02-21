function [phase, pow, maxidx] = getphasetimecourse(D, freq, taper, smoothing)

if ~exist('freq','var'), freq=10; end
if ~exist('taper','var'), taper='hanning'; end
if ~exist('smoothing', 'var'), smoothing=5; end

% get alpha power and phase time course
d = spm2fieldtrip(D);
% mtmconvol
numcycles=2;
cfg=[];
cfg.method = 'mtmconvol';
cfg.taper = taper;
cfg.foi = freq;
if strcmp(taper, 'dpss')
  cfg.tapsmofrq = smoothing;
  cfg.t_ftimwin = repmat(1/smoothing, [1, length(cfg.foi)]);
else
  cfg.t_ftimwin = numcycles*1./cfg.foi;
end
cfg.pad = ceil(d.time{1}(end));
cfg.output = 'fourier';
cfg.toi = 0:1./D.fsample:d.time{1}(end);
freq = ft_freqanalysis(cfg, d);

pow = nanmean(abs(freq.fourierspctrm).^2,4);
for ifoi=1:length(freq.freq)
  [~, maxidx(1,ifoi)] = max(squeeze(pow(:,:,ifoi,:)));
end
pow=squeeze(pow);

phase = angle(squeeze(freq.fourierspctrm))+pi;
end

