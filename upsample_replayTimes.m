function upsampled = upsample_replayTimes(replaytimes, phase, t_window)
% NOTE: not the prettiest solution..
% upsample for state timecourses at higher sampling rate than replay
% scores. 
[~,original] = find(replaytimes);
if size(phase,3)>length(replaytimes)
  Q = size(phase,3) / length(replaytimes);
  upsampled = round(original*Q);
else
  Q = 1;
  upsampled = original;
end
if any(upsampled>size(phase,3)-t_window)
  upsampled(upsampled>size(phase,3)-t_window) = [];
  original(original>size(phase,3)-t_window) = [];
end