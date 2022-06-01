function [itpc, pow] = compute_itpc(spctrm, topidx)

[s1,s2,~] = size(spctrm);
s3=251;
spctrm_tl = zeros(length(topidx), s1, s2, s3);
s0=numel(topidx);
for k=1:s0
  spctrm_tl(k,:,:,:) = spctrm(:,:, topidx(k)-125:topidx(k)+125);
end
amp       = abs(spctrm_tl);
itpc      = spctrm_tl./amp;         % divide by amplitude
itpc      = sum(itpc,1);   % sum angles
itpc      = abs(itpc)/size(spctrm_tl,1);   % take the absolute value and normalize
itpc      = squeeze(itpc);
pow       = squeeze(mean(amp.^2,1));

