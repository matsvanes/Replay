function GoodChannel=find_goodChannel(iSj, is)
% This function finds the good channels in all relevant MEG sessions in the
% Replay data.
fprintf('Find good channels over all sessions of interest \n')
GoodChannel=nan(is.nChan,numel(is.sessionsOfInterest));
iSjGood=is.goodsub(iSj);
cnt=1;
for isession=1:length(is.MEGruns{1,iSjGood})%is.sessionsOfInterest
  tempdir=[is.OPTPath strrep(is.fnDate{iSjGood},'-','_') '/' is.fnMEG{iSjGood} '_' num2str(isession,'%02d') '.ds/highpass_' num2str(is.highpass) '.opt'];
  load(fullfile(tempdir, 'opt.mat'));

  if isempty (opt.results.spm_files_epoched{1,1})
    SessionMEG=spm_eeg_load(fullfile(tempdir,[opt.results.spm_files_basenames{1,1}]));
  else
    SessionMEG=spm_eeg_load(fullfile(tempdir,[opt.results.spm_files_epoched_basenames{1,1}]));
  end

  chan_inds = indchantype(SessionMEG,'meeg','GOOD');
  Channindex  = indchantype(SessionMEG,'meeg','ALL');
  [~,ind]=setdiff(Channindex,chan_inds);

  index=ones(1,is.nChan);
  index(ind)=0;
  GoodChannel(:,cnt)=index;
  cnt=cnt+1;
  clear opt;
end