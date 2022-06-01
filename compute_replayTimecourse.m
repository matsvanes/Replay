function ToRall = compute_replayTimecourse(is,iLCV,Mreplay,iRun, iPerm)
%% Save an inferred Replay Timecourse:

%define peak sequenceness time:
[~,lag] = max(abs(Mreplay));
if iRun==2 && ~is.doPermutation
  usepath = [is.AnalysisPath, 'classifiers/TestResting4Cell/']; % original results from Yunzhe
else
  usepath = [is.matsdir, 'classifiers/TestResting/'];
end

for iSj=1:length(is.goodsub)
  % take that subject's CV parameter selection:
  iL = 5;% hard coded because this is used in the past. No crossvalidation!! %iLCV(iSj);

  % load reactivation timecourse:
  if is.doPermutation==true
    S.Rreds = apply_toRest(iSj, is, [], [], iRun, iPerm);
  else
    S = load([usepath,'Rreds' num2str(iSj) '_' num2str(iRun)]); Rreds = S.Rreds;  % load this subject's preds
  end

  reactivationtimecourse = S.Rreds{1,37,1,iL}; % 37 is 200 ms
  % convert to probability:
  reactivationtimecourse = 1./(1+exp(reactivationtimecourse));

  % and logic operation for replay timecourse:
  T = [0,1,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0; ...
    0,0,0,0,0,1,0,0;0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0];

  ToR = zeros(size(reactivationtimecourse,1),1);
  orig = reactivationtimecourse(1:(end-2*lag),:)*T;
  proj = reactivationtimecourse((1+lag):(end-lag),:);
  ToR(1:end-2*lag,1) = nansum(orig.*proj,2);

  ToRall(iSj,:) = ToR;
end
if ~is.doPermutation
  replaydir = [usepath,'classifiers/Sequence_by_Training_4Cell/'];
  mkdir(replaydir)
  if iRun==1
    fname = 'STUDYII_PreplayOnset.mat';
  else
    fname = 'STUDYII_ReplayOnset.mat';
  end
  save([replaydir,fname],'ToRall');
end