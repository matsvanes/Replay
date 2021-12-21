%% get phase information
% adapted from
% https://github.com/OHBA-analysis/Higgins2020_Neuron/blob/master/utils/NeuronFig2Analyses.m,
% from line 107
% and https://github.com/OHBA-analysis/Higgins2020_Neuron/blob/master/utils/MiscFunctions/plotReplayFig1.m
% from line 40
datadir = [is.studydir,'Neuron2020Analysis/Study2/bfnew_1to45hz/'];
[~,~,triggerpoints,goodsamples]=getSubjectMasks(datadir);

triggerpoints = triggerpoints(2:2:end);
goodsamples = goodsamples(2:2:end);
whichstudy=2;

% load replay times:
if whichstudy==1
  load([wd,'GenericReplayData/STUDYI_ReplayOnset/RwdReplay2ndOnset'],'ToRall');
  replayScores = ToRall;
else
  load([is.AnalysisPath,'classifiers/Sequence_by_Training_4Cell/', 'STUDYII_ReplayOnset'],'ToRall');
  replayScores = ToRall;
  MFStartup_studyII
end
nsub = length(is.goodsub);

for iSj=1:nsub
  RST = spm_eeg_load([datadir, sprintf('sfold_giles_symmetric_f_session%d', iSj*2)]);% use iSj*2 because iSj would refer to Preplay, we want replay
  fsample = RST.fsample;
  t_window = fsample/2;

  % get alpha phase
  [phase{iSj}, pow{iSj}, maxidx(iSj)] = getphasetimecourse(RST);
  phase{iSj} = phase{iSj}(:,triggerpoints{iSj});
  pow{iSj} = pow{iSj}(:,triggerpoints{iSj});
  
  % get the replay times (top 1% replay)
  replaytimes = replayScores(iSj,:) > prctile(replayScores(iSj,:),100-is.topPercentile);
  replaytimes= [0,diff(replaytimes)]==1;% eliminate adjacent points
  replaytimes(1:t_window)=0;replaytimes(end-t_window:end)=0;% eliminate border points:
  % upsample replay times
  replay_idx = upsample_replayTimes(replaytimes, phase{iSj}, t_window);

  % do the same for the control condition (lowest percentile replay)
    % get the replay times (top 1% replay)
  controltimes = replayScores(iSj,:) < prctile(replayScores(iSj,:),is.topPercentile);
  controltimes= [0,diff(controltimes)]==1;% eliminate adjacent points
  controltimes(1:t_window)=0;controltimes(end-t_window:end)=0;% eliminate border points:
  % upsample replay times
  control_idx = upsample_replayTimes(controltimes, phase{iSj}, t_window);

  % get same amount of points for replay and control
  n = min([numel(control_idx), numel(replay_idx)]);
  ix1=randperm(numel(control_idx));
  ix2=randperm(numel(replay_idx));
  control_idx = control_idx(ix1(1:n));
  replay_idx = replay_idx(ix2(1:n));

  % take the phase of these points
  phase_top{iSj} = phase{iSj}(:,replay_idx);
  phase_bottom{iSj} = phase{iSj}(:,control_idx);

  % compare concentration values
  for iParc=1:38
    [~,fval(iSj,iParc)] = circ_ktest(phase_top{iSj}(iParc,:), phase_bottom{iSj}(iParc,:));
  end
end
