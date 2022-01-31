% This script loads in probabilities of the decoding model applied to
% resting state data. Based on source parcellations, the hypothesis is
% tested whether the highest probabilities are phase clustered (e.g. in the
% theta/alpha band), w.r.t. the lowest probabilities (this is a good
% control because classifiers have a trivial phase preference).
% 2021 - Mats van Es

%% Define some variables
whichstudy = 2;
if whichstudy==1
  MFStartup;
else
  MFStartup_studyII;
  load('/Users/matsvanes/Data/YunzheData/StrLearn_MEGexp/BehavData/Exptrial.mat')
  is.dirname='StrLearn_MEGexp/';
  savebase = '/Users/matsvanes/Data/YunzheData/StrLearn_MEGexp/Analysis/phasePreference/';
end
datadir = [is.studydir,'Neuron2020Analysis/Study2/bfnew_1to45hz/'];
if ~isfolder(is.AnalysisPath),   mkdir(is.AnalysisPath), end
is.usePrecomputed = true;
is.iRun = 2;
is.topPercentile=1;
is.lambda = 5; % hardcoded best lambda for all subjects
is.whenInMS = is.tss*is.msPerSample;
is.TOI = 200; % 200ms post stimulus.
is.sessionsOfInterest = [1,2,3,4,8]; % rst, lci, lci, lci, rst
is.nStim = 8;
is.ntrlAll=120; % per session (includes upside down images)
is.ntrlValid=96; % per session
is.nSes = 3;
is.nChan = 273;
is.nTime = 401;
is.iRun = 2;
is.nPerm = 1000;
fsample=250;
is.t_window=fsample/2;
parc=parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
FOI = 'alpha'; % can be 'alpha', 'theta'

%% Get alpha peak frequency
if strcmp(FOI,'alpha')
  for iSj=1:length(is.goodsub)
    sprintf('Determining peak frequency | Subject %d', iSj)
    RST = spm_eeg_load([datadir, sprintf('sfold_giles_symmetric_f_session%d', iSj*2)]);
    doplot=0;
    [~, peakfreq(iSj)] = compute_peakfreq(RST, 2:0.1:20, [],[], doplot);
  end
end

%% Get phases
[~,~,triggerpoints,goodsamples]=getSubjectMasks(datadir);
triggerpoints = triggerpoints(2:2:end); % only need post-task data
goodsamples = goodsamples(2:2:end);

for iSj=1:length(is.goodsub)
  sprintf('Get phase time course | Subject %d', iSj)
  RST = spm_eeg_load([datadir, sprintf('sfold_giles_symmetric_f_session%d', iSj*2)]);

  % get FOI phase
  if strcmp(FOI, 'alpha')
    [phase{iSj}, pow{iSj}, maxidx(iSj)] = getphasetimecourse(RST, peakfreq(iSj));
  elseif strcmp(FOI,'theta')
    [phase{iSj}, pow{iSj}, maxidx(iSj)] = getphasetimecourse(RST, 6,'dpss',2);
  end
  phase{iSj} = phase{iSj}(:,triggerpoints{iSj});
  pow{iSj} = pow{iSj}(:,triggerpoints{iSj});
end

%% phase preference for reactivation (compared to shuffled classifier weights (over stimuli))

for iSj = 1:length(is.goodsub)
  is.doPermutation=false;
  GoodChannel = find_goodChannel(iSj, is); % get good channels
  gnf = train_classifier(iSj,is,Exptrial, GoodChannel); % get decoding models
  iRun=2;
  Rreds = apply_toRest(iSj, is, gnf, GoodChannel, iRun); % get all resting state decoding results
  RredsAll{iSj} = Rreds{1,37,1,is.lambda}; % select the decoding results of the correct model
  probAll{iSj} =  1./(1+exp(RredsAll{iSj})); % transform decoding results to probabilities

  % compute phase non-uniformity of the top percentile probabilities
  fname = [is.AnalysisPath, 'phasePreference/', sprintf('phasePrefReactivation%d_%d', iSj, iRun)];
  [topphase{iSj}, ZtopPhase(iSj,:)] = phase_nonuniformity(probAll{iSj}, phase{iSj}, is, fname);


  % do the same for permutations, where channel weights from the decoding
  % model are permuted over stimuli
  is.doPermutation=true;
  fname = [fname, '_perm'];
  if is.usePrecomputed && isfile([fname, '.mat'])
    tmp=load(fname);
    topphasePermAll{iSj} = tmp.topphasePerm;
    ZtopPhasePermAll(iSj,:,:) = tmp.ZtopPhasePerm;
  else
    for iPerm=1:is.nPerm
      gnfperm=permute_classifierweights(gnf,is);
      Rreds = apply_toRest(iSj, is, gnfperm, GoodChannel, iRun, iPerm); % get all Rreds
      probPerm = 1./(1+exp(Rreds{1,37,1,is.lambda}));
      [topphasePerm{iPerm}, ZtopPhasePerm(:,iPerm)] = phase_nonuniformity(probPerm, phase{iSj}, is, []);
    end
    save(fname, 'topphasePerm', 'ZtopPhasePerm')
    topphasePermAll{iSj} = topphasePerm;
    ZtopPhasePermAll(iSj,:,:) = ZtopPhasePerm;
  end
end

ZtopPhasePerm = mean(ZtopPhasePermAll,3);

phasepref_reactivation = ZtopPhase-ZtopPhasePerm;
GA_phasepref_reactivation = mean(phasepref_reactivation);

% Do signflipping test on the group level
dat=[];
dat.label=parc.labels;
dat.time=0;
dat.dimord = 'rpt_chan_time';
dat.avg = ZtopPhase;

datperm=dat;
datperm.avg = ZtopPhasePerm;

Nsub=length(is.goodsub);
cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfgs.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfgs.correctm = 'bonferroni';
cfgs.tail = 1;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.numrandomization=10000;
stat=ft_timelockstatistics(cfgs,dat,datperm);

% Visualize mean difference
figure
filename = [savebase 'phasepreference_reactivation_', sprintf('%s', FOI)];
parc.savenii(GA_phasepref_reactivation, filename);
osl_render4D([filename, '.nii.gz'], 'visualise', true);

filename = [savebase 'phasepreference_reactivation_', sprintf('%s', FOI), '_tstat'];
parc.savenii(stat.stat, filename);
osl_render4D([filename, '.nii.gz'], 'visualise', true);

%% phase preference for replay (top vs permutation)
% Do the same as above, but now for replay

% observed data
if whichstudy==1
  tmp=load([wd,'GenericReplayData/STUDYI_ReplayOnset/RwdReplay2ndOnset'],'ToRall');
  replayScores = tmp.ToRall;
else
  tmp=load([is.AnalysisPath,'classifiers/Sequence_by_Training_4Cell/', 'STUDYII_ReplayOnset'],'ToRall');
  replayScores = tmp.ToRall;
end

replayScoresPerm = zeros(22,30000,is.nPerm);
for iSj=1:length(is.goodsub)
  is.doPermutation=false;
  % compute phase non-uniformity of the top percentile Replay scores
  fname = [is.AnalysisPath, 'phasePreference/', sprintf('phasePrefReplay%d_%d', iSj, iRun)];
  [topphase{iSj}, ZtopPhase(iSj,:)] = phase_nonuniformity(replayScores(iSj,:)', phase{iSj}, is, fname);
end

% do the same for permutations, where channel weights from the decoding
% model are permuted over stimuli
is.doPermutation=true;
if is.usePrecomputed && isfile([fname, '.mat'])
  tmp=load(fname);
  topphasePermAll{iSj} = tmp.topphasePerm;
  ZtopPhasePermAll(iSj,:,:) = tmp.ZtopPhasePerm;
else
  % We need to use the lag that is used in the observed statistic. In
  % order to find it, we use Mreplay (group mean replay) to feed into
  % compute_sequenceness
  StimAll = load([is.AnalysisPath, 'classifiers/Sequence_by_Training_4Cell/StimAll.mat']);
  sfAll = StimAll.sfAll;
  sbAll = StimAll.sbAll;
  [~, ~, Mreplay, ~] = find_bestLambda(is, sfAll, sbAll);
  is.plot=false;
  for iPerm=1:is.nPerm
    % The permuted Rreds should already have been computed, so no need to
    % do the stuff that comes before. This is done for all subjects at once
    % though, and we want to save the permutation statistic per subject -
    % so that's why we use a seperate loop.
    replayScoresPerm(:,:,iPerm) = compute_replayTimecourse(is,[],Mreplay,iRun,iPerm);
  end
  for iSj=1:length(is.goodsub)
    for iPerm=1:is.nPerm
      [topphasePerm{iPerm}, ZtopPhasePerm(:,iPerm)] = phase_nonuniformity(squeeze(replayScoresPerm(iSj,:,iPerm))', phase{iSj}, is, []);
    end
    fname = [is.AnalysisPath, 'phasePreference/', sprintf('phasePrefReplay%d_%d', iSj, iRun)];
    fname = [fname, '_perm'];
    topphasePermAll{iSj} = topphasePerm;
    ZtopPhasePermAll(iSj,:,:) = ZtopPhasePerm;
    save(fname, 'topphasePerm', 'ZtopPhasePerm')
  end
end

ZtopPhasePerm = mean(ZtopPhasePermAll,3);

phasepref_replay = ZtopPhase-ZtopPhasePerm;
GA_phasepref_replay = mean(phasepref_replay);

% Do signflipping test on the group level
dat=[];
dat.label=parc.labels;
dat.time=0;
dat.dimord = 'rpt_chan_time';
dat.avg = ZtopPhase;

datperm=dat;
datperm.avg = ZtopPhasePerm;

Nsub=length(is.goodsub);
cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfgs.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfgs.correctm = 'bonferroni';
cfgs.tail = 1;
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.numrandomization=10000;
stat=ft_timelockstatistics(cfgs,dat,datperm);

% Visualize mean difference
filename = [savebase 'phasepreference_replay_', sprintf('%s', FOI)];
parc.savenii(GA_phasepref_replay, filename);
osl_render4D([filename, '.nii.gz'], 'visualise', true);

filename = [savebase 'phasepreference_replay_', sprintf('%s', FOI), '_tstat'];
parc.savenii(stat.stat, filename);
osl_render4D([filename, '.nii.gz'], 'visualise', true);
