% This script loads in probabilities of the decoding model applied to
% resting state data. Based on source parcellations, the hypothesis is
% tested whether the highest probabilities are phase clustered (e.g. in the
% theta/alpha band), w.r.t. what could be expected by chance (i.e. phase
% non-uniformity if a strictly positive measure, and is higher when there
% is strong coherence in the data). Chance is estimated by running 1000
% permutations, where in every permutation, the weights of the decoding
% models are permuted over stimuli.
% 2022 - Mats van Es

%% Define some variables
whichstudy = 2;
if whichstudy==1
  MFStartup;
else
  MFStartup_studyII;
  load('/Users/matsvanes/Data/YunzheData/StrLearn_MEGexp/BehavData/Exptrial.mat')
  is.dirname='StrLearn_MEGexp/';
  is.matsdir = '/Users/matsvanes/Data/YunzheData/Mats/Study2/';
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
is.nPerm = 100;
is.K = 12;
fsample=250;
is.plot=false;
is.t_window=fsample/2;
is.t=-.5:1/250:0.5;
is.useMask = 1; % to use a state-mask on the phases
parc=parcellation('fmri_d100_parcellation_with_PCC_reduced_2mm_ss5mm_ds8mm.nii.gz');
FOI = repmat(1:1:30, length(is.goodsub), 1); % can be 'alpha', 'emd'
if ischar(FOI)
  substring=['_', FOI];
else
  substring='';
end
load('Giles39_neighbours.mat')

%% Load HMM and get masks for RSN states 1 and 2
% This section loads in the HMM results, and trims the Gamma to be equal
% length to phase. Next, it can be used to mask the replay events that
% occur when specific states are on (i.e. Parietal Alpha (1), DMN (2)).

% get the HMM results
fname = fullfile(is.studydir, 'Neuron2020Analysis/', sprintf('Study%d',whichstudy), "hmm_1to45hz", "hmm5usingtemplate_parc_giles_symmetric__pcdim80_voxelwise_embed13_K12_big1_dyn_modelhmm.mat");
load(fname)
hmm = hmm_permutestates(hmm,new_state_ordering);
Gamma = hmm.gamma;

% and the timings
fname=fullfile(is.studydir,'Neuron2020Analysis/', sprintf('Study%d',whichstudy), 'hmm_1to45hz/hmm_parc_giles_symmetric__pcdim80_voxelwise_embed13.mat');
load(fname,'hmmT','subj_inds')
scan_T = cell2mat(hmmT);

% load in the timings
datadir = fullfile(is.studydir,'Neuron2020Analysis/', sprintf('Study%d',whichstudy), 'bfnew_1to45hz/');
[maskA,maskB,triggerpoints,goodsamples] = getSubjectMasks(datadir);

nSes = length(goodsamples);
K=size(Gamma,2);

if length(triggerpoints{1})>40000
  Fs = 250;
else
  Fs = 100;
end

t=[-is.t_window:is.t_window]./Fs;
method=1;

R = [[1;1+cumsum(scan_T(1:end-1))'],cumsum(scan_T(1:end))'];
if ~all(R(2:end,1)==R(1:end-1,2)+1)
  % correct any uncorrected offsets in R:
  R(2:end,1) = R(1:end-1,2)+1;
end

for iSes=1:nSes
  % convert Gamma to padded timecourse
  Gam_long = NaN(length(goodsamples{iSes}),K);
  Gam_long(goodsamples{iSes},:) = Gamma(R(iSes,1):R(iSes,2),:);
  Gam_long = Gam_long(triggerpoints{iSes},:);

  G(iSes,:,:) = Gam_long;

  temp1 = zeros(length(goodsamples{iSes}),1);
  temp1(goodsamples{iSes}) = hmm.statepath(R(iSes,1):R(iSes,2));
  temp1 = temp1(triggerpoints{iSes},:);
  vpath(iSes,:) = temp1;
end
vpath = vpath(2:2:end,:);
Gamma=G(is.iRun:2:end,:,:);


%% Get alpha peak frequency
if isstr(FOI) && strcmp(FOI,'alpha')
  for iSj=1:length(is.goodsub)
    sprintf('Determining peak frequency | Subject %d', iSj)
    RST = spm_eeg_load([datadir, sprintf('sfold_giles_symmetric_f_session%d', iSj*2)]);
    doplot=0;
    [~, peakfreq(iSj)] = compute_peakfreq(RST, 2:0.1:20, [],[], doplot);
  end
  FOI=peakfreq';
end

%% Get phases
[~,~,triggerpoints,goodsamples]=getSubjectMasks(datadir);
triggerpoints = triggerpoints(2:2:end); % only need post-task data
goodsamples = goodsamples(2:2:end);

for iSj=1:length(is.goodsub)
  sprintf('Get phase time course | Subject %d', iSj)

  if isstr(FOI) && strcmp(FOI, 'emd')
    tmp=load([is.matsdir, sprintf('phase/emd_session%d', iSj*2)]);
    IF(iSj, :, :) = squeeze(mean(tmp.IF(:,triggerpoints{iSj},:),2));
    phase{iSj} = permute(tmp.IP(:, triggerpoints{iSj},:), [1,3,2]);
  else
    RST = spm_eeg_load([datadir, sprintf('sfold_giles_symmetric_f_session%d', iSj*2)]);
    % get FOI phase
    if exist('peakfreq', 'var') && numel(FOI)==numel(is.goodsub)
      [phase{iSj}, ~, maxidx(iSj,:), tmp] = getphasetimecourse(RST, peakfreq(iSj));
    else
      [phase{iSj}, ~, maxidx(iSj,:), tmp] =  getphasetimecourse(RST, FOI(iSj,:)); % getphasetimecourse(RST, FOI(iSj,:),'dpss',2);
    end
    spctrm{iSj} = squeeze(tmp.fourierspctrm(:,:,:,triggerpoints{iSj}));
    phase{iSj} = phase{iSj}(:,:,triggerpoints{iSj});
  end
end

%% phase preference for reactivation (compared to shuffled classifier weights (over stimuli))
%{
for iSj = 1:length(is.goodsub)
  sprintf('Subject %d', iSj)
  is.doPermutation=false;
  GoodChannel = find_goodChannel(iSj, is); % get good channels
  gnf = train_classifier(iSj,is,Exptrial, GoodChannel); % get decoding models
  Rreds = apply_toRest(iSj, is, gnf, GoodChannel, is.iRun); % get all resting state decoding results
  RredsAll{iSj} = Rreds{1,37,1,is.lambda}; % select the decoding results of the correct model
  probAll{iSj} =  1./(1+exp(RredsAll{iSj})); % transform decoding results to probabilities
  % compute phase non-uniformity of the top percentile probabilities
  fname = [is.matsdir, 'phasePreference/', sprintf('phasePrefReactivation%d_%d', iSj, is.iRun), substring];
  if is.useMask 
    fname = [fname, sprintf('_mask_%dk', is.useMask)];
    temporalMask = vpath(iSj,:)==is.useMask;
    is.topPercentile=5;
  else
    temporalMask=[];
  end
  [topphase{iSj}, ZtopPhase(iSj,:,:)] = phase_nonuniformity(probAll{iSj}, phase{iSj}, is, fname, temporalMask);


  % do the same for permutations, where channel weights from the decoding
  % model are permuted over stimuli
  is.doPermutation=true;
  fname = [fname, '_perm'];
  if is.usePrecomputed && isfile([fname, '.mat'])
    tmp=load(fname);
    ZtopPhasePermAll(iSj,:,:,:) = tmp.ZtopPhasePerm;
  else
    for iPerm=1:is.nPerm
      gnfperm=permute_classifierweights(gnf,is);
      Rreds = apply_toRest(iSj, is, gnfperm, GoodChannel, is.iRun, iPerm); % get all Rreds
      probPerm = 1./(1+exp(Rreds{1,37,1,is.lambda}));
      [~, ZtopPhasePerm(:,:,iPerm)] = phase_nonuniformity(probPerm, phase{iSj}, is, [], temporalMask);
    end
    if ischar(FOI)
      foi=squeeze(IF(iSj,:,:));
    else
      foi=FOI(iSj,:);
    end
    save(fname, 'ZtopPhasePerm', 'foi')
    ZtopPhasePermAll(iSj,:,:,:) = ZtopPhasePerm;
  end
end

ZtopPhasePerm = mean(ZtopPhasePermAll,4);

phasepref_reactivation = ZtopPhase-ZtopPhasePerm;
WGA_phasepref_reactivation = phasepref_reactivation./std(ZtopPhasePermAll,[],4);
GA_phasepref_reactivation = squeeze(mean(phasepref_reactivation));

% Do signflipping test on the group level
dat=[];
dat.label=parc.labels;
if ischar(FOI)
  dat.freq = squeeze(mean(mean(IF)));
else
  dat.freq=FOI(1,:);
end
dat.dimord = 'rpt_chan_freq';
dat.z = ZtopPhase;

datperm=dat;
datperm.z = ZtopPhasePerm;

Nsub=length(is.goodsub);
cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfgs.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfgs.tail = 1;
% cfgs.correctm = 'cluster';%'bonferroni';
cfgs.clustertail = 1;
cfgs.neighbours = [];
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.numrandomization=10000;
cfgs.parameter='z';
stat=ft_freqstatistics(cfgs,dat,datperm);
stat.avg = GA_phasepref_reactivation;
stat.weighted_avg = squeeze(mean(WGA_phasepref_reactivation));

filename = [is.matsdir, 'phasePreference/', 'phasepreference_reactivation', substring];
if is.useMask, filename = [filename, sprintf('_mask_%dk', is.useMask)]; end
save(filename, 'stat', 'cfgs', 'phasepref_reactivation', 'WGA_phasepref_reactivation')

% Visualize mean difference
filename = [is.matsdir 'phasePreference/', 'phasepreference_reactivation', substring, '_GA'];
if is.useMask, filename = [filename, sprintf('_mask_%dk', is.useMask)]; end
parc.savenii(GA_phasepref_reactivation, filename);
plot_surface_movie(parc, stat, 'avg', [], false)

filename = [is.matsdir 'phasePreference/', 'phasepreference_reactivation', substring, '_WGA'];
if is.useMask, filename = [filename, sprintf('_mask_%dk', is.useMask)]; end
parc.savenii(WGA_phasepref_reactivation, filename);
plot_surface_movie(parc, stat, 'weighted_avg', [], false)

filename = [is.matsdir 'phasePreference/', 'phasepreference_reactivation', substring, '_tstat'];
if is.useMask, filename = [filename, sprintf('_mask_%dk', is.useMask)]; end
parc.savenii(stat.stat, filename);
plot_surface_movie(parc, stat, 'stat', [], false)
%}
%% phase preference for replay (top vs permutation)
% Do the same as above, but now for replay
is.doPermutation=false;
% observed data
if whichstudy==1
  tmp=load([wd,'GenericReplayData/STUDYI_ReplayOnset/RwdReplay2ndOnset'],'ToRall');
  replayScores = tmp.ToRall;
else
  StimAll = load([is.AnalysisPath, 'classifiers/Sequence_by_Training_4Cell/StimAll.mat']);
  sfAll = StimAll.sfAll;
  sbAll = StimAll.sbAll;
  [~, ~, Mreplay, ~] = find_bestLambda(is, sfAll, sbAll);
  ToRall = compute_replayTimecourse(is, [], Mreplay, is.iRun);
  replayScores = ToRall;
end

replayScoresPerm = zeros(22,30000,is.nPerm);
for iSj=1:length(is.goodsub)
  try
    fname = [is.matsdir, 'itpc/', sprintf('itpc%d_%d', iSj, is.iRun), substring];
    if is.useMask, fname = [fname, sprintf('_mask_%dk', is.useMask)]; end
    load(fname);
    itpc{iSj} = tl.itpc;
    pow{iSj} = tl.pow;
  catch
    % compute phase non-uniformity of the top percentile Replay scores
    fname = [is.matsdir, 'itpc/', sprintf('itpc%d_%d', iSj, is.iRun), substring];
    if is.useMask
      fname = [fname, sprintf('_mask_%dk', is.useMask)];
      temporalMask = vpath(iSj,:)==is.useMask;
      is.topPercentile=5;
    else
      temporalMask=[];
    end
    [~, topidx{iSj}] = get_percentilePhase(replayScores(iSj,:), phase{iSj}, is.t_window, is.topPercentile, temporalMask);
    %   [topphase{iSj}, ZtopPhase(iSj,:,:)] = phase_nonuniformity(replayScores(iSj,:)', phase{iSj}, is, fname, temporalMask);
    [itpc{iSj}, pow{iSj}] = compute_itpc(spctrm{iSj}, topidx{iSj});
    tl=[];
    tl.itpc=itpc{iSj};
    tl.pow=pow{iSj};
    tl.t=pow;
    if ischar(FOI)
      tl.f=squeeze(IF(iSj,:,:));
    else
      tl.f=FOI(iSj,:);
    end
    save(fname, 'tl')
  end
end

% do the same for permutations, where channel weights from the decoding
% model are permuted over stimuli
is.doPermutation=true;
try
  for iSj=1:length(is.goodsub)
    fname = [is.matsdir, 'itpc/', sprintf('itpc%d_%d', iSj, is.iRun), substring];
    if is.useMask, fname = [fname, sprintf('_mask_%dk', is.useMask)]; end
    fname = [fname, '_perm'];
    load(fname);
    itpcPerm_u(iSj,:,:,:) = tl.itpc_u;
    itpcPerm_std(iSj,:,:,:) = tl.itpc_std;
    powPerm_u(iSj,:,:,:) = tl.pow_std;
    powPerm_std(iSj,:,:,:) = tl.pow_std;
  end
catch
  % We need to use the lag that is used in the observed statistic. In
  % order to find it, we use Mreplay (group mean replay) to feed into
  % compute_sequenceness
  if is.useMask
    temporalMask = vpath(iSj,:)==is.useMask;
    is.topPercentile=5;
  else
    temporalMask=[];
  end
  is.plot=false;
  StimAll = load([is.AnalysisPath, 'classifiers/Sequence_by_Training_4Cell/StimAll.mat']);
  sfAll = StimAll.sfAll;
  sbAll = StimAll.sbAll;
  [~, ~, Mreplay, ~] = find_bestLambda(is, sfAll, sbAll);
  for iSj=1:length(is.goodsub)
    itpcPerm = zeros(size(spctrm{iSj},1), size(spctrm{iSj},2), 251, is.nPerm);
    powPerm = zeros(size(spctrm{iSj},1), size(spctrm{iSj},2), 251, is.nPerm);
    for iPerm=1:is.nPerm
      sprintf('%d | %d', iSj, iPerm)
      % The permuted Rreds should already have been computed, so no need to
      % do the stuff that comes before. This is done for all subjects at once
      % though, and we want to save the permutation statistic per subject -
      % so that's why we use a seperate loop.
      replayScoresPerm = compute_replayTimecourse(is,[],Mreplay,is.iRun,iPerm);
      [~, topidxPerm] = get_percentilePhase(replayScoresPerm, phase{iSj}, is.t_window, is.topPercentile, temporalMask);
      [itpcPerm(:,:,:,iPerm), powPerm(:,:,:,iPerm)] = compute_itpc(spctrm{iSj}, topidxPerm);
    end
    if ischar(FOI)
      foi=squeeze(IF(iSj,:,:));
    else
      foi=FOI(iSj,:);
    end
    fname = [is.matsdir, 'itpc/', sprintf('itpc%d_%d', iSj, is.iRun), substring];
    if is.useMask, fname = [fname, sprintf('_mask_%dk', is.useMask)]; end
    fname = [fname, '_perm'];
    tl=[];
    tl.itpc_u = mean(itpcPerm,4);
    tl.itpc_std = std(itpcPerm,[],4);
    tl.pow_u = mean(powPerm,4);
    tl.pow_std = std(powPerm,[],4);
    tl.f=foi;
    tl.t=is.t;
    save(fname, "tl")
    itpcPerm_u(iSj,:,:,:) = tl.itpc_u;
    itpcPerm_std(iSj,:,:,:) = tl.itpc_std;
    powPerm_u(iSj,:,:,:) = tl.pow_std;
    powPerm_std(iSj,:,:,:) = tl.pow_std;
  end
end
%
itpc = permute(cat(4,itpc{:}), [4,1,2,3]);

itpc_diff = itpc-itpcPerm_u;
WGA_itpc = squeeze(mean(itpc_diff./itpcPerm_std));
GA_itpc = squeeze(mean(itpc_diff));

% Do signflipping test on the group level
dat=[];
dat.label=parc.labels;
dat.time = is.t;
if ischar(FOI)
  dat.freq = squeeze(mean(mean(IF)));
else
  dat.freq=FOI(1,:);
end
dat.dimord = 'rpt_chan_freq_time';
dat.z = itpc; % itpc./itpcPerm_std

datperm=dat;
datperm.z = itpcPerm_u; % ZtopPhasePerm./std(ZtopPhasePermAll,[],4)

if is.useMask
  spatialMask = get_spatialHMMmask;
  spatialMask = spatialMask(:,is.useMask);
  dat.z(:,isnan(spatialMask),:,:) = NaN;
  datperm.z(:,isnan(spatialMask),:,:) = NaN;
  nparc_correct = nansum(spatialMask);
else
  nparc_correct = 38;
  spatialMask = ones(38,1);
end

Nsub=length(is.goodsub);
cfgs=[];
cfgs.method = 'montecarlo';
cfgs.statistic = 'depsamplesT';
cfgs.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfgs.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfgs.tail = 1;
cfgs.correctm = 'cluster';
cfgs.clustertail = 1;
cfgs.neighbours = [];
cfgs.ivar = 1;
cfgs.uvar = 2;
cfgs.numrandomization=1000;
cfgs.parameter='z';
cfgs.alpha=0.05/nparc_correct;
cfgs.clusteralpha = 0.05;
osl_shutdown()
addpath(('/Volumes/T5_OHBA/software/fieldtrip'))
ft_defaults
stat=ft_freqstatistics(cfgs,dat,datperm);
restoredefaultpath
startup
stat.avg = GA_itpc;
stat.weighted_avg = WGA_itpc;

filename = [is.matsdir 'itpc/', 'itpc_replay', substring];
if is.useMask, filename = [filename, sprintf('_mask_%dk', is.useMask)]; end
save(filename, 'stat', 'cfgs', 'itpc_diff')

% Visualize mean difference
filename = [is.matsdir 'itpc/', 'phasepreference_replay', substring, '_GA'];
if is.useMask, filename = [filename, sprintf('_mask_%dk', is.useMask)]; end
parc.savenii(GA_itpc, filename);
plot_surface_movie(parc, stat, 'avg', [], false)

filename = [is.matsdir 'itpc/', 'phasepreference_replay', substring, '_WGA'];
if is.useMask, filename = [filename, sprintf('_mask_%dk', is.useMask)]; end
parc.savenii(WGA_itpc, filename);
plot_surface_movie(parc, stat, 'weighted_avg', [], false)

filename = [is.matsdir 'itpc/', 'phasepreference_replay', substring, '_tstat'];
if is.useMask, filename = [filename, sprintf('_mask_%dk', is.useMask)]; end
parc.savenii(stat.stat, filename);
plot_surface_movie(parc, stat, 'stat', [], false)
%}
%% simple erp analysis


for iSj=1:length(is.goodsub)
[~, topidx{iSj}] = get_percentilePhase(replayScores(iSj,:), phase{iSj}, is.t_window, is.topPercentile, temporalMask);
end

for iSj=1:22
  RST = spm_eeg_load([datadir, sprintf('sfold_giles_symmetric_f_session%d', iSj*2)]);
  RST=RST(:, find(triggerpoints{iSj}));

  tl{iSj} = zeros(38,501,length(topidx{iSj}));
  for k=1:length(topidx{iSj})
    tl{iSj}(:,:,k) = RST(:,topidx{iSj}(k)-250:topidx{iSj}(k)+250);
  end
  tl{iSj} = mean(tl{iSj},3);
end
tl_long = cat(3,tl{:});
tl=tl_long(:,126:end-125,:);
figure; for k=1:38, subplot(4,10,k), shadedErrorBar(t, mean(tl(k,:,:)-mean(tl(k,:,:),2),3), std(tl(k,:,:)-mean(tl(k,:,:),2),[],3)), end

%% spectral power analysis
for iSj=1:22
  X{iSj} = zeros(38,30, 251,length(topidx{iSj}));
  for k=1:length(topidx{iSj})
    X{iSj}(:,:,:,k) = abs(spctrm{iSj}(:,:,topidx{iSj}-125:topidx{iSj}+125)).^2;
  end
  X{iSj} = mean(X{iSj},4);
end
X = cat(4,X{:});

figure; for k=1:38, subplot(4,10,k), imagesc(t,1:30,squeeze(mean(X(k,:,:,:),4))), axis xy, colormap(inferno), end

% power spectrum of ERP
cfg=[];
cfg.method = 'mtmconvol';
cfg.taper='hanning';
cfg.foi=1:30;
cfg.t_ftimwin = 2*1./cfg.foi;
cfg.output = 'pow';
cfg.toi = -1:1/250:1;
cfg.pad=4;

dat=rmfield(dat, {'freq', 'z'});
dat.dimord = 'rpt_chan_time';
dat.trial = permute(tl_long, [3,1,2]);
dat.time = -1:1/250:1;
X_erp = ft_freqanalysis(cfg, dat);

figure; for k=1:38, subplot(4,10,k), imagesc(t, 1:30, squeeze(X_erp.powspctrm(k,:,126:end-125))), axis xy, colormap(inferno), title(sprintf('%d',k)), end


