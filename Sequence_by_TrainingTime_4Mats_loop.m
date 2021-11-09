% These settings were the ones used in Yunzhe's Cell paper. L1 regression,
% training on one time point only.
%% Train decoders - stim code
%clear all;
% data info and directories
whichstudy = 2;
if whichstudy==1
  MFStartup;
else
  MFStartup_studyII;
  load('/Volumes/T5_OHBA/data/replay/StrLearn_MEGexp/BehavData/Exptrial.mat')
  is.dirname='StrLearn_MEGexp/';
end
if ~isfolder(is.AnalysisPath),   mkdir(is.AnalysisPath), end
is.usePrecomputed = true;

% variables
is.whenInMS = is.tss*is.msPerSample;
is.TOI = 200; % 200ms post stimulus.
is.sessionsOfInterest = [1,2,3,4,8]; % rst, lci, lci, lci, rst
is.nStim = 8;
is.ntrlAll=120; % per session (includes upside down images)
is.ntrlValid=96; % per session
is.nSes = 3;
is.nChan = 273;
is.nTime = 401;
is.iRun = [2];


sfAll = nan(length(is.lgncy)+1, is.nShuf, length(is.goodsub), 1, length(is.ENL1));
sbAll = nan(length(is.lgncy)+1, is.nShuf, length(is.goodsub), 1, length(is.ENL1));

%% Analysis pipeline
for iSj=1:length(is.goodsub)
  sprintf('%d',iSj)
  iSjGood=is.goodsub(iSj);

  % Find Good Channels in all sessions of interest
  GoodChannel = find_goodChannel(iSj, is);

  % train classifier
  gnf = train_classifier(iSj,is,Exptrial, GoodChannel);

  % apply classifier to Resting State data
  % compute for pre-task (iRun=1), and post task (iRun=2)
  cnt=1;
  for i=is.iRun
    Rreds{cnt} = apply_toRest(iSj, is, gnf, GoodChannel, i); 
    cnt=cnt+1;
  end

  if is.usePrecomputed
    StimAll = load([is.AnalysisPath, 'classifiers/Sequence_by_Training_4Cell/StimAll.mat']);
    sfAll = StimAll.sfAll;
    sbAll = StimAll.sbAll;
  else
    % compute sequenceness, only on iRun=2
    [sf, sb] = compute_sequenceness(iSj, is, 2);
    sfAll(:,:,iSj,:,:) = sf;
    sbAll(:,:,iSj,:,:) = sb;
  end
end

if ~is.usePrecomputed
  save([is.AnalysisPath,'classifiers/Sequence_by_Training_4Cell/StimAll'],...
    'sfAll','sbAll','is','-v7.3');
end


%% Replay
% Plot replay effect broken by training time
plot_replayEffect

% cross-validate to choose a single lambda parameter value:
[iLCV, sfbCV, Mreplay, Sreplay] = find_bestLambda(is, sfAll, sbAll);

% Compute replay time course
for i=is.iRun
  ToRall{i} = compute_replayTimecourse(is,iLCV,Mreplay,i);
end

% plot for all lambdas
plot_sequenceness
