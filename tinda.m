%% Load HMM and get masks for RSN states 1 and 2
% This section loads in the HMM results, and computes the long term
% asymmetry 

whichstudy = 2;
if whichstudy==1
  MFStartup;
else
  MFStartup_studyII;
  load('/Users/matsvanes/Data/YunzheData/StrLearn_MEGexp/BehavData/Exptrial.mat')
  is.dirname='StrLearn_MEGexp/';
  is.matsdir = '/Users/matsvanes/Data/YunzheData/Mats/Study2/';
end
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
is.K = 12;
is.plot=false;
fsample=250;
is.t_window=fsample/2;

addpath(genpath('/Users/matsvanes/Documents/Werk/scripts/Tinda/'))

% get the HMM results
fname = fullfile(is.studydir, 'Neuron2020Analysis/', sprintf('Study%d',whichstudy), "hmm_1to45hz", "hmm5usingtemplate_parc_giles_symmetric__pcdim80_voxelwise_embed13_K12_big1_dyn_modelhmm.mat");
load(fname)
hmm = hmm_permutestates(hmm,new_state_ordering);
Gamma = hmm.gamma;
K=hmm.K;

% and the timings
fname=fullfile(is.studydir,'Neuron2020Analysis/', sprintf('Study%d',whichstudy), 'hmm_1to45hz/hmm_parc_giles_symmetric__pcdim80_voxelwise_embed13.mat');
load(fname,'hmmT','subj_inds')

% load in the subject masks
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


scan_T = cell2mat(hmmT);
R = [[1;1+cumsum(scan_T(1:end-1))'],cumsum(scan_T(1:end))'];
if ~all(R(2:end,1)==R(1:end-1,2)+1)
    % correct any uncorrected offsets in R:
    R(2:end,1) = R(1:end-1,2)+1;
end

for iSes=1:nSes
  vpath{iSes} = hmm.statepath(R(iSes,1):R(iSes,2));
end

vpath=vpath(is.iRun:2:end);
hmmT=hmmT(is.iRun:2:end);

[FO,pvals,t_intervals] = computeLongTermAsymmetry(vpath,hmmT,K);

bonf_ncomparisons = K.^2-K;
mean_direction = squeeze(mean(FO(:,:,1,:)-FO(:,:,2,:),4));
mean_assym = squeeze(mean((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3),4));

hmm_1stlevel.FO_intervals = FO;
hmm_1stlevel.FO_assym = squeeze((FO(:,:,1,:)-FO(:,:,2,:))./mean(FO,3));

figure; 
x=permute(hmm_1stlevel.FO_assym, [3,1,2]);
x = x(:,~eye(12));
x = mean(abs(x),2);
distributionPlot(x)

%%
% this section determines the optimal state ordering for a circular plot; it
% then determines whether such a sequentially organised network could arise
% by chance by random shuffles of the rows of the transmat

optimalseqfile = [is.matsdir,'tinda/', 'bestseq.mat'];
if ~isfile(optimalseqfile)
    bestsequencemetrics = optimiseSequentialPattern(FO);
    save(optimalseqfile,'bestsequencemetrics');
else
    load(optimalseqfile);
end
bestseq = bestsequencemetrics{2};
% save each subject's rotational strength by this metric:
disttoplot_manual = zeros(12,1);
for i3=1:12
    disttoplot_manual(bestseq(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end
angleplot = exp(sqrt(-1)*(angle(disttoplot_manual.')-angle(disttoplot_manual)));

rotational_momentum = imag(sum(sum(angleplot.*hmm_1stlevel.FO_assym)));
hmm_1stlevel.rotational_momentum = squeeze(rotational_momentum);


% plot as circular diagram:
cyclicalstateplot(bestseq,mean_direction,pvals<(0.05));

%% Vpath Circle analysis
triggerpoints=triggerpoints(2:2:end);
goodsamples=goodsamples(2:2:end);
for i=1:22
  temp1 = zeros(length(goodsamples{i}),1);
  temp1(goodsamples{i}) = getCircleVpath(vpath{i}, bestseq);
  vpathcircle{i} = double(temp1(triggerpoints{i}));
  trialonset = zeros(75000,1);
  trialonset(topidx{i})=1;
  [tmp, vpathcircle_freq, vpathcircle_time] = fourieranalysis_circleVpath(vpathcircle{i}, trialonset);
  circlepow_evoked(i,:,:) = squeeze(nanmean(abs(tmp).^2));
  circlespctrm_evoked(i,:,:) = squeeze(nanmean(tmp));
end

for i=1:22
for k=1:length(topidx{i})
Q{i}(k,:) = vpathcircle{i}(topidx{i}(k)-250:topidx{i}(k)+250);
end
end
for k=1:22
Q2(k,:) = nanmean(Q{k});
end


%% Within a window, sum over the positions in the unit circle and do spectral analysis on this
% get replay scores
% observed data

is.doPermutation=false;
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

cnt=1;
for iSes=is.iRun:2:nSes 
    temp1 = zeros(length(goodsamples{iSes}),1);
    temp1(goodsamples{iSes}) = hmm.statepath(R(iSes,1):R(iSes,2));
    temp1 = temp1(triggerpoints{iSes},:);
    vpath{cnt} = temp1;
    cnt=cnt+1;
end


% get top percentile probabilities
topperc = replayScores > prctile(replayScores(:),100-is.topPercentile);
topperc = [zeros(size(topperc,1),1),diff(topperc,[],2)]==1;% eliminate adjacent points
topperc(:,1:is.t_window)=0;topperc(:,end-is.t_window:end)=0;% eliminate border points:



fs=250;
W = 125; %set arbitrary half second window for averaging
Noverlap = 1; % number of points to overlap over successive windows
t_last = 0;

% get the circle plot positions
optimalseqfile = [is.matsdir,'tinda/', 'bestseq.mat'];
load(optimalseqfile);
bestseq = bestsequencemetrics{2};

% save each subject's rotational strength by this metric:
disttoplot_manual = zeros(12,1);
for i3=1:12
  disttoplot_manual(bestseq(i3)) = exp(sqrt(-1)*i3/12*2*pi);
end

for iSj=1:22
 if 0
   GoodChannel = find_goodChannel(iSj, is); % get good channels
   gnf = train_classifier(iSj,is,Exptrial, GoodChannel); % get decoding models
   Rreds = apply_toRest(iSj, is, gnf, GoodChannel, is.iRun); % get all resting state decoding results
   RredsAll = Rreds{1,37,1,is.lambda}; % select the decoding results of the correct model
   replayScores =  1./(1+exp(RredsAll'));
   topperc = replayScores > prctile(replayScores(:),100-is.topPercentile);
   topperc = [zeros(size(topperc,1),1),diff(topperc,[],2)]==1;% eliminate adjacent points
   topperc(:,1:is.t_window)=0;topperc(:,end-is.t_window:end)=0;% eliminate border points:
   topidx = upsample_replayTimes(topperc, zeros(1,1,75000), is.t_window);
 else

  topidx = upsample_replayTimes(topperc(iSj,:), zeros(1,1,75000), is.t_window);
 end
  vpcircle = vpath{iSj};
  for k=1:12
    vpcircle(vpcircle==k) = disttoplot_manual(k);
  end
  
  % smooth vpcircle with half second window
  vpcircle=smoothdata(vpcircle,'movmean',W);

  % get replay locked phase information
  vpcircle_timelocked{iSj} = zeros(length(topidx), Fs+1);
  for i=1:length(topidx)
    vpcircle_timelocked{iSj}(i,:) = vpcircle(topidx(i)-W:topidx(i)+W);
  end
  vpcircle_timelocked_group(iSj,:) = mean(vpcircle_timelocked{iSj});
end


figure; for k=1:22
subplot(3,8,k), plot(t,angle(mean(vpcircle_timelocked{k}))), ylim([-pi, pi]), title(sprintf('subj %d', k)), xlabel('time (s)')
end
suptitle('cyclical state phase timelocked to replay')
figure; plot(t, angle(mean(vpcircle_timelocked_group)))
title('cyclical state phase timelocked to replay - GROUP AVERAGE'), ylim([-pi,pi]), xlabel('time (s)')


