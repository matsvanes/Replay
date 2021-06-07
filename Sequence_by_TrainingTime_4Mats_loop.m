% These settings were the ones used in Yunzhe's Cell paper. L1 regression,
% training on one time point only.
%% Train decoders - stim code
%clear all;
% data info and directories
whichstudy = 2;
if whichstudy==1
  data_dir = '/Volumes/T5_OHBA/data/replay/Replaydata4Cam/';
  MFStartup;
else
  data_dir = '/Volumes/T5_OHBA/data/replay/StrLearn_MEGexp';
  MFStartup_studyII;
  load('/Volumes/T5_OHBA/data/replay/StrLearn_MEGexp/BehavData/Exptrial.mat')
end
analysisdir = '/Volumes/T5_OHBA/analysis/replay/';
if ~isfolder(analysisdir),   mkdir(analysisdir), end
dirname='L1regression/';

% variables
whenInMS = is.tss*is.msPerSample;
TOI = 0.200; % 200ms post stimulus.
sessionsOfInterest = [1,2,3,4,8]; % rst, lci, lci, lci, rst
nStim = 8;
ntrlAll=120; % per session (includes upside down images)
ntrlValid=96; % per session
nSes = 3;

% Find Good Channels in all sessions of interest
GoodChannel=nan(length(is.goodsub),273,numel(sessionsOfInterest)); %nsubject*nsensors*nruns
for iSj=1:length(is.goodsub)
  iSjGood=is.goodsub(iSj);
  cnt=1;
  for isession=sessionsOfInterest
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
    
    index=ones(1,273);
    index(ind)=0;
    GoodChannel(iSj,:,cnt)=index;
    cnt=cnt+1;
    clear opt;
  end
end

%% Create Linear classifier from localiser trials
%%%%%% this first for loop iterates through each subject's FLI data and
%%%%%% saves an ensemble of linear classifiers as a structure called 'gnf'
mkdir([analysisdir,'classifiers/TrainedStim4Cell/',dirname]);
saveWrap = @(gnf, whenInMS, iSj) save([analysisdir,'classifiers/TrainedStim4Cell/',dirname,'/Lsj' num2str(iSj)], 'gnf', 'whenInMS');

for iSj=1:length(is.goodsub)
  iSjGood=is.goodsub(iSj);
  tic
  % creates gnf, which is a cell array of structs, with dimensions of 1 * 8 states * L1p
  gnf = cell(nStim, length(is.ENL1));
  lciInd = find(ismember(is.MEGruns{iSjGood},'lci')); % 'lci' denotes the functional localiser data
  data= nan(273,401,ntrlValid,length(lciInd)); % channel*time points(-3s - 1s)*trials*sessions
  trialindex=zeros(ntrlValid,length(lciInd));
  
  % Getting the data
  for iSes = 1:nSes
    
    % load the epoched data
    dir_temp=[is.OPTPath strrep(is.fnDate{iSjGood},'-','_') '/' is.fnMEG{iSjGood} '_' num2str(lciInd(iSes),'%02d') '.ds/highpass_' num2str(is.highpass) '.opt'];
    opt=load(fullfile(dir_temp, 'opt.mat'));
    opt=opt.opt;
    DLCI=spm_eeg_load(fullfile(dir_temp,[opt.results.spm_files_epoched_basenames{1,1}]));
    
    % load good channels:
    chan_MEG = indchantype(DLCI,'meeg');
    
    % load good trials
    allconds=[{'S1'},{'S2'},{'S3'},{'S4'},{'S5'},{'S6'},{'S7'},{'S8'}];
    
    all_trls = sort(DLCI.indtrial(allconds(:)));
    good_trls = sort(DLCI.indtrial(allconds(:),'good'));
    good_ind=ismember(all_trls,good_trls);
    
    trialindex(good_ind,iSes)=1;
    
    % get clean data:
    cldata = DLCI(chan_MEG,:,good_trls);
    data(:,:,good_ind,iSes)= cldata;
  end
  time=DLCI.time;
  % reshape;
  Cleandata=reshape(data,[size(data,1),size(data,2),size(data,3)*size(data,4)]); % combine three runs of data;
  Goodtrialindex=reshape(trialindex,[size(trialindex,1)*size(trialindex,2),1]);  % combine three runs of good trials;
  
  % Stimuli -> State mapping!!!
  % get the STIMULI information from the functional localiser
  trainedlist=load([is.rootBehav datestr(datenum(is.fnDate{iSjGood})) '/' is.fnSID{iSjGood} '_' is.fnTrainClass{iSjGood}{1}], 'permutedStateList'); % those are index of Exptrial, NOT ACTUAL STATES!!!
  upsdindex=load([is.rootBehav datestr(datenum(is.fnDate{iSjGood})) '/' is.fnSID{iSjGood} '_' is.fnTrainClass{iSjGood}{1}], 'upsd');
  lcitrials=load([is.rootBehav datestr(datenum(is.fnDate{iSjGood})) '/' is.fnSID{iSjGood} '_' is.fnTrainClass{iSjGood}{1}], 'lciTrials');
  correctindex=zeros(nSes*ntrlAll,1);
  RTindex=zeros(nSes*ntrlAll,1);
  ISI=zeros(nSes*ntrlAll,1);
  
  for n=1:nSes*ntrlAll
    correctindex(n,1)=lcitrials.lciTrials{n}.correct;
    RTindex(n,1)=lcitrials.lciTrials{n}.RT;
    ISI(n,1)=lcitrials.lciTrials{n}.thisISI;
  end
  
  % convert from STIMULUS labels to TASK SEQUENCE labels (stimuli to
  % sequence order was randomised across participants)
  Statemapping=nan(nSes*ntrlAll,1);
  S = load([is.rootBehav datestr(datenum(is.fnDate{iSjGood})) '/' is.fnSID{iSjGood} '_' is.fnTrainClass{iSjGood}{2}], 'data');
  
  for n=1:nSes*ntrlAll
    stimlabel=Exptrial(trainedlist.permutedStateList(n)).label;
    Statemapping(n)=S.data.Subdata.Exptrial(strmatch(stimlabel,{S.data.Subdata.Exptrial(:).label},'exact')).class;
  end
  
  labStim = Statemapping;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Subset data for further analysis - train classifer!
  subsetindex=ones(nSes*ntrlValid,1); % maximum number of actual lci trials
  % get rid of upsside down stimuli
  valididx = upsdindex.upsd==0;
  labStm_rightside =labStim(valididx);
  correctindex_rightside=correctindex(valididx);
  RTindex_rightside=RTindex(valididx);
  ISI_rightside=ISI(valididx);
  
  % get rid of incorrect trials in behavioral response OR bad trials
  % in the MEG scan OR long RT
  for itrial=1:nSes*ntrlValid
    if Goodtrialindex(itrial)==0 || RTindex_rightside(itrial)>1000 || RTindex_rightside(itrial)<250 %|| correctindex_rightside(itrial)==0
      subsetindex(itrial)=0;
    end
  end
  
  labStim_subset=labStm_rightside(logical(subsetindex));
  Cleandata_subset=Cleandata(:,:,logical(subsetindex));
  isis=ISI_rightside(logical(subsetindex));
  
  stimlabel=labStim_subset; % state labels for further analysis - it is now in STATE NOT STIMULI!
  
  % Subset channels that are useable for all sessions of interest
  sjchannel=squeeze(GoodChannel(iSj,:,:));
  goodchannindex= nansum(sjchannel')==length(sessionsOfInterest);
  lcidata = Cleandata_subset(logical(goodchannindex),:,:);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Now we train a classifier on the functional localiser data at 200 MS
  % post stimulus onset
  [nchan, ntim, ntrl] = size(lcidata);
  nullData=zeros(nchan,ntrl);
  trueData=zeros(nchan,ntrl);
  % Note we train classifiers to discriminate the different stimuli
  % AND a null condition defined as the data at the ISI timepoint
  % when fixation cross appears before stimulus
  for ii=1:ntrl
    nullData(:,ii) = squeeze(lcidata(:,nearest(time, 0) - ceil(isis(ii)/is.msPerSample),ii)); % 301 is the onset of lex stimuli in the -3s-1s trials
    trueData(:,ii) = squeeze(lcidata(:,nearest(time, TOI), ii)); % this is 200 ms after onset of the stimulus (onset stimulus @ sample 301)
  end
  
  % May need to turn it off for time lagged regression | MVE: different
  % scaling for trueData and nullData? --> Cam: ideally they would have
  % the same scaling, but it changed the initial results--> legacy
  % stuff
  nullData = scaleFunc(nullData');  % note, there are other options for how to scale the data
  trueData = scaleFunc(trueData');
  
  % Set up the classifer
  for iC=1:nStim   % train classifiers on null and iT sample of data.
    labels = [stimlabel == iC; zeros(size(nullData,1),1)];
    for iL1=1:length(is.ENL1)  % loop over L1 penalty values
      l1p = is.ENL1(iL1); l2p = 0; % alpha = L1/(2*L2+L1) ; lambda = 2*L2+L1
      [beta, fitInfo] = lassoglm([trueData; nullData], labels, 'binomial', ...
        'Alpha', l1p / (2*l2p+l1p), 'Lambda', 2*l2p + l1p, 'Standardize', false);
      gnf{iC,iL1}.beta = beta; gnf{iC,iL1}.Intercept = fitInfo.Intercept;
    end
  end
  saveWrap(gnf, whenInMS, iSj);  % anonymous function to call 'save' in parfor
  disp(['sj' num2str(iSj) ': No. of Good trial = ' num2str(sum(subsetindex))]);
  toc;
end

%% Apply classifiers to resting state data
%%%%%%%%% this second for loop iterates through subjects applying the
%%%%%%%%% ensemble of clasifiers saved above to the resting state data,
%%%%%%%%% saving the ensemble of outputted reactivation timecourses as a
%%%%%%%%% structure called Rreds:

iRun = 2; % select the 2nd resting session
mkdir([analysisdir,'classifiers/TestResting4Cell/',dirname]);
saveWrap = @(Rreds, is, iSj, iRun) save([analysisdir,'classifiers/TestResting4Cell/',dirname,'Rreds' num2str(iSj) '_' num2str(iRun)], 'Rreds', 'is', '-v7.3');

for iSj= 1:length(is.goodsub)
  iSjGood=is.goodsub(iSj);
  % load classifier
  gnf = load([analysisdir,'classifiers/TrainedStim4Cell/',dirname,'Lsj' num2str(iSj)], 'gnf') ; % get the gnf variable with the regression models
  gnf =gnf.gnf;
  
  % Load Resting State Data
  dstInd = find(ismember(is.MEGruns{iSjGood},'rst'));
  dir_temp=[is.OPTPath strrep(is.fnDate{iSjGood},'-','_') '/' is.fnMEG{iSjGood} '_' num2str(dstInd(iRun),'%02d') '.ds/highpass_' num2str(is.highpass) '.opt'];
  opt=load(fullfile(dir_temp, 'opt.mat'));
  opt=opt.opt;
  RST=spm_eeg_load(fullfile(dir_temp,[opt.results.spm_files_basenames{1,1}]));
  
  % Good timepoints/trials
  good_samples = ~all(badsamples(RST,':',':',':'));
  
  % Select MEG channel:
  chan_meg = indchantype(RST,'meeg');
  rRST = RST(chan_meg,:,:);
  
  sjchannel=squeeze(GoodChannel(iSj,:,:));
  goodchannindex= nansum(sjchannel')==length(sessionsOfInterest);
  rRST = rRST(logical(goodchannindex),:,:);
  
  % get the interval between the start and the end
  evtypes = RST.events; evtypes = {evtypes(:).type}';
  evvals = RST.events; evvals = {evvals(:).value}';
  evvals(all(cellfun(@ischar,evvals),2),:) = {0}; % replace string with 0
  evvals=cell2mat(evvals); %convert cell to double
  evtimes = RST.events; evtimes = [evtimes(:).time]';
  
  RestingTriggerVals = [87, 88]; %MVE: is mistake in code or comment??--> % 87 for the resting state BEFORE reward learning, 88 for the resting state AFTER reward learning
  RestStmInds = strcmp('frontpanel trigger', evtypes) & ismember(evvals, RestingTriggerVals(iRun));    % onset of Resting
  RestStmTimes = evtimes(RestStmInds);
  
  EndTriggerVals = [99];
  EndStmInds = strcmp('frontpanel trigger', evtypes) & ismember(evvals, EndTriggerVals);  % end of Resting
  EndStmTimes = evtimes(EndStmInds);
  
  wholetimeindex=zeros(size(RST,2),1);
  
  wholetimeindex(floor(RestStmTimes(end)*100):floor(RestStmTimes(end)*100)+30000-1)=1;
  
  % set the artifacts to zero
  badwithin_index= wholetimeindex==1 & double(good_samples')==0;
  rRST(:,badwithin_index)=0;
  data = rRST(:,logical(wholetimeindex));
  
  data = data'; % transform data format to nsamples*nsensors
  data = scaleFunc(data);
  
  % apply classifier to the clean data
  if isempty(RestStmTimes)
    nTr = 1;
  else
    nTr = length(RestStmTimes);
  end
  
  Rreds = cell(nTr,length(is.ENL1)); % permutation shuffle X classifier training times X ntrials X l1 params
  for iTr=1:nTr
    for iL1=1:length(is.ENL1)
      for iC=1:nStim
        % use regression models
        Rreds{iTr,iL1}(:,iC) = -(data * gnf{iC,iL1}.beta + gnf{iC,iL1}.Intercept);   % work on the X*beta space! MVE: the minus is because of a weird matlab definition. Now high positive values mean high probability
      end
    end
  end
  saveWrap(Rreds, is, iSj, iRun);
  disp(['sj' num2str(iSj) ' finished']);
end


%% Compute sequenceness on resting state activations
%%%%%%%%% This third for loop iterates over each subject's ensemble of
%%%%%%%%% reactivation timecourses, computing the Sequenceness measure for
%%%%%%%%% each timecourse
iRun = 2;
is.nShuf=factorial(4);
[~,~,uniquePerms] = uperms(1:4,is.nShuf);
% number of time lags * number of shuffles * number of subjects * 1 *
% number of lambda values
sfAll = nan(length(is.lgncy)+1, is.nShuf, length(is.goodsub), 1, length(is.ENL1));
sbAll = nan(length(is.lgncy)+1, is.nShuf, length(is.goodsub), 1, length(is.ENL1));
% transition matrix. backward is transpose of this
Tfwd = [0,1,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,1,0,0;0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0];

for iSj= 1:length(is.goodsub)
  
  sf = cell(is.nShuf, 1);
  sb = cell(is.nShuf, 1);
  
  S = load([analysisdir,'classifiers/TestResting4Cell/',dirname,'Rreds' num2str(iSj) '_' num2str(iRun)]); Rreds = S.Rreds;  % load this subject's preds
  
  nTr = size(Rreds,3);
  L1l = length(is.ENL1);
  maxLag = length(is.lgncy);
  tic
  for iTr=1:nTr
    tmp = Rreds(iTr,:);
    prSj = permute(shiftdim(cat(3,tmp{:}),-1),[2,3,1,4]); % this is samples*states*alpha*lambda
    prSj = prSj(:,:, :); % alpha and trainingTimes put into a single dimension for the vectorization of sequenceness
    
    % Sequence
    nP = size(prSj,3);
    sf_temp=nan(maxLag,is.nShuf,nP); % time lag (forward sequence) x permuation x (lambdas*trainingtimes)
    sb_temp=nan(maxLag,is.nShuf,nP); % time lag (backward sequence) x permuation x (lambdas*trainingtimes)
    scf_temp=nan(maxLag,is.nShuf,nP);
    scb_temp=nan(maxLag,is.nShuf,nP);
    
    for vec=1:nP
      
      X=squeeze(prSj(:,:,vec));
      
      if ~any(~isnan(X(:)))
        continue
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Temporally shift the data
      nbins=maxLag+1;
      warning off
      dm=[];
      for kk=1:nStim
        temp=toeplitz(X(:,kk), zeros(nbins,1));
        temp=temp(:,2:end);
        dm=[dm temp]; % this is now samples x (timelags*nStim)
      end
      
      warning on
      
      Y=X;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % constructing alpha cycle regressors
      betas = nan(nStim*maxLag, nStim);
      bins=10; % because 10 Hz = 10 samples @ 100 Hz
      % CAM'S GONNA FIGURE OUT WHAT IS HAPPENING HERE! (only regressing
      % at each alpha cycle lag??)
      for ilag=1:bins%maxLag
        temp_zinds = (1:bins:nStim*maxLag) + ilag - 1; % indices of all alpha cycles away from the lag of interest
        temp = pinv([dm(:,temp_zinds) ones(length(dm(:,temp_zinds)),1)])*Y; % 49 alpha cycles are regressed out
        betas(temp_zinds,:)=temp(1:end-1,:);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Regressing the templates (forward and backward) sequences on the
      % time-lagged regression weights organised into a transition matrix
      % (transition matrix per time lag)
      X=[];
      Y=[];
      betasr=reshape(betas,[maxLag nStim nStim]);
      betasnbins64=reshape(betas,[maxLag nStim^2]);
      
      for iShuf=1:is.nShuf
        if iShuf==1  % the real transition matrix
          T1 = Tfwd; T2 = Tfwd';
        else  % construct random transition matrix
          X = Tfwd(1:4,1:4); % Sequence 1: 1->2->3->4
          Y = Tfwd(5:8,5:8); % Sequence 2: 5->6->7->8
          
          rp = uniquePerms(iShuf,:);  % use the 30 unique permutations (is.nShuf should be set to 29)
          Tfwd_temp=zeros(8,8);
          Tfwd_temp(1:4,5:8)= X(rp,rp);
          Tfwd_temp(5:8,1:4)=Y(rp,rp);
          T1 = Tfwd_temp;%Tfwd(rp,rp);
          T2 = T1'; % backwards transitions
        end
        
        bbb=pinv([T1(:) T2(:) squash(eye(nStim)) squash(ones(nStim))])*(betasnbins64'); % MVE: what's happening here?? Why do we multiply the alpha cycle's betas with the other betas??
        sf_temp(:,iShuf,vec)=bbb(1,:);
        sb_temp(:,iShuf,vec)=bbb(2,:);
      end
    end
    
    for iShuf=1:is.nShuf
      sf{iShuf,iTr} = nan(1, 1, length(is.lgncy)+1, L1l);
      sb{iShuf,iTr} = nan(1, 1, length(is.lgncy)+1, L1l);
      
      sf{iShuf,iTr}(1, 1, 2:end, :, :)=squeeze(sf_temp(:,iShuf,:));
      sb{iShuf,iTr}(1, 1, 2:end, :, :)=squeeze(sb_temp(:,iShuf,:));
    end
  end
  
  % RESHAPE
  sf2 = permute(cell2mat(sf), [3 1 2 4]); % sf2 = latency(61) * shuffles(24) * trials (1) * alpha (L1 regulation = 10)
  sb2 = permute(cell2mat(sb), [3 1 2 4]);
  
  sfAll(:,:,iSj,:,:) = sf2;
  sbAll(:,:,iSj,:,:) = sb2;
  
  disp(['Sub' num2str(iSj) ' Sequence Finished' ])
  toc
end

cd([analysisdir,'classifiers/Sequence_by_Training_4Cell']);
mkdir(dirname);
save([dirname,'StimAll'],'sfAll','sbAll','is','-v7.3');

%% Plot replay effect broken by training time
close all
dir1name='L1regression/';

load([analysisdir,'classifiers/Sequence_by_Training_4Cell/',dir1name,'StimAll.mat']);
nsub=1:length(is.goodsub);
sfAll=squeeze(sfAll(2:end,:,nsub,1,:)); % timelag*nsubs*trainingtime*L1
sbAll=squeeze(sbAll(2:end,:,nsub,1,:));

% Plot sequenceiness over time lag for all training time (100ms: 300ms) -
% without standard error
figure,
x0=40;
y0=40;
width=700;
height=400;
for iL1=1:10
  subplot(2,5,iL1)
  sfb=squeeze(sfAll(:,:,:,iL1)-sbAll(:,:,:,iL1));
  sfb_perms = squeeze(sfb(:,2:end,:));
  sfb=squeeze(sfb(:,1,:));
  Mreplay=squeeze(nanmean(sfb,2));
  Sreplay=squeeze(nanstd(sfb,[],2))/sqrt(size(sfb,2));
  tTimes2=10:10:600;
  
  npThreshAll = max(max(abs(nanmean(sfb_perms,3))));hold on;
  
  set(gcf,'units','points','position',[x0,y0,width,height])
  plot(tTimes2, Mreplay, 'LineWidth', 1.5)
  xlabel('lag (ms)')
  ylabel('L1-reverse/Leftarrow sequenceness /RightarrowL1-forward'),
  grid on;
  title(['Sequenceiness for lambda = ',num2str(iL1)]);
  
  plot([tTimes2(1) tTimes2(end)], -npThreshAll*[1 1], 'k--'), plot([tTimes2(1) tTimes2(end)], npThreshAll*[1 1], 'k--')
end

%% and cross validate to choose a single parameter value:
% MVE: Note that it doens't do cross validation to select the best and same
% parameter value for the whole group. It uses the value that is best in
% the group without the subject of interest.
shuf1=1;
figure,
x0=40;
y0=40;
width=700;
height=400;
for iSjLO = 1:length(nsub)
  iSjLI = setdiff(1:length(nsub),iSjLO);
  sfb=squeeze(nanmean(sfAll(:,shuf1,iSjLI,:)-sbAll(:,shuf1,iSjLI,:),3));
  %CV over all time:
  [~,iLCV(iSjLO)] = max(max(sfb,[],1)); % the best l1 over all other subjects
  sfbCV(iSjLO,:,:) = sfAll(:,:,iSjLO,iLCV(iSjLO))-sbAll(:,:,iSjLO,iLCV(iSjLO));
end

Mreplay=squeeze(nanmean(sfbCV(:,:,1),1));
Sreplay=squeeze(nanstd(sfbCV(:,:,1),[],1))/sqrt(size(sfbCV(:,:,1),12));
tTimes2=10:10:600;

set(gcf,'units','points','position',[x0,y0,width,height])
plot(tTimes2, Mreplay, 'LineWidth', 1.5)
xlabel('lag (ms)')
ylabel('L1-reverse/Leftarrow sequenceness /RightarrowL1-forward'),
title({'Sequenceness, Study II Lasso,',[' criteria maximised across all timepoints']});
npThreshAll = max(max(abs(nanmean(sfbCV(:,:,2:end),1))));hold on;
plot([tTimes2(1) tTimes2(end)], -npThreshAll*[1 1], 'k--'), plot([tTimes2(1) tTimes2(end)], npThreshAll*[1 1], 'k--')
grid on;

%% Save an inferred Replay Timecourse:

%define peak sequenceness time:
[~,lag] = max(abs(Mreplay));
for iSj=1:length(is.goodsub)
  % take that subject's CV parameter selection:
  iL = iLCV(iSj);
  
  % select second resting state session:
  iRun = 2;
  
  % load reactivation timecourse:
  S = load([analysisdir,'classifiers/TestResting4Cell/',dirname,'Rreds' num2str(iSj) '_' num2str(iRun)]); Rreds = S.Rreds;  % load this subject's preds
  
  reactivationtimecourse = S.Rreds{1,iL};
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
replaydir = [analysisdir,'classifiers/Sequence_by_Training_4Cell/',dirname,'STUDYII_ReplayOnset/'];
mkdir(replaydir)
if iRun==1
  fname = 'STUDYII_PreplayOnset.mat';
else
  fname = 'STUDYII_ReplayOnset.mat';
end
save([replaydir,fname],'ToRall');

%% analyse pooled saved values and plot as image:
figure('Position', [1 41 1920 1084]);
dirname='L1regression/';

plotlocs=[1,2,4,5,3,6];

cd(analysisdir);
load([analysisdir,'classifiers/Sequence_by_Training_4Cell/',dirname,'/StimAll.mat']);
nsub=1:length(is.goodsub);
sfAll=squeeze(sfAll(2:end,:,nsub,1,:)); % timelag*nsubs*1*L1
sbAll=squeeze(sbAll(2:end,:,nsub,1,:)); % timelag*nsubs*1*L1

t_points=1:60; % which points to compute sequenceness over - should be all, 1:60
for iL1=1:10
  sfb=squeeze(sfAll(:,:,:,iL1)-sbAll(:,:,:,iL1));
  sfb_perms = squeeze(sfb(:,2:end,:));
  sfb=squeeze(sfb(t_points,1,:));
  Mreplay=squeeze(nanmean(sfb,2));
  Sreplay=squeeze(nanstd(sfb,[],2))/sqrt(size(sfb,2));
  tTimes2=10:10:600;
  npThreshAll = max(max(abs(nanmean(sfb_perms(t_points,:,:),3))));hold on;
  seqmat(iL1,:) = Mreplay ./ npThreshAll;
end
max_replay_points = squeeze(max(seqmat,[],2));
min_replay_points = squeeze(min(seqmat,[],2));
replay_points = max(abs(cat(2,max_replay_points,min_replay_points)),[],2) .* ...
  [[max_replay_points>abs(min_replay_points)] + ...
  [max_replay_points<abs(min_replay_points)]*-1];
subplot(2,3,plotlocs);
imagesc(replay_points);caxis([-2,2]);hold on;
hold on;
pos_sig_points = replay_points>=1;
neg_sig_points = replay_points<=-1;
contour(pos_sig_points,1,'k--');
contour(neg_sig_points,1,'r--');
set(gca,'YTick',1:10);
ylabel('/lambda regularisation penalty');

set(gca,'YTickLabel',0.001:0.001:0.01);
title('L1 regression on single timepoints');

set(gca,'XTick',2:2:20);
set(gca,'XTickLabel',120:20:400);
xlabel('Time of training (msec)');

