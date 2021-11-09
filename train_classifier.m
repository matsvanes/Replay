function gnf = train_classifier(iSj, is, Exptrial, GoodChannel)
%% Create Linear classifier from localiser trials
if is.usePrecomputed
  fprintf('Loading precomputed classifier \n')
  try
    % load classifier
    gnf = load([is.AnalysisPath,'classifiers/TrainedStim4Cell/','Lsj' num2str(iSj)], 'gnf') ; % get the gnf variable with the regression models
    gnf = squeeze(gnf.gnf(1,find(is.whenInMS==is.TOI)+1,:,:));
  catch
    fprintf('Loading in data was not possible - recomputing now \n')
    is.usePrecomputed=false;
    gnf = train_classifier(iSj, is, Exptrial, GoodChannel);
  end
else
  %%%%%% this first for loop iterates through each subject's FLI data and
  %%%%%% saves an ensemble of linear classifiers as a structure called 'gnf'
  if ~isfolder([is.AnalysisPath,'classifiers/TrainedStim4Cell/L1regression']), mkdir([is.AnalysisPath,'classifiers/TrainedStim4Cell/L1regression']), end;
  saveWrap = @(gnf, whenInMS, iSj) save([is.AnalysisPath,'classifiers/TrainedStim4Cell/','/Lsj' num2str(iSj)], 'gnf', 'whenInMS');

  iSjGood=is.goodsub(iSj);
  whenInMS=is.whenInMS;
  tic
  % creates gnf, which is a cell array of structs, with dimensions of 1 * 8 states * L1p
  gnf = cell(is.nStim, length(is.ENL1));
  lciInd = find(ismember(is.MEGruns{iSjGood},'lci')); % 'lci' denotes the functional localiser data
  data= nan(is.nChan,is.nTime,is.ntrlValid,length(lciInd)); % channel*time points(-3s - 1s)*trials*sessions
  trialindex=zeros(is.ntrlValid,length(lciInd));

  % Getting the data
  for iSes = 1:is.nSes

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
  correctindex=zeros(is.nSes*is.ntrlAll,1);
  RTindex=zeros(is.nSes*is.ntrlAll,1);
  ISI=zeros(is.nSes*is.ntrlAll,1);

  for n=1:is.nSes*is.ntrlAll
    correctindex(n,1)=lcitrials.lciTrials{n}.correct;
    RTindex(n,1)=lcitrials.lciTrials{n}.RT;
    ISI(n,1)=lcitrials.lciTrials{n}.thisISI;
  end

  % convert from STIMULUS labels to TASK SEQUENCE labels (stimuli to
  % sequence order was randomised across participants)
  Statemapping=nan(is.nSes*is.ntrlAll,1);
  S = load([is.rootBehav datestr(datenum(is.fnDate{iSjGood})) '/' is.fnSID{iSjGood} '_' is.fnTrainClass{iSjGood}{2}], 'data');

  for n=1:is.nSes*is.ntrlAll
    stimlabel=Exptrial(trainedlist.permutedStateList(n)).label;
    Statemapping(n)=S.data.Subdata.Exptrial(strmatch(stimlabel,{S.data.Subdata.Exptrial(:).label},'exact')).class;
  end

  labStim = Statemapping;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Subset data for further analysis - train classifer!
  subsetindex=ones(is.nSes*is.ntrlValid,1); % maximum number of actual lci trials
  % get rid of upsside down stimuli
  valididx = upsdindex.upsd==0;
  labStm_rightside =labStim(valididx);
  correctindex_rightside=correctindex(valididx);
  RTindex_rightside=RTindex(valididx);
  ISI_rightside=ISI(valididx);

  % get rid of incorrect trials in behavioral response OR bad trials
  % in the MEG scan OR long RT
  for itrial=1:is.nSes*is.ntrlValid
    if Goodtrialindex(itrial)==0 || RTindex_rightside(itrial)>1000 || RTindex_rightside(itrial)<250 %|| correctindex_rightside(itrial)==0
      subsetindex(itrial)=0;
    end
  end

  labStim_subset=labStm_rightside(logical(subsetindex));
  Cleandata_subset=Cleandata(:,:,logical(subsetindex));
  isis=ISI_rightside(logical(subsetindex));

  stimlabel=labStim_subset; % state labels for further analysis - it is now in STATE NOT STIMULI!

  % Subset channels that are useable for all sessions of interest
  goodchannindex= nansum(GoodChannel')==length(is.sessionsOfInterest);
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
    trueData(:,ii) = squeeze(lcidata(:,nearest(time, is.TOI), ii)); % this is 200 ms after onset of the stimulus (onset stimulus @ sample 301)
  end

  % May need to turn it off for time lagged regression | MVE: different
  % scaling for trueData and nullData? --> Cam: ideally they would have
  % the same scaling, but it changed the initial results--> legacy
  % stuff
  nullData = scaleFunc(nullData');  % note, there are other options for how to scale the data
  trueData = scaleFunc(trueData');

  % Set up the classifer
  for iC=1:is.nStim   % train classifiers on null and iT sample of data.
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