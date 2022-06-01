function Rreds = apply_toRest(iSj, is, gnf, GoodChannel, iRun, iPerm)
%% Apply classifiers to resting state data
if ~isfield(is, 'doPermutation')
  is.doPermutation=false;
end
if ~isfield(is, 'doSave')
  is.doSave=true;
end
if iRun==2 && ~is.doPermutation
  usepath = [is.AnalysisPath, 'classifiers/TestResting4Cell/']; % original results from Yunzhe
else
  usepath = [is.matsdir, 'classifiers/TestResting/'];
end
if is.usePrecomputed
  if ~is.doPermutation || iPerm==1
%     fprintf('Loading results - classifier applied to resting state data \n')
  end
  % load classifier
  try
    if is.doPermutation==false
      fname = [usepath, 'Rreds' num2str(iSj),'_', sprintf('%d',iRun)]; 
    else
      fname = [usepath, 'perm/', 'Rreds' num2str(iSj),'_', sprintf('%d',iRun), '_', sprintf('perm%d', iPerm)]; % get the gnf variable with the regression models
    end
    Rreds = load(fname, 'Rreds');% get the gnf variable with the regression models
    Rreds =Rreds.Rreds;
  catch
    if ~is.doPermutation || iPerm==1
      fprintf('Loading in data was not possible - recomputing now \n')
    end
    is.usePrecomputed=false;
    Rreds = apply_toRest(iSj, is, gnf, GoodChannel, iRun, iPerm);
  end
else
  %%%%%%%%% this second for loop iterates through subjects applying the
  %%%%%%%%% ensemble of clasifiers saved above to the resting state data,
  %%%%%%%%% saving the ensemble of outputted reactivation timecourses as a
  %%%%%%%%% structure called Rreds:
  % iRun specifies which resting state session (1=pre, 2=post)

  folder=usepath;
  if is.doPermutation
    folder = [folder, 'perm/'];
  end
  if ~exist(folder, 'dir')
    mkdir(folder);
  end
  fname =  [folder,'Rreds' num2str(iSj) '_' num2str(iRun)];
  if is.doPermutation
    fname = [fname, '_', sprintf('perm%d', iPerm)];
  end
  if is.doSave
    saveWrap = @(Rreds, is, iSj, iRun) save(fname, 'Rreds', 'is', '-v7.3');
  end
  iSjGood=is.goodsub(iSj);

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

  % TODO: NOTE: the goodchannindex is based on all sessions as below. But I
  % don't have the data for the other sessions, so the following will only
  % work if also the goodchannels are redefined.
  % before:
%   sjchannel=squeeze(GoodChannel(iSj,:,:));
  sjchannel=squeeze(GoodChannel);
  goodchannindex= nansum(sjchannel')==length(is.MEGruns{1,iSjGood});
  % now:
%   goodchannindex= nansum(GoodChannel(:,is.sessionsOfInterest)')==length(is.sessionsOfInterest);
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

  Rreds = cell(nTr,37,1,length(is.ENL1)); % permutation shuffle X classifier training times X ntrials X l1 params
  % size is legacy stuff
  for iTr=1:nTr
    for iL1=1:length(is.ENL1)
      for iC=1:is.nStim
        % use regression models
        Rreds{iTr,37,1,iL1}(:,iC) = -(data * gnf{iC,iL1}.beta + gnf{iC,iL1}.Intercept);   % work on the X*beta space! MVE: the minus is because of a weird matlab definition. Now high positive values mean high probability
      end
    end
  end
  if is.doPermutation
    Rreds(1,37,1,setdiff(1:10, is.lambda))=cell(1,9);
    Rreds{1,37,1,is.lambda} = single(Rreds{1,37,1,is.lambda});
  end
  if is.doSave
    saveWrap(Rreds, is, iSj, iRun);
  end
  if ~is.doPermutation || iPerm==is.nPerm
    disp(['sj' num2str(iSj) ' finished']);
  end
end