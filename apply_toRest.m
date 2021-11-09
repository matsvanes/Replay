function Rreds = apply_toRest(iSj, is, gnf, GoodChannel, iRun)
%% Apply classifiers to resting state data
if is.usePrecomputed
  fprintf('Loading results - classifier applied to resting state data \n')
  % load classifier
  try
    Rreds = load([is.AnalysisPath,'classifiers/TestResting4Cell/', 'Rreds' num2str(iSj),'_', sprintf('%d',iRun)], 'Rreds') ; % get the gnf variable with the regression models
    Rreds =Rreds.Rreds;
  catch
    fprintf('Loading in data was not possible - recomputing now \n')
    is.usePrecomputed=false;
    Rreds = apply_toRest(iSj, is, gnf, GoodChannel, iRun);
  end
else
  %%%%%%%%% this second for loop iterates through subjects applying the
  %%%%%%%%% ensemble of clasifiers saved above to the resting state data,
  %%%%%%%%% saving the ensemble of outputted reactivation timecourses as a
  %%%%%%%%% structure called Rreds:
  % iRun specifies which resting state session (1=pre, 2=post)

  mkdir([is.AnalysisPath,'classifiers/TestResting4Cell/',is.dirname]);
  saveWrap = @(Rreds, is, iSj, iRun) save([is.AnalysisPath,'classifiers/TestResting4Cell/','Rreds' num2str(iSj) '_' num2str(iRun)], 'Rreds', 'is', '-v7.3');
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
  saveWrap(Rreds, is, iSj, iRun);
  Rreds = squeeze(Rreds(:,37,1,:));
  disp(['sj' num2str(iSj) ' finished']);
end