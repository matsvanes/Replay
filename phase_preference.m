% This script should load in pre-trained decoding models (trained on
% localiser trials) of the Replay dataset (studyII). It should then apply this to
% resting state data. The probabilities in the post-task resting state data
% should be higher than on the pre-task resting state. The working
% hypothesis is that the phase concentration of the highest percentile
% probabilities is also higher for the post-task compared to the pre-task
% resting state data
% 2021 - Mats van Es

% define some study specific variables.
lambda = 0.006; % used by https://doi.org/10.1016/j.cell.2019.06.012 (VERIFY!)
MFStartup_studyII;
lambda_idx = find(is.ENL1==lambda);
analysisdir = '/Volumes/T5_OHBA/analysis/replay/';
dirname='L1regression/';
whenInMS = is.tss*is.msPerSample;
TOI = 0.200; % 200ms post stimulus.
sessionsOfInterest = [1,2,3,4,8]; % rst, lci, lci, lci, rst
pretask=1;
posttask=2;
nStim = 8;
ntrlAll=120; % per session (includes upside down images)
ntrlValid=96; % per session
nSes = 3;
perc=1; % highest percentile probabilities to use


mkdir([analysisdir,'phase']);
GoodChannel=nan(length(is.goodsub),273,numel(sessionsOfInterest)); %nsubject*nsensors*nruns
for iSj = 1%:length(is.goodsub)
  figure;
  iSjGood=is.goodsub(iSj);
  
  % find good channels in all sessions of interest
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
  
  % load in decoding model
  gnf = load([analysisdir,'classifiers/TrainedStim4Cell/',dirname,'Lsj' num2str(iSj)], 'gnf') ; % get the gnf variable with the regression models
  % select the one corresponding to the correct lambda
  gnf = gnf.gnf(:,lambda_idx);
  
  % load resting state data
  for iRun = [pretask, posttask]
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
    data(data==0)=nan; % remove zeros
    %     % apply classifier to the clean data
    %     if isempty(RestStmTimes)
    %       nTr = 1;
    %     else
    %       nTr = length(RestStmTimes);
    %     end
    
    for iC=1:nStim
      % use regression models
      Rreds{iRun}(:,iC) = -(data * gnf{iC}.beta + gnf{iC}.Intercept);   % work on the X*beta space! MVE: the minus is because of a weird matlab definition. Now high positive values mean high probability
    end
    prob{iRun} = 1./(1+exp(Rreds{iRun}));
    
    % find highest percentile probabilities
    x=prob{iRun}; x(isnan(x))=-Inf;
    [val, idx] = sort(prob{iRun}, 'descend', 'MissingPlacement', 'last');
    n = round((size(val,1)-sum(isnan(data(:,1))))*perc/100);
    idx_top = idx(1:n,:);
    val_top{iRun} = val(1:n,:);
    
    % find lowest percentile probabilities
    [val, idx] = sort(prob{iRun}, 'ascend', 'MissingPlacement', 'last');
    n = round((size(val,1)-sum(isnan(data(:,1))))*perc/100);
    idx_bottom = idx(1:n,:);
    val_bottom{iRun} = val(1:n,:);
    
    % simulate alpha signal --> TODO: replace this by using the alpha
    % signal of a specific region
    time = 0:1/RST.fsample:(size(data,1)-1)/RST.fsample;
    f=10;
    clear i
    a = cos(2*pi*f*time) + i*sin(2*pi*f*time); % raw alpha signal
    phase = angle(a);
    phase = repmat(phase, [8 1])';
    for k=1:8
      % take the phases of the highest/lowest probabilities
      phase_top{iRun}(:,k) = phase(idx_top(:,k));
      phase_bottom{iRun}(:,k) = phase(idx_bottom(:,k));
      % probably we can use the parametric test, but this assumes max 1
      % mode, and a "gaussian" distribution around it.
      % alternatives are the nonparametric Hodges-Ajne test (circ_otest),
      % but this doesn't return a test statistic. Other options is to
      % compare the distribution with a uniform distribution using the
      % Kuiper test (circ_kuipertest) or Winson's U2 test.
      [pval(iRun,k), z(iRun,k)] = circ_rtest(phase_top{iRun}(:,k)); 
%       [pval(iRun,k), m(iRun,k)] = circ_otest(phase_top{iRun}(:,k));
      subplot(2,8,(iRun-1)*8+k);
      rose(phase_top{iRun}(:,k)); hold on
      rose(phase_bottom{iRun}(:,k));
      if (iRun-1)*8+k == 4 
        title('pre-task Resting State')
      elseif (iRun-1)*8+k == 12
        title('post-task Resting State')
      end
    end
    for k=1:8
      [~,fval(iRun,k)] = circ_ktest(phase_top{iRun}(:,k), phase_bottom{iRun}(:,k));
    end
  end
  % optionally - use Watson's U2 statistic to compare the two circular
  % distributions. See https://uk.mathworks.com/matlabcentral/fileexchange/43543-pierremegevand-watsons_u2
  % or possibly use circ_ktest for equal concentration parameter?
end


