%% Train decoders - stim code
%clear all;
whichstudy = 2;
if whichstudy==1
    data_dir = '/Volumes/T5_OHBA/data/replay/Replaydata4Cam/';
    MFStartup;
else
    data_dir = '/Volumes/T5_OHBA/data/replay/StrLearn_MEGexp/';
    MFStartup_studyII;
end
analysisdir = [data_dir,'Analysis\'];
if ~isdir(analysisdir)
    mkdir(analysisdir)
end
%addpath(genpath('D:\Documents\spm12'));  
whenInMS = is.tss*is.msPerSample;
is.whichTimes=28:47; % 10th is 37 -> 200 ms
GoodChannel=nan(length(is.goodsub),273,10); %nsubject*nsensors*nruns
if whichstudy==1
    
else
    load('F:\My Data\Yunzhe data\StrLearn_MEGexp\Analysis\Exptrial.mat')
end


% Good Channels
for iSj=1:length(is.goodsub)
    %iSjGood=is.goodsub(iSj);
    iSjGood = iSj;
    for isession=1:length(is.MEGruns{1,iSjGood})
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
        GoodChannel(iSj,:,isession)=index;
        clear opt;
    end
end

% this loop runs through the different methods options I tested in Chapter
% 4 of my thesis. To run the method Yunzhe used in his Cell paper, select
% iCond=1 and ignore options 2 to 7

for iCond=1
    % note for iCond=5 and 6 need to manually amend from parfor loops below to for, due to
    % how HMM code calls sparsity params
if iCond==1
    L1regress=true;
    HMMclassifier=false;
    windowing=false; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
    dirname='L1regression\';
elseif iCond==2
    L1regress=false;
    HMMclassifier=false;
    windowing=false; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
    dirname='L2regression\';
elseif iCond==3
    L1regress=true;
    HMMclassifier=false;
    windowing=true; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
    dirname='L1regression_window3\';
elseif iCond==4
    L1regress=false;
    HMMclassifier=false;
    windowing=true; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
    dirname='L2regression_window3\';
elseif iCond==5
    L1regress=false;
    HMMclassifier=true;
    windowing=false; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
    dirname='ARDsingletimepoint\';
elseif iCond==6
    L1regress=false;
    HMMclassifier=true;
    windowing=true; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
    dirname='ARDregression_window3\';
elseif iCond==7
    L1regress=false;
    HMMclassifier=true;
    windowing=false; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
    dirname='ARDsingletimepoint_bigprior\';
end


%%%%%% this first for loop iterates through each subject's FLI data and
%%%%%% saves an ensemble of linear classifiers as a structure called 'gnf'
mkdir([analysisdir,'classifiers\TrainedStim4Cell\',dirname]);
saveWrap = @(gnf, whenInMS, iSj) save([analysisdir,'classifiers\TrainedStim4Cell\',dirname,'\Lsj' num2str(iSj)], 'gnf', 'whenInMS');

for iSj=1:length(is.goodsub)
    iSjGood=is.goodsub(iSj);
    tic
    % creates gnf, which is a cell array of structs, with dimensions of 1 * trainTime * 8 states * L1p
    gnf = cell(1, length(is.tss), 8, length(is.ENL1));    
    lciInd = find(ismember(is.MEGruns{iSjGood},'lci')); % 'lci' denotes the functional localiser data
    data= nan(273,401,96,length(lciInd)); % channel*time points(-3s - 1s)*trials*sessions
    trialindex=zeros(96,length(lciInd));

    % Getting the data
    for ii = 1:length(lciInd)

        % load the epoched data
        dir_temp=[is.OPTPath strrep(is.fnDate{iSjGood},'-','_') '\' is.fnMEG{iSjGood} '_' num2str(lciInd(ii),'%02d') '.ds\highpass_' num2str(is.highpass) '.opt'];
        opt=load(fullfile(dir_temp, 'opt.mat'));
        opt=opt.opt;
        DLCI=spm_eeg_load(fullfile(dir_temp,[opt.results.spm_files_epoched_basenames{1,1}]));

        % load good channels:
        chan_MEG = indchantype(DLCI,'meeg');

        % load good trials
        allconds=DLCI.condlist;
        allconds=[{'S1'},{'S2'},{'S3'},{'S4'},{'S5'},{'S6'},{'S7'},{'S8'}];

        all_trls = sort([DLCI.indtrial(allconds(:))]);
        good_trls = sort([DLCI.indtrial(allconds(:),'good')]);
        good_ind=ismember(all_trls,good_trls);

        trialindex(good_ind,ii)=1;

        % get clean data:
        cldata = DLCI(chan_MEG,:,good_trls);
        data(:,:,good_ind,ii)= cldata;
    end
    % reshape;
    Cleandata=reshape(data,[size(data,1),size(data,2),size(data,3)*size(data,4)]); % combine three runs of data;
    Goodtrialindex=reshape(trialindex,[size(trialindex,1)*size(trialindex,2),1]);  % combine three runs of good trials;

    % Stimuli -> State mapping!!!        
    % get the STIMULI information from the functional localiser
    trainedlist=load([is.rootBehav datestr(datenum(is.fnDate{iSjGood})) '/' is.fnSID{iSjGood} '_' is.fnTrainClass{iSjGood}{1}], 'permutedStateList'); % those are index of Exptrial, NOT ACTUAL STATES!!!
    upsdindex=load([is.rootBehav datestr(datenum(is.fnDate{iSjGood})) '/' is.fnSID{iSjGood} '_' is.fnTrainClass{iSjGood}{1}], 'upsd');
    lcitrials=load([is.rootBehav datestr(datenum(is.fnDate{iSjGood})) '/' is.fnSID{iSjGood} '_' is.fnTrainClass{iSjGood}{1}], 'lciTrials');
    correctindex=[];
    RTindex=[];
    labStim=[];
    ISI=[];

    for n=1:360
        correctindex(n,1)=lcitrials.lciTrials{n}.correct;
        RTindex(n,1)=lcitrials.lciTrials{n}.RT;
        ISI(n,1)=lcitrials.lciTrials{n}.thisISI;
    end

    % convert from STIMULUS labels to TASK SEQUENCE labels (stimuli to
    % sequence order was randomised across participants)
    Statemapping=nan(360,1);
    S = load([is.rootBehav datestr(datenum(is.fnDate{iSjGood})) '/' is.fnSID{iSjGood} '_' is.fnTrainClass{iSjGood}{2}], 'data');

    for n=1:360
        stimlabel=Exptrial(trainedlist.permutedStateList(n)).label;
        Statemapping(n)=S.data.Subdata.Exptrial(strmatch(stimlabel,{S.data.Subdata.Exptrial(:).label},'exact')).class;
    end

    labStim = Statemapping;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subset data for further analysis - train classifer!
    subsetindex=ones(288,1);
    % get rid of upsside down stimuli
    labStm_rightside =labStim(upsdindex.upsd==0); 
    correctindex_rightside=correctindex(upsdindex.upsd==0);
    RTindex_rightside=RTindex(upsdindex.upsd==0);
    ISI_rightside=ISI(upsdindex.upsd==0);

    % get rid of incorrect trials in behavioral response OR bad trials
    % in the MEG scan OR long RT

    for itrial=1:288
        if Goodtrialindex(itrial)==0 || RTindex_rightside(itrial)>1000 || RTindex_rightside(itrial)<250 %|| correctindex_rightside(itrial)==0
            subsetindex(itrial)=0;
        end
    end

    labStim_subset=labStm_rightside(logical(subsetindex));
    Cleandata_subset=Cleandata(:,:,logical(subsetindex));
    isis=ISI_rightside(logical(subsetindex));

    stimlabel=labStim_subset; % state labels for further analysis - it is now in STATE NOT STIMULI!

    % if use eyeblink, get rid of channels influenced by eyeblink, mostly in
    % the frontal area, Otherwise, subset channels that are useable for both two runs
    if is.eyeBlinkSens==1
        %commonmisschannel=any(any(isnan(Cleandata_subset),3),2);
        lcidata = Cleandata_subset(~eyeBlinkSens,:,:);
    else
        sjchannel=squeeze(GoodChannel(iSj,:,:));
        goodchannindex= nansum(sjchannel')==length(is.MEGruns{1,iSjGood});
        %lcidata(any(any(isnan(lcidata),3),2),:,:) = [];
        lcidata = Cleandata_subset(logical(goodchannindex),:,:); 
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now we train a classifier on the functional localiser data by
    % iterating over timepoints listed in is.whichtimes:
    for iTT = 1:length(is.whichTimes) % which specific time to be trained on 
        nullData=[];
        trueData=[];
        % Note we train classifiers to discriminate the different stimuli
        % AND a null condition defined as the data at the ISI timepoint
        % when fixation cross appears before stimulus
        for ii=1:size(lcidata,3)
            nullData(:,ii) = squeeze(lcidata(:,301-ceil(isis(ii)/is.msPerSample),ii)); % 301 is the onset of lex stimuli in the -3s-1s trials 
            trueData(:,ii) = squeeze(lcidata(:,301+is.tss(is.whichTimes(iTT)),ii));
            if windowing
                %add samples from square window one point to either side:
                t_null = 301-ceil(isis(ii)/is.msPerSample)-1;if t_null==0;t_null=1;end
                nullData(:,ii+size(lcidata,3)) = squeeze(lcidata(:,t_null,ii)); % 301 is the onset of lex stimuli in the -3s-1s trials 
                trueData(:,ii+size(lcidata,3)) = squeeze(lcidata(:,301+is.tss(is.whichTimes(iTT))-1,ii));
                nullData(:,ii+2*size(lcidata,3)) = squeeze(lcidata(:,301-ceil(isis(ii)/is.msPerSample)+1,ii)); % 301 is the onset of lex stimuli in the -3s-1s trials 
                trueData(:,ii+2*size(lcidata,3)) = squeeze(lcidata(:,301+is.tss(is.whichTimes(iTT))+1,ii));
            end
        end

        %nullData' = nullData' - nullData'*sDst;  % spatial highpass
        %trueData' = trueData' - trueData'*sDst;  % spatial highpass        

        % May need to turn it off for time lagged regression
        nullData = scaleFunc(nullData');  % note, there are other options for how to scale the data  
        trueData = scaleFunc(trueData');    

        % Set up the classifer
        for iShuf=1 % ignore this - previously applied label permutations, now defunct

            for iC=1:8   % train classifiers on null and iT sample of data.
                %disp(['sj' num2str(iSj) ', iShuf=' num2str(iShuf) ', iT=' num2str(is.whichTimes(iTT)) ', iC=' num2str(iC)])
                labels = [stimlabel == iC; zeros(size(nullData,1),1)];
                if windowing
                     labels = [repmat(stimlabel == iC,3,1); zeros(size(nullData,1),1)];
                end
                if iShuf==1
                    labels = labels;
                else
                    labels = labels(randperm(length(labels)));
                end

                for iL1=1:length(is.ENL1)  % loop over L1 penalty values
                    if L1regress
                        l1p = is.ENL1(iL1); l2p = 0; % alpha = L1/(2*L2+L1) ; lambda = 2*L2+L1
                        [beta, fitInfo] = lassoglm([trueData; nullData], labels, 'binomial', ...
                            'Alpha', l1p / (2*l2p+l1p), 'Lambda', 2*l2p + l1p, 'Standardize', false);
                    elseif HMMclassifier
                        [beta,fitInfo] = decodeHMMwrapper([trueData;nullData],labels,iL1);
                    else
                        %implies L2regression
                        l2p = [0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,50];
                        [beta, fitInfo] = lassoglm([trueData; nullData], labels, 'binomial', ...
                            'Alpha', realmin, 'Lambda', l2p(iL1), 'Standardize', false);
                    end
                    gnf{iShuf,is.whichTimes(iTT),iC,iL1}.beta = beta; gnf{iShuf,is.whichTimes(iTT),iC,iL1}.Intercept = fitInfo.Intercept;
                end
            end
        end
    end

    saveWrap(gnf, whenInMS, iSj);  % anonymous function to call 'save' in parfor
    disp(['sj' num2str(iSj) ': No. of Good trial = ' num2str(sum(subsetindex))]);
    toc;
end

%%%%%%%%% this second for loop iterates through subjects applying the
%%%%%%%%% ensemble of clasifiers saved above to the resting state data,
%%%%%%%%% saving the ensemble of outputted reactivation timecourses as a 
%%%%%%%%% structure called Rreds:

iRun = 2; % select the 2nd resting session
is.nShuf=1;
mkdir([analysisdir,'classifiers\TestResting4Cell\',dirname]);
saveWrap = @(Rreds, is, iSj, iRun) save([analysisdir,'classifiers\TestResting4Cell\',dirname,'Rreds' num2str(iSj) '_' num2str(iRun)], 'Rreds', 'is', '-v7.3');

for iSj= 1:length(is.goodsub)
    iSjGood=is.goodsub(iSj);
    % load classifier
    gnf = load([analysisdir,'classifiers\TrainedStim4Cell\',dirname,'Lsj' num2str(iSj)], 'gnf') ; % get the gnf variable with the regression models
    gnf =gnf.gnf;

    % Load Resting State Data
    dstInd = find(ismember(is.MEGruns{iSjGood},'rst'));
    dir_temp=[is.OPTPath strrep(is.fnDate{iSjGood},'-','_') '\' is.fnMEG{iSjGood} '_' num2str(dstInd(iRun),'%02d') '.ds\highpass_' num2str(is.highpass) '.opt'];
    opt=load(fullfile(dir_temp, 'opt.mat'));
    opt=opt.opt;    
    RST=spm_eeg_load(fullfile(dir_temp,[opt.results.spm_files_basenames{1,1}]));

    % Good timepoints/trials
    good_samples = ~all(badsamples(RST,':',':',':'));

    % Select MEG channel:
    chan_meg = indchantype(RST,'meeg');
    rRST = RST(chan_meg,:,:);

    % if use eyeblink, get rid of channels influenced by eyeblink, mostly in the frontal area
    if is.eyeBlinkSens==1
        rRST = rRST(~eyeBlinkSens,:,:);
    else
        sjchannel=squeeze(GoodChannel(iSj,:,:));
        goodchannindex= nansum(sjchannel')==length(is.MEGruns{1,iSjGood});
        rRST = rRST(logical(goodchannindex),:,:);
    end

    % get the interval between the start and the end
    evtypes = RST.events; evtypes = {evtypes(:).type}';
    evvals = RST.events; evvals = {evvals(:).value}';
    evvals(all(cellfun(@ischar,evvals),2),:) = {0}; % replace string with 0
    evvals=cell2mat(evvals); %convert cell to double
    evtimes = RST.events; evtimes = [evtimes(:).time]';

    RestingTriggerVals = [87, 88];  % 77 for the resting state BEFORE reward learning, 88 for the resting state AFTER reward learning
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

%     if is.delete==1
%         is.Zebcorrection =1;
%     else
%         is.Zebcorrection =0;
%     end
% 
%     if is.Zebcorrection ==1
%         [data, nRej] = artToZero(data);  % set artifacts to zeros (should not have this commented.)
%         disp(['sj' num2str(iSj) ', Resting' num2str(iRun) ': rejected ' num2str(100*nRej/prod(size(data))) '% of data'])
%     end

    data = scaleFunc(data);    

    % apply classifier to the clean data
    if isempty(RestStmTimes) 
        nTr = 1;
    else
        nTr = length(RestStmTimes);
    end

    Rreds = cell(is.nShuf,length(is.whichTimes),nTr,length(is.ENL1)); 

    for whichShuf=1:is.nShuf
        for iTT = 1:length(is.whichTimes)
            shuf1=is.whichTimes(iTT);
            for iTr=1:nTr
                for iL1=1:length(is.ENL1)
                    for iC=1:8
                        % use regression models 
                        %Rreds{whichShuf,iT,iTr,iL1}(:,iC) = 1 ./ (1 + exp(-(data * gnf{whichShuf,iT,iC,iL1}.beta + gnf{whichShuf,iT,iC,iL1}.Intercept)));   % matlab's built-in lassoglm
                        Rreds{whichShuf,shuf1,iTr,iL1}(:,iC) = -(data * gnf{whichShuf,shuf1,iC,iL1}.beta + gnf{whichShuf,shuf1,iC,iL1}.Intercept);   % work on the X*beta space!
                    end
                end
            end
        end
    end
    saveWrap(Rreds, is, iSj, iRun);    
    disp(['sj' num2str(iSj) ' finished']);    
end

%%%%%%%%% This third for loop iterates over each subject's ensemble of 
%%%%%%%%% reactivation timecourses, computing the Sequenceness measure for
%%%%%%%%% each timecourse
iRun = 2;
is.nShuf=24;
nsub=is.goodsub;
%load('half_between_permutation.mat');
[~,~,uniquePerms] = uperms([1:4],is.nShuf);
sfAll = nan(length(is.lgncy)+1, is.nShuf, length(is.whichSubj), 1, length(is.whichTimes), length(is.ENL1));
sbAll = nan(length(is.lgncy)+1, is.nShuf, length(is.whichSubj), 1, length(is.whichTimes), length(is.ENL1));

Tfwd = [0,1,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,1,0,0;0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0];

for iSj= 1:length(is.goodsub)

    maxTrials = 1; 

    sf = cell(is.nShuf, maxTrials);
    sb = cell(is.nShuf, maxTrials);

    S = load([analysisdir,'classifiers\TestResting4Cell\',dirname,'Rreds' num2str(iSj) '_' num2str(iRun)]); Rreds = S.Rreds;  % load this subject's preds

    nTr = size(Rreds,3);
    L1l = length(is.ENL1); 
    Ltt = length(is.whichTimes);
    maxLag = length(is.lgncy);
    nstates=8;

    tic

    for iTr=1:nTr
        tmp = squeeze(shiftdim(Rreds(1,:,iTr,:),-1)); 
        mpty = all(all(cellfun(@isempty, tmp),2),3); tmp(mpty,:,:) = [];  % remove empty trainingTimes
        prSj = permute(cell2mat(shiftdim(tmp,-2)),[1 2 4 3]); % this is samples*states*alpha*lambda*trainingTimes
        prSj = prSj(:,:, :); % alpha and trainingTimes put into a single dimension for the vectorization of sequenceness

       % Sequence 
       nP = size(prSj,3);
       sf_temp=nan(maxLag,is.nShuf,nP);
       sb_temp=nan(maxLag,is.nShuf,nP);
       scf_temp=nan(maxLag,is.nShuf,nP);
       scb_temp=nan(maxLag,is.nShuf,nP);

       for vec=1:nP

           X=squeeze(prSj(:,:,vec));

           if ~any(~isnan(X(:))) 
            continue
           end

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           nbins=maxLag+1;

           warning off
           dm=[toeplitz(X(:,1),[zeros(nbins,1)])];
           dm=dm(:,2:end);

           for kk=2:nstates
               temp=toeplitz(X(:,kk),[zeros(nbins,1)]);
               temp=temp(:,2:end);
               dm=[dm temp]; 
           end

           warning on

           Y=X;       

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           betas = nan(nstates*maxLag, nstates);
           bins=10;

           for ilag=1:bins%maxLag
               zinds = (1:maxLag:nstates*maxLag) + ilag - 1; 
               temp_zinds = (1:bins:nstates*maxLag) + ilag - 1; 
    %            betas(zinds,:)=pinv(dm(:,zinds))*Y;      
    %            temp = pinv([dm(:,zinds) ones(length(dm(:,zinds)),1)])*Y;
               temp = pinv([dm(:,temp_zinds) ones(length(dm(:,temp_zinds)),1)])*Y;
               betas(temp_zinds,:)=temp(1:end-1,:);
           end   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           X=[];
           Y=[];
           betasr=reshape(betas,[maxLag nstates nstates]);
           betasnbins64=reshape(betas,[maxLag nstates^2]);

           for iShuf=1:is.nShuf
               if iShuf==1  % the real transition matrix
                   T1 = Tfwd; T2 = Tfwd';
%                    T3 = Vfwd; T4 = Vfwd';
               else  % construct random transition matrix
                   X=[];
                   Y=[];
                   Xv=[];
                   Yv=[];

                   X = Tfwd (1:4,1:4); % Sequence 1: 1->2->3->4
                   Y = Tfwd (5:8,5:8); % Sequence 2: 5->6->7->8

%                    Xv = Vfwd (1:4,1:4);
%                    Yv = Vfwd (5:8,5:8);

                   rp = uniquePerms(iShuf,:);  % use the 30 unique permutations (is.nShuf should be set to 29)
                   Tfwd_temp=zeros(8,8);
                   Tfwd_temp(1:4,5:8)= X(rp,rp);
                   Tfwd_temp(5:8,1:4)=Y(rp,rp);        
                   T1 = Tfwd_temp;%Tfwd(rp,rp);
                   T2 = T1'; % backwards transitions

%                    Vfwd_temp=zeros(8,8);
%                    Vfwd_temp(1:4,5:8)= Xv(rp,rp);
%                    Vfwd_temp(5:8,1:4)=Yv(rp,rp);        
%                    T3 = Vfwd_temp;%Tfwd(rp,rp);
%                    T4 = T3'; % backwards transitions        
               end

               bbb=pinv([T1(:) T2(:) squash(eye(nstates)) squash(ones(nstates))])*(betasnbins64');
               sf_temp(:,iShuf,vec)=bbb(1,:);
               sb_temp(:,iShuf,vec)=bbb(2,:);    
           end
       end

       for iShuf=1:is.nShuf
           sf{iShuf,iTr} = nan(1, 1, length(is.lgncy)+1, L1l, Ltt);
           sb{iShuf,iTr} = nan(1, 1, length(is.lgncy)+1, L1l, Ltt);

           sf{iShuf,iTr}(1, 1, 2:end, :, :)=reshape(squeeze(sf_temp(:,iShuf,:)),[maxLag L1l Ltt]);
           sb{iShuf,iTr}(1, 1, 2:end, :, :)=reshape(squeeze(sb_temp(:,iShuf,:)),[maxLag L1l Ltt]);         
       end

    end

    % RESHAPE
    % aim is -> latency(61) * shuffles(20) * trials (120, padded) * trainTimes (1 OR 2) * alpha (L1 regulation = 10)
    sf2 = permute(cell2mat(sf), [3 1 2 5 4]); % sf2 = latency(61) * shuffles(20) * trials (1) * trainTimes (1/2) * alpha (L1 regulation = 10)
    sf2 = cat(3, sf2, nan(length(is.lgncy)+1, is.nShuf, maxTrials - size(sf2,3), Ltt, L1l));  % [pad with nans for trials after the end]
    sb2 = permute(cell2mat(sb), [3 1 2 5 4]);
    sb2 = cat(3, sb2, nan(length(is.lgncy)+1, is.nShuf, maxTrials - size(sb2,3), Ltt, L1l));  % [pad with nans for trials after the end]

    sfAll(:,:,iSj,:,:,:) = sf2; 
    sbAll(:,:,iSj,:,:,:) = sb2;

    disp(['Sub' num2str(iSj) ' Sequence Finished' ])  
    toc
end

cd([analysisdir,'classifiers\Sequence_by_Training_4Cell']);
mkdir(dirname);
save([dirname,'StimAll'],'sfAll','sbAll','is','-v7.3');
end
%% Plot replay effect broken by training time
close all
%cd('C:\Users\chiggins\My Data\Yunzhe data\StrLearn_MEGexp\Analysis');

dir1name='L1regression\';
dir2name='L2regression\';
dir3name='L1regression_window3\';
dir4name='L2regression_window3\';
dir5name='ARDsingletimepoint\';
dir6name='ARDregression_window3\';
dir7name='ARDsingletimepoint_bigprior';

load([analysisdir,'classifiers\Sequence_by_Training_4Cell\',dir1name,'StimAll.mat']);
%load(['/Volumes/CamsHD2/YunzheData/StrLearn_MEGexp/Analysis/classifiers/Sequence_by_Training_4Cell/',dir5name,'StimAll.mat']);
%nsub=is.goodsub;
nsub=1:length(is.goodsub);
iTT=9;
sfAll=squeeze(sfAll(2:end,:,nsub,1,iTT,:)); % timelag*nsubs*trainingtime*L1
sbAll=squeeze(sbAll(2:end,:,nsub,1,iTT,:)); % timelag*nsubs*trainingtime*L1



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
   % set(groot,'defaultAxesColorOrder',jet(size(sfb,3)));    
    plot(tTimes2, Mreplay, 'LineWidth', 1.5)
    xlabel('lag (ms)')     
    ylabel('L1-reverse\Leftarrow sequenceness \RightarrowL1-forward'),
    %ylim([0 0.15])        
    grid on;
    %legend(cellfun(@(x) ['Train: ' strtrim(x) 'ms'], cellstr(num2str((110:10:300)')), 'UniformOutput', false))
    title(['Sequenceiness for lambda = ',num2str(iL1)]);

    
    plot([tTimes2(1) tTimes2(end)], -npThreshAll*[1 1], 'k--'), plot([tTimes2(1) tTimes2(end)], npThreshAll*[1 1], 'k--')

end

%% and cross validate to choose a single parameter value:
shuf1=1;
figure, 
x0=40;
y0=40;
width=700;
height=400;
%for iTCriteria=1:10
for iSjLO = 1:length(nsub)
    iSjLI = setdiff(1:length(nsub),iSjLO);
    sfb=squeeze(nanmean(sfAll(:,shuf1,iSjLI,:)-sbAll(:,shuf1,iSjLI,:),3));
    %CV over all time:
    [~,iLCV(iSjLO)] = max(max(sfb,[],1));
    %iLCV(iSjLO)=5;
    % or CV over single timepoint:
    %[~,iLCV(iSjLO)] = max(sfb(iTCriteria,:),[],2);

    sfbCV(iSjLO,:,:) = sfAll(:,:,iSjLO,iLCV(iSjLO))-sbAll(:,:,iSjLO,iLCV(iSjLO));
end

Mreplay=squeeze(nanmean(sfbCV(:,:,1),1));
Sreplay=squeeze(nanstd(sfbCV(:,:,1),[],1))/sqrt(size(sfbCV(:,:,1),12));   
tTimes2=10:10:600;

%subplot(2,5,iTCriteria);
set(gcf,'units','points','position',[x0,y0,width,height])
%set(groot,'defaultAxesColorOrder',jet(size(sfb,3)));    
plot(tTimes2, Mreplay, 'LineWidth', 1.5)
xlabel('lag (ms)')     
ylabel('L1-reverse\Leftarrow sequenceness \RightarrowL1-forward'),
%ylim([0 0.15])        

%legend(cellfun(@(x) ['Train: ' strtrim(x) 'ms'], cellstr(num2str((110:10:300)')), 'UniformOutput', false))
if L1regress
    title({'Sequenceness, Study II Lasso,',[' criteria maximised across all timepoints']});
elseif HMMclassifier
    title({'Sequenceness, Study II ARD,',[' criteria maximised across all timepoints']});
else
    title({'Sequenceness, Study II Ridge,',[' criteria maximised across all timepoints']});
end

npThreshAll = max(max(abs(nanmean(sfbCV(:,:,2:end),1))));hold on;
plot([tTimes2(1) tTimes2(end)], -npThreshAll*[1 1], 'k--'), plot([tTimes2(1) tTimes2(end)], npThreshAll*[1 1], 'k--')
grid on;
%end

%% Save an inferred Replay Timecourse:

%define peak sequenceness time:
[~,lag] = max(abs(Mreplay));
for iSj=1:length(is.goodsub)
    % take that subject's CV parameter selection:
    iL = iLCV(iSj);
    
    % select second resting state session:
    iRun = 2;
    
    % load reactivation timecourse:
    S = load([analysisdir,'classifiers\TestResting4Cell\',dirname,'Rreds' num2str(iSj) '_' num2str(iRun)]); Rreds = S.Rreds;  % load this subject's preds
    
    reactivationtimecourse = S.Rreds{1,is.whichTimes(iTT),:,iL};
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
replaydir = [analysisdir,'classifiers\Sequence_by_Training_4Cell\',dirname,'STUDYII_ReplayOnset\'];
mkdir(replaydir)
if iRun==1
    fname = 'STUDYII_PreplayOnset.mat';
else
    fname = 'STUDYII_ReplayOnset.mat';
end
save([replaydir,fname],'ToRall');

%% analyse pooled saved values and plot as image:
figure('Position', [1 41 1920 1084]);
for iCond=1:6
   
    if iCond==1
        L1regress=true;
        HMMclassifier=false;
        windowing=false; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
        dirname='L1regression\';
    elseif iCond==2
        L1regress=false;
        HMMclassifier=false;
        windowing=false; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
        dirname='L2regression\';
    elseif iCond==3
        L1regress=true;
        HMMclassifier=false;
        windowing=true; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
        dirname='L1regression_window3\';
    elseif iCond==4
        L1regress=false;
        HMMclassifier=false;
        windowing=true; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
        dirname='L2regression_window3\';
    elseif iCond==5
        L1regress=false;
        HMMclassifier=true;
        windowing=false; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
        dirname='ARDsingletimepoint\';
    elseif iCond==6
        L1regress=true;
        HMMclassifier=true;
        windowing=true; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
        dirname='ARDregression_window3\';
%     elseif iCond==6
%         L1regress=true;
%         HMMclassifier=true;
%         windowing=true; % windowing: denotes classifiers trained on square window of samples around the given is.whichTimes
%         dirname='ARDsingletimepoint_bigprior\';
     end
     plotlocs=[1,2,4,5,3,6];
         
    cd(analysisdir);
    load([analysisdir,'classifiers\Sequence_by_Training_4Cell\',dirname,'\StimAll.mat']);
    nsub=1:length(is.goodsub);
    sfAll=squeeze(sfAll(2:end,:,nsub,1,:,:)); % timelag*nsubs*trainingtime*L1
    sbAll=squeeze(sbAll(2:end,:,nsub,1,:,:)); % timelag*nsubs*trainingtime*L1

    t_points=[3:5]; % which points to compute sequenceness over - should be all, 1:60
    for iL1=1:10
        for iTT=1:length(is.whichTimes)
        %subplot(2,5,iL1)
            sfb=squeeze(sfAll(:,:,:,iTT,iL1)-sbAll(:,:,:,iTT,iL1));

            sfb_perms = squeeze(sfb(:,2:end,:));

            sfb=squeeze(sfb(t_points,1,:));
            Mreplay=squeeze(nanmean(sfb,2));
            Sreplay=squeeze(nanstd(sfb,[],2))/sqrt(size(sfb,2));   
            tTimes2=10:10:600;

            npThreshAll = max(max(abs(nanmean(sfb_perms(t_points,:,:),3))));hold on;
            %plot([tTimes2(1) tTimes2(end)], -npThreshAll*[1 1], 'k--'), plot([tTimes2(1) tTimes2(end)], npThreshAll*[1 1], 'k--')
            seqmat(iL1,iTT,:) = Mreplay ./ npThreshAll;
        end
    end
    %imagesc(seqmat(:,:));
    max_replay_points{iCond} = squeeze(max(seqmat,[],3));
    min_replay_points{iCond} = squeeze(min(seqmat,[],3));
    replay_points{iCond} = max(abs(cat(3,max_replay_points{iCond},min_replay_points{iCond})),[],3) .* ...
        [[max_replay_points{iCond}>abs(min_replay_points{iCond})] + ...
        [max_replay_points{iCond}<abs(min_replay_points{iCond})]*-1];
    subplot(2,3,plotlocs(iCond));
    imagesc(replay_points{iCond});caxis([-2,2]);hold on;
    hold on;
    pos_sig_points = replay_points{iCond}>=1;
    neg_sig_points = replay_points{iCond}<=-1;
    contour(pos_sig_points,1,'k--');
    contour(neg_sig_points,1,'r--');
    set(gca,'YTick',1:10);
    if iCond<5
        ylabel('\lambda regularisation penalty');
    else
        ylabel('\alpha regularisation prior');
    end
    if iCond==1 || iCond==3
        set(gca,'YTickLabel',0.001:0.001:0.01);
        
        if iCond==1
            title('L1 regression on single timepoints');
        else
            title('L1 regression on windowed timepoints');
        end
    elseif iCond==2 || iCond==4
        l2p = [0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,50];
        set(gca,'YTickLabel',l2p);
        if iCond==2
            title('L2 regression on single timepoints');
        else
            title('L2 regression on windowed timepoints');
        end
    else
        alpha = round(exp([-5:1:5]),2,'significant');
        set(gca,'YTickLabel',alpha);
        if iCond==5
            title('ARD regression on single timepoints');
        else
            title('ARD regression on windowed timepoints');
%             alpha = round(exp([0:2:20]),2);
%             set(gca,'YTickLabel',alpha);
        end
    end
    set(gca,'XTick',2:2:20);
    set(gca,'XTickLabel',120:20:400);
    xlabel('Time of training (msec)');
    if iCond==2
        l1=legend({'Positive significant sequenceness','Reverse significant sequenceness'});
        l1.Position=[0.5967 0.4974 0.1219 0.0337];
    end
    
    
end
%h{1}=contour([NaN,NaN;NaN,NaN],1,'k--');
%h{2}=contour([NaN,NaN;NaN,NaN],1,'r--');


