function [sf2, sb2] = compute_sequenceness(iSj, is, iRun)
%% Compute sequenceness on resting state activations
if is.usePrecomputed
  fprintf('Loading results - sequenceness \n')
  % load sequenceness and select for specific subject
  try
    s = load([is.AnalysisPath,'classifiers/Sequence_by_Training_4Cell/StimAll']);
    sf2 = squeeze(s.sfAll(:,:,iSj,:,:));
    sb2 = squeeze(s.sbAll(:,:,iSj,:,:));
  catch
    fprintf('Loading in data was not possible - recomputing now \n')
    is.usePrecomputed=false;
    [sf2, sb2] = compute_sequenceness(iSj, is, iRun);
  end
else

  %%%%%%%%% This third for loop iterates over each subject's ensemble of
  %%%%%%%%% reactivation timecourses, computing the Sequenceness measure for
  %%%%%%%%% each timecourse
  is.nShuf=factorial(4);
  [~,~,uniquePerms] = uperms(1:4,is.nShuf);
  % number of time lags * number of shuffles * number of subjects * 1 *
  % number of lambda values
  sfAll = nan(length(is.lgncy)+1, is.nShuf, length(is.goodsub), 1, length(is.ENL1));
  sbAll = nan(length(is.lgncy)+1, is.nShuf, length(is.goodsub), 1, length(is.ENL1));
  % transition matrix. backward is transpose of this
  Tfwd = [0,1,0,0,0,0,0,0;0,0,1,0,0,0,0,0;0,0,0,1,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,1,0,0;0,0,0,0,0,0,1,0;0,0,0,0,0,0,0,1;0,0,0,0,0,0,0,0];

  sf = cell(is.nShuf, 1);
  sb = cell(is.nShuf, 1);

  S = load([is.AnalysisPath,'classifiers/TestResting4Cell/','Rreds' num2str(iSj) '_' num2str(iRun)]); Rreds = S.Rreds;  % load this subject's preds

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
      for kk=1:is.nStim
        temp=toeplitz(X(:,kk), zeros(nbins,1));
        temp=temp(:,2:end);
        dm=[dm temp]; % this is now samples x (timelags*is.nStim)
      end

      warning on

      Y=X;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % constructing alpha cycle regressors
      betas = nan(is.nStim*maxLag, is.nStim);
      bins=10; % because 10 Hz = 10 samples @ 100 Hz
      % CAM'S GONNA FIGURE OUT WHAT IS HAPPENING HERE! (only regressing
      % at each alpha cycle lag??)
      for ilag=1:bins%maxLag
        temp_zinds = (1:bins:is.nStim*maxLag) + ilag - 1; % indices of all alpha cycles away from the lag of interest
        temp = pinv([dm(:,temp_zinds) ones(length(dm(:,temp_zinds)),1)])*Y; % 49 alpha cycles are regressed out
        betas(temp_zinds,:)=temp(1:end-1,:);
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Regressing the templates (forward and backward) sequences on the
      % time-lagged regression weights organised into a transition matrix
      % (transition matrix per time lag)
      X=[];
      Y=[];
      betasr=reshape(betas,[maxLag is.nStim is.nStim]);
      betasnbins64=reshape(betas,[maxLag is.nStim^2]);

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

        bbb=pinv([T1(:) T2(:) squash(eye(is.nStim)) squash(ones(is.nStim))])*(betasnbins64'); % MVE: what's happening here?? Why do we multiply the alpha cycle's betas with the other betas??
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
