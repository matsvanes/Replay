% according to Cam, a linear classifier will always have a certain phase
% preference. Meaning even if the classifier is trained on data that does
% not contain a strong alpha signal, when it is then applied to resting
% state data that has a strong alpha phase coherence over channels - we
% will consequentially have phase dependent reactivation.
% We can test this by using the learned weights to decode noise (which
% should be at chance level). If we then add an alpha signal with strong
% coherence over channels (either random channels or according to strength
% of decoding weights) - and we still have chance level decoding, we're
% probably fine.
%close all, clear
clear
rng('shuffle')
decoding_model = 'lassoglm'; % can be 'lassoglm', 'svmrbf'
useReplayModels = false;

doplot=true;

% create "resting state" data
nchan = 12;
fs = 1000;
time = 1/fs:1/fs:10;
plot_t = time(1:1000); % selection of the time axis


% for the alpha signal
cnt=1;
AA = 0.5%logspace(-2, 1, 20);
for A = AA%[0.01, 0.05, 0.1, 0.5, 1, 5] % general cosine amplitude
  for rep=1
    f=10; % Hz
    chanphase = 'rand'; % can be 'beta', 'same', or 'rand';
    nStim=2;
    alpha_amp = A*sort(normrnd(0.5,0.15,[nchan,1])); % different weights for each channel
    switch chanphase
      case 'same' % we can give all channels the same phase
        theta=2*pi*rand(1)*ones(nchan,1);
      case 'rand' % random phase for each chan
        theta=2*pi*rand(nchan,1);
      case 'beta' % different phase for each channel, depending on beta weight
        theta=(2*pi)*rand(nchan,1); theta(idx)=theta;
    end


    if useReplayModels
      % study info
      MFStartup_studyII;
      analysisdir = '/Volumes/T5_OHBA/analysis/replay/';

      % load trained decoding model

      iSj=1;
      iSjGood=is.goodsub(iSj);
      % load classifier
      load(['/Users/matsvanes/Data/YunzheData/StrLearn_MEGexp/Analysis/classifiers/TrainedStim4Cell/Lsj' num2str(iSj)], 'gnf') ; % get the gnf variable with the regression models

      iL1=5; % lambda = 0.003

      gnf=squeeze(gnf(:,37,:,iL1)); % contains the learned decoding model for 8 stimuli

      for ii=1:nStim
        beta{ii} = gnf{ii}.beta;
        intercept{ii} = gnf{ii}.Intercept;
      end
    else
      % simulate decoding model
      nrep=500;
      truedata_labels = [];
      for ii=1:nStim
        X{ii} = rand(nrep,nchan);
        X{ii} = X{ii} + rand(1);%rand(1,nchan);
        truedata_labels = [truedata_labels; ii*ones(nrep,1)];
      end
      nulldata=rand(nrep*nStim,nchan);
      nulldata = scaleFunc(nulldata);
      truedata=cat(1,X{:});
      truedata = scaleFunc(truedata);
      for ii=1:nStim

        labels = [truedata_labels == ii; zeros(size(nulldata,1),1)];
        switch decoding_model
          case 'lassoglm'
            l2p=0;l1p=0.003;
            [beta{ii}, fitInfo{ii}] = lassoglm([truedata; nulldata], labels, 'binomial', ...
              'Alpha', l1p / (2*l2p+l1p), 'Lambda', 2*l2p + l1p, 'Standardize', false);
            intercept{ii} = fitInfo{ii}.Intercept;
            % beta{ii} = rand(12,1)%.*randi([0,1],12,1);
            % beta{ii}=(beta{ii}-(beta{ii}~=0).*mean(beta{ii}(beta{ii}~=0)))./(std(beta{ii}(beta{ii}~=0)));

          case 'svmrbf'
            model{ii} = fitrsvm([truedata; nulldata], labels, 'KernelFunction', 'rbf', 'KernelScale', 'auto', 'Standardize', true);
        end
      end
    end

    data_woa = rand(nchan,numel(time)); % random data

    if useReplayModels
      % find the largest average beta to construct weights for alpha (largest
      % beta weights should also have highest alpha signal).
      for ii=1:nStim
        betas(ii,:) = beta{ii};
      end
      avgbeta = mean(abs(betas),1);

      % give the sensors with largest abs(beta) also the largest alpha amp
      [~, idx] = sort(avgbeta);
      alpha_amp(idx)=alpha_amp;
      figure; plot(avgbeta, alpha_amp, '.'); xlabel('avgbeta'), ylabel('alpha amp'), title('those sensors with high beta weights also have strongest alpha components')
    end
    alpha_signal = alpha_amp.*sawtooth(2*pi*f*time+theta);%cos(2*pi*f*time + theta);
    data_wa = data_woa + alpha_signal;
    data_woa = data_woa ;
    if doplot
      figure;
      subplot(3,1,1); plot(plot_t,[1:10]' + data_woa(end-9:end,1:numel(plot_t))), title('random data')
      subplot(3,1,2); plot(plot_t,[1:10]' + alpha_signal(end-9:end,1:numel(plot_t))), title('raw alpha')
      subplot(3,1,3); plot(plot_t,[1:10]' + data_wa(end-9:end,1:numel(plot_t))); title('random data with alpha')
      suptitle(sprintf('alpha SNR: %d', A/0.5))
    end
    clear i
    c = cos(2*pi*f*time) + i*sin(2*pi*f*time); % raw alpha signal
    phase = angle(c);
    phase = repmat(phase, [nStim 1])';
    if doplot
      figure; subplot(2,1,1), plot(plot_t, real(c(1:numel(plot_t)))), xlabel('time'), title('amplitude')
      subplot(2,1,2), plot(plot_t, phase(1:numel(plot_t),1)), xlabel('time'), title('phase')
    end
    data_woa = scaleFunc(data_woa);
    data_wa = scaleFunc(data_wa);

    if useReplayModels
      %   f1=figure; f1.WindowState='maximized';
      subplot(2,1,1), plot(plot_t, data_woa(idx(end), 1:numel(plot_t))), title('random data without alpha'), xlabel('time')
      subplot(2,1,2), plot(plot_t, data_wa(idx(end), 1:numel(plot_t))), title('random data with alpha'), xlabel('time')
    end

    %% Now apply the decoding models on the simulated data
    % demean spatial map at each time point
    if 1
      data_woa  = (data_woa-mean(data_woa,1));
      data_wa  = (data_wa-mean(data_wa,1));
      %   data_woa  = zscore(data_woa,[],2);
      %   data_wa  = zscore(data_wa,[],2);
    end

    % demean beta maps
    if 0
      for i=1:nStim
        % only the nonzero element
        beta{i}(beta{i}~=0) = (beta{i}(beta{i}~=0) - mean(beta{i}(beta{i}~=0)));%./std(beta{i}(beta{i}~=0));
        % all elements
        % beta{i} = demean(beta{i});
      end
    end


    r_woa = zeros(numel(time), nStim);
    r_wa = zeros(numel(time), nStim);
    for ii=1:nStim
      switch decoding_model
        case 'lassoglm'
          r_woa(:,ii) = (data_woa' * beta{ii} + intercept{ii});
          r_wa(:,ii) = (data_wa' * beta{ii} + intercept{ii});
        case 'svmrbf'
          r_woa(:,ii) = predict(model{ii}, data_woa');
          r_wa(:,ii) = predict(model{ii}, data_wa');
        case 'spatialmap'
          r_woa(:,ii) = abs(beta{ii}'*data_woa);
          r_wa(:,ii) = abs(beta{ii}'*data_wa);
      end
    end

    % apply logistic sigmoid to go to probabilities:
    if strcmp(decoding_model, 'lassoglm') || strcmp(decoding_model, 'svmrbf')
      prob_woa=1./(1+exp(-r_woa));
      prob_wa=1./(1+exp(-r_wa));
    end

    % threshold probabilities
    perc=5;

    % iThr_woa = prob_woa>0.95;
    [~, sortIdx_woa] = sort(prob_woa, 'descend');
    sortIdx_woa_sel = sortIdx_woa(1:size(prob_woa,1)*perc/100,:);

    % iThr_wa = prob_wa>0.95;
    [~, sortIdx_wa] = sort(prob_wa, 'descend');
    sortIdx_wa_sel = sortIdx_wa(1:size(prob_wa,1)*perc/100,:);

    % phases at highest percentile probability
    phase_woa_probthr = phase(sortIdx_woa_sel);
    phase_wa_probthr = phase(sortIdx_wa_sel);

    % plot probabilities and selected probabilities
    % f2=figure; f2.WindowState='maximized';
    if doplot
      subplot(1,2,1); plot(plot_t,prob_woa(1:numel(plot_t),:)), title('probability - no alpha')
      subplot(1,2,2); plot(plot_t,prob_wa(1:numel(plot_t),:)), title('probability - with alpha')
    end
    %{
% do the same thing for random probabilities (permutations)
nperm=1000;
phase_woa_perm = zeros(nStim,numel(time)*perc/100, nperm);
phase_wa_perm = zeros(nStim,numel(time)*perc/100, nperm);
for ii=1:nperm
  for j=1:nStim
    permIdx_woa = randperm(numel(time),numel(time)*perc/100);
    phase_woa_perm(j,:,ii) = phase(permIdx_woa,j);
    
    permIdx_wa = randperm(numel(time),numel(time)*perc/100);
    phase_wa_perm(j,:,ii) = phase(permIdx_wa,j);
  end
end
    %}
    % alpha phase of highest probabilities
    % f3=figure; f3.WindowState='maximized';
    for ii = 1:nStim
      % W/o alpha
      [p_woa, z_woa] = circ_rtest(phase_woa_probthr(:,ii));
      [p_wa, z_wa] = circ_rtest(phase_wa_probthr(:,ii));
      zwoa(cnt,rep) = z_woa;
      zwa(cnt,rep) = z_wa;
      if doplot
        axesHandle(ii,1) = subplot(2,nStim,ii);
        polarAxesHandle(ii,1) = polaraxes('Units',axesHandle(ii,1).Units,'Position',axesHandle(ii,1).Position);
        delete(axesHandle(ii,1));
        %   polarhistogram(polarAxesHandle(ii,1),reshape(phase_woa_perm(ii,:,:), 1, []),20); hold on
        polarhistogram(polarAxesHandle(ii,1),repmat(phase_woa_probthr(:,ii), 500,1),20)
        title(sprintf('w/o alpha \n stimulus %d \n nonuniformity pval=%s', ii, num2str(round(p_woa,2, 'significant'))))

        % With alpha
        axesHandle(ii,2) = subplot(2,nStim,ii+nStim);
        polarAxesHandle(ii,2) = polaraxes('Units',axesHandle(ii,2).Units,'Position',axesHandle(ii,2).Position);
        delete(axesHandle(ii,2));
        %   polarhistogram(polarAxesHandle(ii,2),reshape(phase_wa_perm(ii,:,:), 1, []),20); hold on
        polarhistogram(polarAxesHandle(ii,2),repmat(phase_wa_probthr(:,ii), 500,1),20)
        title(sprintf('w/ alpha \n stimulus %d \n nonuniformity pval=%.1e', ii, p_wa))
      end
    end
    % suptitle(sprintf('alpha phase preference for highest %d percent probabilities',perc))
  end
  cnt=cnt+1;
end

% compare correlations phase with decoding probability (with and without alpha)
% for ii =1:nStim
%   r_prob_phase_woa(ii) = circ_corrcl(phase(:,ii), log10(prob_woa(:,ii)));
%   r_prob_phase_wa(ii) = circ_corrcl(phase(:,ii), log10(prob_wa(:,ii)));
% end
% figure; plot(r_prob_phase_woa, 'O'), hold on, plot(r_prob_phase_wa, 'O'), xlabel('stimulus'), ylabel('R'),title('circular correlation between phase and probability'), legend({'w/o alpha','w/ alpha'})

