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
close all, clear
useReplayModels = false;
A = 0.2; % general cosine amplitude
f=10; % Hz
chanphase = 'rand'; % can be 'beta', 'same', or 'rand';

% create "resting state" data
nchan = 242;
fs = 1000;
time = 1/fs:1/fs:30;
plot_t = time(1:1000); % selection of the time axis

if useReplayModels
  % study info
  MFStartup_studyII;
  analysisdir = '/Volumes/T5_OHBA/analysis/replay/';
  
  % load trained decoding model
  
  iSj=1;
  iSjGood=is.goodsub(iSj);
  % load classifier
  load([analysisdir,'classifiers/TrainedStim4Cell/','L1regression/','Lsj' num2str(iSj)], 'gnf') ; % get the gnf variable with the regression models
  
  iL1=3; % lambda = 0.003
  
  gnf=gnf(:,iL1); % contains the learned decoding model for 8 stimuli
  
  for ii=1:8
    beta{ii} = gnf{ii}.beta;
    intercept(ii) = gnf{ii}.Intercept;
  end
else
  % simulate decoding model
  nrep=500;
  truedata_labels = [];
  for ii=1:8
    X{ii} = rand(nrep,nchan);
    X{ii} = X{ii} + rand(1,nchan);
    truedata_labels = [truedata_labels; ii*ones(nrep,1)];
  end
  nulldata = scaleFunc(rand(nrep*8,nchan));
  truedata = scaleFunc(cat(1,X{:}));
  for ii=1:8
    
    labels = [truedata_labels == ii; zeros(size(nulldata,1),1)];
    l2p=0;l1p=0.003;
    [beta{ii}, fitInfo{ii}] = lassoglm([truedata; nulldata], labels, 'binomial', ...
      'Alpha', l1p / (2*l2p+l1p), 'Lambda', 2*l2p + l1p, 'Standardize', false);
    intercept{ii} = fitInfo{ii}.Intercept;
  end
  
end

data_woa = rand(nchan,numel(time)); % random data

% find the largest average beta to construct weights for alpha (largest
% beta weights should also have highest alpha signal).
for ii=1:8
  betas(ii,:) = beta{ii};
end
avgbeta = mean(abs(betas),1);
alpha_amp = A*sort(normrnd(0.5,0.15,[nchan,1])); % different weights for each channel

% give the sensors with largest abs(beta) also the largest alpha amp
[~, idx] = sort(avgbeta);
alpha_amp(idx)=alpha_amp;
figure; plot(avgbeta, alpha_amp, '.'); xlabel('avgbeta'), ylabel('alpha amp'), title('those sensors with high beta weights also have strongest alpha components')


switch chanphase
  case 'same' % we can give all channels the same phase
    theta=2*pi*rand(1)*ones(nchan,1);
  case 'rand' % random phase for each chan
    theta=2*pi*rand(nchan,1);
  case 'beta' % different phase for each channel, depending on beta weight
    theta=(2*pi)*rand(nchan,1); theta(idx)=theta;
end
alpha_signal = alpha_amp.*cos(2*pi*f*time + theta);
data_wa = data_woa + alpha_signal;
figure; subplot(2,1,1); plot(plot_t,[1:10]' + alpha_signal(1:10,1:numel(plot_t)))
subplot(2,1,2); plot(plot_t,[1:10]' + data_wa(1:10,1:numel(plot_t)))

clear i
c = cos(2*pi*f*time) + i*sin(2*pi*f*time); % raw alpha signal
phase = angle(c);
phase = repmat(phase, [8 1])';
figure; subplot(2,1,1), plot(plot_t, real(c(1:numel(plot_t)))), xlabel('time'), title('amplitude')
subplot(2,1,2), plot(plot_t, phase(1:numel(plot_t),1)), xlabel('time'), title('phase')

data_woa = scaleFunc(data_woa);
data_wa = scaleFunc(data_wa);
f1=figure; f1.WindowState='maximized';
subplot(2,1,1), plot(plot_t, data_woa(idx(end), 1:numel(plot_t))), title('data without alpha'), xlabel('time')
subplot(2,1,2), plot(plot_t, data_wa(idx(end), 1:numel(plot_t))), title('data with alpha'), xlabel('time')

%% Now apply the decoding models on the simulated data
% First the data w/o alpha
r_woa = zeros(numel(time), 8);
r_wa = zeros(numel(time), 8);
for ii=1:8
  r_woa(:,ii) = (data_woa' * beta{ii} + intercept{ii});
  r_wa(:,ii) = (data_wa' * beta{ii} + intercept{ii});
end

% apply logistic sigmoid to go to probabilities:
prob_woa=1./(1+exp(-r_woa));
prob_wa=1./(1+exp(-r_wa));

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
f2=figure; f2.WindowState='maximized';
subplot(2,2,1); plot(plot_t,prob_woa(1:numel(plot_t),:)), title('probability - no alpha')
subplot(2,2,2); pth = prob_woa; pth(pth<prob_woa(sortIdx_woa(perc*100)))=nan; plot(plot_t,pth(1:numel(plot_t),:)), title('probability (thresholded) - no alpha')
subplot(2,2,3); plot(plot_t,prob_wa(1:numel(plot_t),:)), title('probability - with alpha')
subplot(2,2,4); pth = prob_wa; pth(pth<prob_wa(sortIdx_wa(perc*100)))=nan; plot(plot_t,pth(1:numel(plot_t),:)), title('probability (thresholded) - with alpha')
%{
% do the same thing for random probabilities (permutations)
nperm=1000;
phase_woa_perm = zeros(8,numel(time)*perc/100, nperm);
phase_wa_perm = zeros(8,numel(time)*perc/100, nperm);
for ii=1:nperm
  for j=1:8
    permIdx_woa = randperm(numel(time),numel(time)*perc/100);
    phase_woa_perm(j,:,ii) = phase(permIdx_woa,j);
    
    permIdx_wa = randperm(numel(time),numel(time)*perc/100);
    phase_wa_perm(j,:,ii) = phase(permIdx_wa,j);
  end
end
%}
% alpha phase of highest probabilities
f3=figure; f3.WindowState='maximized';
for ii = 1:8
  % W/o alpha
  p_nonuniform = circ_otest(phase_woa_probthr(:,ii));
  axesHandle(ii,1) = subplot(2,8,ii);
  polarAxesHandle(ii,1) = polaraxes('Units',axesHandle(ii,1).Units,'Position',axesHandle(ii,1).Position);
  delete(axesHandle(ii,1));
  %   polarhistogram(polarAxesHandle(ii,1),reshape(phase_woa_perm(ii,:,:), 1, []),20); hold on
  polarhistogram(polarAxesHandle(ii,1),repmat(phase_woa_probthr(:,ii), 500,1),20)
  title(sprintf('w/o alpha \n stimulus %d \n nonuniformity pval=%s', ii, num2str(round(p_nonuniform,2, 'significant'))))
  
  % With alpha
  p_nonuniform = circ_otest(phase_wa_probthr(:,ii));
  axesHandle(ii,2) = subplot(2,8,ii+8);
  polarAxesHandle(ii,2) = polaraxes('Units',axesHandle(ii,2).Units,'Position',axesHandle(ii,2).Position);
  delete(axesHandle(ii,2));
  %   polarhistogram(polarAxesHandle(ii,2),reshape(phase_wa_perm(ii,:,:), 1, []),20); hold on
  polarhistogram(polarAxesHandle(ii,2),repmat(phase_wa_probthr(:,ii), 500,1),20)
  title(sprintf('w/ alpha \n stimulus %d \n nonuniformity pval=%.1e', ii, p_nonuniform))
end
suptitle(sprintf('alpha phase preference for highest %d percent probabilities',perc))

% compare correlations phase with decoding probability (with and without alpha)
for ii =1:8
  r_prob_phase_woa(ii) = circ_corrcl(phase(:,ii), log10(prob_woa(:,ii)));
  r_prob_phase_wa(ii) = circ_corrcl(phase(:,ii), log10(prob_wa(:,ii)));
end
figure; plot(r_prob_phase_woa, 'O'), hold on, plot(r_prob_phase_wa, 'O'), xlabel('stimulus'), ylabel('R'),title('circular correlation between phase and probability'), legend({'w/o alpha','w/ alpha'})








