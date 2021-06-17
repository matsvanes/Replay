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
clear
% study info
MFStartup_studyII;
analysisdir = '/Volumes/T5_OHBA/analysis/replay/';

% load trained decoding model
iSj=1;
iSjGood=is.goodsub(iSj);
% load classifier
load([analysisdir,'classifiers/TrainedStim4Cell/','L1regression/','Lsj' num2str(iSj)], 'gnf') ; % get the gnf variable with the regression models

iL1=3; % lambda = 0.003

gnf=gnf(:,iL1);

% create "resting state" data
nchan = size(gnf{1}.beta,1);
fs = 100;
time = 1/fs:1/fs:300;
data_woa = rand(nchan,numel(time));

% find the largest average beta to construct weights for alpha (largest
% beta weights should also have highest alpha signal).
for i=1:8
  betas(i,:) = gnf{i}.beta;
end
avgbeta = mean(betas,1);
alpha_amp = 1; % general alpha amplitude
alpha_amp = alpha_amp*sort(rand(nchan,1)); % different weights for each channel

% give the sensors with largest abs(beta) also the largest alpha amp
[~, idx] = sort(abs(avgbeta));
alpha_amp(idx)=alpha_amp;

f=10;
% we can give all channels the same phase
phase=2*pi*rand(1)*zeros(nchan,1);
% alternatively, different phase for each channel
phase=(2*pi)*rand(nchan,1);
alpha_signal = alpha_amp.*cos(2*pi*f*time + phase);
data_wa = data_woa + alpha_signal;


f1=figure; f1.WindowState='maximized';
subplot(2,1,1), plot(time(1:100), data_woa(idx(end), 1:100)), title('data without alpha'), xlabel('time')
subplot(2,1,2), plot(time(1:100), data_wa(idx(end), 1:100)), title('data with alpha'), xlabel('time')

%% Now apply the decoding models on the simulated data
% First the data w/o alpha
r_woa = zeros(numel(time), 8);
r_wa = zeros(numel(time), 8);
for i=1:8
  r_woa(:,i) = -(data_woa' * gnf{i}.beta + gnf{i}.Intercept);
  r_wa(:,i) = -(data_wa' * gnf{i}.beta + gnf{i}.Intercept);
end

% apply logistic sigmoid to go to probabilities:
prob_woa=1./(1+exp(-r_woa));
prob_wa=1./(1+exp(-r_wa));

% threshold probabilities
iThr_woa = prob_woa>0.95;
iThr_wa = prob_wa>0.95;

t=time(1:100);
f2=figure; f2.WindowState='maximized';
subplot(2,2,1); plot(t,prob_woa(1:100,:)), ylim([0.5 1]), title('probability - no alpha')
subplot(2,2,2); pth=iThr_woa(1:100,:).*prob_woa(1:100,:); pth(pth==0)=nan; plot(t,pth), ylim([0.5 1]), title('probability (thresholded) - no alpha')
subplot(2,2,3); plot(t,prob_wa(1:100,:)), ylim([0.5 1]), title('probability - with alpha')
subplot(2,2,4); pth=iThr_wa(1:100,:).*prob_wa(1:100,:); pth(pth==0)=nan; plot(t,pth), ylim([0.5 1]), title('probability (thresholded) - with alpha')

c = cos(2*pi*f*time); % raw alpha signal
phase = angle(fft(c));

% alpha phase of highest probabilities
for i=1:8
  X_woa{i} = phase(iThr_woa(:,i));
  X_wa{i} = phase(iThr_wa(:,i));
end

f3=figure; f3.WindowState='maximized';
for ii = 1:8
  axesHandle(ii,1) = subplot(2,8,ii);
  polarAxesHandle(ii,1) = polaraxes('Units',axesHandle(ii,1).Units,'Position',axesHandle(ii,1).Position);
  delete(axesHandle(ii,1));
  polarhistogram(polarAxesHandle(ii,1),X_woa{ii})
  title(sprintf('no alpha \n stimulus %d', ii))
  
  axesHandle(ii,2) = subplot(2,8,ii+8);
  polarAxesHandle(ii,2) = polaraxes('Units',axesHandle(ii,2).Units,'Position',axesHandle(ii,2).Position);
  delete(axesHandle(ii,2));
  polarhistogram(polarAxesHandle(ii,2),X_wa{ii})
  title(sprintf('with alpha \n stimulus %d', ii))
end
suptitle('alpha phase preference for highest 5% probabilities')
figure(2);
figure(1);
