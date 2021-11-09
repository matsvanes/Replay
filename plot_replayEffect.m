%% Plot replay effect broken by training time
close all
dir1name='L1regression/';

load([is.AnalysisPath,'classifiers/Sequence_by_Training_4Cell/','StimAll.mat']);
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