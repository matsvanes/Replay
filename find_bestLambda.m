function [iLCV, sfbCV, Mreplay, Sreplay] = find_bestLambda(is, sfAll, sbAll)
%% and cross validate to choose a single parameter value:
% MVE: Note that it doens't do cross validation to select the best and same
% parameter value for the whole group. It uses the value that is best in
% the group without the subject of interest.
shuf1=1;
nsub=length(is.goodsub);
if is.plot,figure, end
x0=40;
y0=40;
width=700;
height=400;
for iSjLO = 1:nsub
  iSjLI = setdiff(1:length(nsub),iSjLO);
  sfb=squeeze(nanmean(sfAll(:,shuf1,iSjLI,:)-sbAll(:,shuf1,iSjLI,:),3));
  %CV over all time:
  [~,iLCV(iSjLO)] = max(max(sfb,[],1)); % the best l1 over all other subjects
  sfbCV(iSjLO,:,:) = sfAll(:,:,iSjLO,iLCV(iSjLO))-sbAll(:,:,iSjLO,iLCV(iSjLO));
end

Mreplay=squeeze(nanmean(sfbCV(:,:,1),1));
Sreplay=squeeze(nanstd(sfbCV(:,:,1),[],1))/sqrt(size(sfbCV(:,:,1),12));
tTimes2=0:10:600;

if is.plot
  set(gcf,'units','points','position',[x0,y0,width,height])
  plot(tTimes2, Mreplay, 'LineWidth', 1.5)
  xlabel('lag (ms)')
  ylabel('L1-reverse/Leftarrow sequenceness /RightarrowL1-forward'),
  title({'Sequenceness, Study II Lasso,',[' criteria maximised across all timepoints']});
  npThreshAll = max(max(abs(nanmean(sfbCV(:,:,2:end),1))));hold on;
  plot([tTimes2(1) tTimes2(end)], -npThreshAll*[1 1], 'k--'), plot([tTimes2(1) tTimes2(end)], npThreshAll*[1 1], 'k--')
  grid on;
end