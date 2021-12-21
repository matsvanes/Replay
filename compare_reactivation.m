% this scripts will load in the Rreds from study 2 (should be expanded to
% study 1) and compare the pre- and post- task resting state probabilities
% (top 1%).

whichstudy = 2;
if whichstudy==1
  MFStartup;
else
  MFStartup_studyII;
  load('/Volumes/T5_OHBA/data/replay/StrLearn_MEGexp/BehavData/Exptrial.mat')
  is.dirname='StrLearn_MEGexp/';
end
if ~isfolder(is.AnalysisPath),   mkdir(is.AnalysisPath), end
is.usePrecomputed = true;

% variables
is.whenInMS = is.tss*is.msPerSample;
is.TOI = 200; % 200ms post stimulus.
is.sessionsOfInterest = [1,2,3,4,8]; % rst, lci, lci, lci, rst
is.nStim = 8;
is.ntrlAll=120; % per session (includes upside down images)
is.ntrlValid=96; % per session
is.nSes = 3;
is.nChan = 273;
is.nTime = 401;
is.iRun = [1 2];
is.topPercentile=1; % define top 1% probabilities as reactivation

% get subject specific lambda
if is.usePrecomputed
  StimAll = load([is.AnalysisPath, 'classifiers/Sequence_by_Training_4Cell/StimAll.mat']);
  sfAll = StimAll.sfAll;
  sbAll = StimAll.sbAll;
else
  % compute sequenceness, only on iRun=2
  [sf, sb] = compute_sequenceness(iSj, is, 2);
  sfAll(:,:,iSj,:,:) = sf;
  sbAll(:,:,iSj,:,:) = sb;
end
[iLCV, sfbCV, Mreplay, Sreplay] = find_bestLambda(is, sfAll, sbAll);

% get probabilities
for iSj=1:length(is.goodsub)
  for i=is.iRun
    tmp = apply_toRest(iSj, is, [], [], i); % get all Rreds
    RredsAll{iSj,i} = tmp{1,37,1,iLCV(iSj)}; % get the one with the subject specific lambda
    prob =  1./(1+exp(RredsAll{iSj,i})); % the Rreds are the negative of how they usually are defined. That's why we're using +Rreds instead of minus (which is normally in the logsigmoid function)
    probAll{iSj,i} = prob;
    save([is.AnalysisPath, 'classifiers/TestResting4Cell/', sprintf('prob%d_%d',iSj,i)], 'prob')
  end
end


% find highest percentile probabilities (= reactivation)
n = size(probAll{1},1)*is.topPercentile/100;
val_top = zeros(length(is.goodsub),n,is.nStim,2);
for iSj=1:length(is.goodsub)
  for i=is.iRun
    x=probAll{iSj,i}; x(isnan(x))=-Inf;
    [val, idx] = sort(probAll{iSj,i}, 'descend', 'MissingPlacement', 'last');
    val_top(iSj,:,:,i) = val(1:n,:);
  end
end

% test the hypothesis that the top percentile probabilities in the
% post-task resting state data is on average higher than in the pre-task
% resting state.
x1 = reshape(val_top(:,:,:,1), 22,[]);
x2 = reshape(val_top(:,:,:,2), 22,[]);

x=squeeze(mean(val_top(:,:,:,2)-val_top(:,:,:,1),2));
[h,p,ci,stats] = ttest(mean(x2,2), mean(x1,2),'Tail','right');

%% plot
cmap = flipud(brewermap(2,'RdBu'));
f1 = figure;
subplot(1,2,1)
h2 = raincloud_plot(mean(x2,2), 'box_on', 1, 'color', cmap(2,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0);
h1 = raincloud_plot(mean(x1,2), 'box_on', 1, 'color', cmap(1,:), 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,...
    'box_col_match', 0);

legend([h1{1} h2{1}], {'pre-task RS', 'post-task RS'});
set(gca,'Ylim', [-4 8], 'Xlim', [0 0.8]);
xlabel('classification probability'); ylabel('Probability density')
box off

subplot(1,2,2)
plot([1,2],[mean(x1,2), mean(x2,2)],'-o','color', [.6 .6 .6]), 
xlim([0.5 2.5]), hold on, 
plot([1,2], [mean(x1(:)), mean(x2(:))], '-*k', 'LineWidth',2)
xticks([1,2])
xticklabels({'pre-task', 'post-task'})
ylabel('classification probability')

suptitle(sprintf('pre-/post- task difference in top percentile probabilities (reactivation) \n p = %d', p));