%% NOTE: THIS IS NOT WORKING PROPERLY

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
         

    load([is.AnalysisPath,'classifiers/Sequence_by_Training_4Cell/','StimAll.mat']);
    nsub=1:length(is.goodsub);
    sfAll=squeeze(sfAll(2:end,:,nsub,1,:,:)); % timelag*nsubs*trainingtime*L1
    sbAll=squeeze(sbAll(2:end,:,nsub,1,:,:)); % timelag*nsubs*trainingtime*L1

    t_points=1:60%[3:5]; % which points to compute sequenceness over - should be all, 1:60
    for iL1=1:10
        for iTT=1:length(is.whichTimes)
        %subplot(2,5,iL1)
            sfb=squeeze(sfAll(:,:,iTT,iL1)-sbAll(:,:,iTT,iL1));

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
